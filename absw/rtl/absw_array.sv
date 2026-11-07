// Adaptive banded Smith-Waterman (seed extension) array.
//
// W processing elements compute one anti-diagonal of the DP table per cycle.
// PE k holds row i = lo + k, column j = d - i.  After each anti-diagonal the band
// either moves "right" (lo unchanged) or "down" (lo + 1).  In adaptive mode the move
// goes toward whichever band end has the higher H score; in static mode moves
// alternate, which is a fixed diagonal band.
//
// Bit-exact with absw/golden.py extend(policy="adaptive" | "static", W=W).
//
// Recurrences (same as BWA ksw_extend2, see golden.py):
//   M(i,j) = H(i-1,j-1) > 0 ? H(i-1,j-1) + S(t_i, q_j) : 0
//   E(i,j) = i > 1 ? max(E(i-1,j) - e_del, M(i-1,j) - (o_del+e_del), 0) : 0
//   F(i,j) = j > 1 ? max(F(i,j-1) - e_ins, M(i,j-1) - (o_ins+e_ins), 0) : 0
//   H(i,j) = max(M, E, F)
//   row 0 / column 0 hold h0-derived boundary scores; cells outside the band are 0.
//
// Data movement (all neighbour-only):
//   s0 = move into this diagonal (1 = down), s1 = the move before it.
//   up   (i-1, j)   at d-1 : PE k + s0 - 1
//   left (i,   j-1) at d-1 : PE k + s0
//   diag (i-1, j-1) at d-2 : PE k + s0 + s1 - 1
//   query chars enter at PE 0 and shift toward PE W-1 on a right move;
//   target chars enter at PE W-1 and shift toward PE 0 on a down move.
//
// Reset is synchronous (rst_n sampled on clk): ASAP7 has only one async-reset flop.
//
// Usage: while idle, write target and query symbols (0..3 = ACGT, 4 = N) through
// wr_*; then pulse start with tlen, qlen, h0, adaptive valid.  done rises when the
// result is ready and stays high until the next start.

module absw_array #(
    parameter int W        = 16,     // PEs; must be a power of two (score reduction tree)
    parameter int SW       = 16,     // score width (signed)
    parameter int IW       = 16,     // index width (signed)
    parameter int MAXLEN   = 128,    // max sequence length held on chip
    parameter int MATCH    = 1,
    parameter int MISMATCH = 4,
    parameter int AMBIG    = 1,
    parameter int O_DEL    = 6,
    parameter int E_DEL    = 1,
    parameter int O_INS    = 6,
    parameter int E_INS    = 1
) (
    input  logic                 clk,
    input  logic                 rst_n,
    // sequence load (only while idle)
    input  logic                 wr_en,
    input  logic                 wr_sel,          // 0 = target (reference), 1 = query (read)
    input  logic [$clog2(MAXLEN)-1:0] wr_addr,    // 0-based position
    input  logic [2:0]           wr_data,
    // job
    input  logic                 start,
    input  logic [IW-1:0]        tlen,
    input  logic [IW-1:0]        qlen,
    input  logic signed [SW-1:0] h0,
    input  logic                 adaptive,
    output logic                 busy,
    output logic                 done,
    output logic signed [SW-1:0] score,
    output logic signed [SW-1:0] gscore,
    output logic signed [IW-1:0] max_i,
    output logic signed [IW-1:0] max_j,
    output logic [31:0]          cycles
);

    localparam logic [2:0] SYM_NONE = 3'd7;
    localparam int LO0 = 1 - (W - 1) / 2;          // band centred on cell (1,1)
    localparam int LOGW = $clog2(W);

    typedef enum logic [2:0] {S_IDLE, S_INIT, S_RUN, S_FLUSH, S_DONE} state_t;
    state_t state;

    // ---------------- sequence buffers ----------------
    logic [2:0] tbuf [MAXLEN];
    logic [2:0] qbuf [MAXLEN];

    always_ff @(posedge clk) begin
        if (wr_en && state != S_RUN && state != S_INIT) begin
            if (wr_sel) qbuf[wr_addr] <= wr_data;
            else        tbuf[wr_addr] <= wr_data;
        end
    end

    // ---------------- job registers ----------------
    logic signed [IW-1:0] tlen_r, qlen_r, lo, d;
    logic signed [SW-1:0] h0_r;
    logic                 adapt_r, s0, s1, last_dir;
    logic [LOGW:0]        init_cnt;

    // ---------------- per-PE state ----------------
    logic signed [SW-1:0] H1 [W];   // H at d-1
    logic signed [SW-1:0] H2 [W];   // H at d-2
    logic signed [SW-1:0] E1 [W];
    logic signed [SW-1:0] F1 [W];
    logic signed [SW-1:0] M1 [W];
    logic [2:0]           tch [W];
    logic [2:0]           qch [W];

    // ---------------- per-PE combinational results ----------------
    logic signed [SW-1:0] hn [W], en [W], fn [W], mn [W];
    logic                 vld [W];

    function automatic logic signed [SW-1:0] smax(input logic signed [SW-1:0] a,
                                                   input logic signed [SW-1:0] b);
        return (a > b) ? a : b;
    endfunction

    function automatic logic signed [SW-1:0] subst(input logic [2:0] a, input logic [2:0] b);
        if (a == 3'd4 || b == 3'd4) return SW'(-AMBIG);
        return (a == b) ? SW'(MATCH) : SW'(-MISMATCH);
    endfunction

    // boundary score at distance x from the origin along row 0 (ins) or column 0 (del)
    function automatic logic signed [SW-1:0] bnd(input logic signed [IW-1:0] x,
                                                  input logic signed [SW-1:0] hz,
                                                  input logic signed [SW-1:0] o,
                                                  input logic signed [SW-1:0] e);
        logic signed [SW-1:0] v;
        if (x == 0) return hz;
        v = hz - o - e * SW'(x);
        return (v > 0) ? v : '0;
    endfunction

    genvar k;
    generate
        for (k = 0; k < W; k++) begin : g_pe
            localparam int KU = (k > 0)     ? k - 1 : 0;
            localparam int KL = (k < W - 1) ? k + 1 : W - 1;
            logic signed [IW-1:0] i, j;
            logic signed [SW-1:0] e_up, m_up, f_left, m_left, hd, s, m, e, f, h;

            always_comb begin
                i = lo + IW'(k);
                j = d - i;
                vld[k] = (i >= 1) && (i <= tlen_r) && (j >= 1) && (j <= qlen_r);

                // up (i-1, j): PE k-1 after a right move, PE k after a down move
                if (s0)          begin e_up = E1[k];  m_up = M1[k];  end
                else if (k > 0)  begin e_up = E1[KU]; m_up = M1[KU]; end
                else             begin e_up = '0;     m_up = '0;     end

                // left (i, j-1): PE k after a right move, PE k+1 after a down move
                if (!s0)           begin f_left = F1[k];  m_left = M1[k];  end
                else if (k < W-1)  begin f_left = F1[KL]; m_left = M1[KL]; end
                else               begin f_left = '0;     m_left = '0;     end

                // diag (i-1, j-1) at d-2: PE k + s0 + s1 - 1
                case ({s0, s1})
                    2'b00:          hd = (k > 0)     ? H2[KU] : '0;
                    2'b01, 2'b10:   hd = H2[k];
                    default:        hd = (k < W - 1) ? H2[KL] : '0;
                endcase
                if (i == 1)      hd = bnd(j - 1, h0_r, SW'(O_INS), SW'(E_INS));   // row 0
                else if (j == 1) hd = bnd(i - 1, h0_r, SW'(O_DEL), SW'(E_DEL));   // column 0

                s = subst(tch[k], qch[k]);
                m = (hd > 0) ? hd + s : '0;
                e = (i > 1) ? smax(e_up - SW'(E_DEL), smax(m_up - SW'(O_DEL + E_DEL), '0)) : '0;
                f = (j > 1) ? smax(f_left - SW'(E_INS), smax(m_left - SW'(O_INS + E_INS), '0)) : '0;
                h = smax(m, smax(e, f));

                hn[k] = vld[k] ? h : '0;
                en[k] = vld[k] ? e : '0;
                fn[k] = vld[k] ? f : '0;
                mn[k] = vld[k] ? m : '0;
            end
        end
    endgenerate

    // ---------------- band steering ----------------
    logic signed [SW-1:0] h_top, h_bot;
    logic                 down;

    always_comb begin
        h_top = vld[0]     ? hn[0]     : SW'(-1);
        h_bot = vld[W - 1] ? hn[W - 1] : SW'(-1);
        if (adapt_r && h_top != h_bot) down = (h_bot > h_top);
        else                           down = !last_dir;
        if (lo + 1 > tlen_r)                          down = 1'b0;
        if ((d + 1) - (lo + IW'(W - 1)) > qlen_r)     down = 1'b1;
    end

    // ---------------- symbol injection (one read port per buffer) ----------------
    logic signed [IW-1:0] t_idx, q_idx;
    logic [2:0]           t_new, q_new;

    always_comb begin
        if (state == S_INIT) begin
            t_idx = IW'(LO0) + IW'(init_cnt);
            q_idx = IW'(2 - LO0 - (W - 1)) + IW'(init_cnt);
        end else begin
            t_idx = lo + IW'(W);          // new bottom row after a down move
            q_idx = d + 1 - lo;           // new PE 0 column after a right move
        end
        t_new = (t_idx >= 1 && t_idx <= tlen_r) ? tbuf[t_idx - 1] : SYM_NONE;
        q_new = (q_idx >= 1 && q_idx <= qlen_r) ? qbuf[q_idx - 1] : SYM_NONE;
    end

    // ---------------- main sequencer ----------------
    always_ff @(posedge clk) begin
        if (!rst_n) begin
            state <= S_IDLE;
            lo <= '0; d <= '0; s0 <= 1'b0; s1 <= 1'b0; last_dir <= 1'b1;
            init_cnt <= '0;
            tlen_r <= '0; qlen_r <= '0; h0_r <= '0; adapt_r <= 1'b0;
            for (int x = 0; x < W; x++) begin
                H1[x] <= '0; H2[x] <= '0; E1[x] <= '0; F1[x] <= '0; M1[x] <= '0;
                tch[x] <= SYM_NONE; qch[x] <= SYM_NONE;
            end
        end else begin
            case (state)
                S_IDLE, S_DONE: if (start) begin
                    tlen_r <= tlen; qlen_r <= qlen; h0_r <= h0; adapt_r <= adaptive;
                    lo <= IW'(LO0); d <= IW'(2);
                    s0 <= 1'b0; s1 <= 1'b0; last_dir <= 1'b1;
                    init_cnt <= '0;
                    for (int x = 0; x < W; x++) begin
                        H1[x] <= '0; H2[x] <= '0; E1[x] <= '0; F1[x] <= '0; M1[x] <= '0;
                    end
                    state <= S_INIT;
                end
                S_INIT: begin
                    // fill: PE k ends with t[LO0 + k] and q[(2 - LO0) - k]
                    for (int x = 0; x < W - 1; x++) tch[x] <= tch[x + 1];
                    tch[W - 1] <= t_new;
                    for (int x = W - 1; x > 0; x--) qch[x] <= qch[x - 1];
                    qch[0] <= q_new;
                    init_cnt <= init_cnt + 1'b1;
                    if (init_cnt == (LOGW+1)'(W - 1)) state <= S_RUN;
                end
                S_RUN: begin
                    for (int x = 0; x < W; x++) begin
                        H2[x] <= H1[x]; H1[x] <= hn[x];
                        E1[x] <= en[x]; F1[x] <= fn[x]; M1[x] <= mn[x];
                    end
                    if (down) begin
                        for (int x = 0; x < W - 1; x++) tch[x] <= tch[x + 1];
                        tch[W - 1] <= t_new;
                    end else begin
                        for (int x = W - 1; x > 0; x--) qch[x] <= qch[x - 1];
                        qch[0] <= q_new;
                    end
                    s1 <= s0; s0 <= down; last_dir <= down;
                    lo <= lo + IW'(down);
                    d  <= d + 1;
                    if (d == tlen_r + qlen_r) state <= S_FLUSH;
                end
                S_FLUSH: state <= S_DONE;
                default: state <= S_IDLE;
            endcase
        end
    end

    // ---------------- score tracking, one stage behind ----------------
    logic signed [SW-1:0] hv [W];
    logic                 vv [W];
    logic signed [IW-1:0] lo_r, d_r;
    logic                 sv;

    always_ff @(posedge clk) begin
        if (!rst_n) begin
            sv <= 1'b0; lo_r <= '0; d_r <= '0;
            for (int x = 0; x < W; x++) begin hv[x] <= '0; vv[x] <= 1'b0; end
        end else begin
            sv <= (state == S_RUN);
            lo_r <= lo; d_r <= d;
            for (int x = 0; x < W; x++) begin hv[x] <= hn[x]; vv[x] <= vld[x]; end
        end
    end

    // max over valid PEs; ties go to the lowest k (smallest i), as numpy argmax does
    logic signed [SW-1:0] tv [2*W];
    logic                 tf [2*W];
    logic [LOGW-1:0]      tk [2*W];
    logic signed [SW-1:0] g_val;
    logic                 g_hit;

    always_comb begin
        for (int x = 0; x < W; x++) begin
            tv[W + x] = hv[x]; tf[W + x] = vv[x]; tk[W + x] = LOGW'(x);
        end
        for (int x = W - 1; x >= 1; x--) begin
            if (tf[2*x] && (!tf[2*x+1] || tv[2*x] >= tv[2*x+1])) begin
                tv[x] = tv[2*x];   tf[x] = 1'b1;       tk[x] = tk[2*x];
            end else begin
                tv[x] = tv[2*x+1]; tf[x] = tf[2*x+1];  tk[x] = tk[2*x+1];
            end
        end
        tv[0] = '0; tf[0] = 1'b0; tk[0] = '0;

        // the one valid cell (if any) in the last query column, for gscore
        g_val = '0; g_hit = 1'b0;
        for (int x = 0; x < W; x++)
            if (vv[x] && (d_r - (lo_r + IW'(x)) == qlen_r)) begin g_val = hv[x]; g_hit = 1'b1; end
    end

    always_ff @(posedge clk) begin
        if (!rst_n) begin
            score <= '0; gscore <= '0; max_i <= '0; max_j <= '0;
        end else if ((state == S_IDLE || state == S_DONE) && start) begin
            score <= h0; gscore <= SW'(-1); max_i <= '0; max_j <= '0;
        end else if (sv) begin
            if (tf[1] && tv[1] > score) begin
                score <= tv[1];
                max_i <= lo_r + IW'(tk[1]);
                max_j <= d_r - (lo_r + IW'(tk[1]));
            end
            if (g_hit && g_val > gscore) gscore <= g_val;
        end
    end

    // ---------------- status ----------------
    assign busy = (state == S_INIT) || (state == S_RUN) || (state == S_FLUSH);
    assign done = (state == S_DONE);

    always_ff @(posedge clk) begin
        if (!rst_n)                                              cycles <= '0;
        else if ((state == S_IDLE || state == S_DONE) && start)  cycles <= '0;
        else if (busy)                                           cycles <= cycles + 1;
    end

endmodule
