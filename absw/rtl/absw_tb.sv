// Testbench: runs every case in a vector file (golden vectors) through absw_array
// and checks score, gscore, max_i, max_j against the golden model, plus the cycle
// count (W fill + tlen + qlen - 1 anti-diagonals + 1 flush).
//
//   +vec=<file>   vector file (default vectors.txt)
module absw_tb;
    parameter int W      = 16;
    parameter int MAXLEN = 256;
    localparam int SW = 16, IW = 16;

    logic clk = 1'b0, rst_n = 1'b0;
    always #5 clk = ~clk;

    logic                        wr_en = 1'b0, wr_sel = 1'b0, start = 1'b0, adaptive = 1'b0;
    logic [$clog2(MAXLEN)-1:0]   wr_addr = '0;
    logic [2:0]                  wr_data = '0;
    logic [IW-1:0]               tlen = '0, qlen = '0;
    logic signed [SW-1:0]        h0 = '0;
    logic                        busy, done;
    logic signed [SW-1:0]        score, gscore;
    logic signed [IW-1:0]        max_i, max_j;
    logic [31:0]                 cycles;

    // Gate-level: +define+GATE_TOP=<netlist module>; W and MAXLEN must match the synthesis run.
`ifdef GATE_TOP
    `GATE_TOP dut (
`else
    absw_array #(.W(W), .MAXLEN(MAXLEN)) dut (
`endif
        .clk, .rst_n, .wr_en, .wr_sel, .wr_addr, .wr_data,
        .start, .tlen, .qlen, .h0, .adaptive,
        .busy, .done, .score, .gscore, .max_i, .max_j, .cycles
    );

    initial begin
        string fn;
        int fd, n, tl, ql, h0v, ad, es, eg, ei, ej, sym, x, t, pass, fail, rc;
        longint total;
        pass = 0; fail = 0; total = 0;
        if (!$value$plusargs("vec=%s", fn)) fn = "vectors.txt";
        fd = $fopen(fn, "r");
        if (fd == 0) begin $display("cannot open %s", fn); $finish; end
        rc = $fscanf(fd, "%d", n);

        repeat (3) @(negedge clk);
        rst_n = 1'b1;

        for (t = 0; t < n; t++) begin
            rc = $fscanf(fd, "%d %d %d %d %d %d %d %d", tl, ql, h0v, ad, es, eg, ei, ej);
            for (x = 0; x < tl; x++) begin
                rc = $fscanf(fd, "%d", sym);
                @(negedge clk); wr_en = 1'b1; wr_sel = 1'b0; wr_addr = x[$clog2(MAXLEN)-1:0]; wr_data = sym[2:0];
            end
            for (x = 0; x < ql; x++) begin
                rc = $fscanf(fd, "%d", sym);
                @(negedge clk); wr_en = 1'b1; wr_sel = 1'b1; wr_addr = x[$clog2(MAXLEN)-1:0]; wr_data = sym[2:0];
            end
            @(negedge clk);
            wr_en = 1'b0;
            tlen = tl[IW-1:0]; qlen = ql[IW-1:0]; h0 = h0v[SW-1:0]; adaptive = ad[0];
            start = 1'b1;
            @(negedge clk);
            start = 1'b0;
            while (!done) @(negedge clk);

            total += cycles;
            if (score == es[SW-1:0] && gscore == eg[SW-1:0] && max_i == ei[IW-1:0] && max_j == ej[IW-1:0]
                && cycles == 32'(W + tl + ql)) begin
                pass++;
            end else begin
                fail++;
                if (fail <= 10)
                    $display("FAIL test %0d (tlen %0d qlen %0d h0 %0d %s): got score %0d gscore %0d at (%0d,%0d) in %0d cycles; expected %0d %0d at (%0d,%0d) in %0d",
                             t, tl, ql, h0v, ad ? "adaptive" : "static", score, gscore, max_i, max_j,
                             cycles, es, eg, ei, ej, W + tl + ql);
            end
        end
        $display("W=%0d MAXLEN=%0d: %0d/%0d tests pass, %0d total array cycles", W, MAXLEN, pass, n, total);
        if (fail != 0) $display("RESULT: FAIL");
        else           $display("RESULT: PASS");
        $finish;
    end
endmodule
