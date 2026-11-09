// Behavioral stand-ins for the ASAP7 flip-flops, used only for gate-level simulation
// with Verilator, which can't run the UDP tables in asap7sc7p5t_SEQ_*.v (same workaround
// as ORFS's platforms/asap7/verilog/stdcell/dff.v). VCS uses the real cell models.
// A netlist that uses a sequential cell not listed here (scan, latch, clock gate)
// fails to build with a "cannot find module" error.

`define ABSW_DFF(name, edge_, out, inv) \
module name (out, D, CLK); \
    output reg out; \
    input D, CLK; \
    always @(edge_ CLK) out <= inv D; \
endmodule

`ABSW_DFF(DFFHQNx1_ASAP7_75t_R, posedge, QN, ~)
`ABSW_DFF(DFFHQNx2_ASAP7_75t_R, posedge, QN, ~)
`ABSW_DFF(DFFHQNx3_ASAP7_75t_R, posedge, QN, ~)
`ABSW_DFF(DFFHQx4_ASAP7_75t_R,  posedge, Q,   )
`ABSW_DFF(DFFLQNx1_ASAP7_75t_R, negedge, QN, ~)
`ABSW_DFF(DFFLQNx2_ASAP7_75t_R, negedge, QN, ~)
`ABSW_DFF(DFFLQNx3_ASAP7_75t_R, negedge, QN, ~)
`ABSW_DFF(DFFLQx4_ASAP7_75t_R,  negedge, Q,   )

`undef ABSW_DFF
