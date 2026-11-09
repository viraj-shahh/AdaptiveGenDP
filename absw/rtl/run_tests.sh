#!/bin/bash
# Build absw_array and check it against the C++ golden model for several band widths.
# Needs a C++17 compiler, plus a simulator:
#   SIM=verilator (default; found on PATH or in the conda env "hw")
#   SIM=vcs       (Synopsys VCS, e.g. on CAEN)
# Knobs (environment variables):
#   WS="16 8 32"   band widths to test
#   MAXLEN=256     on-chip sequence buffer size (longest sequence tested)
#   TSV=pairs.tsv  test these pairs (query<TAB>target<TAB>h0, e.g. real read segments)
#                  instead of the synthetic set; N=500 keeps the first 500 that fit
#   NETLIST=f.v    gate-level: simulate a synthesized netlist instead of the RTL.
#                  WS and MAXLEN must match the synthesis run. CELLS= points at the
#                  ASAP7 cell Verilog (default: the copy fetched by `make libs`).
set -euo pipefail
SIM=${SIM:-verilator}
if [ "$SIM" = verilator ]; then
    command -v verilator >/dev/null || export PATH="$PATH:/opt/anaconda3/envs/hw/bin"
fi
NETLIST=${NETLIST:+$(cd "$(dirname "$NETLIST")" && pwd)/$(basename "$NETLIST")}
TSV=${TSV:+$(cd "$(dirname "$TSV")" && pwd)/$(basename "$TSV")}
cd "$(dirname "$0")"
MAXLEN=${MAXLEN:-256}
CELLS=${CELLS:-../flow/.asap7_src/flow/platforms/asap7/verilog/stdcell}
# GNU make can't build in a path with spaces ("EECS 570"), so compile outside the repo.
OBJ=${OBJ:-"${TMPDIR:-/tmp}/absw_rtl_obj"}
mkdir -p build "$OBJ"
GOLDEN=../build/golden
mkdir -p ../build
c++ -std=c++17 -O2 -Wall -o "$GOLDEN" ../golden/golden.cpp

TAG=rtl
DESIGN=(absw_array.sv)
VCS_GATE=()
VL_GATE=()
if [ -n "$NETLIST" ]; then
    TAG=gate
    [ -f "$NETLIST" ] || { echo "no netlist at $NETLIST"; exit 1; }
    [ -d "$CELLS" ] || { echo "no ASAP7 cell models at $CELLS (run: make libs)"; exit 1; }
    # Top module = the one no other module instantiates (Genus may emit helper modules).
    TOP=${GATE_TOP:-$(awk '
        $1 == "module" { sub(/\(.*/, "", $2); mods[$2] = 1; next }
        $1 ~ /^[A-Za-z_][A-Za-z0-9_$]*$/ && $2 ~ /^[\\A-Za-z_]/ && $1 != "input" && $1 != "output" \
            && $1 != "wire" && $1 != "inout" && $1 != "assign" && $1 != "reg" { used[$1] = 1 }
        END { for (m in mods) if (!(m in used)) print m }' "$NETLIST" | head -1)}
    echo "gate-level: $NETLIST (top module $TOP)"
    DESIGN=("$NETLIST")
    VCS_GATE=(+define+GATE_TOP="$TOP" +delay_mode_zero +notimingcheck +nospecify)
    for f in "$CELLS"/asap7sc7p5t_*_RVT_TT_*.v; do VCS_GATE+=(-v "$f"); done
    # Verilator: real combinational cells (UDP blocks stripped; no combinational cell uses
    # them), behavioral flip-flops in place of the UDP-based SEQ library.
    mkdir -p "$OBJ/cells"
    VL_GATE=(-DGATE_TOP="$TOP" -Wno-lint -Wno-style
             --x-initial unique asap7_seq_verilator.v)
    for f in "$CELLS"/asap7sc7p5t_*_RVT_TT_*.v; do
        case $f in *SEQ*) continue;; esac
        awk '/^primitive/{skip=1} !skip{print} /^endprimitive/{skip=0}' "$f" > "$OBJ/cells/$(basename "$f")"
        VL_GATE+=("$OBJ/cells/$(basename "$f")")
    done
fi

FAILED=0
for W in ${WS:-16 8 32}; do
    VEC="build/vec_${TAG}_w$W.txt"
    if [ -n "$TSV" ]; then
        "$GOLDEN" vectors -W "$W" --maxlen "$MAXLEN" --tsv "$TSV" ${N:+-n "$N"} -o "$VEC"
    else
        "$GOLDEN" vectors -W "$W" --maxlen "$MAXLEN" -o "$VEC"
    fi
    B="$OBJ/${TAG}_w$W"
    mkdir -p "$B"
    if [ "$SIM" = vcs ]; then
        vcs -full64 -sverilog -timescale=1ns/1ps +lint=none \
            -pvalue+absw_tb.W="$W" -pvalue+absw_tb.MAXLEN="$MAXLEN" -top absw_tb ${VCS_GATE[@]+"${VCS_GATE[@]}"} \
            -Mdir="$B/csrc" -o "$B/sim" "${DESIGN[@]}" absw_tb.sv > "build/build_${TAG}_w$W.log" 2>&1 \
            || { tail -30 "build/build_${TAG}_w$W.log"; exit 1; }
    else
        verilator --binary --timing --timescale 1ns/1ps -Wno-fatal -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND \
            -GW="$W" -GMAXLEN="$MAXLEN" --top-module absw_tb ${VL_GATE[@]+"${VL_GATE[@]}"} \
            -Mdir "$B" -o sim "${DESIGN[@]}" absw_tb.sv > "build/build_${TAG}_w$W.log" 2>&1 \
            || { tail -30 "build/build_${TAG}_w$W.log"; exit 1; }
    fi
    # +verilator+rand+reset+2: flops start random, so a design that relies on power-up
    # values fails here too (VCS starts them at X).
    "$B/sim" +vec="$VEC" +verilator+rand+reset+2 | tee "build/run_${TAG}_w$W.log" | grep -E "FAIL|pass|RESULT"
    grep -q "RESULT: PASS" "build/run_${TAG}_w$W.log" || FAILED=1
done
exit $FAILED
