#!/bin/bash
# Build absw_array and check it against golden.py for several band widths.
# Needs python3 with numpy, plus a simulator:
#   SIM=verilator (default; found on PATH or in the conda env "hw")
#   SIM=vcs       (Synopsys VCS, e.g. on CAEN)
# WS="16 8 32" picks the band widths, MAXLEN the on-chip sequence buffer size.
set -euo pipefail
SIM=${SIM:-verilator}
if [ "$SIM" = verilator ]; then
    command -v verilator >/dev/null || export PATH="$PATH:/opt/anaconda3/envs/hw/bin"
fi
cd "$(dirname "$0")"
MAXLEN=${MAXLEN:-256}
# GNU make can't build in a path with spaces ("EECS 570"), so compile outside the repo.
OBJ=${OBJ:-"${TMPDIR:-/tmp}/absw_rtl_obj"}
mkdir -p build "$OBJ"
for W in ${WS:-16 8 32}; do
    python3 gen_vectors.py -W "$W" --maxlen "$MAXLEN" -o "build/vec_w$W.txt"
    mkdir -p "$OBJ/w$W"
    if [ "$SIM" = vcs ]; then
        vcs -full64 -sverilog -timescale=1ns/1ps +lint=none \
            -pvalue+absw_tb.W="$W" -pvalue+absw_tb.MAXLEN="$MAXLEN" -top absw_tb \
            -Mdir="$OBJ/w$W/csrc" -o "$OBJ/w$W/sim" absw_array.sv absw_tb.sv > "build/build_w$W.log" 2>&1 \
            || { tail -30 "build/build_w$W.log"; exit 1; }
    else
        verilator --binary --timing --timescale 1ns/1ps -Wno-fatal -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND \
            -GW="$W" -GMAXLEN="$MAXLEN" --top-module absw_tb \
            -Mdir "$OBJ/w$W" -o sim absw_array.sv absw_tb.sv > "build/build_w$W.log" 2>&1 \
            || { tail -30 "build/build_w$W.log"; exit 1; }
    fi
    "$OBJ/w$W/sim" +vec="build/vec_w$W.txt" | grep -E "FAIL|pass|RESULT"
done
