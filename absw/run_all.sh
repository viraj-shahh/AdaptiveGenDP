#!/bin/bash
# Reproduce everything in absw/ from a clean checkout (macOS arm64 or Linux).
#   1. build the GenDP simulator and the ksw reference
#   2. run GenDP's BSW kernel on synthetic pairs and check it against ksw
#   3. validate the Python golden model against ksw
#   4. compare static vs adaptive bands on short and long reads
set -euo pipefail
ABSW="$(cd "$(dirname "$0")" && pwd)"
GENDP="$ABSW/../gendp"
D="$ABSW/data"
mkdir -p "$D"

echo "== build simulator and ksw-test"
# Upstream flags use -static-libasan (gcc only) and OpenMP; neither is needed here.
make -C "$GENDP" -s CXXFLAGS="-g -w -fPIC -std=c++11 -O3 -fpermissive" LINKFLAGS="" sim
make -C "$GENDP/kernel/bwa-mem" -s CXX=c++ CXXFLAGS="-O3 -std=c++11 -I$ABSW/compat -I$GENDP -w"

echo "== GenDP BSW simulator vs ksw (256 short pairs, h0=20 keeps scores in 8-bit range)"
python3 "$ABSW/gen_pairs.py" -n 256 --h0 20 --seed 2 \
    --gendp "$D/smoke_gendp.txt" --ksw "$D/smoke_ksw.csv" --tsv "$D/smoke.tsv"
"$GENDP/kernel/bwa-mem/ksw-test" -i "$D/smoke_ksw.csv" -o "$D/smoke_ksw_out.txt" -x -n 256 > /dev/null
(cd "$GENDP" && mkdir -p bsw_sim_results instructions/bsw && python3 scripts/bsw_instruction_generator.py \
    && ./sim -k 1 -i "$D/smoke_gendp.txt" -o "$D/smoke_sim_out.txt" -s -n -1 > "$D/smoke_sim_result.txt" 2>&1)
BAD=$(cd "$GENDP" && python3 scripts/bsw_check_correctness.py "$D/smoke_sim_out.txt" "$D/smoke_ksw_out.txt" | wc -l | tr -d ' ')
echo "simulator mismatches vs ksw: $BAD / 256"

echo "== golden model"
python3 "$ABSW/golden.py" selftest
python3 "$ABSW/golden.py" validate "$D/smoke_ksw.csv" "$D/smoke_ksw_out.txt"

echo "== band policies"
python3 "$ABSW/gen_pairs.py" -n 200 --qlen 100 --h0 20 --seed 3 --tsv "$D/short.tsv" > /dev/null
python3 "$ABSW/gen_pairs.py" -n 60 --qlen 1000 --tlen-extra 200 --sub-rate 0.04 --indel-rate 0.06 \
    --max-indel 2 --del-frac 0.65 --h0 30 --max-len 0 --seed 4 --tsv "$D/long_drift.tsv" > /dev/null
for f in short long_drift; do
    echo "-- $f"
    python3 "$ABSW/golden.py" compare "$D/$f.tsv" -W 8 16 32
done
