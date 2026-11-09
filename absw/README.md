# absw: adaptive banded Smith-Waterman on GenDP

GenDP fixes the DP band before a kernel starts (paper Sec. 7.6.2: "does not support adaptive
or dynamic band choice"). This folder adds an adaptive band: after each anti-diagonal, the
array moves the band right or down toward whichever end scores higher.

Run everything: `bash absw/run_all.sh` (about 1 minute).

## Files

| File | What it is |
|---|---|
| `golden/absw_golden.hpp` | **The golden model (C++17).** An exact port of `ksw_extend2`, plus an anti-diagonal engine with `full`, `rowband`, `static`, `adaptive` band policies. Bit-exact with the RTL |
| `golden/golden.cpp` | Command-line tool: `selftest`, `validate`, `compare`, `run`, `pairs` (synthetic data), `vectors` (RTL tests). Built by `make golden` into `build/golden` |
| `run_all.sh` | Builds the simulator and ksw-test, runs all checks and comparisons |
| `compat/` | Empty `emmintrin.h` and stub `omp.h` so `ksw-test` builds on Apple Silicon (no source edits) |
| `golden.py`, `gen_pairs.py`, `rtl/gen_vectors.py` | The original Python model. No longer used by any flow; kept as an independent cross-check |

The C++ model replaced the Python one after being shown identical to it: 7488/7488 pair-runs
match on score, gscore, best-cell position, and cell count (4 datasets x 4 policies x W in
{4,8,16,32}), and the `compare` tables match character for character. It is about 90x faster
(0.55 s vs 50 s for the long-read comparison). Its random generator is a fixed splitmix64, so
`pairs` and `vectors` produce the same data on macOS (clang) and Linux (gcc).

## What is verified

- GenDP BSW simulator vs ksw-test: 256/256 pairs match (score and gscore).
- C++ `ksw_extend2` port vs ksw-test: 256/256 identical on all 6 output fields.
- Anti-diagonal engine (`rowband`) vs the port with z-drop and pruning off: 768/768 identical.
- The engine keeps only three anti-diagonals (O(length) memory: 3.3 MB at 10 kb, down from
  1.9 GB for the earlier full-table version, and ~3x faster). Checked against the full-table
  version: 66,000 random runs (all policies, W from 1 to 33, empty to 400 bp, h0 from 0 to 100)
  plus 1/5/10 kb pairs match on score, gscore, best cell, cell count, and every band move.

## First results (synthetic data, not yet real reads)

`% optimal` = pairs where the banded score equals the full-DP score.

| Workload | static W=16 | adaptive W=16 | static W=32 |
|---|---|---|---|
| 100bp, 2% sub, 1% indel (GenDP's BSW regime) | 100% | 100% | 100% |
| 1kb, 4% sub, 6% indel, 65% deletions (drift) | 8.3% | 100% | 53.3% |

- These are the numbers `make golden` prints. The earlier Python-generated dataset with the same
  settings gave 15% / 100% / 60%; same conclusion, different random draw.
- Longer segments widen the gap (synthetic, 5% sub, 5% indel, unbiased indels, 30 pairs each):

  | Segment | static W=16 | adaptive W=16 | static W=32 |
  |---|---|---|---|
  | 1 kb | 70.0% | 100% | 100% |
  | 2 kb | 53.3% | 100% | 93.3% |
  | 5 kb | 13.3% | 100% | 66.7% |
  | 10 kb | 0.0% | 100% | 43.3% |

- Adaptive bands don't matter for short reads; they matter for long reads, where small indels
  add up to drift.
- Adaptive bands do **not** recover a single gap wider than the band. Until the new diagonal is
  inside the band, both end scores are the same, so there is no signal to steer on.
- `static` and `adaptive` use the same W-PE array. The only new hardware is the per-anti-diagonal
  compare of the two end scores and the 1-bit shift that follows.

## RTL (`rtl/`)

`absw_array.sv` is a W-PE anti-diagonal array that computes one anti-diagonal per cycle, plus a
band-steering unit. The `adaptive` input picks adaptive or static (alternating) steering, so one
design gives both sides of the comparison. Each PE talks only to its neighbours. The query enters
at PE 0 and the reference at PE W-1, and each shifts only on its own kind of move.

- `bash absw/rtl/run_tests.sh` (or `make sim`) builds with Verilator and checks every case
  against the C++ golden model for W = 16, 8, 32: score, gscore, best-cell position, and the exact cycle count
  (W + tlen + qlen).
- Status: 258/258 tests pass at each W, with sequences up to 256 bases (`make sim`), and also
  up to 1 kb and 4 kb (`MAXLEN=1024` / `4096 bash rtl/run_tests.sh`). Scores and indices are
  16-bit, so lengths up to ~32k fit.
- Real data: `make sim TSV=pairs.tsv SIM_MAXLEN=2048` turns any pairs file (same format as
  `golden compare` takes) into RTL tests. Pairs longer than `SIM_MAXLEN` are skipped and counted.
- Gate level: `make gatesim` (same `W`/`MAXLEN`/`CLK_PS` as the `make absw` run) simulates the
  Genus netlist with the ASAP7 cell models and runs the same tests. It was checked locally on a
  Yosys-synthesized ASAP7 netlist (26.5k cells), with every flip-flop starting at a random value:
  258/258 pass. Swapping 3 of its AND gates for OR gates fails the run. Three deliberately broken RTL variants (flipped steering
  compare, wrong diagonal mux, wrong gap-open penalty) each failed 52-112 of the 178 tests in an
  earlier, smaller vector set.
- In these 256-base tests, adaptive never scores below static. It scores higher in 30/129 cases
  at W=8 and 4/129 at W=16.
- Verilator is in the conda env `hw` (Homebrew here is the Intel build and doesn't work).

## Makefile

Run from `absw/`:

| Command | Does |
|---|---|
| `make absw` | Synthesizes the Smith-Waterman array in ASAP7 with Cadence Genus. Downloads the ASAP7 libraries on first run. Options: `W=16 MAXLEN=128 CLK_PS=1000 VT=RVT CORNER=TT` |
| `make sim` | RTL vs golden model at W = 8, 16, 32. Use `SIM=vcs` on CAEN if Verilator isn't there. `TSV=pairs.tsv` tests your own pairs (`SIM_MAXLEN=2048` for 1 kb reads, `N=200` to cap the count) |
| `make gatesim` | Synthesized netlist (from `make absw`, same knobs) vs golden model. `SIM=vcs` uses the real ASAP7 cell models; Verilator swaps in simple flip-flop models |
| `make golden` | Golden-model selftest plus the static vs adaptive comparison on long reads |
| `make baseline` | Builds the GenDP simulator and BWA `ksw-test` and checks them against each other |
| `make orfs ORFS=<path>` | Full RTL-to-GDSII with OpenROAD-flow-scripts, if you have it installed |
| `make libs` | Re-fetches only the ASAP7 liberty and LEF files |

Results land in `flow/genus/W16_L128_1000ps_RVT_TT/reports/`: area, timing, power, gates, QoR.
To find the maximum clock, lower `CLK_PS` until the slack in `timing.rpt` goes negative.

### On CAEN

```
git clone -b adaptive-dp https://github.com/viraj-shahh/AdaptiveGenDP.git && cd AdaptiveGenDP/absw
module avail 2>&1 | grep -i -E "genus|vcs|verilator"   # find the module names, then module load them
make absw
```

### Not yet verified

- `make absw` (Genus), `SIM=vcs` (RTL and gate level), and `make orfs` haven't been run, because no Cadence, VCS, or OpenROAD
  is available on the Mac. Expect small fixes on the first CAEN run.
- The library download in `make absw`, `make sim` (Verilator, including `TSV=`), `make gatesim`
  (Verilator, on a Yosys netlist), `make golden`, and `make baseline` are tested.
- If Genus picks a cell the ORFS Verilog models don't have, `make gatesim` stops with a missing
  module error; send that output over.
- ASAP7 libraries checked: units are ps / fF / pW, and INVx1 area is 0.044 um², so these are real
  1x areas.
- The RTL uses synchronous reset, because ASAP7's only async-reset flop is a single-size set/reset
  cell that would inflate area.
- GDSII still needs place-and-route: either `make orfs`, or an Innovus script, which can be written
  once `make absw` runs.

## Gotchas found so far

- GenDP BSW runs in 8-bit lanes and saturates at 127. Inputs with `h0 + qlen > 127` mismatch
  ksw. Long reads need 16-bit lanes.
- `gendp/kernel/bwa-mem/ksw.cpp:377` (`ksw_extend2_new`) sets `w = qlen`, so the GenDP cycle
  model of BSW computes the full table for these reads, not a narrow band.
- The repo only implements BSW, PairHMM, POA, and Chain. DTW and Bellman-Ford appear in the paper
  but not in the code.

## Next

1. `make absw` on CAEN: clock, area, power. Then `make gatesim SIM=vcs`.
2. Place and route to GDSII.
3. Add real long reads from GenDP's own dataset (`gendp-datasets.tar.gz`, POA input: ONT reads)
   alongside the synthetic pairs.
