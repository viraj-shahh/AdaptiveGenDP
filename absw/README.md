# absw: adaptive banded Smith-Waterman on GenDP

GenDP fixes the DP band before a kernel starts (paper Sec. 7.6.2: "does not support adaptive
or dynamic band choice"). This folder adds an adaptive band: after each anti-diagonal, the
array moves the band right or down toward whichever end scores higher.

Run everything: `bash absw/run_all.sh` (about 1 minute).

## Files

| File | What it is |
|---|---|
| `gen_pairs.py` | Synthetic read/reference pairs in GenDP-sim, ksw-test, and TSV formats |
| `golden.py` | Golden model: exact port of `ksw_extend2`, plus an anti-diagonal engine with `full`, `rowband`, `static`, `adaptive` band policies |
| `run_all.sh` | Builds the simulator and ksw-test, runs all checks and comparisons |
| `compat/` | Empty `emmintrin.h` and stub `omp.h` so `ksw-test` builds on Apple Silicon (no source edits) |

## What is verified

- GenDP BSW simulator vs ksw-test: 256/256 pairs match (score and gscore).
- `ksw_extend2_port` vs ksw-test: 256/256 identical on all 6 output fields.
- Anti-diagonal engine (`rowband`) vs the port with z-drop and pruning off: 768/768 identical.

## First results (synthetic data, not yet real reads)

`% optimal` = pairs where the banded score equals the full-DP score.

| Workload | static W=16 | adaptive W=16 | static W=32 |
|---|---|---|---|
| 100bp, 2% sub, 1% indel (GenDP's BSW regime) | 100% | 100% | 100% |
| 1kb, 4% sub, 6% indel, 65% deletions (drift) | 15% | 100% | 60% |

- Adaptive bands don't matter for short reads; they matter for long reads, where small indels
  add up to drift.
- Adaptive bands do **not** recover a single gap wider than the band. Until the new diagonal is
  inside the band, both end scores are the same, so there is no signal to steer on.
- `static` and `adaptive` use the same W-PE array. The only new hardware is the per-anti-diagonal
  compare of the two end scores and the 1-bit shift that follows.

## Gotchas found so far

- GenDP BSW runs in 8-bit lanes and saturates at 127. Inputs with `h0 + qlen > 127` mismatch
  ksw. Long reads need 16-bit lanes.
- `gendp/kernel/bwa-mem/ksw.cpp:377` (`ksw_extend2_new`) sets `w = qlen`, so the GenDP cycle
  model of BSW computes the full table for these reads, not a narrow band.
- The repo only implements BSW, PairHMM, POA, and Chain. DTW and Bellman-Ford appear in the paper
  but not in the code.
- Verilator, Icarus, and Yosys are not installed on this machine yet.

## Next

1. Add real long reads from GenDP's own dataset (`gendp-datasets.tar.gz`, POA input: ONT reads)
   alongside the synthetic pairs.
2. Add adaptive BSW as kernel 5 in the GenDP simulator (band-shift decision at the array level),
   to get cycle counts.
3. RTL: a W-PE anti-diagonal array plus the band-steering unit, checked against `golden.py`.
