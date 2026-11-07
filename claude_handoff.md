# Handoff: EECS 570 project (GenDP extension)

Context for continuing this work in a new Claude chat. Written 2026-11-06.
Also attach: GenDP.pdf (ISCA 2023 paper), spec.md (an earlier, over-scoped GenDP-MT design),
"Additional DP Algos for 570 Project.pdf", and 570-1.pdf (original proposal, now outdated).

## Who and when

- Viraj Shah, EECS 570 (Fall 2026), team of 4. RTL due ~Nov 30, testing Dec 5, 6-page paper and
  poster due Dec 8.
- Viraj wants hands-on progress, not plans, timelines, or work splits. His teammates won't start
  until he does.

## Decisions made

- **Project = extend GenDP with new DP kernels.** The original proposal's minimap2 profiling,
  ksw2 baseline, chr20 reads, and k-mer models are **dropped**. Don't bring minimap2 back.
- **Main kernel: adaptive banded Smith-Waterman (ABSW).** GenDP Sec. 7.6.2 says it doesn't
  support adaptive or dynamic bands. The new hardware is one compare of the two band-end scores
  per anti-diagonal, plus a 1-bit right/down shift.
- **Second kernel: Nussinov RNA folding.** An interval DP over a triangular table, a dependency
  pattern GenDP doesn't have. Needs only add and max.
- **Dropped: Value Iteration and Policy Evaluation.**
  - They're really sparse matrix-vector products, not table DPs.
  - They'd need a convergence check, gather hardware, 1-per-cycle MAC, and floating point.
  - spec.md's own model projects only a modest win over one CPU core.
  - Deterministic Value Iteration is just Bellman-Ford, which GenDP already claims.
- **spec.md is a reference, not the plan.** It's far too large for the time available
  (multithreading, AXI, DMA, six kernels, full programmable RTL).

## Repo state (github GenDP artifact, local at ~/Downloads/EECS 570/GenDP)

Untracked new folders, nothing committed:

- `absw/`: `gen_pairs.py` (synthetic pairs), `golden.py` (exact Python port of BWA `ksw_extend2`
  plus an anti-diagonal engine with full / rowband / static / adaptive band policies),
  `run_all.sh` (one-command repro), `compat/` (stub `omp.h` and `emmintrin.h` so `ksw-test`
  builds on Apple Silicon), `README.md`.
- `nussinov/`: `golden.py` (two references, spec.md conformance vectors, strip/column-stream
  PE-array mapping model with cycle counts), `cpu_baseline.c`, `README.md`.

## Verified results

- The GenDP simulator builds on the M2 Mac. Upstream Makefile flags break with clang; the working
  flags are in `absw/run_all.sh`. GenDP's BSW kernel matches BWA `ksw-test` on 256/256 synthetic
  pairs.
- `ksw_extend2` port vs `ksw-test`: 256/256 identical on all 6 fields. The anti-diagonal engine
  matches the port on 768/768 cases.
- **ABSW accuracy** (% of pairs where the banded score equals full DP):

| Workload | static W=16 | adaptive W=16 | static W=32 |
|---|---|---|---|
| 100bp short reads (GenDP's BSW regime) | 100% | 100% | 100% |
| 1kb long reads, 6% indels, deletion-biased drift | 15% | 100% | 60% |

  - Adaptive bands follow gradual drift but can't recover one gap wider than the band.
- **Nussinov**:
  - All 9 spec.md vectors pass.
  - The array model matches the reference on full tables (100 random sequences x 4 PE counts).
    Negative controls fail as they should.
  - Throughput: array is the cycle model at 2 GHz assumed; CPU is one M2 Pro core.

| n | CPU 1 core | 16 PEs | 64 PEs |
|---|---|---|---|
| 256 | 2.9 Gops/s | 23.8 | 66.0 |
| 512 | 4.3 Gops/s | 27.7 | 93.2 |

  - Caveats: the CPU is one core, so a multithreaded baseline is still needed. The real claim
    should be per mm² and per watt. The model uses 1 row per PE and 1 lane; spec.md packs up to
    4x more per PE.

## Gotchas found

- GenDP BSW uses 8-bit lanes and saturates at 127. Long reads need 16-bit lanes.
- `gendp/kernel/bwa-mem/ksw.cpp:377` (`ksw_extend2_new`) sets `w = qlen`, so GenDP's BSW
  baseline effectively computes the full table.
- The repo only implements BSW, PairHMM, POA, and Chain. DTW and Bellman-Ford are in the paper
  but not the code.
- Verilator, Icarus, and Yosys aren't installed locally. Genus is available on CAEN.

## How to measure the value of the new kernels

- New kernels don't speed up GenDP's existing ones. They add capability.
- **ABSW:** compare cycles per alignment against static BSW at equal accuracy.
- **Nussinov:** compare throughput per mm² against a CPU baseline.
- Area comes from GenDP Table 7 (0.035 mm² per integer PE, 28 nm) plus synthesis of new units.
- Use Python cycle models rather than hand-writing GenDP simulator instructions.

## RTL plan

- A full programmable GenDP PE in RTL isn't feasible in the time left.
- Build a fixed-function ABSW anti-diagonal array plus the band-steering unit, about 1-1.5 weeks.
  Verify in Verilator against `absw/golden.py`, synthesize for clock and area.
- Risk: the steering decision may lengthen the critical path. Test a 1-step-delayed steering
  decision in the golden model first.
- Nussinov stays model-level only.

## Open next steps

1. Nussinov: area estimate, a multithreaded CPU baseline, and 2 rows per PE / 2 lanes in the model.
2. ABSW: cycle model, the delayed-steering experiment, and real long reads from GenDP's own
   dataset tarball (POA input has ONT reads).
3. ABSW RTL plus a Verilator testbench against the golden model.
