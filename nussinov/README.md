# nussinov: RNA folding on a GenDP-style PE chain

Nussinov is an interval DP over a triangular table. Each cell reads its whole row and column:
`N(i,j) = max(N(i+1,j-1) + d(i,j), max_k N(i,k) + N(k+1,j))`. That is a dependency pattern
GenDP doesn't cover (it has grid, graph, and 1D), and the reach is far beyond its 128-cell SPM range.
The compute needs only add and max, which GenDP's compute unit already has.

## Files

| File | What it is |
|---|---|
| `golden.py verify` | Two independent references, spec.md Sec. 14.2 vectors, and the PE-array mapping model checked on full tables |
| `golden.py perf` | Cycle model of the mapping: ops per cycle, utilization, Gops/s at an assumed clock |
| `cpu_baseline.c` | Single-thread C baseline with a vectorizable inner loop (transposed column copy) |

```
python3 golden.py verify
python3 golden.py perf -n 64 128 256 -P 4 16 64
cc -O3 -march=native -o cpu_baseline cpu_baseline.c && ./cpu_baseline data/seq_256.txt 3
```

## Mapping (simplified from spec.md Sec. 11.5)

Strips of P rows run bottom to top, with PE p owning row `i_bot - p`. For each column j one stream
flows down the chain: first `N(j,j) .. N(i_bot+1,j)` from memory, then each PE appends its own
finished cell. PE p consumes and forwards every word ahead of its own, so it sees exactly the
column it needs, in order. The model computes each cell only from what that PE can see (its row
store, the passing stream, and the previous column's last word).

## Verified

- All 9 conformance vectors in spec.md Sec. 14.2 pass, and the two references agree on full tables.
- The array model matches the reference on full tables: 100 random sequences x P in {1,2,4,8}.
- Negative controls: removing the pair term, or breaking the previous-column word, fails 100/100 checks.
- C baseline answers match Python (n = 64, 128, 256).

## First numbers

Array: cycle model, T_TAIL = 16, 256 cycles between strips, 2 GHz assumed (GenDP's assumption).
CPU: Apple M2 Pro, one core, `-O3 -march=native`, best of 3.

| n | CPU 1 core | 16 PEs | 64 PEs | Utilization, 16 PEs |
|---|---|---|---|---|
| 64 | 0.87 Gops/s | 9.9 | 13.7 | 31% |
| 128 | 1.40 | 17.4 | 35.1 | 55% |
| 256 | 2.91 | 23.8 | 66.0 | 74% |
| 512 | 4.29 | 27.7 | 93.2 | 87% |

Caveats to state in the paper:

- The CPU number is one core. The M2 Pro has 8 performance cores, and anti-diagonals of the table
  parallelize, so a multithreaded CPU could get near 16-PE throughput. The defensible win is
  throughput per mm² and per watt, GenDP's own metrics, not raw throughput.
- One row per PE, one lane. spec.md's mapping uses 2 rows per PE (both compute units) and 2 x 16-bit
  lanes, up to 4x more per PE. This model doesn't include that yet.
- T_TAIL (per-column overhead) is an estimate. Short sequences are dominated by it.
- No link back-pressure is modelled (deep FIFOs assumed; spec.md Sec. 12.1 shows the depth needed).

## Next

- Area: GenDP Table 7 gives 0.035 mm² per integer PE at 28 nm. Add SRAM for row stores and column
  buffering, then compute Gops/s per mm² against the CPU core area scaled to the same node.
- Multithreaded CPU baseline, so the comparison is honest.
- 2 rows per PE and 2 lanes in the model.
