#!/usr/bin/env python3
"""Nussinov RNA folding: golden model and a cycle model of a linear PE-array mapping.

Recurrence (spec.md Sec. 2.2), 0-based, N(i,j) = 0 for j <= i:

    N(i,j) = max( N(i+1,j-1) + d(i,j),  max_{i <= k < j} N(i,k) + N(k+1,j) )
    d(i,j) = s(x_i, x_j) if j - i > L else 0         s = 1 for AU, GC, GU pairs

Three layers, each checked against the one before:

  ref_bifurcation()  numpy implementation of the recurrence above
  ref_classic()      the textbook 4-case form (independent check)
  ArrayModel         the strip / column-stream mapping on a chain of P PEs
                     (spec.md Sec. 11.5, simplified to one row per PE and one lane).
                     It computes the table only from what each PE can see (its own
                     row store, the words streaming past, and the previous column's
                     last word), so a full-table match shows the dataflow is
                     correct, not just the math.  It also counts cycles.

Cycle model assumptions (state them in the paper):
  * a PE does one add-max per cycle on a stream word (one row per PE, one lane)
  * after its last word in a column a PE spends T_TAIL cycles finishing the cell and
    emitting it; a word hops one PE per cycle
  * the head streams 1 word per cycle from L1; link FIFOs are deep enough that they
    never back-pressure (spec.md Sec. 12.1 shows depth >= T_TAIL + 2 achieves this)
  * strips run back to back with STRIP_GAP cycles between them
"""
import argparse
import sys
import time

import numpy as np

PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
NEG = -(1 << 30)                       # neutral value: loses every max


def pair_matrix(seq, L):
    n = len(seq)
    d = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + L + 1, n):
            if (seq[i], seq[j]) in PAIRS:
                d[i, j] = 1
    return d


def ref_bifurcation(seq, L=0):
    n = len(seq)
    d = pair_matrix(seq, L)
    N = np.zeros((n + 1, n + 1), dtype=np.int64)    # extra row/col so N[i+1, j] exists
    for span in range(1, n):
        for i in range(0, n - span):
            j = i + span
            inner = N[i + 1, j - 1] if i + 1 <= j - 1 else 0
            split = int(np.max(N[i, i:j] + N[i + 1:j + 1, j]))
            N[i, j] = max(inner + d[i, j], split)
    return N[:n, :n]


def ref_classic(seq, L=0):
    n = len(seq)
    d = pair_matrix(seq, L)
    N = [[0] * n for _ in range(n)]
    for span in range(1, n):
        for i in range(0, n - span):
            j = i + span
            best = max(N[i + 1][j], N[i][j - 1])
            best = max(best, (N[i + 1][j - 1] if i + 1 <= j - 1 else 0) + d[i, j])
            for k in range(i + 1, j):
                best = max(best, N[i][k] + N[k + 1][j])
            N[i][j] = best
    return np.array(N, dtype=np.int64)


def useful_ops(n):
    """Add-max steps in the split term: sum over cells of (j - i)."""
    return n * (n * n - 1) // 6


class ArrayModel:
    """Strip / column-stream mapping on P PEs, one row per PE.

    Strips of P rows run bottom to top.  PE p owns row i = i_bot - p, so PE0 holds the
    bottom row of the strip.  For column j a single stream flows down the chain:
      N(j, j), N(j-1, j), ..., N(i_bot+1, j)    from memory (rows below the strip)
      then PE0 appends N(i_bot, j), PE1 appends N(i_bot-1, j), ...
    PE p consumes every word before its own and forwards it, so it sees N(m, j) for
    every m in (i, j] in descending order, which is exactly the split term's column.
    """

    def __init__(self, seq, L, P, t_tail=16, strip_gap=256):
        self.seq, self.L, self.P = seq, L, P
        self.t_tail, self.strip_gap = t_tail, strip_gap
        self.n = len(seq)
        self.d = pair_matrix(seq, L)

    def run(self):
        n, P = self.n, self.P
        store = np.full((n, n), NEG, dtype=np.int64)   # column store in memory, N(m, j)
        cycles = 0
        n_pad = -(-n // P) * P
        for i_top in range(n_pad - P, -1, -P):        # bottom strip first
            i_bot = i_top + P - 1
            rows = [i_bot - p for p in range(P)]
            rowstore = [dict() for _ in range(P)]     # PE-local: k -> N(i, k)
            slast = [NEG] * P
            # timing state
            pe_free = [0] * P                          # cycle a PE can start its next column
            head_t = 0                                 # next cycle the head can push a word
            strip_end = 0
            for j in range(max(i_top, 0), n):
                mem = [store[m, j] if m != j else 0 for m in range(j, i_bot, -1)]  # A_j words
                stream = list(mem)
                # timing: arrival time of each stream word at the current PE
                arrive = []
                for _ in mem:
                    arrive.append(head_t)
                    head_t += 1
                for p, i in enumerate(rows):
                    # ---- function: only stream words, own row store, slast ----
                    acc = NEG
                    last = slast[p]
                    for m_idx, w in enumerate(stream):
                        m = j - m_idx if m_idx < len(mem) else i_bot - (m_idx - len(mem))
                        if m > i:
                            left = rowstore[p].get(m - 1, 0 if m - 1 == i else NEG)
                            if left > NEG and w > NEG:
                                acc = max(acc, left + w)
                        last = w
                    if i >= n or i > j:
                        v = NEG
                    elif i == j:
                        v = 0
                    else:
                        inner = max(slast[p], 0)          # N(i+1, j-1); empty interval = 0
                        v = max(acc, inner + int(self.d[i, j]))
                    slast[p] = last
                    if i < n:
                        rowstore[p][j] = v
                    stream.append(v)
                    # ---- timing ----
                    t = pe_free[p]
                    done = []
                    for a in arrive:
                        t = max(t, a) + 1                 # one add-max per cycle
                        done.append(t)
                    emit = (t if arrive else pe_free[p]) + self.t_tail
                    pe_free[p] = emit
                    arrive = [x + 1 for x in done] + [emit + 1]   # hop to next PE
                strip_end = max(strip_end, arrive[-1])
                for p, i in enumerate(rows):
                    if i < n:
                        store[i, j] = rowstore[p][j]
            cycles += strip_end + self.strip_gap
        N = np.where(store == NEG, 0, store)
        N = np.triu(N)
        return N, cycles


def conformance():
    """spec.md Sec. 14.2 vectors."""
    return [
        ("A", 0, 0), ("AAAA", 0, 0), ("ACGU", 0, 2), ("GGGAAAUCC", 0, 3), ("GGGAAAUCC", 3, 3),
        ("GCAUCUAUGC", 0, 4), ("GCAUCUAUGC", 3, 3), ("GGGGAAAACCCC", 3, 4),
        ("CAGAUUUUCAUAUUAUGCAGAAAAUCUACUUCGCCUGAUACGAGUCGGUUAUCUUCGGAUACUG", 3, 24),
    ]


def rand_seq(n, rng):
    return "".join(rng.choice(list("ACGU"), n))


def cmd_verify(a):
    bad = 0
    for s, L, want in conformance():
        r1 = ref_bifurcation(s, L)
        r2 = ref_classic(s, L)
        got = int(r1[0, -1])
        ok = got == want and np.array_equal(r1, r2)
        bad += not ok
        print(f"{'ok ' if ok else 'BAD'} n={len(s):3d} L={L}  N(0,n-1)={got} (spec says {want})")
    rng = np.random.default_rng(1)
    for t in range(a.random):
        n, L = int(rng.integers(1, 60)), int(rng.integers(0, 4))
        s = rand_seq(n, rng)
        r1 = ref_bifurcation(s, L)
        if not np.array_equal(r1, ref_classic(s, L)):
            bad += 1
            print(f"BAD refs differ: {s} L={L}")
        for P in (1, 2, 4, 8):
            N, _ = ArrayModel(s, L, P).run()
            if not np.array_equal(N, r1):
                bad += 1
                if bad < 5:
                    print(f"BAD array model P={P}: {s} L={L}")
    print(f"{a.random} random sequences x P in {{1,2,4,8}}: "
          f"{'all full tables match' if not bad else f'{bad} failures'}")
    return 1 if bad else 0


def cmd_perf(a):
    rng = np.random.default_rng(7)
    print(f"T_TAIL={a.t_tail}  STRIP_GAP={a.strip_gap}  clock={a.ghz} GHz (assumed)")
    print(f"{'n':>6}{'P':>5}{'cycles':>12}{'ops/cycle':>11}{'util':>8}{'Gops/s':>9}")
    for n in a.n:
        s = rand_seq(n, rng)
        ops = useful_ops(n)
        for P in a.P:
            N, cyc = ArrayModel(s, 3, P, a.t_tail, a.strip_gap).run()
            opc = ops / cyc
            print(f"{n:>6}{P:>5}{cyc:>12}{opc:>11.2f}{100 * opc / P:>7.1f}%{opc * a.ghz:>9.2f}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify", help="conformance vectors, two refs, array model vs refs")
    v.add_argument("--random", type=int, default=100)
    p = sub.add_parser("perf", help="cycle model: ops per cycle and utilization")
    p.add_argument("-n", type=int, nargs="+", default=[64, 128, 256])
    p.add_argument("-P", type=int, nargs="+", default=[4, 16, 64])
    p.add_argument("--t-tail", type=int, default=16)
    p.add_argument("--strip-gap", type=int, default=256)
    p.add_argument("--ghz", type=float, default=2.0, help="GenDP's paper assumes 2 GHz")
    a = ap.parse_args()
    sys.exit({"verify": cmd_verify, "perf": cmd_perf}[a.cmd](a))


if __name__ == "__main__":
    main()
