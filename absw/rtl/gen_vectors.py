#!/usr/bin/env python3
"""Make RTL test vectors with expected results from golden.py.

Output format (read by absw_tb.sv):
    <number of tests>
    tlen qlen h0 adaptive score gscore max_i max_j
    <tlen target symbols>
    <qlen query symbols>
    ...
Symbols are 0..3 (ACGT) and 4 (N).
"""
import argparse
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import golden                      # noqa: E402
from gen_pairs import mutate       # noqa: E402

BASES = "ACGT"


def cases(rng, maxlen, n_random):
    """Yield (target, query, h0) as strings."""
    rs = lambda n: "".join(rng.choice(BASES) for _ in range(n))
    # edge cases
    yield "A", "A", 5
    yield "A", "C", 5
    yield "ACGT", "A", 10
    yield "A", "ACGT", 10
    yield "NNNN", "ACGT", 10
    s = rs(40)
    yield s, s, 20                                   # perfect match
    yield rs(maxlen), rs(8), 20                      # very long target, short query
    yield rs(8), rs(maxlen), 20                      # short target, long query
    yield rs(60), rs(60), 3                          # unrelated, low h0
    # short reads (GenDP's BSW regime)
    for _ in range(n_random):
        q = rng.randint(10, min(100, maxlen - 20))
        ref = rs(q + 20)
        read = mutate(ref, rng, 0.03, 0.03, 3)[:q]
        yield ref, read, rng.randint(1, 40)
    # long, drifting reads: stress the steering
    for _ in range(n_random):
        q = rng.randint(maxlen // 2, maxlen - 30)
        ref = rs(maxlen)
        read = mutate(ref, rng, 0.04, 0.08, 2, del_frac=rng.choice([0.25, 0.5, 0.75]))[:q]
        yield ref, read, rng.randint(10, 40)
    # heavy one-sided drift: a static band falls off, the adaptive one should follow
    for _ in range(n_random):
        q = rng.randint(maxlen // 2, maxlen - 60)
        ref = rs(maxlen)
        read = mutate(ref, rng, 0.03, 0.12, 2, del_frac=rng.choice([0.1, 0.9]))[:q]
        if rng.random() < 0.5:                       # insertions need room in the read, not the ref
            ref, read = ref[:q], read
        yield ref, read, rng.randint(10, 40)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", default="vectors.txt")
    ap.add_argument("-W", type=int, default=16)
    ap.add_argument("--maxlen", type=int, default=256)
    ap.add_argument("-n", type=int, default=40, help="random cases per category")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    out = []
    for t, q, h0 in cases(rng, a.maxlen, a.n):
        assert len(t) <= a.maxlen and len(q) <= a.maxlen and len(q) >= 1 and len(t) >= 1
        tt, qq = golden.encode(t), golden.encode(q)
        for adaptive in (1, 0):
            r = golden.extend(qq, tt, h0, "adaptive" if adaptive else "static", W=a.W)
            out.append(f"{len(t)} {len(q)} {h0} {adaptive} {r['score']} {r['gscore']} "
                       f"{r['max_i']} {r['max_j']}\n"
                       + " ".join(map(str, tt)) + "\n" + " ".join(map(str, qq)) + "\n")
    with open(a.o, "w") as f:
        f.write(f"{len(out)}\n")
        f.writelines(out)
    print(f"wrote {len(out)} tests (W={a.W}, maxlen={a.maxlen}) to {a.o}")


if __name__ == "__main__":
    main()
