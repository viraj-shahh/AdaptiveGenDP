#!/usr/bin/env python3
"""Generate synthetic read/reference pairs for BSW and adaptive-band experiments.

Each pair is a reference window plus a read derived from it with substitutions
and indels.  Writes up to two formats:

  --gendp   GenDP simulator format: 3 lines per pair -> h0, ref digits, query digits
            (A=0 C=1 G=2 T=3), as read by gendp/bsw.cpp:loadPairs.
  --ksw     ksw-test format: "QUERY,TARGET,h0,w" per line (ACGT letters).

Also --tsv: query<TAB>target<TAB>h0 for the Python golden model.

Example:
  python3 gen_pairs.py -n 64 --qlen 100 --indel-rate 0.02 --gendp t.txt --ksw t.csv
"""
import argparse
import random

BASES = "ACGT"


def mutate(seq, rng, sub_rate, indel_rate, max_indel, del_frac=0.5):
    out = []
    i = 0
    while i < len(seq):
        r = rng.random()
        if r < indel_rate * del_frac:               # deletion from read
            i += rng.randint(1, max_indel)
            continue
        if r < indel_rate:                          # insertion into read
            out.extend(rng.choice(BASES) for _ in range(rng.randint(1, max_indel)))
        c = seq[i]
        if rng.random() < sub_rate:
            c = rng.choice([b for b in BASES if b != c])
        out.append(c)
        i += 1
    return "".join(out)


def make_pairs(n, qlen, tlen_extra, sub_rate, indel_rate, max_indel, h0, seed, max_len, del_frac=0.5):
    rng = random.Random(seed)
    pairs = []
    while len(pairs) < n:
        ref = "".join(rng.choice(BASES) for _ in range(qlen + tlen_extra))
        read = mutate(ref, rng, sub_rate, indel_rate, max_indel, del_frac)[:qlen]
        if not read or (max_len and (len(read) >= max_len or len(ref) >= max_len)):
            continue
        pairs.append((read, ref, h0))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, default=64)
    ap.add_argument("--qlen", type=int, default=100)
    ap.add_argument("--tlen-extra", type=int, default=20, help="reference is qlen + this long")
    ap.add_argument("--sub-rate", type=float, default=0.02)
    ap.add_argument("--indel-rate", type=float, default=0.01)
    ap.add_argument("--max-indel", type=int, default=3)
    ap.add_argument("--del-frac", type=float, default=0.5,
                    help="fraction of indels that are deletions from the read; != 0.5 gives drift")
    ap.add_argument("--h0", type=int, default=30, help="seed score carried into the extension")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-len", type=int, default=128,
                    help="GenDP BSW buffers are 128 wide; set 0 to disable the limit")
    ap.add_argument("--gendp")
    ap.add_argument("--ksw")
    ap.add_argument("--tsv")
    a = ap.parse_args()

    pairs = make_pairs(a.n, a.qlen, a.tlen_extra, a.sub_rate, a.indel_rate, a.max_indel,
                       a.h0, a.seed, a.max_len, a.del_frac)
    if a.gendp:
        with open(a.gendp, "w") as f:
            for q, t, h0 in pairs:
                f.write(f"{h0}\n{''.join(str(BASES.index(c)) for c in t)}\n"
                        f"{''.join(str(BASES.index(c)) for c in q)}\n")
    if a.ksw:
        with open(a.ksw, "w") as f:
            for q, t, h0 in pairs:
                f.write(f"{q},{t},{h0},100\n")
    if a.tsv:
        with open(a.tsv, "w") as f:
            for q, t, h0 in pairs:
                f.write(f"{q}\t{t}\t{h0}\n")
    print(f"wrote {len(pairs)} pairs")


if __name__ == "__main__":
    main()
