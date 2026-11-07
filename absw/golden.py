#!/usr/bin/env python3
"""Golden model for banded and adaptive-banded affine-gap seed extension.

Two layers:

1. ksw_extend2_port(): a line-for-line port of ksw_extend2() in
   gendp/kernel/bwa-mem/ksw.cpp (BWA-MEM seed extension, what GenDP's BSW
   kernel computes).  Validated field-for-field against the C binary
   (`python3 golden.py validate`).

2. extend(): an anti-diagonal DP engine using the same recurrences, with
   pluggable band policies:
     full      every cell (ground truth)
     rowband   |i - j| <= w, ksw's static band (no zero-pruning, no z-drop)
     static    W cells per anti-diagonal, band alternates right/down
               (a fixed diagonal band, expressed the way a W-PE array sees it)
     adaptive  W cells per anti-diagonal, band moves right or down each step
               toward the end with the higher H score (Suzuki-Kasahara style;
               the rule used in Liao et al., ASAP 2018 / GenDP ref [44])
   `static` and `adaptive` use identical hardware; only the 1-bit move
   decision per anti-diagonal differs.  That decision is the new hardware.

Conventions (match ksw.cpp): matrix indices are 1-based, i = target (ref)
position, j = query (read) position.  Row 0 and column 0 hold the h0-derived
boundary.  A cell outside the band is treated as H = E = F = M = 0, which is
how ksw's band behaves too.
"""
import argparse
import sys

import numpy as np

MATCH, MISMATCH, AMBIG = 1, 4, 1          # sa, sb, ambiguous penalty
GAP_O, GAP_E = 6, 1
END_BONUS, ZDROP = 5, 100
NT4 = {c: i for i, c in enumerate("ACGT")}


def encode(s):
    return np.array([NT4.get(c, 4) for c in s.upper()], dtype=np.int32)


def score_matrix():
    m = np.full((5, 5), -AMBIG, dtype=np.int32)
    for a in range(4):
        for b in range(4):
            m[a, b] = MATCH if a == b else -MISMATCH
    return m


MAT = score_matrix()


# --------------------------------------------------------------------------
# 1. Exact port of ksw_extend2 (row-major, with BWA's pruning and z-drop)
# --------------------------------------------------------------------------
def ksw_extend2_port(query, target, h0, w, zdrop=ZDROP, end_bonus=END_BONUS,
                     o_del=GAP_O, e_del=GAP_E, o_ins=GAP_O, e_ins=GAP_E):
    """Returns (score, qle, tle, gtle, gscore, max_off, cells)."""
    q, t = list(query), list(target)
    qlen, tlen = len(q), len(t)
    oe_del, oe_ins = o_del + e_del, o_ins + e_ins
    assert h0 > 0
    eh_h = [0] * (qlen + 1)
    eh_e = [0] * (qlen + 1)
    eh_h[0] = h0
    eh_h[1] = h0 - oe_ins if h0 > oe_ins else 0
    j = 2
    while j <= qlen and eh_h[j - 1] > e_ins:
        eh_h[j] = eh_h[j - 1] - e_ins
        j += 1
    mx = int(MAT.max())
    max_ins = max(int((qlen * mx + end_bonus - o_ins) / e_ins + 1.0), 1)
    w = min(w, max_ins)
    max_del = max(int((qlen * mx + end_bonus - o_del) / e_del + 1.0), 1)
    w = min(w, max_del)

    best, max_i, max_j, max_ie, gscore, max_off = h0, -1, -1, -1, -1, 0
    beg, end, cells = 0, qlen, 0
    for i in range(tlen):
        f, m, mj = 0, 0, -1
        row = MAT[t[i]]
        if beg < i - w:
            beg = i - w
        if end > i + w + 1:
            end = i + w + 1
        if end > qlen:
            end = qlen
        if beg == 0:
            h1 = max(h0 - (o_del + e_del * (i + 1)), 0)
        else:
            h1 = 0
        j = beg
        while j < end:
            cells += 1
            M, e = eh_h[j], eh_e[j]
            eh_h[j] = h1
            M = M + int(row[q[j]]) if M else 0
            h = max(M, e, f)
            mj = mj if m > h else j
            m = m if m > h else h
            tt = max(M - oe_del, 0)
            e = max(e - e_del, tt)
            f = max(f - e_ins, tt)          # repo uses oe_del here too (o_ins line is commented out)
            h1 = h
            eh_e[j] = e
            j += 1
        eh_h[end] = h1
        eh_e[end] = 0
        if j == qlen:
            if not gscore > h1:
                max_ie = i
            gscore = max(gscore, h1)
        if m == 0:
            break
        if m > best:
            best, max_i, max_j = m, i, mj
            max_off = max(max_off, abs(mj - i))
        elif zdrop > 0:
            if i - max_i > mj - max_j:
                if best - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop:
                    break
            else:
                if best - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop:
                    break
        j = beg
        while j < end and eh_h[j] == 0 and eh_e[j] == 0:
            j += 1
        beg = j
        j = end
        while j >= beg and eh_h[j] == 0 and eh_e[j] == 0:
            j -= 1
        end = j + 2 if j + 2 < qlen else qlen
    return best, max_j + 1, max_i + 1, max_ie + 1, gscore, max_off, cells


# --------------------------------------------------------------------------
# 2. Anti-diagonal engine with band policies
# --------------------------------------------------------------------------
def _band_rows(policy, d, lo, W, w, tlen, qlen):
    """Rows i (1-based) computed on anti-diagonal d (cells (i, d - i))."""
    i_min, i_max = max(1, d - qlen), min(tlen, d - 1)
    if policy == "full":
        a, b = i_min, i_max
    elif policy == "rowband":                     # |i - j| <= w  ->  |2i - d| <= w
        a = max(i_min, -((w - d) // 2))           # ceil((d - w) / 2)
        b = min(i_max, (d + w) // 2)
    else:                                         # static / adaptive: W-wide window
        a, b = max(i_min, lo), min(i_max, lo + W - 1)
    return a, b


def extend(query, target, h0, policy="full", W=16, w=None,
           o_del=GAP_O, e_del=GAP_E, o_ins=GAP_O, e_ins=GAP_E, trace=False):
    """Seed extension with the given band policy.

    Returns dict(score, max_i, max_j, gscore, cells, moves[, band]).
    No z-drop or zero-pruning, so policies differ only in which cells exist.
    """
    q, t = np.asarray(query), np.asarray(target)
    qlen, tlen = len(q), len(t)
    oe_del, oe_ins = o_del + e_del, o_ins + e_ins
    if w is None:
        w = W - 1
    H = np.zeros((tlen + 1, qlen + 1), dtype=np.int64)
    E = np.zeros_like(H)
    F = np.zeros_like(H)
    Mx = np.zeros_like(H)
    H[0, :] = np.maximum(0, h0 - o_ins - e_ins * np.arange(qlen + 1))
    H[:, 0] = np.maximum(0, h0 - o_del - e_del * np.arange(tlen + 1))
    H[0, 0] = h0
    S = MAT[t[:, None], q[None, :]]               # S[i-1, j-1]

    best, bi, bj = h0, 0, 0
    gscore, cells = -1, 0
    moves = []                                    # 1 = down, 0 = right
    band = [] if trace else None
    lo = 1 - (W - 1) // 2                         # centred on cell (1,1)
    last_dir = 1
    for d in range(2, tlen + qlen + 1):
        a, b = _band_rows(policy, d, lo, W, w, tlen, qlen)
        if a <= b:
            i = np.arange(a, b + 1)
            j = d - i
            hd = H[i - 1, j - 1]
            M = np.where(hd > 0, hd + S[i - 1, j - 1], 0)
            # E(i,j) from (i-1,j); F(i,j) from (i,j-1).  Row 1 / column 1 start at 0.
            e = np.where(i > 1, np.maximum(E[i - 1, j] - e_del, np.maximum(Mx[i - 1, j] - oe_del, 0)), 0)
            f = np.where(j > 1, np.maximum(F[i, j - 1] - e_ins, np.maximum(Mx[i, j - 1] - oe_ins, 0)), 0)
            h = np.maximum(np.maximum(M, e), f)
            Mx[i, j], E[i, j], F[i, j], H[i, j] = M, e, f, h
            cells += len(i)
            k = int(np.argmax(h))
            if h[k] > best:
                best, bi, bj = int(h[k]), int(i[k]), int(j[k])
            if j[0] == qlen:                      # smallest i on the diagonal has the largest j
                gscore = max(gscore, int(h[0]))
            if trace:
                band.append((d, a, b))
        if policy in ("static", "adaptive"):
            # Decide how the W-wide window moves for anti-diagonal d+1.
            top_i, bot_i = lo, lo + W - 1
            def end_score(ii):
                jj = d - ii
                if 1 <= ii <= tlen and 1 <= jj <= qlen:
                    return int(H[ii, jj])
                return -1                         # outside the matrix loses
            ht, hb = end_score(top_i), end_score(bot_i)
            if policy == "adaptive" and ht != hb:
                down = hb > ht
            else:
                down = not last_dir               # alternate: follows the main diagonal
            # Never steer completely off the matrix.
            if lo + 1 > tlen:
                down = False
            if (d + 1) - (lo + W - 1) > qlen:     # right move would leave the bottom end past qlen
                down = True
            last_dir = down
            moves.append(1 if down else 0)
            if down:
                lo += 1
    out = dict(score=best, max_i=bi, max_j=bj, gscore=gscore, cells=cells, moves=moves)
    if trace:
        out["band"] = band
    return out


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def read_ksw_csv(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            if not line[:1].isalpha():
                continue
            qs, ts, h0, w = line.strip().split(",")[:4]
            pairs.append((qs, ts, int(h0), int(w)))
    return pairs


def cmd_validate(a):
    """Port vs ksw-test output, then engine(rowband) vs port without z-drop/pruning."""
    pairs = read_ksw_csv(a.input)
    ref = [list(map(int, l.split(","))) for l in open(a.expected) if l.strip()]
    bad = 0
    for k, ((qs, ts, h0, w), exp) in enumerate(zip(pairs, ref)):
        got = list(ksw_extend2_port(encode(qs), encode(ts), h0, w)[:6])
        if got != exp:
            bad += 1
            if bad <= 5:
                print(f"port mismatch pair {k}: got {got} expected {exp}")
    print(f"port vs ksw-test: {len(pairs) - bad}/{len(pairs)} identical on all 6 fields")

    bad2 = 0
    for k, (qs, ts, h0, w) in enumerate(pairs):
        q, t = encode(qs), encode(ts)
        for ww in (w, 8, 3):
            p = ksw_extend2_port_noprune(q, t, h0, ww)
            r = extend(q, t, h0, "rowband", w=ww)
            if (p[0], p[4]) != (r["score"], r["gscore"]):
                bad2 += 1
                if bad2 <= 5:
                    print(f"engine mismatch pair {k} w={ww}: port {(p[0], p[4])} engine {(r['score'], r['gscore'])}")
    print(f"engine(rowband) vs port (no z-drop, no pruning), w in {{given,8,3}}: "
          f"{3 * len(pairs) - bad2}/{3 * len(pairs)} identical score and gscore")
    return 1 if (bad or bad2) else 0


def ksw_extend2_port_noprune(query, target, h0, w):
    """ksw_extend2 with z-drop off and the 'beg = 0; end = qlen' debug line enabled."""
    q, t = list(query), list(target)
    qlen, tlen = len(q), len(t)
    oe = GAP_O + GAP_E
    eh_h = [0] * (qlen + 1)
    eh_e = [0] * (qlen + 1)
    eh_h[0] = h0
    eh_h[1] = h0 - oe if h0 > oe else 0
    j = 2
    while j <= qlen and eh_h[j - 1] > GAP_E:
        eh_h[j] = eh_h[j - 1] - GAP_E
        j += 1
    best, gscore = h0, -1
    for i in range(tlen):
        f, m = 0, 0
        beg, end = max(0, i - w), min(qlen, i + w + 1)
        if beg >= end:
            continue
        h1 = max(h0 - (GAP_O + GAP_E * (i + 1)), 0) if beg == 0 else 0   # H(i, beg-1) = 0 when out of band
        for j in range(beg, end):
            M, e = eh_h[j], eh_e[j]
            eh_h[j] = h1
            M = M + int(MAT[t[i], q[j]]) if M else 0
            h = max(M, e, f)
            m = max(m, h)
            tt = max(M - oe, 0)
            e = max(e - GAP_E, tt)
            f = max(f - GAP_E, tt)
            h1 = h
            eh_e[j] = e
        eh_h[end] = h1
        eh_e[end] = 0
        if end == qlen:
            gscore = max(gscore, h1)
        best = max(best, m)
    return best, None, None, None, gscore


def cmd_compare(a):
    rows = []
    with open(a.tsv) as fh:
        for line in fh:
            qs, ts, h0 = line.rstrip("\n").split("\t")
            rows.append((encode(qs), encode(ts), int(h0)))
    if a.n:
        rows = rows[:a.n]
    opt = [extend(q, t, h0, "full") for q, t, h0 in rows]
    print(f"{len(rows)} pairs, mean qlen {np.mean([len(r[0]) for r in rows]):.0f}, "
          f"mean full-DP cells {np.mean([o['cells'] for o in opt]):.0f}")
    print(f"{'policy':<10}{'W':>5}{'cells/pair':>12}{'% optimal':>11}{'mean loss':>11}{'worst loss':>12}")
    for W in a.W:
        for pol in ("rowband", "static", "adaptive"):
            res = [extend(q, t, h0, pol, W=W) for q, t, h0 in rows]
            loss = np.array([o["score"] - r["score"] for o, r in zip(opt, res)])
            cells = np.mean([r["cells"] for r in res])
            print(f"{pol:<10}{W:>5}{cells:>12.0f}{100 * np.mean(loss == 0):>10.1f}%"
                  f"{loss.mean():>11.2f}{loss.max():>12d}")
    return 0


def cmd_selftest(_a):
    # Perfect match: score = h0 + qlen, regardless of band.
    s = encode("ACGTACGTACGTACGTACGT")
    for pol in ("full", "rowband", "static", "adaptive"):
        r = extend(s, s, 10, pol, W=4)
        assert r["score"] == 30, (pol, r["score"])
    rng = np.random.default_rng(0)
    ref = rng.integers(0, 4, 400).astype(np.int32)
    # Drift: one deleted base every 25, so the path drifts 12 diagonals over 300bp.
    # A static W=8 band falls off; the adaptive band follows it.
    keep = np.ones(330, dtype=bool)
    keep[25::25] = False
    read = ref[:330][keep][:300]
    full = extend(read, ref, 20, "full")["score"]
    stat = extend(read, ref, 20, "static", W=8)["score"]
    adap = extend(read, ref, 20, "adaptive", W=8)["score"]
    assert adap == full, (adap, full)
    assert stat < full, (stat, full)
    # Known limit: one gap wider than the band is lost by both (no steering signal
    # until the new diagonal is inside the band).
    read = np.concatenate([ref[:60], ref[84:184]])
    big = [extend(read, ref, 20, p, W=8)["score"] for p in ("full", "static", "adaptive")]
    assert big[1] < big[0] and big[2] < big[0], big
    print(f"selftest ok  (drift 12 over 300bp, W=8: full {full}, static {stat}, adaptive {adap}; "
          f"24bp gap, W=8: {big})")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="check the port against ksw-test output")
    v.add_argument("input", help="ksw-test input csv (QUERY,TARGET,h0,w)")
    v.add_argument("expected", help="ksw-test -x output")
    c = sub.add_parser("compare", help="score loss of band policies vs full DP")
    c.add_argument("tsv", help="query<TAB>target<TAB>h0 (gen_pairs.py --tsv)")
    c.add_argument("-W", type=int, nargs="+", default=[8, 16, 32])
    c.add_argument("-n", type=int, default=0)
    sub.add_parser("selftest")
    a = ap.parse_args()
    sys.exit({"validate": cmd_validate, "compare": cmd_compare, "selftest": cmd_selftest}[a.cmd](a))


if __name__ == "__main__":
    main()
