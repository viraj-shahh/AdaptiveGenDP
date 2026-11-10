#!/usr/bin/env python3
"""Cut real ONT read-to-reference alignments into fixed-length segment pairs (read segment + matching genome segment).

Input: minimap2 PAF with CIGAR (run with -c), the reads (FASTQ) and the reference (FASTA).
Output: TSV in the golden model's format (query \\t target \\t h0), one segment pair per line.

Each pair starts at a matched base on minimap2's alignment path (like an extension from a
seed) and covers L read bases. The target is the reference segment that minimap2 aligned to those
read bases, plus SLACK extra bases so the extension has room to run past the true end.
A side file records each pair's indel drift (how far the true path leaves the diagonal).

usage: cut_segment_pairs.py reads.paf reads.fastq ref.fasta [--outdir D] [-L 1000] [--slack 50] [-n 5000] [--seed 1]
       Writes D/seg_<reads>_<ref>_L<L>_slack<S>_n<n>_seed<seed>.tsv, plus .drift (per-pair drift)
       and .info (command, date, input paths + sizes + md5, parameters, counts).
"""
import argparse, datetime, hashlib, os, random, re, shlex, sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def read_fastq(path, wanted):
    out = {}
    with open(path) as f:
        while True:
            h = f.readline()
            if not h: break
            s = f.readline().strip(); f.readline(); f.readline()
            name = h[1:].split()[0]
            if name in wanted: out[name] = s
    return out


def segment_pairs(rec, read, ref, L, slack):
    """Yield (query, target, max |drift|, net drift) segment pairs along one alignment."""
    qname, qlen, qs, qe, strand, tname, tlen, ts, te = rec[:9]
    qlen, qs, qe, ts, te = int(qlen), int(qs), int(qe), int(ts), int(te)
    cigar = next(x[5:] for x in rec[12:] if x.startswith("cg:Z:"))
    q = read if strand == "+" else revcomp(read)
    if strand == "-": qs, qe = qlen - qe, qlen - qs
    # Walk the CIGAR, recording the ref position paired with each read position.
    qi, ti = qs, ts
    q2t = {}             # read pos -> ref pos, at aligned (M) columns only
    path = []            # (read pos, ref pos) after every op step, for drift
    first = {}           # read pos -> index of its first entry in path
    for n, op in re.findall(r"(\d+)([MID])", cigar):
        n = int(n)
        for _ in range(n):
            if op == "M": q2t[qi] = ti; first.setdefault(qi, len(path)); qi += 1; ti += 1
            elif op == "I": first.setdefault(qi, len(path)); qi += 1
            else: ti += 1
            path.append((qi, ti))
    start = qs
    while start + L <= qe:
        if start in q2t:
            t0 = q2t[start]
            end = start + L
            t1 = next((q2t[k] for k in range(end, qe) if k in q2t), te)
            # drift = (ref advance) - (read advance) relative to the segment's start
            i0, i1 = first[start], first.get(end, len(path))
            d = [(tj - t0) - (qj - start) for qj, tj in path[i0:i1]]
            yield (q[start:end], ref[t0:min(t1 + slack, len(ref))],
                   max(abs(x) for x in d) if d else 0, d[-1] if d else 0)
            start = end
        else:
            start += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf"); ap.add_argument("fastq"); ap.add_argument("fasta")
    ap.add_argument("--outdir", default=".",
                    help="output named seg_<reads>_<ref>_L<L>_slack<S>_n<n>_seed<seed>.tsv (+ .drift, .info)")
    ap.add_argument("-L", type=int, default=1000, help="read bases per segment")
    ap.add_argument("--slack", type=int, default=50, help="extra ref bases past the true end")
    ap.add_argument("-n", type=int, default=5000, help="segment pairs to sample (0 = all)")
    ap.add_argument("--mapq", type=int, default=60)
    ap.add_argument("--h0", type=int, default=30, help="seed score, same default as golden pairs")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    stem = lambda p: os.path.basename(p).split(".")[0]          # guppy_hac.fastq -> guppy_hac
    a.out = (f"{a.outdir}/seg_{stem(a.fastq)}_{stem(a.fasta)}"
             f"_L{a.L}_slack{a.slack}_n{a.n}_seed{a.seed}.tsv")

    recs = [l.rstrip("\n").split("\t") for l in open(a.paf)]
    recs = [r for r in recs if int(r[11]) >= a.mapq and "tp:A:P" in r[12:]]
    ref = read_fasta(a.fasta)
    reads = read_fastq(a.fastq, {r[0] for r in recs})
    rng, allw, total = random.Random(a.seed), [], 0
    for r in recs:                       # reservoir sample, so memory stays at n pairs
        for w in segment_pairs(r, reads[r[0]], ref[r[5]], a.L, a.slack):
            total += 1
            if not a.n or len(allw) < a.n: allw.append(w)
            else:
                k = rng.randrange(total)
                if k < a.n: allw[k] = w
    print(f"{len(recs)} primary alignments (MAPQ>={a.mapq}) -> {total} segment pairs of {a.L} bp",
          file=sys.stderr)
    with open(a.out, "w") as f, open(a.out + ".drift", "w") as g:
        g.write("max_abs_drift\tnet_drift\n")
        for qseq, tseq, md, nd in allw:
            f.write(f"{qseq}\t{tseq}\t{a.h0}\n")
            g.write(f"{md}\t{nd}\n")
    with open(a.out + ".info", "w") as f:    # provenance: what this file was made from, and how
        f.write(f"created: {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"command: {' '.join(shlex.quote(x) for x in sys.argv)}\n")
        f.write(f"cwd: {os.getcwd()}\n")
        for label, p in (("paf", a.paf), ("reads", a.fastq), ("reference", a.fasta)):
            f.write(f"{label}: {os.path.realpath(p)}  size={os.path.getsize(p)}  md5={md5(p)}\n")
        f.write(f"params: L={a.L} slack={a.slack} n={a.n} seed={a.seed} mapq>={a.mapq} h0={a.h0}\n")
        f.write(f"primary alignments used: {len(recs)}\n")
        f.write(f"segment pairs available: {total}, written: {len(allw)}\n")
        f.write("columns (.tsv): read segment, genome segment, h0\n")
        f.write("columns (.drift): max |drift| from diagonal, net drift at segment end\n")
    print(f"wrote {len(allw)} segment pairs to {a.out} (+ .drift, .info)", file=sys.stderr)


if __name__ == "__main__":
    main()
