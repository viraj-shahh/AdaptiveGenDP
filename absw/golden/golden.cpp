// absw golden model command-line tool.  Build: make build/golden (from absw/).
//
//   golden selftest
//   golden validate <ksw.csv> <ksw_out.txt>        port vs ksw-test, engine vs port
//   golden compare  <pairs.tsv> [-W 8 16 32] [-n N] score loss of band policies vs full DP
//   golden run      <pairs.tsv> <policy> [-W n]     per-pair: score gscore max_i max_j cells
//   golden pairs    [options] --tsv F --ksw F --gendp F   synthetic read/reference pairs
//   golden vectors  -o F [-W 16] [--maxlen 256] [-n 40] [--seed 1] [--tsv F]   RTL test vectors
//
// TSV format: query<TAB>target<TAB>h0.  ksw csv: QUERY,TARGET,h0,w.
#include "absw_golden.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace absw;

namespace {

struct Pair {
    std::string q, t;
    int h0;
};

std::vector<Pair> read_tsv(const std::string &path) {
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(2); }
    std::vector<Pair> v;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        Pair p;
        std::string h0;
        std::getline(ss, p.q, '\t');
        std::getline(ss, p.t, '\t');
        std::getline(ss, h0, '\t');
        p.h0 = std::stoi(h0);
        v.push_back(p);
    }
    return v;
}

std::vector<int> parse_ints_after(int argc, char **argv, int &k) {
    std::vector<int> v;
    while (k + 1 < argc && argv[k + 1][0] != '-') v.push_back(std::atoi(argv[++k]));
    return v;
}

int cmd_selftest() {
    // Perfect match: score = h0 + qlen for every policy.
    Seq s = encode("ACGTACGTACGTACGTACGT");
    for (Policy p : {Policy::Full, Policy::RowBand, Policy::Static, Policy::Adaptive}) {
        int sc = extend(s, s, 10, p, 4).score;
        if (sc != 30) { std::printf("selftest FAIL: perfect match %s gave %d\n", policy_name(p), sc); return 1; }
    }
    Rng rng(0);
    std::string ref = rng.bases(400);
    // Drift: one deleted base every 25, so the path drifts 12 diagonals over 300bp.
    std::string read;
    for (int i = 0; i < 330; ++i) if (i == 0 || i % 25 != 0) read.push_back(ref[i]);
    read = read.substr(0, 300);
    Seq rq = encode(read), rt = encode(ref);
    int full = extend(rq, rt, 20, Policy::Full).score;
    int stat = extend(rq, rt, 20, Policy::Static, 8).score;
    int adap = extend(rq, rt, 20, Policy::Adaptive, 8).score;
    if (adap != full || stat >= full) {
        std::printf("selftest FAIL: drift full %d static %d adaptive %d\n", full, stat, adap);
        return 1;
    }
    // Known limit: one gap wider than the band is lost by both.
    Seq gq = encode(ref.substr(0, 60) + ref.substr(84, 100));
    int gf = extend(gq, rt, 20, Policy::Full).score;
    int gs = extend(gq, rt, 20, Policy::Static, 8).score;
    int ga = extend(gq, rt, 20, Policy::Adaptive, 8).score;
    if (!(gs < gf && ga < gf)) {
        std::printf("selftest FAIL: 24bp gap full %d static %d adaptive %d\n", gf, gs, ga);
        return 1;
    }
    std::printf("selftest ok  (drift 12 over 300bp, W=8: full %d, static %d, adaptive %d; "
                "24bp gap, W=8: [%d, %d, %d])\n", full, stat, adap, gf, gs, ga);
    return 0;
}

int cmd_validate(const std::string &csv, const std::string &expected) {
    std::ifstream fi(csv), fe(expected);
    if (!fi || !fe) { std::fprintf(stderr, "cannot open inputs\n"); return 2; }
    std::vector<std::array<std::string, 2>> seqs;
    std::vector<std::pair<int, int>> hw;
    std::string line;
    while (std::getline(fi, line)) {
        if (line.empty() || !std::isalpha((unsigned char)line[0])) continue;
        std::istringstream ss(line);
        std::string q, t, h0, w;
        std::getline(ss, q, ','); std::getline(ss, t, ','); std::getline(ss, h0, ','); std::getline(ss, w, ',');
        seqs.push_back({q, t});
        hw.push_back({std::stoi(h0), std::stoi(w)});
    }
    std::vector<std::array<int, 6>> ref;
    while (std::getline(fe, line)) {
        if (line.empty()) continue;
        std::array<int, 6> r{};
        std::sscanf(line.c_str(), "%d,%d,%d,%d,%d,%d", &r[0], &r[1], &r[2], &r[3], &r[4], &r[5]);
        ref.push_back(r);
    }
    const size_t n = std::min(seqs.size(), ref.size());
    int bad = 0;
    for (size_t k = 0; k < n; ++k) {
        KswResult g = ksw_extend2(encode(seqs[k][0]), encode(seqs[k][1]), hw[k].first, hw[k].second);
        std::array<int, 6> got{g.score, g.qle, g.tle, g.gtle, g.gscore, g.max_off};
        if (got != ref[k] && ++bad <= 5)
            std::printf("port mismatch pair %zu: got %d,%d,%d,%d,%d,%d\n", k, got[0], got[1], got[2], got[3], got[4], got[5]);
    }
    std::printf("port vs ksw-test: %zu/%zu identical on all 6 fields\n", n - bad, n);

    int bad2 = 0;
    for (size_t k = 0; k < seqs.size(); ++k) {
        Seq q = encode(seqs[k][0]), t = encode(seqs[k][1]);
        for (int ww : {hw[k].second, 8, 3}) {
            auto p = ksw_extend2_noprune(q, t, hw[k].first, ww);
            ExtendResult r = extend(q, t, hw[k].first, Policy::RowBand, 16, ww);
            if ((p.first != r.score || p.second != r.gscore) && ++bad2 <= 5)
                std::printf("engine mismatch pair %zu w=%d\n", k, ww);
        }
    }
    std::printf("engine(rowband) vs port (no z-drop, no pruning), w in {given,8,3}: "
                "%zu/%zu identical score and gscore\n", 3 * seqs.size() - bad2, 3 * seqs.size());
    return (bad || bad2) ? 1 : 0;
}

int cmd_compare(const std::string &tsv, std::vector<int> Ws, int limit) {
    auto rows = read_tsv(tsv);
    if (limit > 0 && (size_t)limit < rows.size()) rows.resize(limit);
    std::vector<Seq> qs, ts;
    std::vector<int> opt;
    double qsum = 0, csum = 0;
    for (auto &p : rows) {
        qs.push_back(encode(p.q)); ts.push_back(encode(p.t));
        ExtendResult o = extend(qs.back(), ts.back(), p.h0, Policy::Full);
        opt.push_back(o.score);
        qsum += p.q.size(); csum += o.cells;
    }
    const size_t n = rows.size();
    std::printf("%zu pairs, mean qlen %.0f, mean full-DP cells %.0f\n", n, qsum / n, csum / n);
    std::printf("%-10s%5s%12s%11s%11s%12s\n", "policy", "W", "cells/pair", "% optimal", "mean loss", "worst loss");
    for (int W : Ws) {
        for (Policy pol : {Policy::RowBand, Policy::Static, Policy::Adaptive}) {
            double cells = 0, loss_sum = 0;
            int optimal = 0, worst = 0;
            for (size_t k = 0; k < n; ++k) {
                ExtendResult r = extend(qs[k], ts[k], rows[k].h0, pol, W);
                const int loss = opt[k] - r.score;
                cells += r.cells; loss_sum += loss;
                optimal += loss == 0; worst = std::max(worst, loss);
            }
            std::printf("%-10s%5d%12.0f%10.1f%%%11.2f%12d\n", policy_name(pol), W, cells / n,
                        100.0 * optimal / n, loss_sum / n, worst);
        }
    }
    return 0;
}

int cmd_run(const std::string &tsv, Policy pol, int W) {
    for (auto &p : read_tsv(tsv)) {
        ExtendResult r = extend(encode(p.q), encode(p.t), p.h0, pol, W);
        std::printf("%d %d %d %d %ld\n", r.score, r.gscore, r.max_i, r.max_j, r.cells);
    }
    return 0;
}

int cmd_pairs(int argc, char **argv) {
    int n = 64, qlen = 100, tlen_extra = 20, max_indel = 3, h0 = 30, max_len = 128;
    double sub_rate = 0.02, indel_rate = 0.01, del_frac = 0.5;
    uint64_t seed = 1;
    std::string tsv, ksw, gendp;
    for (int k = 2; k < argc; ++k) {
        std::string a = argv[k];
        auto val = [&]() { return std::string(argv[++k]); };
        if (a == "-n") n = std::stoi(val());
        else if (a == "--qlen") qlen = std::stoi(val());
        else if (a == "--tlen-extra") tlen_extra = std::stoi(val());
        else if (a == "--sub-rate") sub_rate = std::stod(val());
        else if (a == "--indel-rate") indel_rate = std::stod(val());
        else if (a == "--max-indel") max_indel = std::stoi(val());
        else if (a == "--del-frac") del_frac = std::stod(val());
        else if (a == "--h0") h0 = std::stoi(val());
        else if (a == "--seed") seed = std::stoull(val());
        else if (a == "--max-len") max_len = std::stoi(val());
        else if (a == "--tsv") tsv = val();
        else if (a == "--ksw") ksw = val();
        else if (a == "--gendp") gendp = val();
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    Rng rng(seed);
    std::vector<std::pair<std::string, std::string>> pairs;   // (read, ref)
    while ((int)pairs.size() < n) {
        std::string ref = rng.bases(qlen + tlen_extra);
        std::string read = mutate(ref, rng, sub_rate, indel_rate, max_indel, del_frac).substr(0, qlen);
        if (read.empty() || (max_len && ((int)read.size() >= max_len || (int)ref.size() >= max_len))) continue;
        pairs.push_back({read, ref});
    }
    auto digits = [](const std::string &s) { std::string d; for (char c : s) d.push_back('0' + encode(std::string(1, c))[0]); return d; };
    if (!tsv.empty())   { std::ofstream f(tsv);   for (auto &p : pairs) f << p.first << '\t' << p.second << '\t' << h0 << '\n'; }
    if (!ksw.empty())   { std::ofstream f(ksw);   for (auto &p : pairs) f << p.first << ',' << p.second << ',' << h0 << ",100\n"; }
    if (!gendp.empty()) { std::ofstream f(gendp); for (auto &p : pairs) f << h0 << '\n' << digits(p.second) << '\n' << digits(p.first) << '\n'; }
    std::printf("wrote %zu pairs\n", pairs.size());
    return 0;
}

// RTL test vectors (read by rtl/absw_tb.sv):
//   <number of tests>
//   tlen qlen h0 adaptive score gscore max_i max_j
//   <tlen target symbols>
//   <qlen query symbols>
// With --tsv, the tests are the pairs in that file (e.g. real read segments) instead of
// the synthetic set; pairs longer than --maxlen are skipped and counted.
int cmd_vectors(int argc, char **argv) {
    std::string out = "vectors.txt", tsv;
    int W = 16, maxlen = 256, nr = -1;
    uint64_t seed = 1;
    for (int k = 2; k < argc; ++k) {
        std::string a = argv[k];
        auto val = [&]() { return std::string(argv[++k]); };
        if (a == "-o") out = val();
        else if (a == "-W") W = std::stoi(val());
        else if (a == "--maxlen") maxlen = std::stoi(val());
        else if (a == "-n") nr = std::stoi(val());
        else if (a == "--seed") seed = std::stoull(val());
        else if (a == "--tsv") tsv = val();
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    Rng rng(seed);
    struct Case { std::string t, q; int h0; };
    std::vector<Case> cs;
    if (!tsv.empty()) {
        int skipped = 0;
        for (auto &p : read_tsv(tsv)) {
            if (nr >= 0 && (int)cs.size() >= nr) break;
            if (p.q.empty() || p.t.empty() || (int)p.q.size() > maxlen || (int)p.t.size() > maxlen) { ++skipped; continue; }
            cs.push_back({p.t, p.q, p.h0});
        }
        if (skipped) std::printf("skipped %d pairs longer than maxlen %d\n", skipped, maxlen);
        if (cs.empty()) { std::fprintf(stderr, "no usable pairs in %s\n", tsv.c_str()); return 1; }
    } else {
    if (nr < 0) nr = 40;
    // edge cases
    cs.push_back({"A", "A", 5});
    cs.push_back({"A", "C", 5});
    cs.push_back({"ACGT", "A", 10});
    cs.push_back({"A", "ACGT", 10});
    cs.push_back({"NNNN", "ACGT", 10});
    { std::string s = rng.bases(40); cs.push_back({s, s, 20}); }               // perfect match
    { std::string t = rng.bases(maxlen); cs.push_back({t, rng.bases(8), 20}); } // long target, short query
    { std::string t = rng.bases(8); cs.push_back({t, rng.bases(maxlen), 20}); } // short target, long query
    { std::string t = rng.bases(60); cs.push_back({t, rng.bases(60), 3}); }     // unrelated, low h0
    // short reads (GenDP's BSW regime)
    for (int k = 0; k < nr; ++k) {
        int q = rng.randint(10, std::min(100, maxlen - 20));
        std::string ref = rng.bases(q + 20);
        std::string read = mutate(ref, rng, 0.03, 0.03, 3).substr(0, q);
        cs.push_back({ref, read, rng.randint(1, 40)});
    }
    // long, drifting reads: stress the steering
    const double fracs[3] = {0.25, 0.5, 0.75};
    for (int k = 0; k < nr; ++k) {
        int q = rng.randint(maxlen / 2, maxlen - 30);
        std::string ref = rng.bases(maxlen);
        std::string read = mutate(ref, rng, 0.04, 0.08, 2, fracs[rng.randint(0, 2)]).substr(0, q);
        cs.push_back({ref, read, rng.randint(10, 40)});
    }
    // heavy one-sided drift: a static band falls off, the adaptive one should follow
    for (int k = 0; k < nr; ++k) {
        int q = rng.randint(maxlen / 2, maxlen - 60);
        std::string ref = rng.bases(maxlen);
        std::string read = mutate(ref, rng, 0.03, 0.12, 2, rng.randint(0, 1) ? 0.9 : 0.1).substr(0, q);
        if (rng.uniform() < 0.5) ref = ref.substr(0, q);   // insertions need room in the read, not the ref
        cs.push_back({ref, read, rng.randint(10, 40)});
    }
    }

    std::ofstream f(out);
    f << cs.size() * 2 << '\n';
    for (auto &c : cs) {
        if (c.t.empty() || c.q.empty() || (int)c.t.size() > maxlen || (int)c.q.size() > maxlen) {
            std::fprintf(stderr, "bad case length\n");
            return 1;
        }
        Seq tt = encode(c.t), qq = encode(c.q);
        for (int adaptive : {1, 0}) {
            ExtendResult r = extend(qq, tt, c.h0, adaptive ? Policy::Adaptive : Policy::Static, W);
            f << c.t.size() << ' ' << c.q.size() << ' ' << c.h0 << ' ' << adaptive << ' ' << r.score << ' '
              << r.gscore << ' ' << r.max_i << ' ' << r.max_j << '\n';
            for (size_t i = 0; i < tt.size(); ++i) f << (i ? " " : "") << (int)tt[i];
            f << '\n';
            for (size_t i = 0; i < qq.size(); ++i) f << (i ? " " : "") << (int)qq[i];
            f << '\n';
        }
    }
    std::printf("wrote %zu tests (W=%d, maxlen=%d) to %s\n", cs.size() * 2, W, maxlen, out.c_str());
    return 0;
}

void usage() {
    std::fprintf(stderr,
        "usage: golden selftest\n"
        "       golden validate <ksw.csv> <ksw_out.txt>\n"
        "       golden compare <pairs.tsv> [-W 8 16 32] [-n N]\n"
        "       golden run <pairs.tsv> <full|rowband|static|adaptive> [-W n]\n"
        "       golden pairs [-n N --qlen L --tlen-extra E --sub-rate R --indel-rate R --max-indel K\n"
        "                     --del-frac F --h0 H --seed S --max-len M] [--tsv F] [--ksw F] [--gendp F]\n"
        "       golden vectors [-o F] [-W n] [--maxlen L] [-n N] [--seed S] [--tsv pairs.tsv]\n");
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];
    if (cmd == "selftest") return cmd_selftest();
    if (cmd == "validate" && argc >= 4) return cmd_validate(argv[2], argv[3]);
    if (cmd == "compare" && argc >= 3) {
        std::vector<int> Ws{8, 16, 32};
        int limit = 0;
        for (int k = 3; k < argc; ++k) {
            if (!std::strcmp(argv[k], "-W")) Ws = parse_ints_after(argc, argv, k);
            else if (!std::strcmp(argv[k], "-n") && k + 1 < argc) limit = std::atoi(argv[++k]);
        }
        return cmd_compare(argv[2], Ws, limit);
    }
    if (cmd == "run" && argc >= 4) {
        Policy p;
        if (!parse_policy(argv[3], p)) { usage(); return 2; }
        int W = 16;
        for (int k = 4; k + 1 < argc; ++k) if (!std::strcmp(argv[k], "-W")) W = std::atoi(argv[k + 1]);
        return cmd_run(argv[2], p, W);
    }
    if (cmd == "pairs") return cmd_pairs(argc, argv);
    if (cmd == "vectors") return cmd_vectors(argc, argv);
    usage();
    return 2;
}
