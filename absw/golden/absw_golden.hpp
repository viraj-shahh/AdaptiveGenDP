// Golden model for banded and adaptive-banded affine-gap seed extension (C++17).
//
// Two layers:
//
// 1. ksw_extend2(): line-for-line port of ksw_extend2() in
//    gendp/kernel/bwa-mem/ksw.cpp (BWA-MEM seed extension, what GenDP's BSW kernel
//    computes).  Validated field-for-field against the C binary (`golden validate`).
//
// 2. extend(): anti-diagonal DP engine with the same recurrences and pluggable band
//    policies:
//      Full      every cell (ground truth)
//      RowBand   |i - j| <= w, ksw's static band (no zero-pruning, no z-drop)
//      Static    W cells per anti-diagonal, band alternates right/down
//      Adaptive  W cells per anti-diagonal, band moves toward the end with the higher
//                H score (Suzuki-Kasahara style; Liao et al., ASAP 2018)
//    Static and Adaptive use identical hardware; only the 1-bit move per
//    anti-diagonal differs.  extend() is bit-exact with rtl/absw_array.sv.
//
// Conventions (match ksw.cpp): 1-based matrix indices, i = target (reference)
// position, j = query (read) position.  Row 0 / column 0 hold the h0-derived boundary.
// Cells outside the band count as H = E = F = M = 0, as in ksw's band.
// Symbols: A=0 C=1 G=2 T=3, anything else 4 (N).
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace absw {

using Seq = std::vector<uint8_t>;

struct Scoring {
    int match = 1, mismatch = 4, ambig = 1;
    int o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
    int end_bonus = 5, zdrop = 100;
};

inline int subst(int a, int b, const Scoring &s) {
    if (a == 4 || b == 4) return -s.ambig;
    return a == b ? s.match : -s.mismatch;
}

inline Seq encode(const std::string &str) {
    Seq out;
    out.reserve(str.size());
    for (char c : str) {
        switch (c) {
            case 'A': case 'a': out.push_back(0); break;
            case 'C': case 'c': out.push_back(1); break;
            case 'G': case 'g': out.push_back(2); break;
            case 'T': case 't': case 'U': case 'u': out.push_back(3); break;
            default: out.push_back(4);
        }
    }
    return out;
}

// --------------------------------------------------------------------------
// 1. Exact port of ksw_extend2 (row-major, with BWA's pruning and z-drop)
// --------------------------------------------------------------------------
struct KswResult {
    int score, qle, tle, gtle, gscore, max_off;
    long cells;
};

inline KswResult ksw_extend2(const Seq &q, const Seq &t, int h0, int w, const Scoring &s = {},
                             bool zdrop_on = true) {
    const int qlen = (int)q.size(), tlen = (int)t.size();
    const int oe_del = s.o_del + s.e_del, oe_ins = s.o_ins + s.e_ins;
    const int zdrop = zdrop_on ? s.zdrop : 0;
    std::vector<int> eh_h(qlen + 1, 0), eh_e(qlen + 1, 0);
    eh_h[0] = h0;
    eh_h[1] = h0 > oe_ins ? h0 - oe_ins : 0;
    for (int j = 2; j <= qlen && eh_h[j - 1] > s.e_ins; ++j) eh_h[j] = eh_h[j - 1] - s.e_ins;
    const int mx = s.match;                    // largest entry of the scoring matrix
    int max_ins = (int)((double)(qlen * mx + s.end_bonus - s.o_ins) / s.e_ins + 1.);
    max_ins = max_ins > 1 ? max_ins : 1;
    w = w < max_ins ? w : max_ins;
    int max_del = (int)((double)(qlen * mx + s.end_bonus - s.o_del) / s.e_del + 1.);
    max_del = max_del > 1 ? max_del : 1;
    w = w < max_del ? w : max_del;

    int best = h0, max_i = -1, max_j = -1, max_ie = -1, gscore = -1, max_off = 0;
    int beg = 0, end = qlen;
    long cells = 0;
    for (int i = 0; i < tlen; ++i) {
        int f = 0, m = 0, mj = -1, h1, j;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > qlen) end = qlen;
        h1 = beg == 0 ? std::max(h0 - (s.o_del + s.e_del * (i + 1)), 0) : 0;
        for (j = beg; j < end; ++j) {
            ++cells;
            int M = eh_h[j], e = eh_e[j];
            eh_h[j] = h1;
            M = M ? M + subst(t[i], q[j], s) : 0;
            int h = std::max(std::max(M, e), f);
            mj = m > h ? mj : j;
            m = m > h ? m : h;
            int tt = std::max(M - oe_del, 0);
            e = std::max(e - s.e_del, tt);
            f = std::max(f - s.e_ins, tt);       // repo uses oe_del here too (o_ins line commented out)
            h1 = h;
            eh_e[j] = e;
        }
        eh_h[end] = h1;
        eh_e[end] = 0;
        if (j == qlen) {
            max_ie = gscore > h1 ? max_ie : i;
            gscore = gscore > h1 ? gscore : h1;
        }
        if (m == 0) break;
        if (m > best) {
            best = m; max_i = i; max_j = mj;
            max_off = std::max(max_off, std::abs(mj - i));
        } else if (zdrop > 0) {
            if (i - max_i > mj - max_j) {
                if (best - m - ((i - max_i) - (mj - max_j)) * s.e_del > zdrop) break;
            } else {
                if (best - m - ((mj - max_j) - (i - max_i)) * s.e_ins > zdrop) break;
            }
        }
        for (j = beg; j < end && eh_h[j] == 0 && eh_e[j] == 0; ++j) {}
        beg = j;
        for (j = end; j >= beg && eh_h[j] == 0 && eh_e[j] == 0; --j) {}
        end = j + 2 < qlen ? j + 2 : qlen;
    }
    return {best, max_j + 1, max_i + 1, max_ie + 1, gscore, max_off, cells};
}

// ksw_extend2 with z-drop off and zero-pruning off ("beg = 0; end = qlen" debug line).
// Returns {score, gscore}.  Used to validate extend(RowBand).
inline std::pair<int, int> ksw_extend2_noprune(const Seq &q, const Seq &t, int h0, int w,
                                                const Scoring &s = {}) {
    const int qlen = (int)q.size(), tlen = (int)t.size();
    const int oe = s.o_del + s.e_del;
    std::vector<int> eh_h(qlen + 1, 0), eh_e(qlen + 1, 0);
    eh_h[0] = h0;
    eh_h[1] = h0 > oe ? h0 - oe : 0;
    for (int j = 2; j <= qlen && eh_h[j - 1] > s.e_ins; ++j) eh_h[j] = eh_h[j - 1] - s.e_ins;
    int best = h0, gscore = -1;
    for (int i = 0; i < tlen; ++i) {
        int beg = std::max(0, i - w), end = std::min(qlen, i + w + 1);
        if (beg >= end) continue;
        int f = 0, m = 0;
        int h1 = beg == 0 ? std::max(h0 - (s.o_del + s.e_del * (i + 1)), 0) : 0;
        for (int j = beg; j < end; ++j) {
            int M = eh_h[j], e = eh_e[j];
            eh_h[j] = h1;
            M = M ? M + subst(t[i], q[j], s) : 0;
            int h = std::max(std::max(M, e), f);
            m = std::max(m, h);
            int tt = std::max(M - oe, 0);
            e = std::max(e - s.e_del, tt);
            f = std::max(f - s.e_ins, tt);
            h1 = h;
            eh_e[j] = e;
        }
        eh_h[end] = h1;
        eh_e[end] = 0;
        if (end == qlen) gscore = std::max(gscore, h1);
        best = std::max(best, m);
    }
    return {best, gscore};
}

// --------------------------------------------------------------------------
// 2. Anti-diagonal engine with band policies
// --------------------------------------------------------------------------
enum class Policy { Full, RowBand, Static, Adaptive };

inline const char *policy_name(Policy p) {
    switch (p) {
        case Policy::Full: return "full";
        case Policy::RowBand: return "rowband";
        case Policy::Static: return "static";
        default: return "adaptive";
    }
}

inline bool parse_policy(const std::string &s, Policy &p) {
    if (s == "full") p = Policy::Full;
    else if (s == "rowband") p = Policy::RowBand;
    else if (s == "static") p = Policy::Static;
    else if (s == "adaptive") p = Policy::Adaptive;
    else return false;
    return true;
}

struct ExtendResult {
    int score, max_i, max_j, gscore;
    long cells;
    std::vector<uint8_t> moves;     // per anti-diagonal: 1 = down, 0 = right (Static/Adaptive)
};

inline int floor_div(int a, int b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }

// W: PEs (cells per anti-diagonal) for Static/Adaptive.  w: half-width for RowBand
// (defaults to W - 1, the same cell budget as W cells per anti-diagonal).
inline ExtendResult extend(const Seq &q, const Seq &t, int h0, Policy policy, int W = 16,
                           int w = -1, const Scoring &s = {}) {
    const int qlen = (int)q.size(), tlen = (int)t.size();
    const int oe_del = s.o_del + s.e_del, oe_ins = s.o_ins + s.e_ins;
    if (w < 0) w = W - 1;
    const int cols = qlen + 1;
    std::vector<int> H((size_t)(tlen + 1) * cols, 0), E(H.size(), 0), F(H.size(), 0), Mx(H.size(), 0);
    auto at = [cols](int i, int j) { return (size_t)i * cols + j; };
    for (int j = 0; j <= qlen; ++j) H[at(0, j)] = std::max(0, h0 - s.o_ins - s.e_ins * j);
    for (int i = 0; i <= tlen; ++i) H[at(i, 0)] = std::max(0, h0 - s.o_del - s.e_del * i);
    H[at(0, 0)] = h0;

    ExtendResult r{h0, 0, 0, -1, 0, {}};
    int lo = 1 - (W - 1) / 2;                 // centred on cell (1,1)
    bool last_dir = true;
    for (int d = 2; d <= tlen + qlen; ++d) {
        const int i_min = std::max(1, d - qlen), i_max = std::min(tlen, d - 1);
        int a, b;
        if (policy == Policy::Full) {
            a = i_min; b = i_max;
        } else if (policy == Policy::RowBand) {   // |2i - d| <= w
            a = std::max(i_min, -floor_div(w - d, 2));
            b = std::min(i_max, floor_div(d + w, 2));
        } else {
            a = std::max(i_min, lo); b = std::min(i_max, lo + W - 1);
        }
        if (a <= b) {
            int diag_best = 0, diag_bi = -1;
            for (int i = a; i <= b; ++i) {
                const int j = d - i;
                const int hd = H[at(i - 1, j - 1)];
                const int M = hd > 0 ? hd + subst(t[i - 1], q[j - 1], s) : 0;
                const int e = i > 1 ? std::max(E[at(i - 1, j)] - s.e_del, std::max(Mx[at(i - 1, j)] - oe_del, 0)) : 0;
                const int f = j > 1 ? std::max(F[at(i, j - 1)] - s.e_ins, std::max(Mx[at(i, j - 1)] - oe_ins, 0)) : 0;
                const int h = std::max(std::max(M, e), f);
                Mx[at(i, j)] = M; E[at(i, j)] = e; F[at(i, j)] = f; H[at(i, j)] = h;
                if (diag_bi < 0 || h > diag_best) { diag_best = h; diag_bi = i; }  // first max, like argmax
            }
            r.cells += b - a + 1;
            if (diag_best > r.score) { r.score = diag_best; r.max_i = diag_bi; r.max_j = d - diag_bi; }
            if (d - a == qlen) r.gscore = std::max(r.gscore, H[at(a, qlen)]);
        }
        if (policy == Policy::Static || policy == Policy::Adaptive) {
            auto end_score = [&](int ii) {
                const int jj = d - ii;
                return (ii >= 1 && ii <= tlen && jj >= 1 && jj <= qlen) ? H[at(ii, jj)] : -1;
            };
            const int ht = end_score(lo), hb = end_score(lo + W - 1);
            bool down = (policy == Policy::Adaptive && ht != hb) ? (hb > ht) : !last_dir;
            if (lo + 1 > tlen) down = false;                       // never steer off the matrix
            if ((d + 1) - (lo + W - 1) > qlen) down = true;
            last_dir = down;
            r.moves.push_back(down ? 1 : 0);
            if (down) ++lo;
        }
    }
    return r;
}

// --------------------------------------------------------------------------
// Portable RNG (splitmix64): same streams with gcc and clang, unlike <random>
// distributions, so test vectors are identical on every machine.
// --------------------------------------------------------------------------
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double uniform() { return (double)(next() >> 11) * 0x1.0p-53; }
    int randint(int lo, int hi) { return lo + (int)(next() % (uint64_t)(hi - lo + 1)); }  // inclusive
    char base() { return "ACGT"[next() & 3]; }
    std::string bases(int n) { std::string r(n, 'A'); for (auto &c : r) c = base(); return r; }
};

// Read derived from seq with substitutions and indels (same algorithm as gen_pairs.py).
inline std::string mutate(const std::string &seq, Rng &rng, double sub_rate, double indel_rate,
                          int max_indel, double del_frac = 0.5) {
    std::string out;
    size_t i = 0;
    while (i < seq.size()) {
        const double r = rng.uniform();
        if (r < indel_rate * del_frac) {           // deletion from the read
            i += rng.randint(1, max_indel);
            continue;
        }
        if (r < indel_rate) {                       // insertion into the read
            const int n = rng.randint(1, max_indel);
            for (int k = 0; k < n; ++k) out.push_back(rng.base());
        }
        char c = seq[i];
        if (rng.uniform() < sub_rate) {
            char nc;
            do { nc = rng.base(); } while (nc == c);
            c = nc;
        }
        out.push_back(c);
        ++i;
    }
    return out;
}

}  // namespace absw
