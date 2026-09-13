// Exact enumeration of the single-direction (horizontal) wrapping probability
// pi({1,0}) for nearest-neighbour square SITE percolation on an L x L torus.
//
// Independent ground-truth control for issue #576 Part 1.  Brute force over all
// 2^(L^2) configurations; for each we test whether a cluster winds once around
// the horizontal (period-1) cycle.  Method: double the grid horizontally to
// columns 0..2L-1 (columns L..2L-1 are a copy of 0..L-1), keep rows periodic.
// A cluster winds horizontally iff some component contains both (r,0) and (r,L)
// for one row r.
//
// Output: exact count / 2^(L^2) as a reduced fraction, plus decimal.
// L = 2,3,4 are fast; L = 5 (2^25 ~= 33.6M configs) runs in well under a minute
// in C++; L = 6 (2^36) is intentionally NOT attempted.

#include <cstdio>
#include <cstdint>
#include <vector>
#include <numeric>

using u64 = uint64_t;

static int L;
static int N;            // L*L sites in fundamental domain
static std::vector<int> nbr[4]; // up to 4 neighbours, -1 if none

// visited bitmask over the doubled grid (rows L, cols 2L): index = r*(2L)+c
static inline int didx(int r, int c) { return r * (2 * L) + c; }

// BFS queue (indices into doubled grid; size L*2L)
static std::vector<int> queue;
static u64 visited; // bitmask, needs <= L*2L bits; for L<=6 that is <=72 bits -> use two u64

static inline void set_vis(int i) { if (i < 64) visited |= (u64(1) << i); else visited |= (u64(1) << (i - 64)); }
static inline bool get_vis(int i) { return i < 64 ? (visited >> i) & 1 : (visited >> (i - 64)) & 1; }

// config occupancy for fundamental domain: bit b set means (r,c) occupied where
// idx = r*L + c.  Doubled grid occupancy for col>=L is copy of col-L.

static inline bool occ(int r, int c, u64 cfg) {
    if (c < L) return (cfg >> (r * L + c)) & 1;
    return (cfg >> (r * L + (c - L))) & 1;
}

static bool wraps(u64 cfg) {
    visited = 0;
    int total = L * (2 * L);
    for (int start = 0; start < total; ++start) {
        int sr = start / (2 * L);
        int sc = start % (2 * L);
        if (!occ(sr, sc, cfg)) continue;
        if (get_vis(start)) continue;
        // BFS this component
        queue.clear();
        queue.push_back(start);
        set_vis(start);
        bool found_wrap = false;
        // also record whether component contains a (r,0) and a (r,L)
        bool has0[1] = {false}; // placeholder, handled below
        // We detect wrap by checking, for the component, whether it contains
        // both (r,0) and (r,L) for some row r.  Do it after BFS via a second pass
        // over the collected component cells.
        std::vector<int> comp;
        while (!queue.empty()) {
            int cur = queue.back(); queue.pop_back();
            comp.push_back(cur);
            int r = cur / (2 * L);
            int c = cur % (2 * L);
            // neighbours
            // up
            int up = (r == 0 ? L - 1 : r - 1);
            if (occ(up, c, cfg) && !get_vis(didx(up, c))) { set_vis(didx(up, c)); queue.push_back(didx(up, c)); }
            int dn = (r == L - 1 ? 0 : r + 1);
            if (occ(dn, c, cfg) && !get_vis(didx(dn, c))) { set_vis(didx(dn, c)); queue.push_back(didx(dn, c)); }
            // left
            int lc = c - 1;
            if (lc >= 0 && occ(r, lc, cfg) && !get_vis(didx(r, lc))) { set_vis(didx(r, lc)); queue.push_back(didx(r, lc)); }
            // right
            int rc = c + 1;
            if (rc < 2 * L && occ(r, rc, cfg) && !get_vis(didx(r, rc))) { set_vis(didx(r, rc)); queue.push_back(didx(r, rc)); }
        }
        // check wrap: any row r with both (r,0) and (r,L) in component
        for (int r = 0; r < L; ++r) {
            bool a0 = false, aL = false;
            for (int cell : comp) {
                if (cell / (2 * L) != r) continue;
                int cc = cell % (2 * L);
                if (cc == 0) a0 = true;
                else if (cc == L) aL = true;
            }
            if (a0 && aL) { found_wrap = true; break; }
        }
        if (found_wrap) return true;
    }
    return false;
}

int main(int argc, char** argv) {
    std::vector<int> sizes;
    for (int i = 2; i <= 5; ++i) sizes.push_back(i);
    for (int sz : sizes) {
        L = sz; N = L * L;
        u64 total_cfg = u64(1) << N;
        u64 cnt = 0;
        for (u64 cfg = 0; cfg < total_cfg; ++cfg) {
            if (wraps(cfg)) ++cnt;
        }
        // reduce fraction
        u64 g = std::gcd(cnt, total_cfg);
        u64 num = cnt / g;
        u64 den = total_cfg / g;
        double p = double(cnt) / double(total_cfg);
        printf("L=%d  count=%llu / 2^%d = %llu/%llu  p=%.12f\n",
               L, (unsigned long long)cnt, N,
               (unsigned long long)num, (unsigned long long)den, p);
        fflush(stdout);
    }
    return 0;
}
