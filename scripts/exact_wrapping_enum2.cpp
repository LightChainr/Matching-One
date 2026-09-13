// Independent exact enumeration of horizontal-wrapping counts for NN square
// SITE percolation on the L x L torus, per number of occupied sites.
//
// Second independent method for issue #576 Part 1 (the committed
// exact_wrapping_enum.cpp uses the doubled-grid criterion; this file uses a
// union-find with lattice displacement potentials on the torus itself, in the
// spirit of scripts/matched_torus_reference.py: a component winds in x iff a
// cycle closes with nonzero accumulated x-displacement).
//
// For each k = 0..L^2 we output c_k = #{configs with k occupied sites whose
// occupied subgraph has a component with x-winding != 0} (including configs
// that also wind in y).  Sum_k c_k = total wrapping configurations.
//
// Build:  g++ -O2 -o exact_wrapping_enum2 exact_wrapping_enum2.cpp
// Run:    ./exact_wrapping_enum2 5
//
// L <= 5 is 2^25 ~= 3.4e7 configurations; runs in minutes.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <numeric>
#include <algorithm>

using u64 = uint64_t;
using i64 = int64_t;

struct DSU {
    // nodes 0..N-1: torus sites; displacement potentials relative to root
    int n;
    std::vector<int> parent;
    std::vector<int> dx, dy; // position(child) - position(parent) accumulated
    std::vector<uint8_t> wrapx;

    void init(int n_) {
        n = n_;
        parent.resize(n);
        dx.assign(n, 0);
        dy.assign(n, 0);
        wrapx.assign(n, 0);
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x, int &ox, int &oy) {
        if (parent[x] == x) { ox = 0; oy = 0; return x; }
        int px, py;
        int r = find(parent[x], px, py);
        ox = dx[x] + px;
        oy = dy[x] + py;
        parent[x] = r;
        dx[x] = ox;
        dy[x] = oy;
        return r;
    }
    // edge: pos(j) = pos(i) + (ex, ey)
    void add_edge(int i, int j, int ex, int ey) {
        int ix, iy, jx, jy;
        int ri = find(i, ix, iy);
        int rj = find(j, jx, jy);
        if (ri == rj) {
            // cycle: potentials are unwrapped integers; nonzero mismatch means
            // the cycle has x-winding = mismatch/L != 0
            if (ix + ex - jx != 0) wrapx[ri] = 1;
            return;
        }
        // attach rj under ri: position(rj) should be pos(i)+(ex,ey) - pos(j)
        int rdx = ix + ex - jx;
        int rdy = iy + ey - jy;
        parent[rj] = ri;
        dx[rj] = rdx;
        dy[rj] = rdy;
        wrapx[ri] = wrapx[ri] | wrapx[rj];
    }
    static int Lmod;
};
int DSU::Lmod = 0;

int main(int argc, char **argv) {
    int Lmax = argc > 1 ? atoi(argv[1]) : 5;
    for (int L = 2; L <= Lmax; ++L) {
        int N = L * L;
        DSU::Lmod = L;
        auto id = [&](int r, int c) { return r * L + c; };
        // precompute edge list: horizontal (r,c)-(r,c+1) with (ex,ey)=(1,0);
        // vertical (r,c)-(r+1,c) with (0,1); periodic indices.
        std::vector<int> ea, eb, ex, ey;
        for (int r = 0; r < L; ++r)
            for (int c = 0; c < L; ++c) {
                ea.push_back(id(r, c)); eb.push_back(id(r, (c + 1) % L)); ex.push_back(1); ey.push_back(0);
                ea.push_back(id(r, c)); eb.push_back(id((r + 1) % L, c)); ex.push_back(0); ey.push_back(1);
            }
        int E = (int)ea.size();
        u64 total = u64(1) << N;
        std::vector<i64> cnt(N + 1, 0);
        DSU dsu;
        std::vector<uint8_t> occ(N);
        for (u64 mask = 0; mask < total; ++mask) {
            int k = 0;
            for (int i = 0; i < N; ++i) { occ[i] = (mask >> i) & 1; k += occ[i]; }
            bool wrapped = false;
            // early exit impossible in general; just run
            dsu.init(N);
            for (int e = 0; e < E && !wrapped; ++e) {
                if (occ[ea[e]] && occ[eb[e]]) {
                    dsu.add_edge(ea[e], eb[e], ex[e], ey[e]);
                    int rx, ry;
                    int r = dsu.find(ea[e], rx, ry);
                    if (dsu.wrapx[r]) wrapped = true;
                }
            }
            if (wrapped) cnt[k]++;
        }
        i64 sum = 0;
        for (int k = 0; k <= N; ++k) sum += cnt[k];
        // print polynomial
        printf("L=%d total_wrapping=%lld / 2^%d\n  c_k:", L, (long long)sum, N);
        for (int k = 0; k <= N; ++k) printf(" %lld", (long long)cnt[k]);
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
