// Find the configurations at L=5 where the doubled-grid "same row, columns 0
// and L" criterion (as implemented in exact_wrapping_enum.cpp) disagrees with
// the displacement-union-find ground truth.
//
// Result of the full 2^25 scan: exactly 10 configurations disagree, all of the
// form dsu=1 / doubled=0.  These are spiral configurations (the occupied
// cluster winds in x AND in y simultaneously and has no pure-x winding cycle),
// so its lift to the doubled strip contains no same-row (r,0),(r,L) pair.
// This resolves the L=5 discrepancy of issue #576: the published coefficient
// sum 8853301 is correct; the same-row enumerator undercounts by these 10.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <queue>

using u64 = uint64_t;
static const int L = 5;
static const int N = L * L;

static inline bool occ(int r, int c, u64 cfg) { r %= L;
    return (cfg >> (((r % L) * L + (c % L)))) & 1;
}

// ground truth: DSU with displacement
static int par[N], gx[N], gy[N];
static int find(int x, int &ox, int &oy) {
    if (par[x] == x) { ox = 0; oy = 0; return x; }
    int px, py; int r = find(par[x], px, py);
    ox = gx[x] + px; oy = gy[x] + py;
    par[x] = r; gx[x] = ox; gy[x] = oy;
    return r;
}
static bool wraps_dsu(u64 cfg) {
    for (int i = 0; i < N; ++i) { par[i] = i; gx[i] = gy[i] = 0; }
    for (int r = 0; r < L; ++r)
        for (int c = 0; c < L; ++c) {
            int i = r * L + c;
            if (!occ(r, c, cfg)) continue;
            // horizontal neighbour (r, c+1)
            if (occ(r, c + 1, cfg)) {
                int j = r * L + (c + 1) % L;
                int ix, iy, jx, jy;
                int ri = find(i, ix, iy), rj = find(j, jx, jy);
                if (ri == rj) { if (ix + 1 - jx != 0) return true; }
                else { par[rj] = ri; gx[rj] = ix + 1 - jx; gy[rj] = iy - jy; }
            }
            if (occ(r + 1, c, cfg)) {
                int j = ((r + 1) % L) * L + c;
                int ix, iy, jx, jy;
                int ri = find(i, ix, iy), rj = find(j, jx, jy);
                if (ri == rj) { if (ix - jx != 0) return true; } // dy winding ignored
                else { par[rj] = ri; gx[rj] = ix - jx; gy[rj] = iy + 1 - jy; }
            }
        }
    return false;
}

// doubled-grid criterion exactly as in exact_wrapping_enum.cpp:
// some component contains both (r,0) and (r,L) for one row r
static bool wraps_double(u64 cfg) {
    int W = 2 * L;
    std::vector<uint8_t> vis(L * W, 0);
    for (int start = 0; start < L * W; ++start) {
        int sr = start / W, sc = start % W;
        if (!occ(sr, sc, cfg) || vis[start]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(start); vis[start] = 1;
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            comp.push_back(cur);
            int r = cur / W, c = cur % W;
            int ups[4][2] = {{(r == 0 ? L - 1 : r - 1), c}, {(r == L - 1 ? 0 : r + 1), c},
                             {r, c - 1}, {r, c + 1}};
            for (auto &nb : ups) {
                int nr = nb[0], nc = nb[1];
                if (nc < 0 || nc >= W) continue;
                if (!occ(nr, nc, cfg)) continue;
                int idx = nr * W + nc;
                if (!vis[idx]) { vis[idx] = 1; q.push(idx); }
            }
        }
        for (int r = 0; r < L; ++r) {
            bool a0 = false, aL = false;
            for (int cell : comp) {
                if (cell / W != r) continue;
                if (cell % W == 0) a0 = true;
                else if (cell % W == L) aL = true;
            }
            if (a0 && aL) return true;
        }
    }
    return false;
}

int main() {
    u64 total = u64(1) << N;
    long long diff = 0;
    for (u64 cfg = 0; cfg < total; ++cfg) {
        bool w1 = wraps_dsu(cfg), w2 = wraps_double(cfg);
        if (w1 != w2) {
            ++diff;
            if (diff <= 12) {
                printf("cfg #%lld: dsu=%d double=%d occupied:", (long long)cfg, w1, w2);
                for (int i = 0; i < N; ++i) if ((cfg >> i) & 1) printf(" (%d,%d)", i / L, i % L);
                printf("\n");
            }
        }
    }
    printf("total differing configs: %lld\n", diff);
    return 0;
}
