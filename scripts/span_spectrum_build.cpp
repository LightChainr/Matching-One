// span_build.cpp — one-frontier winding-component transfer with VERTICAL SPAN tracking.
//
// Port of the validated winding_build.cpp engine (same advance()/DSU/gain semantics,
// including the floordiv fix for the matching diagonal), extended with a per-component
// min-row offset so that each retiring winding component reports its vertical span.
//
// State = frontier labels + horizontal-lift gains + winding flags + per-component
// DEPTH = (own row) - (min row of the component), clamped at D_MAX+1.
//
// Transition r -> r+1, per old component (depth d, all its sites at rows <= r):
//   mo := d + 1   (its min row relative to the NEW row r+1)
//   merge rule:   mo_merged = max(mo_a, mo_b)   [min_row of the union]
//   retiring:     span = mo                     [last row r, first row r+1-mo]
//   persisting:   new depth = mo
//   newborn:      depth = 0
// Retiring WINDING components are binned by span: bins 1..D_MAX exact, bin D_MAX+1
// collects span >= D_MAX+1 (the tail). The manuscript's bound (9),
//   0 <= nu_w - nu_(w,<=H) <= w [1-(1-p)^w]^H,
// bounds the omitted tail rigorously.
//
// Output (binary): header, then one fixed record per (state, mask):
//   uint32 next_state; uint8 nret; uint8 bins[8]; uint8 cnts[8]  (24 bytes)
// Record index = state * 2^W + mask, so popcount(mask) is known to the reader.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;

static inline int floordiv(int a, int b) {   // Python // semantics
    return (a >= 0) ? (a / b) : -(((-a) + b - 1) / b);
}

static const int MAXW = 12;
static const int MAXN = 4 * MAXW;
static int D_MAX = 40;                       // runtime

struct DSU {
    int parent[MAXN], delta[MAXN];
    uint8_t wind[MAXN];
    int16_t mo[MAXN];
    void init(int n) {
        for (int i = 0; i < n; ++i) { parent[i] = i; delta[i] = 0; wind[i] = 0; mo[i] = 0; }
    }
    inline int find(int a, int &pot) {
        int root = a, acc = 0;
        while (parent[root] != root) { acc += delta[root]; root = parent[root]; }
        int cur = a, run = 0;
        while (parent[cur] != cur) {
            int nxt = parent[cur];
            int d = delta[cur];
            parent[cur] = root;
            delta[cur] = acc - run;
            run += d;
            cur = nxt;
        }
        pot = (a == root) ? 0 : delta[a];
        return root;
    }
    inline void join(int a, int b, int gain) {
        int da, db;
        int ra = find(a, da), rb = find(b, db);
        if (ra == rb) { wind[ra] |= (uint8_t)(db - da != gain); }
        else {
            parent[rb] = ra; delta[rb] = gain + da - db;
            wind[ra] |= wind[rb];
            if (mo[rb] > mo[ra]) mo[ra] = mo[rb];      // min_row of the union
        }
    }
};

struct State {
    int16_t labels[MAXW];
    int16_t gains[MAXW];
    uint8_t flags[MAXW];
    uint8_t depths[MAXW];
    int nflags;
};

static inline void pack_state(const State &s, int W, string &out) {
    out.resize((size_t)W * 2 + s.nflags * 2);
    size_t p = 0;
    for (int i = 0; i < W; ++i) out[p++] = (char)(s.labels[i] & 0xFF);
    for (int i = 0; i < W; ++i) out[p++] = (char)(s.gains[i] & 0xFF);
    for (int i = 0; i < s.nflags; ++i) { out[p++] = (char)s.flags[i]; out[p++] = (char)s.depths[i]; }
}

struct KeyHash {
    inline size_t operator()(const string &k) const {
        size_t h = 1469598103934665603ULL;
        const unsigned char *p = (const unsigned char *)k.data();
        for (size_t i = 0; i < k.size(); ++i) { h ^= p[i]; h *= 1099511628211ULL; }
        return h;
    }
};

static inline int advance(const State &st, int mask, int W, bool matching, State &nx,
                          uint8_t *ret_bins, uint8_t *ret_cnts, int &nret) {
    static DSU dsu;
    static int oldv[MAXW], newv[MAXW];
    static int rep_of[MAXW];
    static int allroots[MAXN], keptroots[MAXN];
    int nold = 0, nnew = 0;
    for (int i = 0; i < W; ++i) if (st.labels[i] >= 0) oldv[nold++] = i;
    for (int i = 0; i < W; ++i) if ((mask >> i) & 1) newv[nnew++] = i;

    dsu.init(2 * W);
    for (int t = 0; t < W; ++t) rep_of[t] = -1;
    for (int t = 0; t < nold; ++t) {
        int i = oldv[t], k = st.labels[i];
        if (rep_of[k] < 0) rep_of[k] = i;
        else dsu.join(rep_of[k], i, st.gains[i]);
    }
    for (int k = 0; k < W; ++k) if (rep_of[k] >= 0) {
        int pot; int r = dsu.find(rep_of[k], pot);
        dsu.wind[r] = (uint8_t)(st.flags[k] != 0);
        dsu.mo[r] = (int16_t)(st.depths[k] + 1);      // min row relative to the NEW row
    }
    for (int t = 0; t < nnew; ++t) {
        int i = newv[t];
        int j = i + 1; if (j == W) j = 0;
        if ((mask >> j) & 1) dsu.join(W + i, W + j, (i + 1) / W);
        int lo = matching ? -1 : 0, hi = matching ? 1 : 0;
        for (int dx = lo; dx <= hi; ++dx) {
            int jj = (i + dx) % W; if (jj < 0) jj += W;
            if (st.labels[jj] >= 0) dsu.join(W + i, jj, floordiv(i + dx, W));
        }
    }
    int na = 0, nk = 0;
    for (int t = 0; t < nold; ++t) { int p2; allroots[na++] = dsu.find(oldv[t], p2); }
    for (int t = 0; t < nnew; ++t) { int p2; keptroots[nk++] = dsu.find(W + newv[t], p2); }
    sort(allroots, allroots + na); na = (int)(unique(allroots, allroots + na) - allroots);
    sort(keptroots, keptroots + nk); nk = (int)(unique(keptroots, keptroots + nk) - keptroots);

    nret = 0;
    for (int a = 0, b = 0; a < na; ++a) {
        while (b < nk && keptroots[b] < allroots[a]) ++b;
        if (b >= nk || keptroots[b] != allroots[a]) {
            if (dsu.wind[allroots[a]]) {
                int bin = dsu.mo[allroots[a]];
                if (bin > D_MAX + 1) bin = D_MAX + 1;   // tail bin
                int found = -1;
                for (int q = 0; q < nret; ++q) if (ret_bins[q] == bin) { found = q; break; }
                if (found < 0 && nret < 8) { ret_bins[nret] = (uint8_t)bin; ret_cnts[nret] = 0; found = nret++; }
                if (found >= 0) ret_cnts[found] += 1;
            }
        }
    }
    // successor
    nx.nflags = 0;
    for (int i = 0; i < W; ++i) { nx.labels[i] = -1; nx.gains[i] = 0; }
    static int tag[MAXN], origin[MAXN];
    for (int i = 0; i < 2 * W; ++i) tag[i] = -1;
    for (int t = 0; t < nnew; ++t) {
        int i = newv[t];
        int pot; int r = dsu.find(W + i, pot);
        if (tag[r] < 0) {
            tag[r] = nx.nflags; origin[r] = pot;
            nx.flags[nx.nflags] = dsu.wind[r] ? 1 : 0;
            nx.depths[nx.nflags] = (uint8_t)min(dsu.mo[r], (int16_t)(D_MAX + 1));
            nx.nflags++;
        }
        nx.labels[i] = (int16_t)tag[r];
        nx.gains[i] = (int16_t)(dsu.wind[r] ? 0 : pot - origin[r]);
    }
    return nret;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s WIDTH 0|1(matching) D_MAX out.bin\n", argv[0]); return 2; }
    int W = atoi(argv[1]);
    bool matching = atoi(argv[2]) != 0;
    D_MAX = atoi(argv[3]);
    const char *outpath = argv[4];
    if (W < 2 || W > 8) { fprintf(stderr, "width must be 2..8 (8 retire-bin slots)\n"); return 2; }
    auto t0 = chrono::steady_clock::now();

    State empty; for (int i = 0; i < W; ++i) { empty.labels[i] = -1; empty.gains[i] = 0; }
    empty.nflags = 0;

    vector<State> states; states.push_back(empty);
    unordered_map<string, int, KeyHash> index;
    index.reserve(1 << 22);
    { string k; pack_state(empty, W, k); index[k] = 0; }

    string kb;
    State nxt, cur;
    size_t cursor = 0;
    while (cursor < states.size()) {
        cur = states[cursor];
        for (int mask = 0; mask < (1 << W); ++mask) {
            uint8_t rb[8], rc[8]; int nr;
            advance(cur, mask, W, matching, nxt, rb, rc, nr);
            pack_state(nxt, W, kb);
            auto it = index.find(kb);
            if (it == index.end()) {
                index.emplace(kb, (int)states.size());
                states.push_back(nxt);
            }
        }
        ++cursor;
        if ((cursor & 16383) == 0)
            fprintf(stderr, "  bfs states=%zu t=%.1fs\n", states.size(),
                    chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    }
    size_t nstates = states.size();
    fprintf(stderr, "BFS W=%d matching=%d D_MAX=%d states=%zu t=%.1fs\n", W, (int)matching, D_MAX,
            nstates, chrono::duration<double>(chrono::steady_clock::now() - t0).count());

    // transition table + binary output in one pass
    FILE *f = fopen(outpath, "wb");
    if (!f) { fprintf(stderr, "cannot open output\n"); return 3; }
    int32_t hdr[4] = { W, matching ? 1 : 0, D_MAX, (int32_t)nstates };
    fwrite(hdr, 4, 4, f);
    const size_t MASKS = (size_t)1 << W;
    vector<unsigned char> rec(24);
    for (size_t s = 0; s < nstates; ++s) {
        cur = states[s];
        for (int mask = 0; mask < (int)MASKS; ++mask) {
            uint8_t rb[8], rc[8]; int nr;
            advance(cur, mask, W, matching, nxt, rb, rc, nr);
            pack_state(nxt, W, kb);
            uint32_t nj = (uint32_t)index.find(kb)->second;
            memset(rec.data(), 0, 24);
            memcpy(rec.data(), &nj, 4);
            rec[4] = (uint8_t)nr;
            for (int q = 0; q < nr && q < 8; ++q) { rec[5 + q] = rb[q]; rec[13 + q] = rc[q]; }
            fwrite(rec.data(), 1, 24, f);
        }
        if ((s & 8191) == 0)
            fprintf(stderr, "  table %zu/%zu t=%.1fs\n", s, nstates,
                    chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    }
    fclose(f);
    fprintf(stderr, "DONE W=%d matching=%d D_MAX=%d states=%zu total=%.1fs -> %s\n", W, (int)matching,
            D_MAX, nstates, chrono::duration<double>(chrono::steady_clock::now() - t0).count(), outpath);
    return 0;
}
