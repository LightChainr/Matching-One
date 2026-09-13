// Faithful, allocation-free C++ port of the supplied one-frontier winding-component
// transfer (scripts/cylinder_winding_intensity.py). advance()/empty_state()/
// reward_lump() semantics are reproduced exactly; only the width limit and the
// dense Fraction linear algebra of the Python reference are changed.
//
// Output: the reward-preserving lumped chain aggregated by row popcount, i.e. all
// that is needed to evaluate nu_w(p) as an exact rational downstream.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <chrono>

using namespace std;

static inline int floordiv(int a, int b) {   // Python // semantics; C++ / truncates
    return (a >= 0) ? (a / b) : -(((-a) + b - 1) / b);
}

static const int MAXW = 16;
static const int MAXN = 4 * MAXW;      // 2*width DSU nodes

struct DSU {
    int parent[MAXN], delta[MAXN];
    uint8_t wind[MAXN];
    void init(int n) {
        for (int i = 0; i < n; ++i) { parent[i] = i; delta[i] = 0; wind[i] = 0; }
    }
    inline int find(int a, int &pot) {
        // iterative path compression
        int root = a, acc = 0;
        while (parent[root] != root) { acc += delta[root]; root = parent[root]; }
        // second pass: compress, accumulating potential relative to root
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
        else { parent[rb] = ra; delta[rb] = gain + da - db; wind[ra] |= wind[rb]; }
    }
};

struct State {
    int16_t labels[MAXW];   // -1 when unoccupied
    int16_t gains[MAXW];
    uint8_t flags[MAXW];
    int nflags;
};

static inline void pack_state(const State &s, int W, string &out) {
    out.resize((size_t)W * 2 + s.nflags);
    size_t p = 0;
    for (int i = 0; i < W; ++i) out[p++] = (char)(s.labels[i] & 0xFF);
    for (int i = 0; i < W; ++i) out[p++] = (char)(s.gains[i] & 0xFF);
    for (int i = 0; i < s.nflags; ++i) out[p++] = (char)s.flags[i];
}

struct KeyHash {
    inline size_t operator()(const string &k) const {
        size_t h = 1469598103934665603ULL;
        const unsigned char *p = (const unsigned char *)k.data();
        size_t n = k.size();
        for (size_t i = 0; i < n; ++i) { h ^= p[i]; h *= 1099511628211ULL; }
        return h;
    }
};

// Returns reward; writes the successor into `nx`.
static inline int advance(const State &st, int mask, int W, bool matching, State &nx) {
    static DSU dsu;
    static int oldv[MAXW], newv[MAXW];
    static int rep_of[MAXW];        // label -> first site
    static int root_of[MAXW];       // label -> root
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
        root_of[k] = r;
    }
    (void)root_of;
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
    int reward = 0;
    for (int a = 0, b = 0; a < na; ++a) {
        while (b < nk && keptroots[b] < allroots[a]) ++b;
        if (b >= nk || keptroots[b] != allroots[a]) reward += dsu.wind[allroots[a]] ? 1 : 0;
    }
    // build successor
    nx.nflags = 0;
    for (int i = 0; i < W; ++i) { nx.labels[i] = -1; nx.gains[i] = 0; }
    static int tag[MAXN], origin[MAXN];
    for (int i = 0; i < 2 * W; ++i) tag[i] = -1;
    for (int t = 0; t < nnew; ++t) {
        int i = newv[t];
        int pot; int r = dsu.find(W + i, pot);
        if (tag[r] < 0) { tag[r] = nx.nflags; origin[r] = pot; nx.flags[nx.nflags++] = dsu.wind[r] ? 1 : 0; }
        nx.labels[i] = (int16_t)tag[r];
        nx.gains[i] = (int16_t)(dsu.wind[r] ? 0 : pot - origin[r]);
    }
    return reward;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s WIDTH 0|1(matching) [out.json]\n", argv[0]); return 2; }
    int W = atoi(argv[1]);
    bool matching = atoi(argv[2]) != 0;
    const char *outpath = (argc > 3) ? argv[3] : nullptr;
    if (W < 2 || W > MAXW) { fprintf(stderr, "width out of range\n"); return 2; }
    auto t0 = chrono::steady_clock::now();

    State empty; for (int i = 0; i < W; ++i) { empty.labels[i] = -1; empty.gains[i] = 0; }
    empty.nflags = 0;

    vector<State> states; states.push_back(empty);
    unordered_map<string, int, KeyHash> index;
    index.reserve(1 << 20);
    { string k; pack_state(empty, W, k); index[k] = 0; }

    string kb;
    State nxt, cur;
    size_t cursor = 0;
    while (cursor < states.size()) {
        cur = states[cursor];
        for (int mask = 0; mask < (1 << W); ++mask) {
            advance(cur, mask, W, matching, nxt);
            pack_state(nxt, W, kb);
            auto it = index.find(kb);
            if (it == index.end()) {
                index.emplace(kb, (int)states.size());
                states.push_back(nxt);
            }
        }
        ++cursor;
        if ((cursor & 32767) == 0)
            fprintf(stderr, "  bfs states=%zu t=%.1fs\n", states.size(),
                    chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    }
    size_t nstates = states.size();
    double tBFS = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "BFS width=%d matching=%d states=%zu transitions=%zu t=%.1fs\n",
            W, (int)matching, nstates, nstates * (size_t)(1 << W), tBFS);

    // ---- cache the transition table so refinement/aggregation do not recompute ----
    const size_t MASKS = (size_t)1 << W;
    fprintf(stderr, "caching transition table: %zu entries (%.2f GB)\n",
            nstates * MASKS, (double)nstates * MASKS * 4 / 1073741824.0);
    vector<uint32_t> tab(nstates * MASKS);
    for (size_t i = 0; i < nstates; ++i) {
        cur = states[i];
        uint32_t *rowp = &tab[i * MASKS];
        for (int mask = 0; mask < (int)MASKS; ++mask) {
            int rw = advance(cur, mask, W, matching, nxt);
            pack_state(nxt, W, kb);
            int j = index.find(kb)->second;
            rowp[mask] = ((uint32_t)j << 3) | (uint32_t)(rw & 7);
        }
    }
    fprintf(stderr, "table cached t=%.1fs\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());

    vector<int> blocks(nstates, 0);
    static vector<int> cnt, stamp, touched;
    int curstamp = 0;
    int iter = 0;
    for (;;) {
        map<string, int> classes;
        vector<int> refined(nstates);
        int maxblock = *max_element(blocks.begin(), blocks.end()) + 1;
        size_t sz = (size_t)(W + 1) * 8 * maxblock;
        if ((int)cnt.size() < (int)sz) { cnt.assign(sz, 0); stamp.assign(sz, 0); curstamp = 0; }
        for (size_t i = 0; i < nstates; ++i) {
            ++curstamp;
            touched.clear();
            cur = states[i];
            const uint32_t *rowp = &tab[i * MASKS];
            for (int mask = 0; mask < (int)MASKS; ++mask) {
                uint32_t v = rowp[mask];
                int j = (int)(v >> 3), rw = (int)(v & 7);
                int pc = __builtin_popcount((unsigned)mask);
                size_t idx = ((size_t)pc * 8 + (size_t)rw) * maxblock + blocks[j];
                if (stamp[idx] != curstamp) { stamp[idx] = curstamp; cnt[idx] = 0; touched.push_back((int)idx); }
                cnt[idx] += 1;
            }
            sort(touched.begin(), touched.end());
            string sig;
            char buf[48];
            for (int idx : touched) {
                int b = idx % maxblock; int rest = idx / maxblock;
                int rw = rest % 8; int pc = rest / 8;
                snprintf(buf, sizeof buf, "%d.%d.%d=%d;", pc, rw, b, cnt[idx]);
                sig += buf;
            }
            auto it = classes.find(sig);
            if (it == classes.end()) { int id = (int)classes.size(); classes.emplace(sig, id); refined[i] = id; }
            else refined[i] = it->second;
        }
        ++iter;
        bool same = (refined == blocks);
        blocks = refined;
        int nb = *max_element(blocks.begin(), blocks.end()) + 1;
        fprintf(stderr, "  lump iter %d -> %d classes t=%.1fs\n", iter, nb,
                chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        if (same) break;
        if (iter > 40) { fprintf(stderr, "no stabilisation\n"); return 3; }
    }
    int nblocks = *max_element(blocks.begin(), blocks.end()) + 1;

    // one REPRESENTATIVE state per block, matching reward_lump's blocks.index(k)
    vector<int> first_of_block(nblocks, -1);
    for (size_t i = 0; i < nstates; ++i) { int b = blocks[i]; if (first_of_block[b] < 0) first_of_block[b] = (int)i; }
    vector<map<pair<int,int>, long long>> flat((size_t)nblocks * (W + 1));
    for (int b = 0; b < nblocks; ++b) {
        const uint32_t *rowp = &tab[(size_t)first_of_block[b] * MASKS];
        for (int mask = 0; mask < (int)MASKS; ++mask) {
            uint32_t v = rowp[mask];
            int j = (int)(v >> 3), rw = (int)(v & 7);
            int pc = __builtin_popcount((unsigned)mask);
            flat[(size_t)b * (W + 1) + pc][make_pair(blocks[j], rw)] += 1;
        }
    }

    string js;
    char buf[256];
    snprintf(buf, sizeof buf,
             "{\"width\":%d,\"matching\":%s,\"states\":%zu,\"lump_blocks\":%d,\"bfs_seconds\":%.3f,\"rows\":[",
             W, matching ? "true" : "false", nstates, nblocks, tBFS);
    js += buf;
    for (int b = 0; b < nblocks; ++b) {
        js += "[";
        for (int m = 0; m <= W; ++m) {
            js += "[";
            bool first = true;
            for (auto &kv : flat[(size_t)b * (W + 1) + m]) {
                if (!first) js += ",";
                first = false;
                snprintf(buf, sizeof buf, "[%d,%d,%lld]", kv.first.first, kv.first.second, kv.second);
                js += buf;
            }
            js += "]";
            if (m < W) js += ",";
        }
        js += "]";
        if (b + 1 < nblocks) js += ",";
    }
    js += "]}";
    if (outpath) { FILE *f = fopen(outpath, "w"); fputs(js.c_str(), f); fclose(f); }
    else fputs(js.c_str(), stdout);
    fprintf(stderr, "DONE width=%d matching=%d states=%zu blocks=%d total=%.1fs\n", W, (int)matching,
            nstates, nblocks, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    return 0;
}
