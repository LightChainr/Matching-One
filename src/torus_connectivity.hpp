#pragma once

// Displacement-aware connectivity on an LxL torus for square-site percolation
// and its NN+NNN matching lattice.
//
// Wrapping convention (Mertens-Ziff / notes/aggressive-research-program.md):
//   horizontal (h): some occupied cluster has nonzero x-winding
//   vertical   (c): some occupied cluster has nonzero y-winding
//   either     (e): horizontal or vertical
//   both       (b): horizontal and vertical
// Closed-path displacement is the covering-space displacement of a cycle;
// a nonzero multiple of the period L in x (resp. y) is H (resp. V) wrapping.
// Periodic edges store the original lattice displacement, not the folded one.
//
// Site index i = y*L + x. Configuration bit i = 1 occupies the site on G
// (square NN) and bit i = 0 occupies it on G* (NN+NNN matching).

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace matching {

struct Edge {
    int to;
    int dx;
    int dy;
};

struct Graph {
    int L = 0;
    int n = 0;
    std::vector<std::vector<Edge>> adj;
};

inline Graph make_graph(int L, bool matching) {
    Graph g;
    g.L = L;
    g.n = L * L;
    g.adj.assign(static_cast<std::size_t>(g.n), {});
    const std::array<std::pair<int, int>, 4> offsets = {{
        {1, 0},
        {0, 1},
        {1, 1},
        {1, -1},
    }};
    const int n_off = matching ? 4 : 2;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = y * L + x;
            for (int k = 0; k < n_off; ++k) {
                const int dx = offsets[static_cast<std::size_t>(k)].first;
                const int dy = offsets[static_cast<std::size_t>(k)].second;
                int nx = x + dx;
                int ny = y + dy;
                // Fold coordinates but KEEP original (dx, dy).
                nx %= L;
                if (nx < 0) {
                    nx += L;
                }
                ny %= L;
                if (ny < 0) {
                    ny += L;
                }
                const int j = ny * L + nx;
                g.adj[static_cast<std::size_t>(i)].push_back(Edge{j, dx, dy});
            }
        }
    }
    return g;
}

// parent-relative displacement: pos(x) - pos(parent(x)) in covering space.
struct DisplacementDSU {
    int n = 0;
    std::vector<int> parent;
    std::vector<int> sz;
    std::vector<int> dx;
    std::vector<int> dy;
    std::vector<uint8_t> wrap_h;
    std::vector<uint8_t> wrap_v;

    explicit DisplacementDSU(int n_sites = 0) { reset(n_sites); }

    void reset(int n_sites) {
        n = n_sites;
        parent.resize(static_cast<std::size_t>(n));
        sz.resize(static_cast<std::size_t>(n));
        dx.resize(static_cast<std::size_t>(n));
        dy.resize(static_cast<std::size_t>(n));
        wrap_h.resize(static_cast<std::size_t>(n));
        wrap_v.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            parent[static_cast<std::size_t>(i)] = i;
            sz[static_cast<std::size_t>(i)] = 1;
            dx[static_cast<std::size_t>(i)] = 0;
            dy[static_cast<std::size_t>(i)] = 0;
            wrap_h[static_cast<std::size_t>(i)] = 0;
            wrap_v[static_cast<std::size_t>(i)] = 0;
        }
    }

    struct FindResult {
        int root;
        int dx;
        int dy;
    };

    FindResult find(int x) {
        if (parent[static_cast<std::size_t>(x)] == x) {
            return FindResult{x, 0, 0};
        }
        const FindResult p = find(parent[static_cast<std::size_t>(x)]);
        const int ndx = dx[static_cast<std::size_t>(x)] + p.dx;
        const int ndy = dy[static_cast<std::size_t>(x)] + p.dy;
        parent[static_cast<std::size_t>(x)] = p.root;
        dx[static_cast<std::size_t>(x)] = ndx;
        dy[static_cast<std::size_t>(x)] = ndy;
        return FindResult{p.root, ndx, ndy};
    }

    // Edge satisfies pos(j) = pos(i) + (edge_dx, edge_dy) in covering space.
    void unite(int i, int j, int edge_dx, int edge_dy) {
        const FindResult a = find(i);
        const FindResult b = find(j);
        const int root_dx = a.dx + edge_dx - b.dx;  // pos(root_j) - pos(root_i)
        const int root_dy = a.dy + edge_dy - b.dy;
        if (a.root == b.root) {
            if (root_dx != 0) {
                wrap_h[static_cast<std::size_t>(a.root)] = 1;
            }
            if (root_dy != 0) {
                wrap_v[static_cast<std::size_t>(a.root)] = 1;
            }
            return;
        }
        if (sz[static_cast<std::size_t>(a.root)] >= sz[static_cast<std::size_t>(b.root)]) {
            parent[static_cast<std::size_t>(b.root)] = a.root;
            dx[static_cast<std::size_t>(b.root)] = root_dx;
            dy[static_cast<std::size_t>(b.root)] = root_dy;
            sz[static_cast<std::size_t>(a.root)] += sz[static_cast<std::size_t>(b.root)];
            wrap_h[static_cast<std::size_t>(a.root)] = static_cast<uint8_t>(
                wrap_h[static_cast<std::size_t>(a.root)] | wrap_h[static_cast<std::size_t>(b.root)]);
            wrap_v[static_cast<std::size_t>(a.root)] = static_cast<uint8_t>(
                wrap_v[static_cast<std::size_t>(a.root)] | wrap_v[static_cast<std::size_t>(b.root)]);
        } else {
            parent[static_cast<std::size_t>(a.root)] = b.root;
            dx[static_cast<std::size_t>(a.root)] = -root_dx;
            dy[static_cast<std::size_t>(a.root)] = -root_dy;
            sz[static_cast<std::size_t>(b.root)] += sz[static_cast<std::size_t>(a.root)];
            wrap_h[static_cast<std::size_t>(b.root)] = static_cast<uint8_t>(
                wrap_h[static_cast<std::size_t>(a.root)] | wrap_h[static_cast<std::size_t>(b.root)]);
            wrap_v[static_cast<std::size_t>(b.root)] = static_cast<uint8_t>(
                wrap_v[static_cast<std::size_t>(a.root)] | wrap_v[static_cast<std::size_t>(b.root)]);
        }
    }
};

struct Observables {
    int k = 0;
    int clusters_G = 0;
    int clusters_Gstar = 0;
    uint8_t H_G = 0;
    uint8_t V_G = 0;
    uint8_t E_G = 0;
    uint8_t B_G = 0;
    uint8_t H_Gstar = 0;
    uint8_t V_Gstar = 0;
    uint8_t E_Gstar = 0;
    uint8_t B_Gstar = 0;
};

inline int popcount_mask(uint64_t mask) {
    return static_cast<int>(__builtin_popcountll(mask));
}

inline bool bit_set(uint64_t mask, int i) {
    return ((mask >> i) & 1ull) != 0ull;
}

// Occupied predicate: for G, occupied iff bit=1; for G*, occupied iff bit=0.
inline void analyze_one_graph(int n, uint64_t mask, bool bit_is_occupied, const Graph& g,
                              DisplacementDSU& dsu, int& n_clusters, uint8_t& wrap_h,
                              uint8_t& wrap_v) {
    dsu.reset(n);
    for (int i = 0; i < n; ++i) {
        const bool occ = bit_set(mask, i) == bit_is_occupied;
        if (!occ) {
            continue;
        }
        for (const Edge& e : g.adj[static_cast<std::size_t>(i)]) {
            const bool occ_to = bit_set(mask, e.to) == bit_is_occupied;
            if (occ_to) {
                dsu.unite(i, e.to, e.dx, e.dy);
            }
        }
    }
    n_clusters = 0;
    wrap_h = 0;
    wrap_v = 0;
    for (int i = 0; i < n; ++i) {
        const bool occ = bit_set(mask, i) == bit_is_occupied;
        if (!occ) {
            continue;
        }
        const auto r = dsu.find(i);
        if (r.root == i) {
            ++n_clusters;
            wrap_h = static_cast<uint8_t>(wrap_h | dsu.wrap_h[static_cast<std::size_t>(i)]);
            wrap_v = static_cast<uint8_t>(wrap_v | dsu.wrap_v[static_cast<std::size_t>(i)]);
        }
    }
}

inline Observables analyze_config(uint64_t mask, const Graph& g, const Graph& gstar,
                                  DisplacementDSU& dsu) {
    Observables o;
    o.k = popcount_mask(mask);
    uint8_t h = 0;
    uint8_t v = 0;
    analyze_one_graph(g.n, mask, true, g, dsu, o.clusters_G, h, v);
    o.H_G = h;
    o.V_G = v;
    o.E_G = static_cast<uint8_t>(h | v);
    o.B_G = static_cast<uint8_t>(h & v);
    analyze_one_graph(gstar.n, mask, false, gstar, dsu, o.clusters_Gstar, h, v);
    o.H_Gstar = h;
    o.V_Gstar = v;
    o.E_Gstar = static_cast<uint8_t>(h | v);
    o.B_Gstar = static_cast<uint8_t>(h & v);
    return o;
}

// V - E + F0 on the occupied NN subgraph of G (black sites). E counts every
// one-sided NN adjacency among occupied sites, which for L=2 includes both
// covering-space images of a periodic pair.
inline int euler_vef0(int L, uint64_t mask, const Graph& g_nn) {
    const int n = L * L;
    const int V = popcount_mask(mask);
    int E = 0;
    for (int i = 0; i < n; ++i) {
        if (!bit_set(mask, i)) {
            continue;
        }
        for (const Edge& e : g_nn.adj[static_cast<std::size_t>(i)]) {
            if (bit_set(mask, e.to)) {
                ++E;
            }
        }
    }
    int F0 = 0;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int a = y * L + x;
            const int b = y * L + ((x + 1) % L);
            const int c = ((y + 1) % L) * L + x;
            const int d = ((y + 1) % L) * L + ((x + 1) % L);
            if (bit_set(mask, a) && bit_set(mask, b) && bit_set(mask, c) && bit_set(mask, d)) {
                ++F0;
            }
        }
    }
    return V - E + F0;
}

struct Microcanonical {
    int L = 0;
    int n = 0;
    std::vector<uint64_t> configuration_count;
    std::vector<uint64_t> sum_clusters_G;
    std::vector<uint64_t> sum_clusters_Gstar_complement;
    std::vector<uint64_t> wrap_H_G;
    std::vector<uint64_t> wrap_V_G;
    std::vector<uint64_t> wrap_E_G;
    std::vector<uint64_t> wrap_B_G;
    std::vector<uint64_t> wrap_H_Gstar;
    std::vector<uint64_t> wrap_V_Gstar;
    std::vector<uint64_t> wrap_E_Gstar;
    std::vector<uint64_t> wrap_B_Gstar;
    uint64_t euler_mismatches = 0;
    uint64_t first_euler_mask = ~0ull;

    explicit Microcanonical(int L_) : L(L_), n(L_ * L_) {
        const std::size_t m = static_cast<std::size_t>(n + 1);
        configuration_count.assign(m, 0);
        sum_clusters_G.assign(m, 0);
        sum_clusters_Gstar_complement.assign(m, 0);
        wrap_H_G.assign(m, 0);
        wrap_V_G.assign(m, 0);
        wrap_E_G.assign(m, 0);
        wrap_B_G.assign(m, 0);
        wrap_H_Gstar.assign(m, 0);
        wrap_V_Gstar.assign(m, 0);
        wrap_E_Gstar.assign(m, 0);
        wrap_B_Gstar.assign(m, 0);
    }

    void add(const Observables& o) {
        const std::size_t k = static_cast<std::size_t>(o.k);
        configuration_count[k] += 1;
        sum_clusters_G[k] += static_cast<uint64_t>(o.clusters_G);
        sum_clusters_Gstar_complement[k] += static_cast<uint64_t>(o.clusters_Gstar);
        wrap_H_G[k] += o.H_G;
        wrap_V_G[k] += o.V_G;
        wrap_E_G[k] += o.E_G;
        wrap_B_G[k] += o.B_G;
        wrap_H_Gstar[k] += o.H_Gstar;
        wrap_V_Gstar[k] += o.V_Gstar;
        wrap_E_Gstar[k] += o.E_Gstar;
        wrap_B_Gstar[k] += o.B_Gstar;
    }

    void absorb(const Microcanonical& o) {
        for (std::size_t k = 0; k < configuration_count.size(); ++k) {
            configuration_count[k] += o.configuration_count[k];
            sum_clusters_G[k] += o.sum_clusters_G[k];
            sum_clusters_Gstar_complement[k] += o.sum_clusters_Gstar_complement[k];
            wrap_H_G[k] += o.wrap_H_G[k];
            wrap_V_G[k] += o.wrap_V_G[k];
            wrap_E_G[k] += o.wrap_E_G[k];
            wrap_B_G[k] += o.wrap_B_G[k];
            wrap_H_Gstar[k] += o.wrap_H_Gstar[k];
            wrap_V_Gstar[k] += o.wrap_V_Gstar[k];
            wrap_E_Gstar[k] += o.wrap_E_Gstar[k];
            wrap_B_Gstar[k] += o.wrap_B_Gstar[k];
        }
        euler_mismatches += o.euler_mismatches;
        if (o.first_euler_mask < first_euler_mask) {
            first_euler_mask = o.first_euler_mask;
        }
    }
};

inline bool euler_identity_holds(int L, uint64_t mask, const Observables& o, const Graph& g_nn) {
    const int euler = euler_vef0(L, mask, g_nn);
    const int lhs = o.clusters_G - o.clusters_Gstar - euler;
    const int h = static_cast<int>(o.H_G) - static_cast<int>(o.H_Gstar);
    const int v = static_cast<int>(o.V_G) - static_cast<int>(o.V_Gstar);
    const int e = static_cast<int>(o.E_G) - static_cast<int>(o.E_Gstar);
    const int b = static_cast<int>(o.B_G) - static_cast<int>(o.B_Gstar);
    return lhs == h && h == v && v == e && e == b;
}

}  // namespace matching
