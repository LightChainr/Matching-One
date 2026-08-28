#include "torus_connectivity.hpp"

#include <cstdint>
#include <iostream>
#include <string>

using matching::DisplacementDSU;
using matching::Graph;
using matching::Observables;
using matching::analyze_config;
using matching::bit_set;
using matching::euler_identity_holds;
using matching::euler_vef0;
using matching::make_graph;
using matching::popcount_mask;

static int failures = 0;

#define CHECK(cond, msg)                                                                       \
    do {                                                                                       \
        if (!(cond)) {                                                                         \
            std::cerr << "FAIL " << msg << "\n";                                               \
            ++failures;                                                                        \
        }                                                                                      \
    } while (0)

static void test_graph_degrees() {
    for (int L = 2; L <= 5; ++L) {
        const Graph g = make_graph(L, false);
        const Graph gs = make_graph(L, true);
        CHECK(g.n == L * L, "G n");
        CHECK(gs.n == L * L, "G* n");
        for (int i = 0; i < g.n; ++i) {
            CHECK(g.adj[static_cast<std::size_t>(i)].size() == 2u, "G degree 2");
            CHECK(gs.adj[static_cast<std::size_t>(i)].size() == 4u, "G* degree 4");
            for (const auto& e : g.adj[static_cast<std::size_t>(i)]) {
                CHECK(e.dx != 0 || e.dy != 0, "nonzero displacement");
            }
        }
        // Original displacements, not folded: +1 never stored as 1-L.
        bool saw_plus_x = false;
        bool saw_wrap_plus_x = false;
        for (int y = 0; y < L; ++y) {
            const int i = y * L + (L - 1);
            for (const auto& e : g.adj[static_cast<std::size_t>(i)]) {
                if (e.dy == 0 && e.dx == 1 && e.to == y * L + 0) {
                    saw_wrap_plus_x = true;
                }
            }
            const int j = y * L + 0;
            for (const auto& e : g.adj[static_cast<std::size_t>(j)]) {
                if (e.dy == 0 && e.dx == 1 && e.to == y * L + 1) {
                    saw_plus_x = true;
                }
            }
        }
        CHECK(saw_plus_x, "interior +x displacement preserved");
        CHECK(saw_wrap_plus_x, "periodic +x displacement preserved (not folded to 1-L)");
    }
}

static void test_hand_configs_L2() {
    const int L = 2;
    const Graph g = make_graph(L, false);
    const Graph gs = make_graph(L, true);
    DisplacementDSU dsu(L * L);

    // Empty: 0 G clusters, 1 G* cluster that wraps both ways on the 2x2 matching torus.
    {
        const Observables o = analyze_config(0, g, gs, dsu);
        CHECK(o.k == 0, "empty k");
        CHECK(o.clusters_G == 0, "empty G clusters");
        CHECK(o.clusters_Gstar == 1, "empty G* one cluster");
        CHECK(o.H_G == 0 && o.V_G == 0, "empty G no wrap");
        CHECK(o.H_Gstar == 1 && o.V_Gstar == 1, "full matching wraps H and V");
        CHECK(euler_identity_holds(L, 0, o, g), "empty Euler");
    }

    // Full occupation: 1 G cluster wrapping both, 0 G*.
    {
        const uint64_t mask = 0xF;
        const Observables o = analyze_config(mask, g, gs, dsu);
        CHECK(o.k == 4, "full k");
        CHECK(o.clusters_G == 1, "full G clusters");
        CHECK(o.clusters_Gstar == 0, "full G* clusters");
        CHECK(o.H_G == 1 && o.V_G == 1 && o.B_G == 1, "full G wraps both");
        CHECK(o.E_Gstar == 0, "full G* no wrap");
        CHECK(euler_identity_holds(L, mask, o, g), "full Euler");
    }

    // Single occupied site, index 0 = (0,0).
    {
        const uint64_t mask = 1ull;
        const Observables o = analyze_config(mask, g, gs, dsu);
        CHECK(o.k == 1, "single k");
        CHECK(o.clusters_G == 1, "single G cluster");
        CHECK(o.H_G == 0 && o.V_G == 0, "single site cannot wrap");
        CHECK(euler_identity_holds(L, mask, o, g), "single Euler");
    }

    // Two adjacent in x: sites 0=(0,0) and 1=(1,0). On L=2 this is a wrapping pair.
    {
        const uint64_t mask = 0x3;
        const Observables o = analyze_config(mask, g, gs, dsu);
        CHECK(o.k == 2, "pair k");
        CHECK(o.clusters_G == 1, "pair one G cluster");
        CHECK(o.H_G == 1, "L=2 adjacent row wraps horizontally");
        CHECK(o.V_G == 0, "pair does not wrap vertically");
        CHECK(euler_identity_holds(L, mask, o, g), "pair Euler");
    }
}

static void test_exhaustive_euler(int L) {
    const Graph g = make_graph(L, false);
    const Graph gs = make_graph(L, true);
    DisplacementDSU dsu(L * L);
    const int n = L * L;
    const uint64_t total = 1ull << n;
    uint64_t bad = 0;
    uint64_t first = ~0ull;
    for (uint64_t mask = 0; mask < total; ++mask) {
        const Observables o = analyze_config(mask, g, gs, dsu);
        CHECK(o.k == popcount_mask(mask), "k popcount");
        CHECK(o.E_G == (o.H_G | o.V_G), "E derived");
        CHECK(o.B_G == (o.H_G & o.V_G), "B derived");
        CHECK(o.E_Gstar == (o.H_Gstar | o.V_Gstar), "E* derived");
        CHECK(o.B_Gstar == (o.H_Gstar & o.V_Gstar), "B* derived");
        if (!euler_identity_holds(L, mask, o, g)) {
            ++bad;
            if (mask < first) {
                first = mask;
            }
        }
    }
    CHECK(bad == 0, (std::string("euler identity L=") + std::to_string(L) +
                     " failures=" + std::to_string(bad) + " first=" + std::to_string(first))
                        .c_str());
}

static void test_l3_full_row() {
    const int L = 3;
    const Graph g = make_graph(L, false);
    const Graph gs = make_graph(L, true);
    DisplacementDSU dsu(L * L);
    // Bottom row y=0: sites 0,1,2.
    const uint64_t mask = 0x7;
    const Observables o = analyze_config(mask, g, gs, dsu);
    CHECK(o.k == 3, "row k");
    CHECK(o.clusters_G == 1, "row one cluster");
    CHECK(o.H_G == 1, "full row wraps H");
    CHECK(o.V_G == 0, "full row no V");
    CHECK(euler_identity_holds(L, mask, o, g), "row Euler");
}

int main() {
    test_graph_degrees();
    test_hand_configs_L2();
    test_l3_full_row();
    test_exhaustive_euler(2);
    test_exhaustive_euler(3);
    if (failures == 0) {
        std::cout << "test_torus_connectivity: PASS\n";
        return 0;
    }
    std::cout << "test_torus_connectivity: FAIL count=" << failures << "\n";
    return 1;
}
