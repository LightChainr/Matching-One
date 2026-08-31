// Two saved N425 checkpoints only: exact pair-clique versus safe-triple census.
// Reuse the frozen production geometry, uint64 RNG, and potential union-find.
// No production runner code, random domain, or Monte Carlo stream is changed.
#define main p334_included_production_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main

struct SavedSnapshot {
    const Geometry& geometry;
    std::vector<std::uint8_t> active;
    HomologyUnionFind uf;
    int rank = 0;
    Vector line{0, 0};

    explicit SavedSnapshot(const Geometry& g)
        : geometry(g), active(g.n, 0), uf(g.quotient) {}

    void insert(int site) {
        if (active[site]) throw std::logic_error("duplicate saved insertion");
        active[site] = 1;
        for (int edge_index : geometry.primal_incident[site]) {
            const Edge& edge = geometry.primal_edges[edge_index];
            if (active[edge.i] && active[edge.j]) uf.add_edge(edge);
        }
        const auto mark = uf.component_mark(site);
        if (mark.rank == 2 || (rank == 1 && mark.rank == 1 &&
                              !same_vector(line, mark.line))) {
            rank = 2;
        } else if (rank == 0 && mark.rank == 1) {
            rank = 1;
            line = mark.line;
        }
    }
};

void write_checkpoint(std::ostream& output, std::uint64_t counter,
                      std::uint64_t expected_squares) {
    constexpr int n = 425, k0 = 252;
    constexpr std::uint64_t seed = 20260831430425ULL;
    const Geometry geometry = make_geometry({425, 268, 0, 1});
    std::vector<int> permutation;
    counter_permutation(n, seed, counter, permutation);
    SavedSnapshot checkpoint(geometry);
    int k1 = 0;
    for (int k = 0; k < k0; ++k) {
        checkpoint.insert(permutation[k]);
        if (k1 == 0 && checkpoint.rank == 1) k1 = k+1;
    }
    if (checkpoint.rank != 1 || !same_vector(checkpoint.line, {12, -19}) ||
        k0-k1 != 10) throw std::logic_error("saved checkpoint signature mismatch");
    std::vector<int> vacant;
    for (int v = 0; v < n; ++v) if (!checkpoint.active[v]) vacant.push_back(v);
    const int d = static_cast<int>(vacant.size());
    std::vector<std::uint8_t> singleton(n, 0);
    std::vector<std::vector<std::uint8_t>> pair(n, std::vector<std::uint8_t>(n, 0));
    std::vector<std::uint64_t> degree(n, 0);
    std::uint64_t b1 = 0, b2 = 0, squares = 0;
    for (int v : vacant) {
        SavedSnapshot child = checkpoint;
        child.insert(v);
        singleton[v] = child.rank == 1;
        b1 += singleton[v];
        if (!singleton[v]) continue;
        for (int w : vacant) {
            if (w <= v) continue;
            SavedSnapshot two = child;
            two.insert(w);
            if (two.rank == 1) {
                pair[v][w] = pair[w][v] = 1;
                ++degree[v]; ++degree[w]; ++b2;
            }
        }
    }
    for (int v : vacant) squares += degree[v]*degree[v];
    if (b1 != 173 || b2 != 14770 || squares != expected_squares) {
        throw std::logic_error("snapshot replay differs from archived exact pair counts");
    }

    std::uint64_t triangles = 0, safe_triples = 0;
    std::vector<std::array<int, 3>> first_nonfaces, all_nonfaces;
    for (int i = 0; i < d; ++i) {
        const int v = vacant[i];
        if (!singleton[v]) continue;
        SavedSnapshot one = checkpoint;
        one.insert(v);
        for (int j = i+1; j < d; ++j) {
            const int w = vacant[j];
            if (!pair[v][w]) continue;
            SavedSnapshot two = one;
            two.insert(w);
            for (int k = j+1; k < d; ++k) {
                const int z = vacant[k];
                if (!pair[v][z] || !pair[w][z]) continue;
                ++triangles;
                SavedSnapshot three = two;
                three.insert(z);
                if (three.rank == 1) ++safe_triples;
                else {
                    all_nonfaces.push_back({v, w, z});
                    if (first_nonfaces.size() < 8) first_nonfaces.push_back({v, w, z});
                }
            }
        }
    }

    output << "{\"seed\":" << seed << ",\"replica_counter\":" << counter
           << ",\"N\":425,\"orientation\":\"second\",\"k0\":252,\"age\":10,"
              "\"period_matrix\":[[425,268],[0,1]],\"ell\":[12,-19],\"d\":" << d
           << ",\"H2\":" << d-b1 << ",\"b1_safe\":" << b1
           << ",\"b2_safe_pairs\":" << b2 << ",\"sum_degree_squared\":" << squares
           << ",\"vacant_triples\":" << static_cast<std::uint64_t>(d)*(d-1)*(d-2)/6
           << ",\"safe_graph_triangles\":" << triangles
           << ",\"actual_safe_triples\":" << safe_triples
           << ",\"minimal_nonfaces_size3\":" << triangles-safe_triples
           << ",\"occupied_prefix_labels\":[";
    for (int i = 0; i < k0; ++i) { if (i) output << ','; output << permutation[i]; }
    output << "],\"first_minimal_nonfaces\":[";
    for (std::size_t i = 0; i < first_nonfaces.size(); ++i) {
        if (i) output << ',';
        const auto triple = first_nonfaces[i];
        output << '[' << triple[0] << ',' << triple[1] << ',' << triple[2] << ']';
    }
    output << "],\"all_minimal_nonfaces\":[";
    for (std::size_t i = 0; i < all_nonfaces.size(); ++i) {
        if (i) output << ',';
        const auto triple = all_nonfaces[i];
        output << '[' << triple[0] << ',' << triple[1] << ',' << triple[2] << ']';
    }
    output << "],\"first_witness_subset_ranks\":[";
    if (!first_nonfaces.empty()) {
        for (int bits = 0; bits < 8; ++bits) {
            SavedSnapshot subset = checkpoint;
            if (bits) output << ',';
            output << "{\"subset_mask\":" << bits << ",\"added_labels\":[";
            bool first = true;
            for (int j = 0; j < 3; ++j) if (bits & (1 << j)) {
                const int v = first_nonfaces[0][j];
                if (!first) output << ',';
                first = false;
                output << v;
                subset.insert(v);
            }
            output << "],\"ambient_rank\":" << subset.rank << '}';
            if (subset.rank != (bits == 7 ? 2 : 1)) {
                throw std::logic_error("minimal nonface subset-rank certificate failed");
            }
        }
    }
    output << "]}";
}

#ifndef P334_SNAPSHOT_LIBRARY
int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "usage: p334_checkpoint_safe_triples OUTPUT.json\n"; return 2; }
    try {
        std::ofstream output(argv[1]);
        if (!output) throw std::runtime_error("cannot open output");
        output << "{\"schema\":\"matching-one/p334-real-safe-triple-census/v2\","
                  "\"new_samples\":0,\"source_commit\":\"6147e22f53902a94e5f133739f2c1d423691d0b8\","
                  "\"checkpoints\":[";
        write_checkpoint(output, 43042514269ULL, 5045796);
        output << ',';
        write_checkpoint(output, 43042505280ULL, 5046876);
        output << "]}\n";
        std::cout << "completed exact two-checkpoint N425 safe-triple census\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 2;
    }
    return 0;
}
#endif
