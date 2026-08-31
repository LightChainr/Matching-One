// Two fixed real N425 snapshots only. Filter frozen minimal pairs/triples,
// then query the existing rank oracle for every remaining four-set.
#define P334_SNAPSHOT_LIBRARY
#include "p334_checkpoint_safe_triples.cpp"

void quartic_checkpoint(std::istream& input, std::ostream& output) {
    constexpr int n = 425, k0 = 252, d = 173;
    constexpr std::uint64_t seed = 20260831430425ULL;
    std::uint64_t counter, expected_four;
    int pair_count, triple_count;
    if (!(input >> counter >> pair_count >> triple_count >> expected_four))
        throw std::runtime_error("missing bounded constraint header");
    const Geometry geometry = make_geometry({425, 268, 0, 1});
    std::vector<int> permutation;
    counter_permutation(n, seed, counter, permutation);
    SavedSnapshot checkpoint(geometry);
    for (int i = 0; i < k0; ++i) checkpoint.insert(permutation[i]);
    if (checkpoint.rank != 1 || !same_vector(checkpoint.line, {12, -19}))
        throw std::logic_error("saved rank-one state mismatch");
    std::vector<int> vacant, index(n, -1);
    for (int v = 0; v < n; ++v) if (!checkpoint.active[v]) {
        index[v] = vacant.size(); vacant.push_back(v);
    }
    if (vacant.size() != d) throw std::logic_error("wrong vacant count");
    std::vector<std::uint8_t> pair(d*d, 0), triple(d*d*d, 0);
    for (int e = 0; e < pair_count; ++e) {
        int u, v; input >> u >> v;
        const int i = index.at(u), j = index.at(v);
        if (i < 0 || j < 0 || i >= j) throw std::logic_error("bad pair labels");
        pair[i*d+j] = 1;
    }
    for (int e = 0; e < triple_count; ++e) {
        int u, v, w; input >> u >> v >> w;
        const int i = index.at(u), j = index.at(v), k = index.at(w);
        if (i < 0 || j < 0 || k < 0 || i >= j || j >= k)
            throw std::logic_error("bad triple labels");
        triple[(i*d+j)*d+k] = 1;
    }
    std::uint64_t candidates = 0, actual_safe = 0;
    std::vector<std::array<int, 4>> quartics;
    const auto started = std::chrono::steady_clock::now();
    for (int i = 0; i < d; ++i) {
        SavedSnapshot one = checkpoint;
        one.insert(vacant[i]);
        for (int j = i+1; j < d; ++j) {
            if (pair[i*d+j]) continue;
            SavedSnapshot two = one;
            two.insert(vacant[j]);
            for (int k = j+1; k < d; ++k) {
                if (pair[i*d+k] || pair[j*d+k] || triple[(i*d+j)*d+k]) continue;
                SavedSnapshot three = two;
                three.insert(vacant[k]);
                for (int l = k+1; l < d; ++l) {
                    if (pair[i*d+l] || pair[j*d+l] || pair[k*d+l] ||
                        triple[(i*d+j)*d+l] || triple[(i*d+k)*d+l] ||
                        triple[(j*d+k)*d+l]) continue;
                    ++candidates;
                    SavedSnapshot four = three;
                    four.insert(vacant[l]);
                    if (four.rank == 1) ++actual_safe;
                    else quartics.push_back({vacant[i], vacant[j], vacant[k], vacant[l]});
                }
            }
        }
        if (i % 40 == 0) std::cerr << counter << " prefix_index=" << i
                                  << " candidates=" << candidates
                                  << " quartics=" << quartics.size() << '\n';
    }
    if (candidates != expected_four) throw std::logic_error("frozen truncated count mismatch");
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    output << "{\"replica_counter\":" << counter << ",\"seed\":" << seed
           << ",\"N\":425,\"orientation\":\"second\",\"k0\":252,\"d\":173,"
              "\"period_matrix\":[[425,268],[0,1]],\"ell\":[12,-19],"
              "\"pair_triple_safe_four_sets\":" << candidates
           << ",\"actual_safe_four_sets\":" << actual_safe
           << ",\"minimal_nonfaces_size4\":" << quartics.size()
           << ",\"single_thread_seconds\":" << std::setprecision(10) << seconds
           << ",\"occupied_prefix_labels\":[";
    for (int i = 0; i < k0; ++i) { if (i) output << ','; output << permutation[i]; }
    output << "],\"all_minimal_quartics\":[";
    for (std::size_t i = 0; i < quartics.size(); ++i) {
        if (i) output << ',';
        const auto e = quartics[i];
        output << '[' << e[0] << ',' << e[1] << ',' << e[2] << ',' << e[3] << ']';
    }
    output << "],\"first_quartic_subset_ranks\":[";
    if (!quartics.empty()) for (int bits = 0; bits < 16; ++bits) {
        SavedSnapshot subset = checkpoint;
        for (int j = 0; j < 4; ++j) if (bits & (1 << j)) subset.insert(quartics[0][j]);
        if (bits) output << ',';
        output << subset.rank;
    }
    output << "]}";
    std::cerr << counter << " DONE candidates=" << candidates << " quartics=" << quartics.size()
              << " seconds=" << seconds << '\n';
}

int main(int argc, char** argv) {
    if (argc != 3) { std::cerr << "usage: p334_safe_quartics CONSTRAINTS.txt OUTPUT.json\n"; return 2; }
    try {
        std::ifstream input(argv[1]);
        std::ofstream output(argv[2]);
        if (!input || !output) throw std::runtime_error("cannot open bounded input/output");
        output << "{\"schema\":\"matching-one/p334-real-safe-quartic-census/v1\","
                  "\"source_commit\":\"d5d2cc89e77ebb2ec6252df75dc858e9c240e6ce\","
                  "\"new_samples\":0,\"threads\":1,\"checkpoints\":[";
        quartic_checkpoint(input, output);
        output << ',';
        quartic_checkpoint(input, output);
        output << "]}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
    return 0;
}
