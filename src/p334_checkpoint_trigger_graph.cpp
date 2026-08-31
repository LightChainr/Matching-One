// Export only the minimal-trigger-pair graph of one archived checkpoint.
#define P334_SNAPSHOT_LIBRARY
#include "p334_checkpoint_safe_triples.cpp"

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "usage: trigger_graph N h12 k0 seed counter OUTPUT.json\n";
        return 2;
    }
    try {
        const int n = std::stoi(argv[1]), h12 = std::stoi(argv[2]), k0 = std::stoi(argv[3]);
        const std::uint64_t seed = std::stoull(argv[4]), counter = std::stoull(argv[5]);
        const Geometry geometry = make_geometry({n, h12, 0, 1});
        std::vector<int> permutation;
        counter_permutation(n, seed, counter, permutation);
        SavedSnapshot snapshot(geometry);
        for (int i = 0; i < k0; ++i) snapshot.insert(permutation[i]);
        if (snapshot.rank != 1) throw std::logic_error("saved checkpoint is not rank one");
        std::vector<int> safe;
        std::vector<std::array<int, 2>> trigger_edges;
        for (int v = 0; v < n; ++v) if (!snapshot.active[v]) {
            SavedSnapshot one = snapshot;
            one.insert(v);
            if (one.rank == 1) safe.push_back(v);
        }
        for (std::size_t i = 0; i < safe.size(); ++i) {
            SavedSnapshot one = snapshot;
            one.insert(safe[i]);
            for (std::size_t j = i+1; j < safe.size(); ++j) {
                SavedSnapshot two = one;
                two.insert(safe[j]);
                if (two.rank == 2) trigger_edges.push_back({safe[i], safe[j]});
            }
        }
        std::ofstream output(argv[6]);
        if (!output) throw std::runtime_error("cannot open graph output");
        output << "{\"N\":" << n << ",\"h12\":" << h12 << ",\"k0\":" << k0
               << ",\"seed\":" << seed << ",\"replica_counter\":" << counter
               << ",\"ell\":[" << snapshot.line.x << ',' << snapshot.line.y
               << "],\"safe_sites\":[";
        for (std::size_t i = 0; i < safe.size(); ++i) {
            if (i) output << ','; output << safe[i];
        }
        output << "],\"minimal_trigger_pairs\":[";
        for (std::size_t i = 0; i < trigger_edges.size(); ++i) {
            if (i) output << ',';
            output << '[' << trigger_edges[i][0] << ',' << trigger_edges[i][1] << ']';
        }
        output << "]}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 2;
    }
    return 0;
}
