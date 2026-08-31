// Additional full-observable readout of the existing e81dd59 counter domain.
// Reuse the original engine/RNG. No geometry-pilot or conditional DP is rerun.
#define main p334_included_threshold_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("usage: p334_full_birth_archive N output.csv");
        const int n = std::stoi(argv[1]);
        if (n != 325 && n != 425) throw std::invalid_argument("only the two original archive sizes");
        const int k0 = n == 325 ? 193 : 252;
        const std::uint64_t seed = n == 325 ? 20260831430325ULL : 20260831430425ULL;
        const std::uint64_t first_counter = n == 325 ? 43032500000ULL : 43042500000ULL;
        const int first_shear = n == 325 ? 57 : 132;
        const int second_shear = n == 325 ? 18 : 268;
        const std::filesystem::path path(argv[2]);
        if (std::filesystem::exists(path) || std::filesystem::exists(path.string()+".metadata.json"))
            throw std::invalid_argument("refuse to overwrite an existing readout");
        if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot open output");
        const Geometry first_geometry = make_geometry({n, first_shear, 0, 1});
        const Geometry second_geometry = make_geometry({n, second_shear, 0, 1});
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        std::vector<int> permutation;
        const auto start = std::chrono::steady_clock::now();
        output << "N,batch,counter,k0,first_k1,first_k2,first_rank,second_k1,second_k2,second_rank\n";
        std::array<std::uint64_t, 9> rank_pairs{};
        for (int replica = 0; replica < 20000; ++replica) {
            const std::uint64_t counter = first_counter + replica;
            counter_permutation(n, seed, counter, permutation);
            const auto first = first_engine.ranks(permutation);
            const auto second = second_engine.ranks(permutation);
            const auto rank_at_checkpoint = [k0](const std::pair<int,int>& births) {
                return static_cast<int>(births.first <= k0) + static_cast<int>(births.second <= k0);
            };
            const int rf = rank_at_checkpoint(first), rs = rank_at_checkpoint(second);
            ++rank_pairs[3*rf+rs];
            output << n << ',' << replica/1000 << ',' << counter << ',' << k0 << ','
                   << first.first << ',' << first.second << ',' << rf << ','
                   << second.first << ',' << second.second << ',' << rs << '\n';
        }
        output.close();
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream metadata(path.string()+".metadata.json");
        metadata << std::setprecision(15)
                 << "{\n  \"source_commit\":\"e81dd59ff6be69056e504e0e81cfeccf73dc5e97\",\n"
                 << "  \"engine\":\"unchanged ThresholdEngine::ranks and counter_permutation\",\n"
                 << "  \"N\":" << n << ",\n  \"k0\":" << k0
                 << ",\n  \"seed\":" << seed << ",\n  \"counter_first\":" << first_counter
                 << ",\n  \"counter_last_exclusive\":" << first_counter+20000
                 << ",\n  \"samples\":20000,\n  \"batches\":20,\n  \"new_random_samples\":0,\n"
                 << "  \"first_period_matrix\":[[" << n << ',' << first_shear << "],[0,1]],\n"
                 << "  \"second_period_matrix\":[[" << n << ',' << second_shear << "],[0,1]],\n"
                 << "  \"rank_pair_counts_row_major\":[";
        for (std::size_t i=0; i<rank_pairs.size(); ++i) { if(i) metadata << ','; metadata << rank_pairs[i]; }
        metadata << "],\n  \"elapsed_seconds\":" << elapsed
                 << ",\n  \"compiler\":\"" << __VERSION__ << "\",\n"
                 << "  \"semantics\":\"K1 first rank>=1, K2 first rank2; original reverse-matching Alexander implementation; same pair permutation and original batch IDs\"\n}\n";
        std::cout << "N" << n << " completed 20000 counters in " << elapsed << " seconds; rank pairs";
        for (const auto count:rank_pairs) std::cout << ' ' << count;
        std::cout << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
