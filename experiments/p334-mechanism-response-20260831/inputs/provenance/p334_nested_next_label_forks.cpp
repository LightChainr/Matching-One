// New conditional suffix sampling on the original common ordered prefixes.
// The rank engine and original-prefix RNG are used unchanged.
#define main p334_included_threshold_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main
#include <zlib.h>

namespace {
constexpr int quartets_per_prefix = 8;
constexpr std::uint64_t fork_seed = 202608311638334ULL;

struct FrozenPrefix {
    std::uint64_t counter;
    int batch;
    int k0;
    int first_rank;
    int second_rank;
};

std::vector<FrozenPrefix> read_prefixes(int n) {
    std::ifstream input("results/p334-full-birth-archive/N"+std::to_string(n)+".csv");
    if (!input) throw std::runtime_error("missing original complete birth archive");
    std::string line;
    std::getline(input, line);
    if (line != "N,batch,counter,k0,first_k1,first_k2,first_rank,second_k1,second_k2,second_rank")
        throw std::runtime_error("unexpected complete birth schema");
    std::vector<FrozenPrefix> rows;
    while (std::getline(input, line)) {
        std::stringstream cells(line);
        std::vector<std::uint64_t> values;
        std::string value;
        while (std::getline(cells, value, ',')) values.push_back(std::stoull(value));
        if (values.size() != 10 || values[0] != static_cast<std::uint64_t>(n))
            throw std::runtime_error("bad prefix source row");
        rows.push_back({values[2], static_cast<int>(values[1]), static_cast<int>(values[3]),
                        static_cast<int>(values[6]), static_cast<int>(values[9])});
    }
    if (rows.size() != 20000) throw std::runtime_error("source must contain all20000 prefixes");
    return rows;
}

std::uint64_t new_stream_key(int n, int prefix_index, int quartet, int group, int slot) {
    const std::uint64_t local_id =
        (static_cast<std::uint64_t>(prefix_index)*quartets_per_prefix+quartet)*8+group*4+slot;
    return splitmix64(fork_seed ^ splitmix64((static_cast<std::uint64_t>(n)<<32)|local_id));
}

void make_new_tail(const std::vector<int>& old_permutation, const std::vector<int>& vacant,
                   int k0, int selected, std::uint64_t stream_key, std::vector<int>& output) {
    output.assign(old_permutation.begin(), old_permutation.begin()+k0);
    output.push_back(selected);
    for (const int label: vacant) if (label != selected) output.push_back(label);
    SplitMixStream rng(stream_key);
    // Only the remaining suffix is shuffled; the next insertion and old prefix stay fixed.
    for (int stop = static_cast<int>(output.size())-1; stop > k0+1; --stop) {
        const int other = k0+1+static_cast<int>(rng.below(static_cast<std::uint64_t>(stop-k0)));
        std::swap(output[stop], output[other]);
    }
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::invalid_argument("usage: p334_nested_next_label_forks N output_dir code_commit");
        const int n = std::stoi(argv[1]);
        if (n != 325 && n != 425) throw std::invalid_argument("only original archive sizes");
        const std::filesystem::path out(argv[2]);
        if (std::filesystem::exists(out)) throw std::runtime_error("refuse to overwrite existing fork directory");
        const auto rows = read_prefixes(n);
        std::filesystem::create_directories(out);
        const std::uint64_t original_seed = n == 325 ? 20260831430325ULL : 20260831430425ULL;
        const int k0 = n == 325 ? 193 : 252;
        const Geometry first_geometry = make_geometry({n, n == 325 ? 57 : 132, 0, 1});
        const Geometry second_geometry = make_geometry({n, n == 325 ? 18 : 268, 0, 1});
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        std::vector<int> old_permutation, vacant, next_permutation;
        const auto start = std::chrono::steady_clock::now();
        std::uint64_t generated = 0, same_next_labels = 0;
        for (int batch = 0; batch < 20; ++batch) {
            std::ostringstream filename;
            filename << "N" << n << ".batch" << std::setw(2) << std::setfill('0') << batch << ".csv.gz";
            const auto path = out/filename.str();
            gzFile output = gzopen(path.c_str(), "wb1");
            if (!output) throw std::runtime_error("cannot open compressed batch output");
            gzprintf(output, "N,batch,counter,k0,first_rank,second_rank,quartet,group,replica,next_label,first_next_rank,second_next_rank,first_k1,first_k2,second_k1,second_k2\n");
            for (int offset = 0; offset < 1000; ++offset) {
                const int prefix_index = batch*1000+offset;
                const auto& prefix = rows[prefix_index];
                if (prefix.batch != batch || prefix.k0 != k0) throw std::runtime_error("original prefix order mismatch");
                counter_permutation(n, original_seed, prefix.counter, old_permutation);
                vacant.assign(old_permutation.begin()+k0, old_permutation.end());
                std::sort(vacant.begin(), vacant.end());
                for (int q = 0; q < quartets_per_prefix; ++q) {
                    int labels[2];
                    for (int group = 0; group < 2; ++group) {
                        SplitMixStream choose(new_stream_key(n, prefix_index, q, group, 0));
                        labels[group] = vacant[choose.below(vacant.size())];
                        for (int replica = 0; replica < 2; ++replica) {
                            make_new_tail(old_permutation, vacant, k0, labels[group],
                                          new_stream_key(n, prefix_index, q, group, replica+1), next_permutation);
                            const auto first = first_engine.ranks(next_permutation);
                            const auto second = second_engine.ranks(next_permutation);
                            const auto rank_after = [k0](const std::pair<int,int>& births) {
                                return static_cast<int>(births.first <= k0+1)+static_cast<int>(births.second <= k0+1);
                            };
                            if (gzprintf(output, "%d,%d,%llu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                                n,batch,static_cast<unsigned long long>(prefix.counter),k0,prefix.first_rank,prefix.second_rank,
                                q,group,replica,labels[group],rank_after(first),rank_after(second),
                                first.first,first.second,second.first,second.second) <= 0)
                                throw std::runtime_error("failed writing compressed fork row");
                            ++generated;
                        }
                    }
                    same_next_labels += labels[0] == labels[1];
                }
            }
            if (gzclose(output) != Z_OK) throw std::runtime_error("failed closing compressed batch");
            std::cout << "N" << n << " batch " << batch << " complete; new tails " << generated << std::endl;
        }
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream metadata(out/"metadata.json");
        metadata << std::setprecision(15)
                 << "{\n  \"code_commit\":\"" << argv[3] << "\",\n"
                 << "  \"prefix_source_commit\":\"9c495ab13e65f2bc93dc0849ee3b73f88724c4b1\",\n"
                 << "  \"N\":" << n << ",\n  \"prefix_count\":20000,\n  \"batch_count\":20,\n"
                 << "  \"quartets_per_prefix\":8,\n  \"new_prefix_samples\":0,\n"
                 << "  \"new_tail_paths\":" << generated << ",\n  \"fork_seed\":" << fork_seed << ",\n"
                 << "  \"original_seed\":" << original_seed << ",\n  \"k0\":" << k0 << ",\n"
                 << "  \"equal_next_label_quartets\":" << same_next_labels << ",\n"
                 << "  \"elapsed_seconds\":" << elapsed << ",\n"
                 << "  \"compiler\":\"" << __VERSION__ << "\",\n"
                 << "  \"host\":\"local Mac CPU\",\n"
                 << "  \"semantics\":\"New conditional suffix MC on old ordered paired prefixes; iid U,V allow equality; four independent suffixes conditional on prefix/U/V; no DP or new prefix population\"\n}\n";
        std::cout << "N" << n << " COMPLETED " << generated << " new tails in " << elapsed << "s" << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
