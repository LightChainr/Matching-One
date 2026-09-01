// Fixed-budget independent suffix extension at existing paired R0/R0 prefixes.
// The unchanged original rank engine and exact contact code are included here.
#define P334_CONTACT_LIBRARY_ONLY
#include "../census/provenance/p334_next_label_contact_coordinates.cpp"

namespace {
constexpr std::uint64_t fork_seed = 202608311638334ULL;
constexpr int extra_quartets = 64;
constexpr std::uint64_t extension_bit = 1ULL << 31;

std::uint64_t extension_stream(int n, int prefix_index, int q, int group, int slot) {
    if (q < 8 || q >= 72 || prefix_index < 0 || prefix_index >= 20000 || group < 0 || group > 1 || slot < 0 || slot > 2)
        throw std::runtime_error("stream coordinate out of fixed extension domain");
    const std::uint64_t low = ((static_cast<std::uint64_t>(prefix_index)*64+(q-8))*8+group*4+slot);
    if (low >= extension_bit) throw std::runtime_error("extension stream domain overflow");
    const std::uint64_t local_id = extension_bit | low;
    return splitmix64(fork_seed ^ splitmix64((static_cast<std::uint64_t>(n)<<32)|local_id));
}

void extension_tail(const std::vector<int>& old, const std::vector<int>& vacant, int k0, int label,
                    std::uint64_t key, std::vector<int>& permutation) {
    permutation.assign(old.begin(), old.begin()+k0);
    permutation.push_back(label);
    for (int v: vacant) if (v != label) permutation.push_back(v);
    SplitMixStream rng(key);
    for (int stop = static_cast<int>(permutation.size())-1; stop > k0+1; --stop) {
        int other = k0+1+static_cast<int>(rng.below(static_cast<std::uint64_t>(stop-k0)));
        std::swap(permutation[stop], permutation[other]);
    }
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: cell00_extension N batch output_dir");
        const int n = std::stoi(argv[1]), batch = std::stoi(argv[2]);
        if ((n != 325 && n != 425) || batch < 0 || batch >= 20) throw std::runtime_error("invalid fixed batch");
        std::filesystem::path out(argv[3]);
        std::filesystem::create_directories(out);
        std::ostringstream name;
        name << "N" << n << ".batch" << std::setw(2) << std::setfill('0') << batch;
        const auto target = out/(name.str()+".csv.gz");
        if (std::filesystem::exists(target)) throw std::runtime_error("refuse to overwrite generated tails");
        std::ifstream source("inputs/prefix_archive/N"+std::to_string(n)+".csv");
        std::string line;
        std::getline(source, line);
        if (line != "N,batch,counter,k0,first_k1,first_k2,first_rank,second_k1,second_k2,second_rank")
            throw std::runtime_error("wrong immutable prefix archive schema");
        const auto first_geometry = make_geometry({n,n==325?57:132,0,1});
        const auto second_geometry = make_geometry({n,n==325?18:268,0,1});
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        Checkpoint first_checkpoint(first_geometry), second_checkpoint(second_geometry);
        const std::uint64_t old_seed = n == 325 ? 20260831430325ULL : 20260831430425ULL;
        gzFile output = gzopen(target.c_str(), "wb1");
        if (!output) throw std::runtime_error("gzip open failed");
        gzprintf(output,"N,batch,counter,k0,first_rank,second_rank,quartet,group,replica,next_label,first_next_rank,second_next_rank,first_k1,first_k2,second_k1,second_k2,first_e,first_c,second_e,second_c\n");
        std::vector<int> old, vacant, permutation;
        std::uint64_t prefix_count = 0, tails = 0, equal_labels = 0;
        const auto start = std::chrono::steady_clock::now();
        int prefix_index = -1;
        while (std::getline(source, line)) {
            ++prefix_index;
            const auto v = parse_csv(line);
            if (v.size() != 10 || v[0] != static_cast<std::uint64_t>(n) || v[1] != static_cast<std::uint64_t>(prefix_index/1000))
                throw std::runtime_error("original prefix ordering changed");
            if (v[1] != static_cast<std::uint64_t>(batch) || v[6] != 0 || v[9] != 0) continue;
            const int k0 = static_cast<int>(v[3]);
            const std::uint64_t counter = v[2];
            counter_permutation(n, old_seed, counter, old);
            first_checkpoint.load(old, k0); second_checkpoint.load(old, k0);
            for (int k = 0; k < k0; ++k) {
                if (first_checkpoint.uf.component_mark(old[k]).rank != 0 || second_checkpoint.uf.component_mark(old[k]).rank != 0)
                    throw std::runtime_error("original prefix reconstruction is not cell00");
            }
            vacant.assign(old.begin()+k0, old.end());
            std::sort(vacant.begin(), vacant.end());
            std::map<int,std::array<ContactMark,2>> cache;
            for (int q = 8; q < 72; ++q) {
                int labels[2];
                for (int group = 0; group < 2; ++group) {
                    SplitMixStream choose(extension_stream(n,prefix_index,q,group,0));
                    labels[group] = vacant[choose.below(vacant.size())];
                    for (int replica = 0; replica < 2; ++replica) {
                        extension_tail(old,vacant,k0,labels[group],extension_stream(n,prefix_index,q,group,replica+1),permutation);
                        const auto first = first_engine.ranks(permutation), second = second_engine.ranks(permutation);
                        const int nr0 = (first.first <= k0+1)+(first.second <= k0+1);
                        const int nr1 = (second.first <= k0+1)+(second.second <= k0+1);
                        auto it = cache.find(labels[group]);
                        if (it == cache.end()) it = cache.emplace(labels[group],std::array<ContactMark,2>{first_checkpoint.read(labels[group],0,nr0),second_checkpoint.read(labels[group],0,nr1)}).first;
                        const auto& a = it->second[0]; const auto& b = it->second[1];
                        if (a.contact_rank != nr0 || b.contact_rank != nr1) throw std::runtime_error("replica immediate-rank disagreement");
                        if (gzprintf(output,"%d,%d,%llu,%d,0,0,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
                            n,batch,static_cast<unsigned long long>(counter),k0,q,group,replica,labels[group],nr0,nr1,
                            first.first,first.second,second.first,second.second,a.e,a.c,b.e,b.c) <= 0) throw std::runtime_error("gzip write failed");
                        ++tails;
                    }
                }
                equal_labels += labels[0] == labels[1];
            }
            ++prefix_count;
        }
        if (prefix_index != 19999 || tails != prefix_count*64*4) throw std::runtime_error("incomplete fixed batch");
        if (gzclose(output) != Z_OK) throw std::runtime_error("gzip close failed");
        double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream metadata(out/(name.str()+".metadata.json"));
        metadata << std::setprecision(15) << "{\"N\":" << n << ",\"batch\":" << batch << ",\"prefixes\":" << prefix_count
                 << ",\"new_quartets_per_prefix\":64,\"quartet_begin\":8,\"quartet_end_inclusive\":71,\"new_tail_paths\":" << tails
                 << ",\"new_prefix_samples\":0,\"original_seed\":" << old_seed << ",\"fork_seed\":" << fork_seed
                 << ",\"new_stream_domain_bit\":31,\"equal_UV_labels\":" << equal_labels << ",\"elapsed_seconds\":" << seconds
                 << ",\"compiler\":\"" << __VERSION__ << "\"}\n";
        std::cout << "N" << n << " batch" << batch << " prefixes=" << prefix_count << " new_tails=" << tails << " seconds=" << seconds << std::endl;
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << std::endl; return 1; }
}
