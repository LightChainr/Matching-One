// Fresh F4 fixed-face/ordinary, two-orientation paired Newman--Ziff histograms.
// Compile on Linux: g++ -O3 -std=c++17 -fopenmp producer.cpp -o f4_producer
// The original Geometry, UF and counter-permutation backend is unchanged.
#define main p337_f4_unused_backend_main
#include "../../src/threshold_rank_orientation_mc.cpp"
#undef main
#include <thread>

#ifndef P337_F4_BUILD_GIT_COMMIT
#define P337_F4_BUILD_GIT_COMMIT "not_embedded; freeze_commit identifies dispatched source"
#endif

namespace {
constexpr std::uint64_t f4_domain_tag = 0x5033333746345631ULL;
constexpr const char* f4_parent = "89f4383d376a53121a14aee725cd9da5d8167674";
constexpr const char* f4_backend_blob = "544ac429257ea9c229eeebb3ab329a62b28bf432";

struct F4Options {
    int n = 0, begin = -1, end = -1, threads = 0;
    std::uint64_t samples = 0, seed = 0;
    std::string prefix, freeze;
    bool has_seed = false;
};

F4Options f4_options(int argc, char** argv) {
    F4Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (++i == argc) throw std::invalid_argument("missing value for " + key);
        const std::string value = argv[i];
        if (key == "--n") o.n = std::stoi(value);
        else if (key == "--samples-per-batch") o.samples = std::stoull(value);
        else if (key == "--batch-begin") o.begin = std::stoi(value);
        else if (key == "--batch-end") o.end = std::stoi(value);
        else if (key == "--seed") { o.seed = std::stoull(value); o.has_seed = true; }
        else if (key == "--threads") o.threads = std::stoi(value);
        else if (key == "--output-prefix") o.prefix = value;
        else if (key == "--freeze-commit") o.freeze = value;
        else throw std::invalid_argument("unknown option " + key);
    }
    if ((o.n != 65 && o.n != 85 && o.n != 130 && o.n != 170) || !o.samples ||
        !o.has_seed || o.begin < 0 || o.begin >= o.end || o.end > 100 ||
        o.threads < 1 || o.prefix.empty() || o.freeze.size() != 40 ||
        o.freeze.find_first_not_of("0123456789abcdef") != std::string::npos)
        throw std::invalid_argument("required: --n {65,85,130,170} --samples-per-batch UINT --batch-begin B --batch-end E (0<=B<E<=100) --seed UINT64 --threads T --output-prefix NEW_PREFIX --freeze-commit FULL_SHA");
    if (o.samples > std::numeric_limits<std::uint64_t>::max()/static_cast<std::uint64_t>(o.end))
        throw std::invalid_argument("replica counter overflow");
    return o;
}

std::array<int,4> f4_face(const Geometry& g) {
    std::array<int,4> face{0, positive_mod(g.a,g.n), positive_mod(g.b,g.n), positive_mod(g.a+g.b,g.n)};
    auto sorted = face;
    std::sort(sorted.begin(),sorted.end());
    if (std::adjacent_find(sorted.begin(),sorted.end()) != sorted.end())
        throw std::logic_error("fixed face must contain four distinct sites");
    return face;
}

class F4Sweep {
    const Geometry& g_;
    HomologyUnionFind uf_;
    std::vector<std::uint8_t> active_;
    std::array<int,4> face_;
    int rank_ = 0;

    void activate(int vertex) {
        active_[vertex] = 1;
        for (const int ei: g_.primal_incident[vertex]) {
            const auto& edge = g_.primal_edges[ei];
            if (active_[edge.i] && active_[edge.j]) uf_.add_edge(edge);
        }
        // Disjoint essential components of an embedded torus have parallel lines;
        // the ambient rank is therefore the maximum component rank here.
        rank_ = std::max(rank_,static_cast<int>(uf_.component_rank(vertex)));
    }

  public:
    explicit F4Sweep(const Geometry& g): g_(g), uf_(g.n,g.a,g.b), active_(g.n), face_(f4_face(g)) {}

    std::pair<int,int> births(const std::vector<int>& permutation, bool forced) {
        uf_.reset(); std::fill(active_.begin(),active_.end(),0); rank_ = 0;
        if (forced) for (int vertex: face_) activate(vertex);
        // Inspect the whole forced face before assigning free-site count zero.
        int first = rank_ >= 1 ? 0 : -1, second = rank_ >= 2 ? 0 : -1;
        int k = 0;
        for (int vertex: permutation) {
            if (second >= 0) break;
            if (active_[vertex]) continue;
            ++k; activate(vertex);
            if (first < 0 && rank_ >= 1) first = k;
            if (rank_ >= 2) second = k;
        }
        const int degree = g_.n-(forced ? 4 : 0);
        if (first < 0 || second < first || second > degree || (!forced && first == 0))
            throw std::logic_error("invalid forward first/second birth");
        return {first,second};
    }
};

struct F4Histogram {
    int degree;
    std::array<std::vector<std::uint64_t>,2> count;
    explicit F4Histogram(int d): degree(d), count{std::vector<std::uint64_t>(d+1),std::vector<std::uint64_t>(d+1)} {}
    void add(std::pair<int,int> k) { ++count[0][k.first]; ++count[1][k.second]; }
};

struct F4Batch {
    // first ordinary, first forced, second ordinary, second forced.
    std::array<F4Histogram,4> mode;
    explicit F4Batch(int n): mode{F4Histogram(n),F4Histogram(n-4),F4Histogram(n),F4Histogram(n-4)} {}
};
}

int main(int argc, char** argv) {
    try {
        const auto o = f4_options(argc,argv);
        const std::filesystem::path csv_path(o.prefix+".hist.csv"), meta_path(o.prefix+".metadata.json");
        if (std::filesystem::exists(csv_path) || std::filesystem::exists(meta_path))
            throw std::runtime_error("refuse to overwrite F4 output");
        if (!csv_path.parent_path().empty()) std::filesystem::create_directories(csv_path.parent_path());
        const auto found = std::find_if(kDesigns.begin(),kDesigns.end(),[&](const PairDesign& d){return d.n == o.n;});
        if (found == kDesigns.end()) throw std::logic_error("missing frozen Gaussian pair");
        const auto& pair = *found;
        const Geometry g1 = make_geometry(pair.a1,pair.b1), g2 = make_geometry(pair.a2,pair.b2);
        const std::uint64_t n_seed = splitmix64(o.seed ^ splitmix64(static_cast<std::uint64_t>(o.n)^f4_domain_tag));
        const int batches = o.end-o.begin;
        std::vector<F4Batch> result;
        result.reserve(batches);
        for (int b = 0; b < batches; ++b) result.emplace_back(o.n);
        std::vector<std::string> errors(batches);
        int actual_threads = 1;
        const auto wall_start = std::chrono::steady_clock::now();
        const std::clock_t cpu_start = std::clock();
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(o.threads);
#pragma omp parallel
        {
#pragma omp single
            actual_threads = omp_get_num_threads();
#pragma omp for schedule(static)
#endif
            for (int local_batch = 0; local_batch < batches; ++local_batch) {
                try {
                    const int batch = o.begin+local_batch;
                    F4Sweep first(g1), second(g2);
                    std::vector<int> permutation;
                    auto& hist = result[local_batch].mode;
                    for (std::uint64_t local = 0; local < o.samples; ++local) {
                        const std::uint64_t replica = static_cast<std::uint64_t>(batch)*o.samples+local;
                        counter_permutation(o.n,n_seed,replica,permutation);
                        hist[0].add(first.births(permutation,false));
                        hist[1].add(first.births(permutation,true));
                        hist[2].add(second.births(permutation,false));
                        hist[3].add(second.births(permutation,true));
                    }
                } catch (const std::exception& error) { errors[local_batch] = error.what(); }
            }
#ifdef _OPENMP
        }
#endif
        for (int b = 0; b < batches; ++b)
            if (!errors[b].empty()) throw std::runtime_error("batch "+std::to_string(o.begin+b)+": "+errors[b]);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_start).count();
        const double cpu = static_cast<double>(std::clock()-cpu_start)/CLOCKS_PER_SEC;
        std::ofstream csv(csv_path);
        if (!csv) throw std::runtime_error("cannot open histogram output");
        csv << "n,batch,orientation,mode,degree,replicas,birth,k,count\n";
        std::uint64_t rows = 0;
        for (int b = 0; b < batches; ++b) for (int mode = 0; mode < 4; ++mode) {
            const auto& hist = result[b].mode[mode];
            for (int birth = 0; birth < 2; ++birth) {
                const auto total = std::accumulate(hist.count[birth].begin(),hist.count[birth].end(),std::uint64_t(0));
                if (total != o.samples) throw std::logic_error("histogram total differs from paired batch samples");
                for (int k = 0; k <= hist.degree; ++k) {
                    csv << o.n << ',' << o.begin+b << ',' << (mode < 2 ? "first" : "second") << ','
                        << (mode%2 ? "forced" : "ordinary") << ',' << hist.degree << ',' << o.samples << ','
                        << (birth ? "second" : "first") << ',' << k << ',' << hist.count[birth][k] << '\n';
                    ++rows;
                }
            }
        }
        csv.close();
        if (!csv) throw std::runtime_error("histogram write failed");
        const auto face1 = f4_face(g1), face2 = f4_face(g2);
        std::ofstream meta(meta_path);
        meta << std::setprecision(17)
             << "{\n\"status\":\"completed\",\n\"N\":" << o.n << ",\n\"samples_per_batch\":" << o.samples
             << ",\n\"batch_begin\":" << o.begin << ",\n\"batch_end_exclusive\":" << o.end
             << ",\n\"total_batch_domain\":100,\n\"paired_permutations\":" << o.samples*batches
             << ",\n\"master_seed\":" << o.seed << ",\n\"N_seed\":" << n_seed
             << ",\n\"N_seed_rule\":\"splitmix64(master_seed XOR splitmix64(N XOR 0x5033333746345631))\""
             << ",\n\"replica_begin\":" << static_cast<std::uint64_t>(o.begin)*o.samples
             << ",\n\"replica_end_exclusive\":" << static_cast<std::uint64_t>(o.end)*o.samples
             << ",\n\"counter_rule\":\"batch*samples_per_batch+local; independent of threads and shard boundaries\""
             << ",\n\"representations\":[[" << g1.a << ',' << g1.b << "],[" << g2.a << ',' << g2.b << "]]"
             << ",\n\"fixed_face_sites\":[[" << face1[0] << ',' << face1[1] << ',' << face1[2] << ',' << face1[3]
             << "],[" << face2[0] << ',' << face2[1] << ',' << face2[2] << ',' << face2[3] << "]]"
             << ",\n\"ordinary_degree\":" << o.n << ",\n\"forced_degree\":" << o.n-4
             << ",\n\"freeze_commit\":\"" << o.freeze << "\",\n\"source_parent_commit\":\"" << f4_parent
             << "\",\n\"backend_blob\":\"" << f4_backend_blob << "\",\n\"build_git_commit\":\"" << P337_F4_BUILD_GIT_COMMIT
             << "\",\n\"compiler\":\"" << __VERSION__ << "\",\n\"requested_threads\":" << o.threads
             << ",\n\"actual_threads\":" << actual_threads << ",\n\"hardware_threads\":" << std::thread::hardware_concurrency()
             << ",\n\"openmp_version\":"
#ifdef _OPENMP
             << _OPENMP
#else
             << 0
#endif
             << ",\n\"elapsed_seconds\":" << elapsed << ",\n\"process_cpu_seconds\":" << cpu
             << ",\n\"histogram_rows\":" << rows
             << ",\n\"semantics\":\"One uniform permutation shared by both geometries and ordinary/forced modes. Forced occupies {0,a,b,a+b} first, then induces the free-site order from that permutation. Forward ambient rank>=1/rank2 birth counts use degree N or N-4; initial forced births may be0. Dense zero bins retained.\""
             << ",\n\"readout\":\"q=r-1=-1+CDF(first)+CDF(second); E=q^2=1-CDF(first)+CDF(second). Each mode uses its own degree and count coordinate.\"\n}\n";
        meta.close();
        if (!meta) throw std::runtime_error("metadata write failed");
        std::cout << "COMPLETE N=" << o.n << " batches=[" << o.begin << ',' << o.end << ") paired_permutations="
                  << o.samples*batches << " histogram_rows=" << rows << " threads=" << actual_threads
                  << " elapsed_seconds=" << elapsed << " process_cpu_seconds=" << cpu << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
