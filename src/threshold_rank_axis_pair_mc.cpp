// Coupled adjacent-size axis threshold-rank engine for Issue #47.
//
// A uniform permutation of the L^2 upper torus is generated once per replica.
// The order induced on the fixed coordinate subset
//     {(x,y): 0<=x,y<L-1}
// is exactly uniform on the (L-1)^2 lower sites.  Thus filtering the upper
// permutation gives an exact-marginal lower permutation while inducing strong
// common-random-number covariance between the two finite tori.
//
// This is tailored to the Mertens-Ziff adjacent-size annihilator
//     L^(13/4) M_L - (L-1)^(13/4) M_(L-1).
//
// Reuse the audited axis geometry/homology/RNG implementation without copying
// it into a second source file.  Rename the included standalone main locally.
#define main matching_one_axis_single_main_disabled_here
#include "threshold_rank_axis_mc.cpp"
#undef main

namespace {

struct PairOutput {
    RankCounts upper;
    RankCounts lower;
    PairOutput(int upper_n, int lower_n) : upper(upper_n), lower(lower_n) {}
};

void restrict_to_lower(const std::vector<int>& upper, int L, std::vector<int>& lower) {
    const int lower_L = L - 1;
    lower.clear();
    lower.reserve(lower_L * lower_L);
    for (const int vertex : upper) {
        const int x = vertex % L;
        const int y = vertex / L;
        if (x < lower_L && y < lower_L) lower.push_back(x + lower_L * y);
    }
    if (static_cast<int>(lower.size()) != lower_L * lower_L) {
        throw std::logic_error("restriction did not produce every lower vertex");
    }
}

struct PairOptions {
    int L = 0;
    std::uint64_t samples = 1000000;
    int batches = 100;
    std::uint64_t seed = 20260828;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void pair_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --L L                upper axis size; lower is L-1 (required)\n"
        << "  --samples N          coupled replicas (default 1000000)\n"
        << "  --batches B          equal aligned batches (default 100)\n"
        << "  --seed S             unsigned 64-bit seed\n"
        << "  --replica-offset K   first sample counter\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix P    writes .hist.csv/.moments.csv/.metadata.json\n"
        << "  --self-test           exact restriction/oracle checks and exit\n";
    std::exit(status);
}

PairOptions parse_pair_options(int argc, char** argv) {
    PairOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string { if (++i >= argc) pair_usage(argv[0], 2); return argv[i]; };
        if (arg == "--L") options.L = parse_number<int>(next(), arg);
        else if (arg == "--samples") options.samples = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(next(), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(next(), arg);
        else if (arg == "--git-commit") options.git_commit = next();
        else if (arg == "--output-prefix") options.output_prefix = next();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") pair_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if (options.L < 3) throw std::invalid_argument("--L must be at least 3 for adjacent production");
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be positive and divisible by batches>=2");
    }
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("replica counter range overflows uint64");
    }
    return options;
}

void pair_self_test() {
    // Keep the exact threshold oracle from the standalone axis engine.
    self_test();

    // Restriction of all 4! permutations to a fixed 3-element subset must
    // produce each 3! ordering exactly four times.
    std::vector<int> upper = {0, 1, 2, 3};
    std::array<int, 6> multiplicity{};
    auto rank3 = [](const std::vector<int>& order) {
        static const std::array<std::array<int, 3>, 6> all = {{{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}};
        for (int i = 0; i < 6; ++i) {
            if (std::equal(order.begin(), order.end(), all[i].begin())) return i;
        }
        throw std::logic_error("unexpected restricted order");
    };
    do {
        std::vector<int> restricted;
        for (const int value : upper) if (value != 3) restricted.push_back(value);
        ++multiplicity[rank3(restricted)];
    } while (std::next_permutation(upper.begin(), upper.end()));
    for (const int count : multiplicity) {
        if (count != 4) throw std::runtime_error("permutation restriction is not uniform");
    }

    // Geometric restriction mapping for L=3 must contain each L=2 label once.
    std::vector<int> large(9), lower;
    std::iota(large.begin(), large.end(), 0);
    restrict_to_lower(large, 3, lower);
    std::sort(lower.begin(), lower.end());
    if (lower != std::vector<int>({0, 1, 2, 3})) {
        throw std::runtime_error("axis coordinate restriction regression failed");
    }
    std::cout << "pair self-test passed: exact permutation restriction is uniform\n";
}

int run_pair(int argc, char** argv) {
    const PairOptions options = parse_pair_options(argc, argv);
    if (options.self_test) { pair_self_test(); return 0; }
    const Geometry upper_geometry = make_axis_geometry(options.L);
    const Geometry lower_geometry = make_axis_geometry(options.L - 1);
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<PairOutput> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) {
        output.emplace_back(upper_geometry.n, lower_geometry.n);
    }
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        PairOutput local(upper_geometry.n, lower_geometry.n);
        ThresholdEngine upper_engine(upper_geometry);
        ThresholdEngine lower_engine(lower_geometry);
        std::vector<int> upper_permutation, lower_permutation;
        const std::uint64_t begin = options.replica_offset + static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(upper_geometry.n, options.seed, replica, upper_permutation);
            restrict_to_lower(upper_permutation, options.L, lower_permutation);
            const auto upper_ranks = upper_engine.ranks(upper_permutation);
            const auto lower_ranks = lower_engine.ranks(lower_permutation);
            local.upper.add(upper_ranks.first, upper_ranks.second);
            local.lower.add(lower_ranks.first, lower_ranks.second);
        }
        output[batch] = std::move(local);
    }

    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const auto hist_path = std::filesystem::path(options.output_prefix.string() + ".hist.csv");
    const auto moments_path = std::filesystem::path(options.output_prefix.string() + ".moments.csv");
    const auto meta_path = std::filesystem::path(options.output_prefix.string() + ".metadata.json");
    std::ofstream hist(hist_path), moments(moments_path);
    if (!hist || !moments) throw std::runtime_error("cannot open output files");
    hist << "pair_L,n,L,role,batch,samples,kind,k,count\n";
    moments << "pair_L,n,L,role,batch,samples,sum_kminus,sum_kplus,sum_kminus2,sum_kplus2,sum_product,sum_gap,sum_gap2\n";

    auto write_role = [&](int batch, const char* role, int L, const RankCounts& counts) {
        const int n = L * L;
        for (int rank = 1; rank <= n; ++rank) {
            if (counts.minus[rank]) hist << options.L << ',' << n << ',' << L << ',' << role << ',' << batch << ',' << counts.samples << ",minus," << rank << ',' << counts.minus[rank] << '\n';
            if (counts.plus[rank]) hist << options.L << ',' << n << ',' << L << ',' << role << ',' << batch << ',' << counts.samples << ",plus," << rank << ',' << counts.plus[rank] << '\n';
        }
        moments << options.L << ',' << n << ',' << L << ',' << role << ',' << batch << ',' << counts.samples << ','
                << counts.sum_minus << ',' << counts.sum_plus << ',' << counts.sum_minus2 << ','
                << counts.sum_plus2 << ',' << counts.sum_product << ',' << counts.sum_gap << ',' << counts.sum_gap2 << '\n';
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write_role(batch, "upper", options.L, output[batch].upper);
        write_role(batch, "lower", options.L - 1, output[batch].lower);
    }
    hist.close(); moments.close();

    std::ostringstream command;
    for (int i = 0; i < argc; ++i) { if (i) command << ' '; command << argv[i]; }
    std::ofstream meta(meta_path);
    meta << "{\n"
         << "  \"engine\": \"coupled adjacent-axis threshold-rank Newman-Ziff\",\n"
         << "  \"coupling\": \"uniform upper permutation restricted to common (L-1)^2 coordinate subset\",\n"
         << "  \"generated_utc\": \"" << utc_now() << "\",\n"
         << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
         << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
         << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
         << "  \"upper_L\": " << options.L << ",\n"
         << "  \"lower_L\": " << (options.L - 1) << ",\n"
         << "  \"samples\": " << options.samples << ",\n"
         << "  \"batches\": " << options.batches << ",\n"
         << "  \"seed\": " << options.seed << ",\n"
         << "  \"replica_offset\": " << options.replica_offset << ",\n"
         << "  \"threads_requested\": " << options.threads << ",\n"
#ifdef _OPENMP
         << "  \"openmp\": true\n";
#else
         << "  \"openmp\": false\n";
#endif
    meta << "}\n";
    std::cout << "completed coupled axis pair L=" << options.L << "/" << (options.L - 1)
              << " samples=" << options.samples << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try { return run_pair(argc, argv); }
    catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
