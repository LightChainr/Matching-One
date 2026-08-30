// Issue #154 Phase-E pilot: a fixed-p local connectivity singlet measured on
// the same counter/permutation stream as the ambient-H1 birth clocks.
//
// This deliberately reuses the production integer-period backend in one
// translation unit.  The local rows are not functions of K: at fixed K they
// depend on which vertices occupy the permutation prefix.

#define main threshold_rank_integer_period_hidden_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main

#include <cmath>

namespace {

struct LocalCounts {
    std::uint64_t black_axis_pairs = 0;
    std::uint64_t white_matching_axis_pairs = 0;
};

LocalCounts local_counts(const Geometry& geometry,
                         const std::vector<int>& permutation, int occupied_count) {
    if (occupied_count < 0 || occupied_count > geometry.n) {
        throw std::invalid_argument("occupied count is outside [0,N]");
    }
    std::vector<std::uint8_t> occupied(geometry.n, 0);
    for (int index = 0; index < occupied_count; ++index) occupied[permutation[index]] = 1;

    HomologyUnionFind black(geometry.quotient), white(geometry.quotient);
    for (const Edge& edge : geometry.primal_edges) {
        if (occupied[edge.i] && occupied[edge.j]) black.add_edge(edge);
    }
    for (const Edge& edge : geometry.matching_edges) {
        if (!occupied[edge.i] && !occupied[edge.j]) white.add_edge(edge);
    }

    LocalCounts result;
    const std::array<Vector, 2> axes = {{{1, 0}, {0, 1}}};
    for (int center = 0; center < geometry.n; ++center) {
        const Vector point = geometry.quotient.representative(center);
        for (const Vector axis : axes) {
            const int left = geometry.quotient.label({point.x - axis.x, point.y - axis.y});
            const int right = geometry.quotient.label({point.x + axis.x, point.y + axis.y});
            if (occupied[left] && occupied[right] &&
                black.find(left).root == black.find(right).root) {
                ++result.black_axis_pairs;
            }
            if (!occupied[left] && !occupied[right] &&
                white.find(left).root == white.find(right).root) {
                ++result.white_matching_axis_pairs;
            }
        }
    }
    return result;
}

int fixed_p_count(int n, double probability, std::uint64_t seed, std::uint64_t replica) {
    SplitMixStream stream(splitmix64(
        seed ^ splitmix64(replica + 0x8cb92ba72f3d8dd7ULL)));
    int count = 0;
    constexpr double inverse_53 = 1.0 / 9007199254740992.0;
    for (int vertex = 0; vertex < n; ++vertex) {
        const double uniform = static_cast<double>(stream.next() >> 11) * inverse_53;
        count += uniform < probability;
    }
    return count;
}

struct PilotOptions {
    int n = 0;
    std::uint64_t samples = 0;
    int batches = 0;
    std::uint64_t seed = 0;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    double p_ref = 0.59274605079;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void pilot_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " --n 65|130 --samples S --batches B"
        << " --seed S --replica-offset K --output-prefix PATH [--threads T]"
        << " [--p-ref P] [--git-commit SHA] [--self-test]\n";
    std::exit(status);
}

PilotOptions parse_pilot_options(int argc, char** argv) {
    PilotOptions options;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        auto value = [&]() -> std::string {
            if (++index >= argc) pilot_usage(argv[0], 2);
            return argv[index];
        };
        if (arg == "--n") options.n = parse_number<int>(value(), arg);
        else if (arg == "--samples") options.samples = parse_number<std::uint64_t>(value(), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(value(), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(value(), arg);
        else if (arg == "--replica-offset") {
            options.replica_offset = parse_number<std::uint64_t>(value(), arg);
        } else if (arg == "--threads") options.threads = parse_number<int>(value(), arg);
        else if (arg == "--p-ref") options.p_ref = std::stod(value());
        else if (arg == "--git-commit") options.git_commit = value();
        else if (arg == "--output-prefix") options.output_prefix = value();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help" || arg == "-h") pilot_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if ((options.n != 65 && options.n != 130) || options.samples == 0 ||
        options.batches <= 0 || options.samples % options.batches != 0 ||
        options.output_prefix.empty() || !(options.p_ref > 0.0 && options.p_ref < 1.0)) {
        throw std::invalid_argument("invalid or incomplete frozen pilot options");
    }
    return options;
}

PairDesign pilot_design(int n) {
    if (n == 65) {
        return {65, 8, 1, {8, -1, 1, 8}, 7, 4, {7, -4, 4, 7}, "N65_q2_parent"};
    }
    if (n == 130) {
        return {130, 11, 3, {11, -3, 3, 11}, 9, 7, {9, -7, 7, 9}, "N130_q2_child"};
    }
    throw std::invalid_argument("pilot supports only N=65 or N=130");
}

struct PilotStats {
    std::uint64_t samples = 0;
    std::uint64_t sum_k1 = 0;
    std::uint64_t sum_k2 = 0;
    std::uint64_t sum_i0 = 0;
    std::uint64_t sum_i1 = 0;
    std::uint64_t sum_i2 = 0;
    std::uint64_t sum_black = 0;
    std::uint64_t sum_white = 0;

    void add(int k1, int k2, int k, const LocalCounts& local) {
        ++samples;
        sum_k1 += k1;
        sum_k2 += k2;
        sum_i0 += k < k1;
        sum_i1 += k1 <= k && k < k2;
        sum_i2 += k2 <= k;
        sum_black += local.black_axis_pairs;
        sum_white += local.white_matching_axis_pairs;
    }
};

struct PilotBatch { PilotStats first, second; };

void pilot_self_test() {
    self_test();
    const Geometry geometry = make_geometry({2, -1, 1, 2});
    std::vector<int> permutation(geometry.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    const LocalCounts empty = local_counts(geometry, permutation, 0);
    const LocalCounts full = local_counts(geometry, permutation, geometry.n);
    if (empty.black_axis_pairs != 0 || empty.white_matching_axis_pairs != 2ULL * geometry.n ||
        full.black_axis_pairs != 2ULL * geometry.n || full.white_matching_axis_pairs != 0) {
        throw std::runtime_error("local connectivity empty/full regression failed");
    }
    const int first = fixed_p_count(65, 0.59274605079, 17, 23);
    const int second = fixed_p_count(65, 0.59274605079, 17, 23);
    if (first != second || first < 0 || first > 65) {
        throw std::runtime_error("fixed-p counter stream regression failed");
    }
    std::cout << "p154 local-singlet self-test passed\n";
}

int pilot_run(int argc, char** argv) {
    const PilotOptions options = parse_pilot_options(argc, argv);
    if (options.self_test) { pilot_self_test(); return 0; }
    const PairDesign design = pilot_design(options.n);
    const Geometry first_geometry = make_geometry(design.first);
    const Geometry second_geometry = make_geometry(design.second);
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<PilotBatch> output(options.batches);
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        PilotBatch local;
        std::vector<int> permutation;
        const std::uint64_t begin = options.replica_offset +
            static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(options.n, options.seed, replica, permutation);
            const int k = fixed_p_count(options.n, options.p_ref, options.seed, replica);
            const auto first_rank = first_engine.ranks(permutation);
            const auto second_rank = second_engine.ranks(permutation);
            local.first.add(first_rank.first, first_rank.second, k,
                            local_counts(first_geometry, permutation, k));
            local.second.add(second_rank.first, second_rank.second, k,
                             local_counts(second_geometry, permutation, k));
        }
        output[batch] = local;
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const auto csv_path = std::filesystem::path(options.output_prefix.string() + ".batches.csv");
    const auto metadata_path = std::filesystem::path(options.output_prefix.string() + ".metadata.json");
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open batch output");
    csv << "n,a,b,orientation,batch,samples,sum_k1,sum_k2,sum_i0,sum_i1,sum_i2,"
           "sum_black_axis_pairs,sum_white_matching_axis_pairs\n";
    auto write = [&](int batch, const char* orientation, int a, int b, const PilotStats& row) {
        csv << options.n << ',' << a << ',' << b << ',' << orientation << ',' << batch << ','
            << row.samples << ',' << row.sum_k1 << ',' << row.sum_k2 << ',' << row.sum_i0 << ','
            << row.sum_i1 << ',' << row.sum_i2 << ',' << row.sum_black << ',' << row.sum_white << '\n';
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write(batch, "first", design.a1, design.b1, output[batch].first);
        write(batch, "second", design.a2, design.b2, output[batch].second);
    }
    csv.close();
    std::ofstream metadata(metadata_path);
    metadata << "{\n"
             << "  \"engine\": \"P154 fixed-p local connectivity singlet pilot\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"N\": " << options.n << ",\n"
             << "  \"first\": [" << design.a1 << ',' << design.b1 << "],\n"
             << "  \"second\": [" << design.a2 << ',' << design.b2 << "],\n"
             << "  \"first_period_matrix\": " << matrix_json(design.first) << ",\n"
             << "  \"second_period_matrix\": " << matrix_json(design.second) << ",\n"
             << "  \"samples\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": " << options.replica_offset + options.samples << ",\n"
             << "  \"p_ref\": " << std::setprecision(17) << options.p_ref << ",\n"
             << "  \"fixed_p_stream\": \"independent counter-derived SplitMix64 Bernoulli prefix length; same K shared by orientation pair\",\n"
             << "  \"black_row\": \"mean over 2N axis placements of NN cluster connectivity between z-e and z+e\",\n"
             << "  \"white_row\": \"mean over 2N axis placements of white matching-cluster connectivity between z-e and z+e\",\n"
             << "  \"elapsed_seconds\": " << elapsed << ",\n"
             << "  \"batches_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed P154 N=" << options.n << " samples=" << options.samples
              << " elapsed=" << elapsed << " s\nwrote " << csv_path << "\nwrote "
              << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return pilot_run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
