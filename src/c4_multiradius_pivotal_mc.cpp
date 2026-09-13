// Multi-radius extension of the fixed-root C4 pivotal-H4 score stream.
//
// One Bernoulli field supplies every radius.  Global crossing/pivotal status is
// evaluated once per colour and then the nested landing marks are read from
// the root-deleted environment.  The long batch format preserves aligned
// covariance across radii, plus/minus channels and N=130/170.

#define C4_LOCAL_ODD_NO_MAIN
#include "c4_local_odd_pivotal_mc.cpp"
#undef C4_LOCAL_ODD_NO_MAIN

#include <set>

namespace {

struct MultiOptions {
    std::uint64_t samples = 20000;
    int batches = 100;
    std::uint64_t seed = 22520260829ULL;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    std::vector<int> radii = {1, 2, 4};
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

std::vector<int> parse_radii(const std::string& text) {
    std::vector<int> radii;
    std::set<int> seen;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) throw std::invalid_argument("--radii contains an empty item");
        const int radius = parse_number<int>(token, "--radii");
        if (radius <= 0) throw std::invalid_argument("radii must be positive");
        if (!seen.insert(radius).second) throw std::invalid_argument("radii must be unique");
        radii.push_back(radius);
    }
    if (radii.empty()) throw std::invalid_argument("--radii requires at least one radius");
    if (!std::is_sorted(radii.begin(), radii.end())) {
        throw std::invalid_argument("radii must be strictly increasing");
    }
    return radii;
}

[[noreturn]] void multi_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N --batches B --seed S --replica-offset K\n"
        << "  --threads T --radii 1,2,4 --git-commit SHA\n"
        << "  --output-prefix PATH writes .batches.csv and .metadata.json\n"
        << "  --self-test runs the inherited exact N=10,R=1 oracle\n";
    std::exit(status);
}

MultiOptions parse_multi_options(int argc, char** argv) {
    MultiOptions options;
    auto need = [&](int& index, const std::string& option) -> std::string {
        if (++index >= argc) throw std::invalid_argument(option + " needs a value");
        return argv[index];
    };
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--samples") options.samples = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--batches") options.batches = parse_number<int>(need(index, argument), argument);
        else if (argument == "--seed") options.seed = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--threads") options.threads = parse_number<int>(need(index, argument), argument);
        else if (argument == "--radii") options.radii = parse_radii(need(index, argument));
        else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") multi_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (options.replica_offset >
        std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

int pivotal_flag(const Geometry& geometry, CrossClassifier& classifier,
                 std::vector<std::uint8_t>& active, bool original_cross) {
    const int root = 0;
    const bool original_root = active[root];
    bool without = false;
    bool with_root = false;
    if (original_root) {
        with_root = original_cross;
        active[root] = 0;
        without = classifier.cross(active);
    } else {
        without = original_cross;
        active[root] = 1;
        with_root = classifier.cross(active);
    }
    active[root] = original_root;
    const int pivotal = static_cast<int>(with_root) - static_cast<int>(without);
    if (pivotal != 0 && pivotal != 1) {
        throw std::logic_error("cross event failed monotonicity");
    }
    return pivotal;
}

struct MultiSample {
    int score_t = 0;
    int score_lambda = 0;
    int global_twice = 0;
    int black_pivotal = 0;
    int white_pivotal = 0;
    std::vector<int> black_h4;
    std::vector<int> white_h4;
};

MultiSample evaluate_multi(const Geometry& geometry, CrossClassifier& classifier,
                           const std::vector<LocalLanding>& landings,
                           std::vector<std::uint8_t>& black) {
    MultiSample value;
    int occupied_even = 0;
    int occupied_odd = 0;
    std::vector<std::uint8_t> white(geometry.n);
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        if (black[vertex]) {
            if (vertex % 2 == 0) ++occupied_even;
            else ++occupied_odd;
        }
        white[vertex] = !black[vertex];
    }
    const bool black_cross = classifier.cross(black);
    const bool white_cross = classifier.cross(white);
    value.black_pivotal = pivotal_flag(geometry, classifier, black, black_cross);
    value.white_pivotal = pivotal_flag(geometry, classifier, white, white_cross);
    value.score_t = 4 * (occupied_even + occupied_odd) - 2 * geometry.n;
    value.score_lambda = 4 * (occupied_even - occupied_odd);
    value.global_twice = static_cast<int>(black_cross) - static_cast<int>(white_cross);
    const bool black_root = black[0];
    const bool white_root = white[0];
    black[0] = 0;
    white[0] = 0;
    value.black_h4.reserve(landings.size());
    value.white_h4.reserve(landings.size());
    for (const LocalLanding& landing : landings) {
        value.black_h4.push_back(value.black_pivotal * landing.h4(black));
        value.white_h4.push_back(value.white_pivotal * landing.h4(white));
    }
    black[0] = black_root;
    white[0] = white_root;
    return value;
}

struct CommonStats {
    std::uint64_t samples = 0;
    std::int64_t sum_score_t = 0;
    std::int64_t sum_score_lambda = 0;
    std::int64_t sum_global_twice = 0;
    std::int64_t global_t = 0;
    std::int64_t global_lambda = 0;
    std::int64_t black_pivotal = 0;
    std::int64_t white_pivotal = 0;
};

struct RadiusStats {
    std::int64_t black_h4 = 0;
    std::int64_t white_h4 = 0;
    std::int64_t h4_plus = 0;
    std::int64_t h4_minus = 0;
    std::int64_t plus_t = 0;
    std::int64_t plus_lambda = 0;
    std::int64_t minus_t = 0;
    std::int64_t minus_lambda = 0;
};

struct MultiBatch {
    CommonStats common;
    std::vector<RadiusStats> radius;

    explicit MultiBatch(std::size_t count = 0) : radius(count) {}

    void add(const MultiSample& value) {
        ++common.samples;
        common.sum_score_t += value.score_t;
        common.sum_score_lambda += value.score_lambda;
        common.sum_global_twice += value.global_twice;
        common.global_t += static_cast<std::int64_t>(value.global_twice) * value.score_t;
        common.global_lambda += static_cast<std::int64_t>(value.global_twice) * value.score_lambda;
        common.black_pivotal += value.black_pivotal;
        common.white_pivotal += value.white_pivotal;
        for (std::size_t index = 0; index < radius.size(); ++index) {
            const int black = value.black_h4[index];
            const int white = value.white_h4[index];
            const int plus = black + white;
            const int minus = black - white;
            radius[index].black_h4 += black;
            radius[index].white_h4 += white;
            radius[index].h4_plus += plus;
            radius[index].h4_minus += minus;
            radius[index].plus_t += static_cast<std::int64_t>(plus) * value.score_t;
            radius[index].plus_lambda += static_cast<std::int64_t>(plus) * value.score_lambda;
            radius[index].minus_t += static_cast<std::int64_t>(minus) * value.score_t;
            radius[index].minus_lambda += static_cast<std::int64_t>(minus) * value.score_lambda;
        }
    }
};

int run_multi(int argc, char** argv) {
    const MultiOptions options = parse_multi_options(argc, argv);
    if (options.self_test) {
        self_test_local();
        return 0;
    }
    const std::array<std::pair<int, int>, 2> designs = {{{11, 3}, {13, 1}}};
    const std::array<Geometry, 2> geometries = {{
        make_c4_geometry(11, 3), make_c4_geometry(13, 1),
    }};
    // Construct once before sampling so an aliased radius fails immediately.
    for (const Geometry& geometry : geometries) {
        for (const int radius : options.radii) {
            (void)LocalLanding(geometry, radius);
        }
    }
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<std::array<MultiBatch, 2>> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) {
        output.push_back({MultiBatch(options.radii.size()), MultiBatch(options.radii.size())});
    }
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        std::array<MultiBatch, 2> local = {{
            MultiBatch(options.radii.size()), MultiBatch(options.radii.size()),
        }};
        std::vector<std::uint8_t> active;
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (int design = 0; design < 2; ++design) {
            CrossClassifier classifier(geometries[design]);
            std::vector<LocalLanding> landings;
            landings.reserve(options.radii.size());
            for (const int radius : options.radii) {
                landings.emplace_back(geometries[design], radius);
            }
            for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                const std::uint64_t replica = first + offset;
                counter_configuration(geometries[design].n, options.seed, replica, active);
                local[design].add(evaluate_multi(
                    geometries[design], classifier, landings, active));
            }
        }
        output[batch] = std::move(local);
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const std::filesystem::path parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string() + ".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open batch CSV");
    csv << "n,a,b,radius,batch,counter_first,counter_last_exclusive,samples,"
           "sum_score_t,sum_score_lambda,sum_global_twice,global_twice_score_t,"
           "global_twice_score_lambda,black_pivotal,white_pivotal,black_h4,"
           "white_h4,h4_plus,h4_minus,h4_plus_score_t,h4_plus_score_lambda,"
           "h4_minus_score_t,h4_minus_score_lambda\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (int design = 0; design < 2; ++design) {
            const CommonStats& common = output[batch][design].common;
            for (std::size_t index = 0; index < options.radii.size(); ++index) {
                const RadiusStats& row = output[batch][design].radius[index];
                csv << geometries[design].n << ',' << designs[design].first << ','
                    << designs[design].second << ',' << options.radii[index] << ','
                    << batch << ',' << first << ',' << first + per_batch << ','
                    << common.samples << ',' << common.sum_score_t << ','
                    << common.sum_score_lambda << ',' << common.sum_global_twice << ','
                    << common.global_t << ',' << common.global_lambda << ','
                    << common.black_pivotal << ',' << common.white_pivotal << ','
                    << row.black_h4 << ',' << row.white_h4 << ',' << row.h4_plus
                    << ',' << row.h4_minus << ',' << row.plus_t << ','
                    << row.plus_lambda << ',' << row.minus_t << ','
                    << row.minus_lambda << '\n';
            }
        }
    }
    csv.close();
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    std::ofstream metadata(metadata_path);
    if (!metadata) throw std::runtime_error("cannot open metadata JSON");
    metadata << "{\n"
             << "  \"schema\": \"matching-one/c4-multiradius-pivotal/v1\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
             << "  \"samples_per_size\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"radii\": [";
    for (std::size_t index = 0; index < options.radii.size(); ++index) {
        if (index) metadata << ',';
        metadata << options.radii[index];
    }
    metadata << "],\n"
             << "  \"radius_semantics\": \"fixed lattice R; physical delta=R/sqrt(N)\",\n"
             << "  \"cross_radius_coupling\": \"same configuration and pivotal flag\",\n"
             << "  \"cross_size_coupling\": \"same seed/counter and prefix-coupled site bits\",\n"
             << "  \"channels\": [\"black_h4\",\"white_h4\",\"h4_plus\",\"h4_minus\"],\n"
             << "  \"designs\": [{\"N\":130,\"a\":11,\"b\":3},"
                "{\"N\":170,\"a\":13,\"b\":1}],\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"batch_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed N=130,N=170 radii=";
    for (std::size_t index = 0; index < options.radii.size(); ++index) {
        if (index) std::cout << ',';
        std::cout << options.radii[index];
    }
    std::cout << " samples=" << options.samples << " elapsed=" << elapsed
              << "\nwrote " << csv_path << "\nwrote " << metadata_path << '\n';
    return 0;
}

}  // namespace

#ifndef C4_MULTIRADIUS_NO_MAIN
int main(int argc, char** argv) {
    try {
        return run_multi(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
#endif
