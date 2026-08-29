// Direct p=1/2 score-function Monte Carlo for the C4 self-matching family.
//
// The implementation reuses the exact Gaussian-quotient homology kernel from
// threshold_rank_orientation_mc.cpp.  No finite differences in t or lambda
// are taken: every configuration supplies both Bernoulli likelihood scores.

#define main threshold_rank_orientation_embedded_main
#include "threshold_rank_orientation_mc.cpp"
#undef main

namespace {

constexpr std::array<const char*, 5> kChannelNames = {
    "cross", "both", "either", "direction_0", "direction_1"
};

Geometry make_self_matching_geometry(int a, int b) {
    if (a % 2 != 1 || b % 2 != 1) {
        throw std::invalid_argument("checkerboard Gaussian periods require odd a,b");
    }
    Geometry geometry = make_geometry(a, b);
    geometry.primal_edges.clear();
    geometry.matching_edges.clear();
    geometry.primal_edges.reserve(3 * geometry.n);
    const std::array<std::tuple<int, int, int>, 4> steps = {{
        {a, 1, 0}, {b, 0, 1}, {a + b, 1, 1}, {a - b, 1, -1},
    }};
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        for (std::size_t index = 0; index < steps.size(); ++index) {
            if (index >= 2 && vertex % 2 != 0) continue;
            const auto [residue, dx, dy] = steps[index];
            geometry.primal_edges.push_back(
                {vertex, positive_mod(vertex + residue, geometry.n), dx, dy});
        }
    }
    if (static_cast<int>(geometry.primal_edges.size()) != 3 * geometry.n) {
        throw std::logic_error("self-matching checkerboard must have 3N edges");
    }
    geometry.matching_edges = geometry.primal_edges;
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = geometry.primal_incident;
    return geometry;
}

class ConfigurationClassifier {
  public:
    explicit ConfigurationClassifier(const Geometry& geometry)
        : geometry_(geometry), union_find_(geometry.n, geometry.a, geometry.b) {}

    std::array<bool, 5> classify(const std::vector<std::uint8_t>& active) {
        union_find_.reset();
        for (const Edge& edge : geometry_.primal_edges) {
            if (active[edge.i] && active[edge.j]) union_find_.add_edge(edge);
        }
        bool cross = false;
        bool direction_0 = false;
        bool direction_1 = false;
        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (!active[vertex]) continue;
            cross = cross || union_find_.component_rank(vertex) == 2;
            direction_0 = direction_0 || union_find_.component_direction_0(vertex);
            direction_1 = direction_1 || union_find_.component_direction_1(vertex);
        }
        const bool either = direction_0 || direction_1;
        const bool both = direction_0 && direction_1;
        return {cross, both, either, direction_0, direction_1};
    }

  private:
    const Geometry& geometry_;
    HomologyUnionFind union_find_;
};

struct BatchStats {
    std::uint64_t samples = 0;
    std::int64_t sum_score_t = 0;
    std::int64_t sum_score_lambda = 0;
    std::uint64_t sum_score_t2 = 0;
    std::uint64_t sum_score_lambda2 = 0;
    std::int64_t sum_score_cross = 0;
    std::array<std::uint64_t, 5> wraps{};
    std::array<std::int64_t, 5> response_t{};
    std::array<std::int64_t, 5> response_lambda{};

    void add(int score_t, int score_lambda, const std::array<bool, 5>& channels) {
        ++samples;
        sum_score_t += score_t;
        sum_score_lambda += score_lambda;
        sum_score_t2 += static_cast<std::uint64_t>(score_t * score_t);
        sum_score_lambda2 += static_cast<std::uint64_t>(score_lambda * score_lambda);
        sum_score_cross += static_cast<std::int64_t>(score_t) * score_lambda;
        for (std::size_t index = 0; index < channels.size(); ++index) {
            if (!channels[index]) continue;
            ++wraps[index];
            response_t[index] += score_t;
            response_lambda[index] += score_lambda;
        }
    }
};

void counter_configuration(const Geometry& geometry, std::uint64_t seed,
                           std::uint64_t replica, std::vector<std::uint8_t>& active,
                           int& score_t, int& score_lambda) {
    active.resize(geometry.n);
    const std::uint64_t stream_key = splitmix64(
        seed ^ splitmix64(replica + 0x8cb92baa3f3d8dd7ULL));
    SplitMixStream generator(stream_key);
    int occupied_even = 0;
    int occupied_odd = 0;
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const bool occupied = (generator.next() >> 63) != 0;
        active[vertex] = occupied;
        if (occupied) {
            if (vertex % 2 == 0) ++occupied_even;
            else ++occupied_odd;
        }
    }
    score_t = 4 * (occupied_even + occupied_odd) - 2 * geometry.n;
    score_lambda = 4 * (occupied_even - occupied_odd);
}

void self_test_tangent() {
    const Geometry geometry = make_self_matching_geometry(3, 1);
    ConfigurationClassifier classifier(geometry);
    BatchStats exact;
    std::vector<std::uint8_t> active(geometry.n);
    for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << geometry.n); ++mask) {
        int occupied_even = 0;
        int occupied_odd = 0;
        for (int vertex = 0; vertex < geometry.n; ++vertex) {
            active[vertex] = (mask >> vertex) & 1U;
            if (active[vertex]) {
                if (vertex % 2 == 0) ++occupied_even;
                else ++occupied_odd;
            }
        }
        const int score_t = 4 * (occupied_even + occupied_odd) - 2 * geometry.n;
        const int score_lambda = 4 * (occupied_even - occupied_odd);
        exact.add(score_t, score_lambda, classifier.classify(active));
    }
    const std::uint64_t configurations = std::uint64_t{1} << geometry.n;
    for (std::size_t index = 0; index < kChannelNames.size(); ++index) {
        if (exact.response_t[index] * 8 != static_cast<std::int64_t>(15 * configurations) ||
            exact.response_lambda[index] * 4 !=
                static_cast<std::int64_t>(5 * configurations)) {
            throw std::runtime_error(std::string("N=10 score oracle failed for ") +
                                     kChannelNames[index]);
        }
    }
    if (exact.sum_score_t != 0 || exact.sum_score_lambda != 0 ||
        exact.sum_score_t2 != 4 * geometry.n * configurations ||
        exact.sum_score_lambda2 != 4 * geometry.n * configurations ||
        exact.sum_score_cross != 0) {
        throw std::runtime_error("N=10 Fisher oracle failed");
    }
    std::cout << "self-test passed: N=10 responses 15/8,5/4 and Fisher=4N I2\n";
}

struct TangentOptions {
    int a = 11;
    int b = 3;
    std::uint64_t samples = 1000000;
    int batches = 100;
    std::uint64_t seed = 2026105501;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void tangent_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --a A --b B          odd Gaussian period a+ib (default 11,3)\n"
        << "  --samples N          p=1/2 configurations (default 1000000)\n"
        << "  --batches B          aligned batches (default 100)\n"
        << "  --seed S             counter RNG seed\n"
        << "  --replica-offset K   first counter\n"
        << "  --threads T          OpenMP threads\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix PATH writes .responses.csv and .metadata.json\n"
        << "  --self-test          exact N=10 oracle then exit\n";
    std::exit(status);
}

TangentOptions parse_tangent_options(int argc, char** argv) {
    TangentOptions options;
    auto need = [&](int& index, const std::string& option) -> std::string {
        if (++index >= argc) throw std::invalid_argument(option + " needs a value");
        return argv[index];
    };
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--a") options.a = parse_number<int>(need(index, argument), argument);
        else if (argument == "--b") options.b = parse_number<int>(need(index, argument), argument);
        else if (argument == "--samples") options.samples = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--batches") options.batches = parse_number<int>(need(index, argument), argument);
        else if (argument == "--seed") options.seed = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--threads") options.threads = parse_number<int>(need(index, argument), argument);
        else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") tangent_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

int run_tangent(int argc, char** argv) {
    const TangentOptions options = parse_tangent_options(argc, argv);
    if (options.self_test) {
        self_test_tangent();
        return 0;
    }
    const Geometry geometry = make_self_matching_geometry(options.a, options.b);
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<BatchStats> batches(options.batches);
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        ConfigurationClassifier classifier(geometry);
        std::vector<std::uint8_t> active;
        BatchStats local;
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
            int score_t = 0;
            int score_lambda = 0;
            counter_configuration(
                geometry, options.seed, first + offset, active, score_t, score_lambda);
            local.add(score_t, score_lambda, classifier.classify(active));
        }
        batches[batch] = local;
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const std::filesystem::path parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string() + ".responses.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open response CSV");
    csv << "n,a,b,batch,samples,sum_score_t,sum_score_lambda,sum_score_t2,"
           "sum_score_lambda2,sum_score_cross";
    for (const char* name : kChannelNames) {
        csv << ',' << name << "_wraps," << name << "_score_t," << name
            << "_score_lambda";
    }
    csv << '\n';
    for (int batch = 0; batch < options.batches; ++batch) {
        const BatchStats& row = batches[batch];
        csv << geometry.n << ',' << options.a << ',' << options.b << ',' << batch
            << ',' << row.samples << ',' << row.sum_score_t << ','
            << row.sum_score_lambda << ',' << row.sum_score_t2 << ','
            << row.sum_score_lambda2 << ',' << row.sum_score_cross;
        for (std::size_t index = 0; index < kChannelNames.size(); ++index) {
            csv << ',' << row.wraps[index] << ',' << row.response_t[index]
                << ',' << row.response_lambda[index];
        }
        csv << '\n';
    }
    csv.close();
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    std::ofstream meta(metadata_path);
    if (!meta) throw std::runtime_error("cannot open metadata JSON");
    meta << "{\n"
         << "  \"engine\": \"C4 self-matching p=1/2 score-function tangent\",\n"
         << "  \"generated_utc\": \"" << utc_now() << "\",\n"
         << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
         << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
         << "  \"N\": " << geometry.n << ",\n"
         << "  \"a\": " << options.a << ", \"b\": " << options.b << ",\n"
         << "  \"samples\": " << options.samples << ", \"batches\": "
         << options.batches << ",\n"
         << "  \"seed\": " << options.seed << ",\n"
         << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
         << "  \"replica_counter_last_exclusive\": "
         << options.replica_offset + options.samples << ",\n"
         << "  \"threads_requested\": " << options.threads << ",\n"
         << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
         << "  \"score_order\": [\"t\",\"lambda\"],\n"
         << "  \"exact_Fisher_per_sample\": [[" << 4 * geometry.n
         << ",0],[0," << 4 * geometry.n << "]],\n"
         << "  \"response_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
         << "}\n";
    std::cout << "completed N=" << geometry.n << " samples=" << options.samples
              << " elapsed=" << elapsed << "\nwrote " << csv_path << "\nwrote "
              << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tangent(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
