// Direct p=1/2 score stream for the local complement-odd pivotal H4 readout.
//
// N=130 (11,3) and N=170 (13,1) use the same seed/counter stream.  Each
// configuration yields the exact Bernoulli scores S_t,S_lambda and two odd
// readouts: global cross half-difference and the fixed-root R=3 pivotal-H4
// half-difference frozen by the N=10 exact oracle.

#define main threshold_rank_orientation_embedded_main
#include "threshold_rank_orientation_mc.cpp"
#undef main

#include <cmath>

namespace {

struct LocalPoint {
    int x;
    int y;
    int vertex;
    int boundary_mask;
    bool root_neighbour;
};

Geometry make_c4_geometry(int a, int b) {
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
        throw std::logic_error("checkerboard triangulation must have 3N edges");
    }
    geometry.matching_edges = geometry.primal_edges;
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = geometry.primal_incident;
    return geometry;
}

class CrossClassifier {
  public:
    explicit CrossClassifier(const Geometry& geometry)
        : geometry_(geometry), union_find_(geometry.n, geometry.a, geometry.b) {}

    bool cross(const std::vector<std::uint8_t>& active) {
        union_find_.reset();
        for (const Edge& edge : geometry_.primal_edges) {
            if (active[edge.i] && active[edge.j]) union_find_.add_edge(edge);
        }
        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (active[vertex] && union_find_.component_crosses(vertex)) return true;
        }
        return false;
    }

  private:
    const Geometry& geometry_;
    HomologyUnionFind union_find_;
};

bool local_adjacent(int x1, int y1, int x2, int y2) {
    const int dx = x2 - x1;
    const int dy = y2 - y1;
    if (std::abs(dx) + std::abs(dy) == 1) return true;
    return std::abs(dx) == 1 && std::abs(dy) == 1 && positive_mod(x1 + y1, 2) == 0;
}

int landing_sector(int x, int y) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    int sector = static_cast<int>(std::floor(
        (std::atan2(static_cast<double>(y), static_cast<double>(x)) + pi / 8.0) /
        (pi / 4.0)));
    return positive_mod(sector, 8);
}

class LocalLanding {
  public:
    LocalLanding(const Geometry& geometry, int radius) : radius_(radius) {
        if (radius <= 0) throw std::invalid_argument("local radius must be positive");
        std::vector<int> seen(geometry.n, 0);
        for (int y = -radius; y <= radius; ++y) {
            for (int x = -radius; x <= radius; ++x) {
                if (x == 0 && y == 0) continue;
                const int vertex = positive_mod(geometry.a * x + geometry.b * y, geometry.n);
                if (seen[vertex]) {
                    throw std::invalid_argument("local annulus is not injective in quotient");
                }
                seen[vertex] = 1;
                const bool boundary = std::max(std::abs(x), std::abs(y)) == radius;
                points_.push_back({
                    x, y, vertex, boundary ? 1 << landing_sector(x, y) : 0,
                    local_adjacent(0, 0, x, y),
                });
            }
        }
        adjacency_.resize(points_.size());
        for (int first = 0; first < static_cast<int>(points_.size()); ++first) {
            for (int second = first + 1; second < static_cast<int>(points_.size()); ++second) {
                if (!local_adjacent(points_[first].x, points_[first].y,
                                    points_[second].x, points_[second].y)) continue;
                adjacency_[first].push_back(second);
                adjacency_[second].push_back(first);
            }
        }
    }

    int h4(const std::vector<std::uint8_t>& active) const {
        const std::vector<int> opened = component_masks(active, true);
        const std::vector<int> closed = component_masks(active, false);
        const bool axis =
            (distinct_pair(opened, 0, 4) && distinct_pair(closed, 2, 6)) ||
            (distinct_pair(opened, 2, 6) && distinct_pair(closed, 0, 4));
        const bool diagonal =
            (distinct_pair(opened, 1, 5) && distinct_pair(closed, 3, 7)) ||
            (distinct_pair(opened, 3, 7) && distinct_pair(closed, 1, 5));
        return static_cast<int>(axis) - static_cast<int>(diagonal);
    }

    int radius() const { return radius_; }

  private:
    std::vector<int> component_masks(const std::vector<std::uint8_t>& active,
                                     bool enabled) const {
        std::vector<std::uint8_t> unseen(points_.size(), 0);
        for (int index = 0; index < static_cast<int>(points_.size()); ++index) {
            unseen[index] = static_cast<bool>(active[points_[index].vertex]) == enabled;
        }
        std::vector<int> masks;
        std::vector<int> stack;
        for (int start = 0; start < static_cast<int>(points_.size()); ++start) {
            if (!unseen[start]) continue;
            unseen[start] = 0;
            stack.assign(1, start);
            bool touches_root = false;
            int mask = 0;
            while (!stack.empty()) {
                const int point = stack.back();
                stack.pop_back();
                touches_root = touches_root || points_[point].root_neighbour;
                mask |= points_[point].boundary_mask;
                for (const int neighbour : adjacency_[point]) {
                    if (!unseen[neighbour]) continue;
                    unseen[neighbour] = 0;
                    stack.push_back(neighbour);
                }
            }
            if (touches_root && mask) masks.push_back(mask);
        }
        return masks;
    }

    static bool distinct_pair(const std::vector<int>& masks, int first, int second) {
        for (int i = 0; i < static_cast<int>(masks.size()); ++i) {
            for (int j = 0; j < static_cast<int>(masks.size()); ++j) {
                if (i != j && (masks[i] & (1 << first)) &&
                    (masks[j] & (1 << second))) return true;
            }
        }
        return false;
    }

    int radius_;
    std::vector<LocalPoint> points_;
    std::vector<std::vector<int>> adjacency_;
};

int pivotal_h4(const Geometry& geometry, CrossClassifier& classifier,
               const LocalLanding& landing, std::vector<std::uint8_t>& active,
               bool original_cross) {
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
    active[root] = 0;
    const int mark = landing.h4(active);
    active[root] = original_root;
    const int pivotal = static_cast<int>(with_root) - static_cast<int>(without);
    if (pivotal != 0 && pivotal != 1) {
        throw std::logic_error("cross event failed monotonicity");
    }
    return pivotal * mark;
}

struct SampleValue {
    int score_t;
    int score_lambda;
    int global_twice;
    int local_twice;
};

SampleValue evaluate(const Geometry& geometry, CrossClassifier& classifier,
                     const LocalLanding& landing, std::vector<std::uint8_t>& black) {
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
    const int black_local = pivotal_h4(
        geometry, classifier, landing, black, black_cross);
    const int white_local = pivotal_h4(
        geometry, classifier, landing, white, white_cross);
    return {
        4 * (occupied_even + occupied_odd) - 2 * geometry.n,
        4 * (occupied_even - occupied_odd),
        static_cast<int>(black_cross) - static_cast<int>(white_cross),
        black_local - white_local,
    };
}

void counter_configuration(int n, std::uint64_t seed, std::uint64_t replica,
                           std::vector<std::uint8_t>& active) {
    active.resize(n);
    const std::uint64_t stream_key = splitmix64(
        seed ^ splitmix64(replica + 0x8cb92baa3f3d8dd7ULL));
    SplitMixStream generator(stream_key);
    for (int vertex = 0; vertex < n; ++vertex) {
        active[vertex] = (generator.next() >> 63) != 0;
    }
}

struct BatchStats {
    std::uint64_t samples = 0;
    std::int64_t sum_score_t = 0;
    std::int64_t sum_score_lambda = 0;
    std::uint64_t sum_score_t2 = 0;
    std::uint64_t sum_score_lambda2 = 0;
    std::int64_t sum_score_cross = 0;
    std::int64_t sum_global_twice = 0;
    std::int64_t sum_local_twice = 0;
    std::int64_t global_t = 0;
    std::int64_t global_lambda = 0;
    std::int64_t local_t = 0;
    std::int64_t local_lambda = 0;

    void add(const SampleValue& value) {
        ++samples;
        sum_score_t += value.score_t;
        sum_score_lambda += value.score_lambda;
        sum_score_t2 += static_cast<std::uint64_t>(value.score_t * value.score_t);
        sum_score_lambda2 +=
            static_cast<std::uint64_t>(value.score_lambda * value.score_lambda);
        sum_score_cross += static_cast<std::int64_t>(value.score_t) * value.score_lambda;
        sum_global_twice += value.global_twice;
        sum_local_twice += value.local_twice;
        global_t += static_cast<std::int64_t>(value.global_twice) * value.score_t;
        global_lambda +=
            static_cast<std::int64_t>(value.global_twice) * value.score_lambda;
        local_t += static_cast<std::int64_t>(value.local_twice) * value.score_t;
        local_lambda += static_cast<std::int64_t>(value.local_twice) * value.score_lambda;
    }
};

void self_test_local() {
    const Geometry geometry = make_c4_geometry(3, 1);
    CrossClassifier classifier(geometry);
    const LocalLanding landing(geometry, 1);
    BatchStats exact;
    std::vector<std::uint8_t> active(geometry.n);
    std::array<int, 5> local_counts{};
    for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << geometry.n); ++mask) {
        for (int vertex = 0; vertex < geometry.n; ++vertex) {
            active[vertex] = (mask >> vertex) & 1U;
        }
        const SampleValue value = evaluate(geometry, classifier, landing, active);
        exact.add(value);
        if (value.local_twice < -2 || value.local_twice > 2) {
            throw std::logic_error("local twice-observable outside exact range");
        }
        ++local_counts[value.local_twice + 2];
    }
    const std::int64_t configurations = std::int64_t{1} << geometry.n;
    if (exact.global_t * 8 != 15 * 2 * configurations ||
        exact.global_lambda * 4 != 5 * 2 * configurations ||
        exact.local_t * 64 != -3 * 2 * configurations ||
        exact.local_lambda * 64 != 11 * 2 * configurations ||
        local_counts != std::array<int, 5>{{0, 88, 848, 88, 0}}) {
        throw std::runtime_error("N=10 local pivotal response oracle failed");
    }
    std::cout << "self-test passed: exact N=10 local/global response matrix\n";
}

struct LocalOptions {
    std::uint64_t samples = 200000;
    int batches = 100;
    std::uint64_t seed = 15520260829ULL;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    int radius = 3;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void local_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N --batches B --seed S --replica-offset K\n"
        << "  --threads T --radius R --git-commit SHA\n"
        << "  --output-prefix PATH writes .batches.csv and .metadata.json\n"
        << "  --self-test exact N=10 R=1 oracle then exit\n";
    std::exit(status);
}

LocalOptions parse_local_options(int argc, char** argv) {
    LocalOptions options;
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
        else if (argument == "--radius") options.radius = parse_number<int>(need(index, argument), argument);
        else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") local_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (options.radius <= 0) throw std::invalid_argument("radius must be positive");
    if (options.replica_offset >
        std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

int run_local(int argc, char** argv) {
    const LocalOptions options = parse_local_options(argc, argv);
    if (options.self_test) {
        self_test_local();
        return 0;
    }
    const std::array<std::pair<int, int>, 2> designs = {{{11, 3}, {13, 1}}};
    const std::array<Geometry, 2> geometries = {{
        make_c4_geometry(11, 3), make_c4_geometry(13, 1),
    }};
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<std::array<BatchStats, 2>> output(options.batches);
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        CrossClassifier classifier130(geometries[0]);
        CrossClassifier classifier170(geometries[1]);
        const LocalLanding landing130(geometries[0], options.radius);
        const LocalLanding landing170(geometries[1], options.radius);
        std::array<CrossClassifier*, 2> classifiers = {{&classifier130, &classifier170}};
        std::array<const LocalLanding*, 2> landings = {{&landing130, &landing170}};
        std::array<BatchStats, 2> local;
        std::vector<std::uint8_t> active;
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
            const std::uint64_t replica = first + offset;
            for (int design = 0; design < 2; ++design) {
                counter_configuration(
                    geometries[design].n, options.seed, replica, active);
                local[design].add(evaluate(
                    geometries[design], *classifiers[design], *landings[design], active));
            }
        }
        output[batch] = local;
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const std::filesystem::path parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string() + ".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open batch CSV");
    csv << "n,a,b,batch,counter_first,counter_last_exclusive,samples,"
           "sum_score_t,sum_score_lambda,sum_score_t2,sum_score_lambda2,"
           "sum_score_cross,sum_global_twice,sum_local_twice,"
           "global_twice_score_t,global_twice_score_lambda,"
           "local_twice_score_t,local_twice_score_lambda\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (int design = 0; design < 2; ++design) {
            const BatchStats& row = output[batch][design];
            csv << geometries[design].n << ',' << designs[design].first << ','
                << designs[design].second << ',' << batch << ',' << first << ','
                << first + per_batch << ',' << row.samples << ',' << row.sum_score_t
                << ',' << row.sum_score_lambda << ',' << row.sum_score_t2 << ','
                << row.sum_score_lambda2 << ',' << row.sum_score_cross << ','
                << row.sum_global_twice << ',' << row.sum_local_twice << ','
                << row.global_t << ',' << row.global_lambda << ',' << row.local_t
                << ',' << row.local_lambda << '\n';
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
             << "  \"schema\": \"matching-one/c4-local-odd-pivotal-score-stream/v1\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
             << "  \"openmp\": "
#ifdef _OPENMP
             << "true,\n"
#else
             << "false,\n"
#endif
             << "  \"threads_requested\": " << options.threads << ",\n"
             << "  \"samples_per_size\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"cross_size_coupling\": \"same seed/counter and prefix-coupled site bits\",\n"
             << "  \"center\": \"p_even=p_odd=1/2\",\n"
             << "  \"radius\": " << options.radius << ",\n"
             << "  \"observable_rows\": [\"global_cross_half_difference\","
                "\"local_pivotal_h4_half_difference\"],\n"
             << "  \"score_columns\": [\"S_t\",\"S_lambda\"],\n"
             << "  \"designs\": [{\"N\":130,\"a\":11,\"b\":3},"
                "{\"N\":170,\"a\":13,\"b\":1}],\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"batch_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed N=130,N=170 samples=" << options.samples
              << " per size elapsed=" << elapsed << "\nwrote " << csv_path
              << "\nwrote " << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_local(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
