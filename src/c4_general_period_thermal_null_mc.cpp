// General-period checkerboard C4 score stream for Issue #155.
//
// The quotient/HNF/homology machinery is the frozen arbitrary-period backend.
// The graph is the self-matching checkerboard triangulation: every NN edge and
// both diagonals from even x+y sites.  Period columns must preserve parity.

#define THRESHOLD_RANK_INTEGER_PERIOD_NO_MAIN
#include "threshold_rank_integer_period_mc.cpp"
#undef THRESHOLD_RANK_INTEGER_PERIOD_NO_MAIN

#include <cmath>
#include <set>

namespace {

struct TNDesign {
    std::string label;
    Matrix periods;
};

struct TNOptions {
    std::uint64_t samples = 20000;
    int batches = 100;
    std::uint64_t seed = 15583020260830ULL;
    std::uint64_t replica_offset = 15500000000ULL;
    int threads = 0;
    std::vector<int> radii{2, 4, 8};
    std::string git_commit = "unknown";
    std::string binary_sha256 = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
    bool validate_only = false;
};

const std::array<TNDesign, 2> kTNDesigns{{
    {"N260_16_2", {16, -2, 2, 16}},
    {"N340_18_4", {18, -4, 4, 18}},
}};

[[noreturn]] void tn_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N --batches B --seed S --replica-offset K\n"
        << "  --threads T --radii 2,4,8 --git-commit SHA\n"
        << "  --binary-sha256 SHA --output-prefix PATH\n"
        << "  --validate-only checks R=8 injectivity; --self-test runs N=10 oracle\n";
    std::exit(status);
}

std::vector<int> tn_parse_radii(const std::string& text) {
    std::vector<int> result;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) throw std::invalid_argument("empty radius");
        result.push_back(parse_number<int>(token, "--radii"));
    }
    if (result != std::vector<int>({2, 4, 8})) {
        throw std::invalid_argument("Issue #155 frozen radii are exactly 2,4,8");
    }
    return result;
}

TNOptions tn_parse_options(int argc, char** argv) {
    TNOptions options;
    auto need = [&](int& index, const std::string& option) {
        if (++index >= argc) throw std::invalid_argument(option + " needs a value");
        return std::string(argv[index]);
    };
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--samples") options.samples = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--batches") options.batches = parse_number<int>(need(index, argument), argument);
        else if (argument == "--seed") options.seed = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(need(index, argument), argument);
        else if (argument == "--threads") options.threads = parse_number<int>(need(index, argument), argument);
        else if (argument == "--radii") options.radii = tn_parse_radii(need(index, argument));
        else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--binary-sha256") options.binary_sha256 = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--validate-only") options.validate_only = true;
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") tn_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test || options.validate_only) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples != 20000 && options.samples != 100000) {
        throw std::invalid_argument("frozen pilot allows only 20000 or one 100000 expansion");
    }
    if (options.batches != 100 || options.samples % 100 != 0) {
        throw std::invalid_argument("frozen pilot requires 100 equal batches");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

int parity(Vector point) {
    return static_cast<int>(positive_mod(point.x + point.y, 2));
}

Geometry make_checkerboard_geometry(Matrix periods) {
    if (positive_mod(periods.a + periods.c, 2) != 0 ||
        positive_mod(periods.b + periods.d, 2) != 0) {
        throw std::invalid_argument("period columns must preserve checkerboard parity");
    }
    Geometry geometry(periods);
    geometry.primal_edges.clear();
    geometry.matching_edges.clear();
    geometry.primal_edges.reserve(3 * geometry.n);
    const std::array<Vector, 4> steps{{{1, 0}, {0, 1}, {1, 1}, {1, -1}}};
    int even_vertices = 0;
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const Vector source = geometry.quotient.representative(vertex);
        const bool even = parity(source) == 0;
        even_vertices += even;
        for (std::size_t index = 0; index < steps.size(); ++index) {
            if (index >= 2 && !even) continue;
            const Vector step = steps[index];
            geometry.primal_edges.push_back({
                vertex,
                geometry.quotient.label({source.x + step.x, source.y + step.y}),
                static_cast<int>(step.x), static_cast<int>(step.y),
            });
        }
    }
    if (2 * even_vertices != geometry.n ||
        static_cast<int>(geometry.primal_edges.size()) != 3 * geometry.n) {
        throw std::logic_error("checkerboard quotient edge/parity count failed");
    }
    geometry.matching_edges = geometry.primal_edges;
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = geometry.primal_incident;
    return geometry;
}

class CrossClassifier {
  public:
    explicit CrossClassifier(const Geometry& geometry)
        : geometry_(geometry), union_find_(geometry.quotient) {}

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

bool checkerboard_adjacent(int x1, int y1, int x2, int y2) {
    const int dx = std::abs(x2 - x1);
    const int dy = std::abs(y2 - y1);
    if (dx + dy == 1) return true;
    return dx == 1 && dy == 1 && positive_mod(x1 + y1, 2) == 0;
}

int landing_sector(int x, int y) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const int sector = static_cast<int>(std::floor(
        (std::atan2(static_cast<double>(y), static_cast<double>(x)) + pi / 8.0) /
        (pi / 4.0)));
    return static_cast<int>(positive_mod(sector, 8));
}

struct Point {
    int x;
    int y;
    int vertex;
    int boundary_mask;
    bool root_neighbour;
};

class LocalLanding {
  public:
    LocalLanding(const Geometry& geometry, int radius, bool euclidean = true)
        : radius_(radius), euclidean_(euclidean) {
        std::vector<int> seen(geometry.n, 0);
        for (int y = -radius; y <= radius; ++y) {
            for (int x = -radius; x <= radius; ++x) {
                if (x == 0 && y == 0) continue;
                const int norm2 = x * x + y * y;
                if (euclidean && norm2 > radius * radius) continue;
                const int vertex = geometry.quotient.label({x, y});
                if (seen[vertex]) {
                    throw std::invalid_argument(
                        "Euclidean landing is not injective for N=" +
                        std::to_string(geometry.n) + ",R=" + std::to_string(radius));
                }
                seen[vertex] = 1;
                const bool boundary = euclidean
                    ? norm2 > (radius - 1) * (radius - 1)
                    : std::max(std::abs(x), std::abs(y)) == radius;
                points_.push_back({x, y, vertex,
                    boundary ? 1 << landing_sector(x, y) : 0,
                    checkerboard_adjacent(0, 0, x, y)});
            }
        }
        adjacency_.resize(points_.size());
        for (int first = 0; first < static_cast<int>(points_.size()); ++first) {
            for (int second = first + 1; second < static_cast<int>(points_.size()); ++second) {
                if (!checkerboard_adjacent(points_[first].x, points_[first].y,
                                           points_[second].x, points_[second].y)) continue;
                adjacency_[first].push_back(second);
                adjacency_[second].push_back(first);
            }
        }
    }

    int h4(const std::vector<std::uint8_t>& active) const {
        const auto opened = component_masks(active, true);
        const auto closed = component_masks(active, false);
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
    bool euclidean_;
    std::vector<Point> points_;
    std::vector<std::vector<int>> adjacency_;
};

int pivotal_h4(CrossClassifier& classifier, const LocalLanding& landing,
               std::vector<std::uint8_t>& active, bool original_cross) {
    const bool original_root = active[0];
    bool without;
    bool with_root;
    if (original_root) {
        with_root = original_cross;
        active[0] = 0;
        without = classifier.cross(active);
    } else {
        without = original_cross;
        active[0] = 1;
        with_root = classifier.cross(active);
    }
    active[0] = 0;
    const int mark = landing.h4(active);
    active[0] = original_root;
    const int pivotal = static_cast<int>(with_root) - static_cast<int>(without);
    if (pivotal != 0 && pivotal != 1) throw std::logic_error("cross event nonmonotone");
    return pivotal * mark;
}

std::vector<std::uint8_t> configuration(int n, std::uint64_t seed,
                                        std::uint64_t replica) {
    std::vector<std::uint8_t> active(n, 0);
    const std::uint64_t counter_key = splitmix64(replica + 0xd1b54a32d192ed03ULL);
    for (int vertex = 0; vertex < n; ++vertex) {
        active[vertex] = (splitmix64(seed ^ counter_key ^
            splitmix64(static_cast<std::uint64_t>(vertex))) >> 63) != 0;
    }
    return active;
}

struct Sample {
    int score_t = 0;
    int score_lambda = 0;
    int global_twice = 0;
    int epsilon_sign_sum = 0;
    std::vector<int> local_twice;
};

Sample evaluate(const Geometry& geometry, CrossClassifier& classifier,
                const std::vector<LocalLanding>& landings,
                std::vector<std::uint8_t>& black) {
    Sample value;
    std::vector<std::uint8_t> white(geometry.n);
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const int sign = black[vertex] ? 1 : -1;
        value.score_t += 2 * sign;
        const Vector representative = geometry.quotient.representative(vertex);
        value.score_lambda += 2 * (parity(representative) == 0 ? sign : -sign);
        white[vertex] = !black[vertex];
    }
    const int root_even = geometry.quotient.label({0, 0});
    const int root_odd = geometry.quotient.label({1, 0});
    value.epsilon_sign_sum = (black[root_even] ? 1 : -1) +
                             (black[root_odd] ? 1 : -1);
    const bool black_cross = classifier.cross(black);
    const bool white_cross = classifier.cross(white);
    value.global_twice = static_cast<int>(black_cross) - static_cast<int>(white_cross);
    value.local_twice.reserve(landings.size());
    for (const auto& landing : landings) {
        const int black_local = pivotal_h4(classifier, landing, black, black_cross);
        const int white_local = pivotal_h4(classifier, landing, white, white_cross);
        value.local_twice.push_back(black_local - white_local);
    }
    return value;
}

struct RadiusStats {
    std::int64_t sum_local_twice = 0;
    std::int64_t local_twice_score_t = 0;
    std::int64_t local_twice_score_lambda = 0;
};

struct BatchStats {
    std::uint64_t samples = 0;
    std::int64_t sum_score_t = 0;
    std::int64_t sum_score_lambda = 0;
    std::uint64_t sum_score_t2 = 0;
    std::uint64_t sum_score_lambda2 = 0;
    std::int64_t sum_score_cross = 0;
    std::int64_t sum_global_twice = 0;
    std::int64_t global_twice_score_t = 0;
    std::int64_t global_twice_score_lambda = 0;
    std::int64_t sum_epsilon_sign = 0;
    std::int64_t epsilon_sign_score_t = 0;
    std::int64_t epsilon_sign_score_lambda = 0;
    std::vector<RadiusStats> radius;

    explicit BatchStats(std::size_t count = 0) : radius(count) {}

    void add(const Sample& value) {
        ++samples;
        sum_score_t += value.score_t;
        sum_score_lambda += value.score_lambda;
        sum_score_t2 += static_cast<std::uint64_t>(value.score_t * value.score_t);
        sum_score_lambda2 += static_cast<std::uint64_t>(value.score_lambda * value.score_lambda);
        sum_score_cross += static_cast<std::int64_t>(value.score_t) * value.score_lambda;
        sum_global_twice += value.global_twice;
        global_twice_score_t += static_cast<std::int64_t>(value.global_twice) * value.score_t;
        global_twice_score_lambda += static_cast<std::int64_t>(value.global_twice) * value.score_lambda;
        sum_epsilon_sign += value.epsilon_sign_sum;
        epsilon_sign_score_t += static_cast<std::int64_t>(value.epsilon_sign_sum) * value.score_t;
        epsilon_sign_score_lambda += static_cast<std::int64_t>(value.epsilon_sign_sum) * value.score_lambda;
        for (std::size_t index = 0; index < radius.size(); ++index) {
            radius[index].sum_local_twice += value.local_twice[index];
            radius[index].local_twice_score_t +=
                static_cast<std::int64_t>(value.local_twice[index]) * value.score_t;
            radius[index].local_twice_score_lambda +=
                static_cast<std::int64_t>(value.local_twice[index]) * value.score_lambda;
        }
    }
};

void tn_self_test() {
    Geometry geometry = make_checkerboard_geometry({3, -1, 1, 3});
    CrossClassifier classifier(geometry);
    // The exact N=10 oracle used the original Chebyshev R=1 registry.  The
    // production R=2,4,8 family below uses the general-period Euclidean
    // registry so R=8 remains injective on N260.
    const std::vector<LocalLanding> landings{LocalLanding(geometry, 1, false)};
    BatchStats totals(1);
    for (std::uint64_t mask = 0; mask < (1ULL << geometry.n); ++mask) {
        std::vector<std::uint8_t> active(geometry.n);
        for (int vertex = 0; vertex < geometry.n; ++vertex) active[vertex] = (mask >> vertex) & 1U;
        totals.add(evaluate(geometry, classifier, landings, active));
    }
    const std::int64_t count = 1LL << geometry.n;
    if (totals.global_twice_score_t * 8 != 15 * 2 * count ||
        totals.global_twice_score_lambda * 4 != 5 * 2 * count ||
        totals.radius[0].local_twice_score_t * 64 != -3 * 2 * count ||
        totals.radius[0].local_twice_score_lambda * 64 != 11 * 2 * count ||
        totals.epsilon_sign_score_t != 4 * count ||
        totals.epsilon_sign_score_lambda != 0) {
        std::ostringstream detail;
        detail << "N=10 exact thermal-null response oracle failed: global="
               << totals.global_twice_score_t << ',' << totals.global_twice_score_lambda
               << " local=" << totals.radius[0].local_twice_score_t << ','
               << totals.radius[0].local_twice_score_lambda << " epsilon="
               << totals.epsilon_sign_score_t << ',' << totals.epsilon_sign_score_lambda;
        throw std::runtime_error(detail.str());
    }
    std::cout << "self-test passed: N=10 exact rows global=(15/8,5/4), "
                 "local=(-3/64,11/64), epsilon=(1,0), alpha*=3/64\n";
}

int tn_run(int argc, char** argv) {
    const TNOptions options = tn_parse_options(argc, argv);
    if (options.self_test) {
        tn_self_test();
        return 0;
    }
    std::array<Geometry, 2> geometries{{
        make_checkerboard_geometry(kTNDesigns[0].periods),
        make_checkerboard_geometry(kTNDesigns[1].periods),
    }};
    for (const auto& geometry : geometries) {
        for (const int radius : options.radii) (void)LocalLanding(geometry, radius);
    }
    if (options.validate_only) {
        std::cout << "validated N260,N340 checkerboard parity and injective Euclidean R=2,4,8\n";
        return 0;
    }
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<std::array<BatchStats, 2>> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) {
        output.push_back({BatchStats(options.radii.size()), BatchStats(options.radii.size())});
    }
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (int design = 0; design < 2; ++design) {
            CrossClassifier classifier(geometries[design]);
            std::vector<LocalLanding> landings;
            for (const int radius : options.radii) landings.emplace_back(geometries[design], radius);
            for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                auto active = configuration(geometries[design].n, options.seed, first + offset);
                output[batch][design].add(evaluate(
                    geometries[design], classifier, landings, active));
            }
        }
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string() + ".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream csv(csv_path);
    csv << "label,n,m00,m01,m10,m11,radius,batch,counter_first,counter_last_exclusive,samples,"
           "sum_score_t,sum_score_lambda,sum_score_t2,sum_score_lambda2,sum_score_cross,"
           "sum_global_twice,global_twice_score_t,global_twice_score_lambda,"
           "sum_epsilon_sign,epsilon_sign_score_t,epsilon_sign_score_lambda,"
           "sum_local_twice,local_twice_score_t,local_twice_score_lambda\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (int design = 0; design < 2; ++design) {
            const auto& common = output[batch][design];
            const Matrix matrix = kTNDesigns[design].periods;
            for (std::size_t index = 0; index < options.radii.size(); ++index) {
                const auto& row = common.radius[index];
                csv << kTNDesigns[design].label << ',' << geometries[design].n << ','
                    << matrix.a << ',' << matrix.b << ',' << matrix.c << ',' << matrix.d << ','
                    << options.radii[index] << ',' << batch << ',' << first << ','
                    << first + per_batch << ',' << common.samples << ','
                    << common.sum_score_t << ',' << common.sum_score_lambda << ','
                    << common.sum_score_t2 << ',' << common.sum_score_lambda2 << ','
                    << common.sum_score_cross << ',' << common.sum_global_twice << ','
                    << common.global_twice_score_t << ',' << common.global_twice_score_lambda << ','
                    << common.sum_epsilon_sign << ',' << common.epsilon_sign_score_t << ','
                    << common.epsilon_sign_score_lambda << ',' << row.sum_local_twice << ','
                    << row.local_twice_score_t << ',' << row.local_twice_score_lambda << '\n';
            }
        }
    }
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    std::ofstream metadata(metadata_path);
    metadata << "{\n"
             << "  \"schema\": \"matching-one/c4-general-period-thermal-null-score-stream/v1\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"binary_sha256\": \"" << json_escape(options.binary_sha256) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"samples_per_design\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": " << options.replica_offset + options.samples << ",\n"
             << "  \"radii\": [2,4,8],\n"
             << "  \"cutoff\": \"euclidean\",\n"
             << "  \"alpha_star\": \"3/64\",\n"
             << "  \"designs\": [{\"label\":\"N260_16_2\",\"N\":260,\"period_matrix\":[[16,-2],[2,16]]},"
                "{\"label\":\"N340_18_4\",\"N\":340,\"period_matrix\":[[18,-4],[4,18]]}],\n"
             << "  \"graph_semantics\": \"self-matching checkerboard triangulation on parity-preserving general-period quotient\",\n"
             << "  \"pairing\": \"configuration/complement half-differences; same seed/counter prefix field across sizes and radii\",\n"
             << "  \"observable_rows\": [\"global_cross_half_difference\",\"local_pivotal_h4_half_difference_plus_3_over_64_epsilon_cell\"],\n"
             << "  \"score_columns\": [\"S_t\",\"S_lambda\"],\n"
             << "  \"epsilon_cell\": \"[(n_(0,0)-1/2)+(n_(1,0)-1/2)]/2\",\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"batch_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed N260,N340 samples=" << options.samples
              << " radii=2,4,8 elapsed=" << elapsed << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return tn_run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
