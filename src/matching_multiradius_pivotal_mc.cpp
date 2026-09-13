// Runtime-Gaussian, multi-radius square-site/matching pivotal H4 stream.
//
// Unlike the checkerboard C4 control, odd-order norm-5 quotients use the
// ordinary NN primal graph and NN+NNN matching graph.  One configuration and
// one pivotal decision per graph feed every requested landing radius.

#define C4_MULTIRADIUS_NO_MAIN
#include "c4_multiradius_pivotal_mc.cpp"
#undef C4_MULTIRADIUS_NO_MAIN

#include <set>

namespace {

struct MatchingDesign {
    std::string label;
    int a;
    int b;
};

struct MatchingOptions {
    std::uint64_t samples = 200000;
    int batches = 200;
    std::uint64_t seed = 22550260829ULL;
    std::uint64_t replica_offset = 15000000000ULL;
    int threads = 0;
    double p = 0.592746050790;
    std::string cutoff = "euclidean";
    std::vector<int> radii = {2, 4, 7, 8};
    std::vector<MatchingDesign> designs = {
        {"n325_first", 17, 6}, {"n325_second", 18, 1},
        {"n425_first", 16, 13}, {"n425_second", 19, 8},
    };
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
    bool validate_only = false;
};

MatchingDesign parse_design(const std::string& text) {
    std::stringstream stream(text);
    std::array<std::string, 3> fields;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (!std::getline(stream, fields[index], ',')) {
            throw std::invalid_argument("--design requires label,a,b");
        }
    }
    std::string extra;
    if (std::getline(stream, extra)) throw std::invalid_argument("--design has too many fields");
    if (fields[0].empty() || fields[0].find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos) {
        throw std::invalid_argument("design label must be nonempty ASCII alnum/_/-");
    }
    const int a = parse_number<int>(fields[1], "--design a");
    const int b = parse_number<int>(fields[2], "--design b");
    if (a <= 0 || b < 0 || std::gcd(a, b) != 1) {
        throw std::invalid_argument("Gaussian reps require a>0,b>=0,gcd(a,b)=1");
    }
    return {fields[0], a, b};
}

[[noreturn]] void matching_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N --batches B --seed S --replica-offset K --p P\n"
        << "  --threads T --radii 2,4,7,8 --design label,a,b [--design ...]\n"
        << "  --cutoff euclidean|chebyshev\n"
        << "  --git-commit SHA --output-prefix PATH\n"
        << "  --validate-only checks every annulus before sampling\n"
        << "  --self-test runs the inherited exact checkerboard oracle\n";
    std::exit(status);
}

MatchingOptions parse_matching_options(int argc, char** argv) {
    MatchingOptions options;
    bool custom_designs = false;
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
        else if (argument == "--p") options.p = std::stod(need(index, argument));
        else if (argument == "--cutoff") options.cutoff = need(index, argument);
        else if (argument == "--radii") options.radii = parse_radii(need(index, argument));
        else if (argument == "--design") {
            if (!custom_designs) {
                options.designs.clear();
                custom_designs = true;
            }
            options.designs.push_back(parse_design(need(index, argument)));
        } else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--validate-only") options.validate_only = true;
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") matching_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test) return options;
    if (!options.validate_only && options.output_prefix.empty()) {
        throw std::invalid_argument("--output-prefix required");
    }
    if (options.designs.empty()) throw std::invalid_argument("at least one design required");
    std::set<std::string> labels;
    for (const auto& design : options.designs) {
        if (!labels.insert(design.label).second) throw std::invalid_argument("design labels must be unique");
    }
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (!(options.p > 0.0 && options.p < 1.0)) throw std::invalid_argument("p must be in (0,1)");
    if (options.cutoff != "euclidean" && options.cutoff != "chebyshev") {
        throw std::invalid_argument("cutoff must be euclidean or chebyshev");
    }
    if (options.replica_offset >
        std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

class MatchingCrossClassifier {
  public:
    MatchingCrossClassifier(const Geometry& geometry, bool matching)
        : geometry_(geometry), edges_(matching ? geometry.matching_edges : geometry.primal_edges),
          union_find_(geometry.n, geometry.a, geometry.b) {}

    bool cross(const std::vector<std::uint8_t>& active) {
        union_find_.reset();
        for (const Edge& edge : edges_) {
            if (active[edge.i] && active[edge.j]) union_find_.add_edge(edge);
        }
        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (active[vertex] && union_find_.component_crosses(vertex)) return true;
        }
        return false;
    }

  private:
    const Geometry& geometry_;
    const std::vector<Edge>& edges_;
    HomologyUnionFind union_find_;
};

struct MatchingPoint {
    int x;
    int y;
    int vertex;
    int boundary_mask;
    bool primal_root_neighbour;
    bool matching_root_neighbour;
};

bool graph_adjacent(int x1, int y1, int x2, int y2, bool matching) {
    const int dx = std::abs(x2 - x1);
    const int dy = std::abs(y2 - y1);
    if (dx + dy == 1) return true;
    return matching && dx == 1 && dy == 1;
}

class MatchingLanding {
  public:
    MatchingLanding(const Geometry& geometry, int radius, bool euclidean)
        : radius_(radius), euclidean_(euclidean) {
        if (radius <= 0) throw std::invalid_argument("local radius must be positive");
        std::vector<int> seen(geometry.n, 0);
        for (int y = -radius; y <= radius; ++y) {
            for (int x = -radius; x <= radius; ++x) {
                if (x == 0 && y == 0) continue;
                const int norm2 = x * x + y * y;
                if (euclidean && norm2 > radius * radius) continue;
                const int vertex = positive_mod(geometry.a * x + geometry.b * y, geometry.n);
                if (seen[vertex]) throw std::invalid_argument(
                    "local annulus is not injective for N=" + std::to_string(geometry.n) +
                    ",R=" + std::to_string(radius));
                seen[vertex] = 1;
                const bool boundary = euclidean
                    ? norm2 > (radius - 1) * (radius - 1)
                    : std::max(std::abs(x), std::abs(y)) == radius;
                points_.push_back({
                    x, y, vertex, boundary ? 1 << landing_sector(x, y) : 0,
                    graph_adjacent(0, 0, x, y, false),
                    graph_adjacent(0, 0, x, y, true),
                });
            }
        }
        primal_adjacency_.resize(points_.size());
        matching_adjacency_.resize(points_.size());
        for (int first = 0; first < static_cast<int>(points_.size()); ++first) {
            for (int second = first + 1; second < static_cast<int>(points_.size()); ++second) {
                if (graph_adjacent(points_[first].x, points_[first].y,
                                   points_[second].x, points_[second].y, false)) {
                    primal_adjacency_[first].push_back(second);
                    primal_adjacency_[second].push_back(first);
                }
                if (graph_adjacent(points_[first].x, points_[first].y,
                                   points_[second].x, points_[second].y, true)) {
                    matching_adjacency_[first].push_back(second);
                    matching_adjacency_[second].push_back(first);
                }
            }
        }
    }

    int h4(const std::vector<std::uint8_t>& active, bool open_matching) const {
        const auto opened = component_masks(active, true, open_matching);
        const auto closed = component_masks(active, false, !open_matching);
        const bool axis =
            (distinct_pair(opened, 0, 4) && distinct_pair(closed, 2, 6)) ||
            (distinct_pair(opened, 2, 6) && distinct_pair(closed, 0, 4));
        const bool diagonal =
            (distinct_pair(opened, 1, 5) && distinct_pair(closed, 3, 7)) ||
            (distinct_pair(opened, 3, 7) && distinct_pair(closed, 1, 5));
        return static_cast<int>(axis) - static_cast<int>(diagonal);
    }

  private:
    std::vector<int> component_masks(const std::vector<std::uint8_t>& active,
                                     bool enabled, bool matching) const {
        const auto& adjacency = matching ? matching_adjacency_ : primal_adjacency_;
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
                touches_root = touches_root || (matching
                    ? points_[point].matching_root_neighbour
                    : points_[point].primal_root_neighbour);
                mask |= points_[point].boundary_mask;
                for (const int neighbour : adjacency[point]) {
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
    std::vector<MatchingPoint> points_;
    std::vector<std::vector<int>> primal_adjacency_;
    std::vector<std::vector<int>> matching_adjacency_;
};

int matching_pivotal_flag(MatchingCrossClassifier& classifier,
                          std::vector<std::uint8_t>& active) {
    const bool original_root = active[0];
    active[0] = 0;
    const bool without = classifier.cross(active);
    active[0] = 1;
    const bool with_root = classifier.cross(active);
    active[0] = original_root;
    const int pivotal = static_cast<int>(with_root) - static_cast<int>(without);
    if (pivotal != 0 && pivotal != 1) throw std::logic_error("nonmonotone crossing event");
    return pivotal;
}

std::vector<std::uint8_t> bernoulli_configuration(
        int n, std::uint64_t seed, std::uint64_t replica, double p) {
    std::vector<std::uint8_t> active(n, 0);
    const std::uint64_t counter_key = splitmix64(replica + 0xd1b54a32d192ed03ULL);
    constexpr double inverse_53 = 1.0 / 9007199254740992.0;
    for (int vertex = 0; vertex < n; ++vertex) {
        const std::uint64_t value = splitmix64(
            seed ^ counter_key ^ splitmix64(static_cast<std::uint64_t>(vertex)));
        active[vertex] = static_cast<double>(value >> 11) * inverse_53 < p;
    }
    return active;
}

struct MatchingSample {
    int primal_pivotal = 0;
    int matching_pivotal = 0;
    std::vector<int> primal_h4;
    std::vector<int> matching_h4;
};

MatchingSample evaluate_matching(
        const Geometry& geometry, MatchingCrossClassifier& primal_classifier,
        MatchingCrossClassifier& matching_classifier,
        const std::vector<MatchingLanding>& landings,
        std::vector<std::uint8_t>& black) {
    std::vector<std::uint8_t> white(geometry.n);
    for (int vertex = 0; vertex < geometry.n; ++vertex) white[vertex] = !black[vertex];
    MatchingSample value;
    value.primal_pivotal = matching_pivotal_flag(primal_classifier, black);
    value.matching_pivotal = matching_pivotal_flag(matching_classifier, white);
    const bool black_root = black[0];
    const bool white_root = white[0];
    black[0] = 0;
    white[0] = 0;
    value.primal_h4.reserve(landings.size());
    value.matching_h4.reserve(landings.size());
    for (const auto& landing : landings) {
        value.primal_h4.push_back(value.primal_pivotal * landing.h4(black, false));
        value.matching_h4.push_back(value.matching_pivotal * landing.h4(white, true));
    }
    black[0] = black_root;
    white[0] = white_root;
    return value;
}

struct MatchingRadiusStats {
    std::int64_t primal_h4 = 0;
    std::int64_t matching_h4 = 0;
    std::int64_t h4_plus = 0;
    std::int64_t h4_minus = 0;
};

struct MatchingBatchStats {
    std::uint64_t samples = 0;
    std::uint64_t primal_pivotal = 0;
    std::uint64_t matching_pivotal = 0;
    std::vector<MatchingRadiusStats> radius;

    explicit MatchingBatchStats(std::size_t count = 0) : radius(count) {}

    void add(const MatchingSample& value) {
        ++samples;
        primal_pivotal += value.primal_pivotal;
        matching_pivotal += value.matching_pivotal;
        for (std::size_t index = 0; index < radius.size(); ++index) {
            const int primal = value.primal_h4[index];
            const int matching = value.matching_h4[index];
            radius[index].primal_h4 += primal;
            radius[index].matching_h4 += matching;
            radius[index].h4_plus += primal + matching;
            radius[index].h4_minus += primal - matching;
        }
    }
};

std::vector<Geometry> matching_geometries(const MatchingOptions& options) {
    std::vector<Geometry> geometries;
    geometries.reserve(options.designs.size());
    for (const auto& design : options.designs) {
        geometries.push_back(make_geometry(design.a, design.b));
    }
    return geometries;
}

void validate_matching_landings(const MatchingOptions& options,
                                const std::vector<Geometry>& geometries) {
    for (std::size_t design = 0; design < geometries.size(); ++design) {
        for (const int radius : options.radii) {
            (void)MatchingLanding(geometries[design], radius, options.cutoff == "euclidean");
        }
    }
}

int run_matching(int argc, char** argv) {
    const MatchingOptions options = parse_matching_options(argc, argv);
    if (options.self_test) {
        self_test_local();
        return 0;
    }
    const std::vector<Geometry> geometries = matching_geometries(options);
    validate_matching_landings(options, geometries);
    if (options.validate_only) {
        std::cout << "validated " << geometries.size() << " designs at radii";
        for (const int radius : options.radii) std::cout << ' ' << radius;
        std::cout << '\n';
        return 0;
    }
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<std::vector<MatchingBatchStats>> output(
        options.batches,
        std::vector<MatchingBatchStats>(
            options.designs.size(), MatchingBatchStats(options.radii.size())));
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::size_t design = 0; design < geometries.size(); ++design) {
            MatchingCrossClassifier primal_classifier(geometries[design], false);
            MatchingCrossClassifier matching_classifier(geometries[design], true);
            std::vector<MatchingLanding> landings;
            landings.reserve(options.radii.size());
            for (const int radius : options.radii) {
                landings.emplace_back(
                    geometries[design], radius, options.cutoff == "euclidean");
            }
            for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                const std::uint64_t replica = first + offset;
                auto black = bernoulli_configuration(
                    geometries[design].n, options.seed, replica, options.p);
                output[batch][design].add(evaluate_matching(
                    geometries[design], primal_classifier, matching_classifier,
                    landings, black));
            }
        }
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const std::filesystem::path parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string() + ".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open batch CSV");
    csv << "label,n,a,b,radius,batch,counter_first,counter_last_exclusive,samples,"
           "primal_pivotal,matching_pivotal,primal_h4,matching_h4,h4_plus,h4_minus\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::size_t design = 0; design < options.designs.size(); ++design) {
            const auto& common = output[batch][design];
            for (std::size_t index = 0; index < options.radii.size(); ++index) {
                const auto& row = common.radius[index];
                csv << options.designs[design].label << ',' << geometries[design].n << ','
                    << options.designs[design].a << ',' << options.designs[design].b << ','
                    << options.radii[index] << ',' << batch << ',' << first << ','
                    << first + per_batch << ',' << common.samples << ','
                    << common.primal_pivotal << ',' << common.matching_pivotal << ','
                    << row.primal_h4 << ',' << row.matching_h4 << ',' << row.h4_plus
                    << ',' << row.h4_minus << '\n';
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
             << "  \"schema\": \"matching-one/matching-multiradius-pivotal/v1\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"samples_per_design\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"p\": \"" << std::setprecision(12) << options.p << "\",\n"
             << "  \"cutoff\": \"" << options.cutoff << "\",\n"
             << "  \"radii\": [";
    for (std::size_t index = 0; index < options.radii.size(); ++index) {
        if (index) metadata << ',';
        metadata << options.radii[index];
    }
    metadata << "],\n  \"designs\": [";
    for (std::size_t index = 0; index < options.designs.size(); ++index) {
        if (index) metadata << ',';
        metadata << "{\"label\":\"" << json_escape(options.designs[index].label)
                 << "\",\"N\":" << geometries[index].n << ",\"a\":"
                 << options.designs[index].a << ",\"b\":" << options.designs[index].b
                 << "}";
    }
    metadata << "],\n"
             << "  \"graph_semantics\": \"black NN primal; complemented white NN+NNN matching\",\n"
             << "  \"landing_registry\": \"absolute lattice axis minus diagonal; no per-design sign relabeling\",\n"
             << "  \"boundary_rule\": \"euclidean: (R-1)^2 < x^2+y^2 <= R^2; chebyshev: max(|x|,|y|)=R\",\n"
             << "  \"cross_radius_coupling\": \"same configuration and pivotal flag\",\n"
             << "  \"cross_design_coupling\": \"same seed/counter; identical fields within equal N and prefix-coupled across N\",\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"batch_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed " << options.designs.size() << " designs samples="
              << options.samples << " radii=";
    for (std::size_t index = 0; index < options.radii.size(); ++index) {
        if (index) std::cout << ',';
        std::cout << options.radii[index];
    }
    std::cout << " elapsed=" << elapsed << "\nwrote " << csv_path
              << "\nwrote " << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_matching(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
