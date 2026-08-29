// Minimal general-period adapter for the P253 fixed-root pivotal-H4 stream.
//
// Observable, Bernoulli counter domain, root toggle, landing registry and
// integer sufficient statistics are the P253 definitions.  Only the cyclic
// Gaussian quotient is replaced by the already-tested arbitrary integer-period
// quotient from threshold_rank_integer_period_mc.cpp.

// A design matrix is row-major and its columns are the lifted period vectors.

//     --design label,m00,m01,m10,m11


#define THRESHOLD_RANK_INTEGER_PERIOD_NO_MAIN
#include "threshold_rank_integer_period_mc.cpp"
#undef THRESHOLD_RANK_INTEGER_PERIOD_NO_MAIN

#include <cmath>
#include <functional>
#include <set>

namespace {

struct GPDesign {
    std::string label;
    Matrix periods;
};

struct GPOptions {
    std::uint64_t samples = 200000;
    int batches = 200;
    std::uint64_t seed = 26725360829ULL;
    std::uint64_t replica_offset = 26725300000ULL;
    int threads = 0;
    double p = 0.592746050790;
    std::string cutoff = "euclidean";
    std::vector<int> radii = {2};
    std::vector<GPDesign> designs = {
        {"n65_first", {8, -1, 1, 8}},
        {"n65_second", {7, -4, 4, 7}},
        {"n130_first", {11, -3, 3, 11}},
        {"n130_second", {9, -7, 7, 9}},
        {"n260_first", {16, -2, 2, 16}},
        {"n260_second", {14, -8, 8, 14}},
        {"n520_first", {22, -6, 6, 22}},
        {"n520_second", {18, -14, 14, 18}},
        {"n85_first", {9, -2, 2, 9}},
        {"n85_second", {7, -6, 6, 7}},
        {"n170_first", {13, -1, 1, 13}},
        {"n170_second", {11, -7, 7, 11}},
        {"n340_first", {18, -4, 4, 18}},
        {"n340_second", {14, -12, 12, 14}},
        {"n680_first", {26, -2, 2, 26}},
        {"n680_second", {22, -14, 14, 22}},
    };
    std::string git_commit = "unknown";
    std::string binary_sha256 = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
    bool validate_only = false;
};

std::vector<int> gp_parse_radii(const std::string& text) {
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

GPDesign gp_parse_design(const std::string& text) {
    std::stringstream stream(text);
    std::array<std::string, 5> fields;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (!std::getline(stream, fields[index], ',')) {
            throw std::invalid_argument(
                "--design requires label,m00,m01,m10,m11");
        }
    }
    std::string extra;
    if (std::getline(stream, extra)) {
        throw std::invalid_argument("--design has too many fields");
    }
    if (fields[0].empty() || fields[0].find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") !=
            std::string::npos) {
        throw std::invalid_argument("design label must be nonempty ASCII alnum/_/-");
    }
    Matrix periods{
        parse_number<Int>(fields[1], "--design m00"),
        parse_number<Int>(fields[2], "--design m01"),
        parse_number<Int>(fields[3], "--design m10"),
        parse_number<Int>(fields[4], "--design m11"),
    };
    (void)QuotientCoordinates(periods);
    return {fields[0], periods};
}

[[noreturn]] void gp_usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N --batches B --seed S --replica-offset K --p P\n"
        << "  --threads T --radii 2 --design label,m00,m01,m10,m11\n"
        << "  --cutoff euclidean|chebyshev\n"
        << "  --git-commit SHA --binary-sha256 SHA256 --output-prefix PATH\n"
        << "  --validate-only checks every local landing before sampling\n"
        << "  --self-test runs primitive cyclic and nonprimitive exact oracles\n";
    std::exit(status);
}

GPOptions gp_parse_options(int argc, char** argv) {
    GPOptions options;
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
        else if (argument == "--radii") options.radii = gp_parse_radii(need(index, argument));
        else if (argument == "--design") {
            if (!custom_designs) {
                options.designs.clear();
                custom_designs = true;
            }
            options.designs.push_back(gp_parse_design(need(index, argument)));
        } else if (argument == "--git-commit") options.git_commit = need(index, argument);
        else if (argument == "--binary-sha256") options.binary_sha256 = need(index, argument);
        else if (argument == "--output-prefix") options.output_prefix = need(index, argument);
        else if (argument == "--validate-only") options.validate_only = true;
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--help") gp_usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.self_test) return options;
    if (!options.validate_only && options.output_prefix.empty()) {
        throw std::invalid_argument("--output-prefix required");
    }
    if (options.designs.empty()) throw std::invalid_argument("at least one design required");
    std::set<std::string> labels;
    for (const auto& design : options.designs) {
        if (!labels.insert(design.label).second) {
            throw std::invalid_argument("design labels must be unique");
        }
    }
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (!(options.p > 0.0 && options.p < 1.0)) {
        throw std::invalid_argument("p must be in (0,1)");
    }
    if (options.cutoff != "euclidean" && options.cutoff != "chebyshev") {
        throw std::invalid_argument("cutoff must be euclidean or chebyshev");
    }
    if (options.replica_offset >
        std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("counter range overflows uint64");
    }
    return options;
}

class GPCrossClassifier {
  public:
    GPCrossClassifier(const Geometry& geometry, bool matching)
        : geometry_(geometry),
          edges_(matching ? geometry.matching_edges : geometry.primal_edges),
          union_find_(geometry.quotient) {}

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

struct GPPoint {
    int x;
    int y;
    int vertex;
    int boundary_mask;
    bool primal_root_neighbour;
    bool matching_root_neighbour;
};

bool gp_graph_adjacent(int x1, int y1, int x2, int y2, bool matching) {
    const int dx = std::abs(x2 - x1);
    const int dy = std::abs(y2 - y1);
    if (dx + dy == 1) return true;
    return matching && dx == 1 && dy == 1;
}

int gp_landing_sector(int x, int y) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    int sector = static_cast<int>(std::floor(
        (std::atan2(static_cast<double>(y), static_cast<double>(x)) + pi / 8.0) /
        (pi / 4.0)));
    return static_cast<int>(positive_mod(sector, 8));
}

class GPLanding {
  public:
    using Labeler = std::function<int(Vector)>;

    GPLanding(const Geometry& geometry, int radius, bool euclidean,
              Labeler labeler = Labeler()) {
        if (radius <= 0) throw std::invalid_argument("local radius must be positive");
        if (!labeler) {
            labeler = [&](Vector point) { return geometry.quotient.label(point); };
        }
        std::vector<int> seen(geometry.n, 0);
        for (int y = -radius; y <= radius; ++y) {
            for (int x = -radius; x <= radius; ++x) {
                if (x == 0 && y == 0) continue;
                const int norm2 = x * x + y * y;
                if (euclidean && norm2 > radius * radius) continue;
                const int vertex = labeler({x, y});
                if (seen[vertex]) {
                    throw std::invalid_argument(
                        "local annulus is not injective for N=" +
                        std::to_string(geometry.n) + ",R=" + std::to_string(radius));
                }
                seen[vertex] = 1;
                const bool boundary = euclidean
                    ? norm2 > (radius - 1) * (radius - 1)
                    : std::max(std::abs(x), std::abs(y)) == radius;
                if (boundary) {
                    const int sector = gp_landing_sector(x, y);
                    boundary_orbits_ |= (sector % 2 == 0) ? 1 : 2;
                }
                points_.push_back({
                    x, y, vertex, boundary ? 1 << gp_landing_sector(x, y) : 0,
                    gp_graph_adjacent(0, 0, x, y, false),
                    gp_graph_adjacent(0, 0, x, y, true),
                });
            }
        }
        primal_adjacency_.resize(points_.size());
        matching_adjacency_.resize(points_.size());
        for (int first = 0; first < static_cast<int>(points_.size()); ++first) {
            for (int second = first + 1; second < static_cast<int>(points_.size()); ++second) {
                if (gp_graph_adjacent(points_[first].x, points_[first].y,
                                      points_[second].x, points_[second].y, false)) {
                    primal_adjacency_[first].push_back(second);
                    primal_adjacency_[second].push_back(first);
                }
                if (gp_graph_adjacent(points_[first].x, points_[first].y,
                                      points_[second].x, points_[second].y, true)) {
                    matching_adjacency_[first].push_back(second);
                    matching_adjacency_[second].push_back(first);
                }
            }
        }
    }

    int boundary_orbits() const { return boundary_orbits_; }

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

    std::vector<GPPoint> points_;
    std::vector<std::vector<int>> primal_adjacency_;
    std::vector<std::vector<int>> matching_adjacency_;
    int boundary_orbits_ = 0;  // bit 0 axis, bit 1 diagonal
};

int gp_pivotal_flag(GPCrossClassifier& classifier,
                    std::vector<std::uint8_t>& active) {
    const bool original_root = active[0];
    active[0] = 0;
    const bool without = classifier.cross(active);
    active[0] = 1;
    const bool with_root = classifier.cross(active);
    active[0] = original_root;
    const int pivotal = static_cast<int>(with_root) - static_cast<int>(without);
    if (pivotal != 0 && pivotal != 1) {
        throw std::logic_error("nonmonotone crossing event");
    }
    return pivotal;
}

std::vector<std::uint8_t> gp_bernoulli_configuration(
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

std::uint64_t gp_field_digest(const std::vector<std::uint8_t>& active,
                              std::uint64_t replica) {
    std::uint64_t value = splitmix64(replica ^ 0x44cdd30c28b8f16dULL);
    for (std::size_t vertex = 0; vertex < active.size(); ++vertex) {
        if (active[vertex]) {
            value ^= splitmix64(static_cast<std::uint64_t>(vertex) +
                                0x9e3779b97f4a7c15ULL);
        }
    }
    return splitmix64(value ^ static_cast<std::uint64_t>(active.size()));
}

struct GPSample {
    int primal_pivotal = 0;
    int matching_pivotal = 0;
    std::vector<int> primal_h4;
    std::vector<int> matching_h4;
};

GPSample gp_evaluate(const Geometry& geometry, GPCrossClassifier& primal_classifier,
                     GPCrossClassifier& matching_classifier,
                     const std::vector<GPLanding>& landings,
                     std::vector<std::uint8_t>& black) {
    std::vector<std::uint8_t> white(geometry.n);
    for (int vertex = 0; vertex < geometry.n; ++vertex) white[vertex] = !black[vertex];
    GPSample value;
    value.primal_pivotal = gp_pivotal_flag(primal_classifier, black);
    value.matching_pivotal = gp_pivotal_flag(matching_classifier, white);
    const bool black_root = black[0];
    const bool white_root = white[0];
    black[0] = 0;
    white[0] = 0;
    value.primal_h4.reserve(landings.size());
    value.matching_h4.reserve(landings.size());
    for (const auto& landing : landings) {
        value.primal_h4.push_back(
            value.primal_pivotal * landing.h4(black, false));
        value.matching_h4.push_back(
            value.matching_pivotal * landing.h4(white, true));
    }
    black[0] = black_root;
    white[0] = white_root;
    return value;
}

bool gp_sample_equal(const GPSample& first, const GPSample& second) {
    return first.primal_pivotal == second.primal_pivotal &&
           first.matching_pivotal == second.matching_pivotal &&
           first.primal_h4 == second.primal_h4 &&
           first.matching_h4 == second.matching_h4;
}

struct GPRadiusStats {
    std::int64_t primal_h4 = 0;
    std::int64_t matching_h4 = 0;
    std::int64_t h4_plus = 0;
    std::int64_t h4_minus = 0;
};

struct GPBatchStats {
    std::uint64_t samples = 0;
    std::uint64_t primal_pivotal = 0;
    std::uint64_t matching_pivotal = 0;
    std::uint64_t common_field_digest = 0;
    std::vector<GPRadiusStats> radius;

    explicit GPBatchStats(std::size_t count = 0) : radius(count) {}

    void add(const GPSample& sample, std::uint64_t digest) {
        ++samples;
        primal_pivotal += sample.primal_pivotal;
        matching_pivotal += sample.matching_pivotal;
        common_field_digest ^= digest;
        for (std::size_t index = 0; index < radius.size(); ++index) {
            const int primal = sample.primal_h4[index];
            const int matching = sample.matching_h4[index];
            radius[index].primal_h4 += primal;
            radius[index].matching_h4 += matching;
            radius[index].h4_plus += primal + matching;
            radius[index].h4_minus += primal - matching;
        }
    }
};

// Independent tiny quotient labelling: breadth-first representatives and a
// direct adj(P) divisibility test.  It intentionally does not use the HNF
// label method exercised by production.
class DirectLabels {
  public:
    explicit DirectLabels(Matrix periods) : periods_(periods), order_(safe_abs(
            determinant(periods), "direct quotient determinant")) {
        representatives_.push_back({0, 0});
        for (std::size_t cursor = 0;
             cursor < representatives_.size() && representatives_.size() <
                 static_cast<std::size_t>(order_);
             ++cursor) {
            const Vector source = representatives_[cursor];
            for (const Vector step : std::array<Vector, 4>{{
                     {1, 0}, {-1, 0}, {0, 1}, {0, -1}}}) {
                const Vector candidate{source.x + step.x, source.y + step.y};
                if (find(candidate) < 0) representatives_.push_back(candidate);
            }
        }
        if (representatives_.size() != static_cast<std::size_t>(order_)) {
            throw std::logic_error("direct quotient BFS did not find every coset");
        }
    }

    int label(Vector point) const {
        const int value = find(point);
        if (value < 0) throw std::logic_error("direct quotient label not found");
        return value;
    }

    Vector representative(int value) const { return representatives_.at(value); }

  private:
    bool equivalent(Vector first, Vector second) const {
        const __int128 dx = static_cast<__int128>(first.x) - second.x;
        const __int128 dy = static_cast<__int128>(first.y) - second.y;
        const __int128 det = determinant(periods_);
        const __int128 numerator0 = static_cast<__int128>(periods_.d) * dx -
                                    static_cast<__int128>(periods_.b) * dy;
        const __int128 numerator1 = -static_cast<__int128>(periods_.c) * dx +
                                    static_cast<__int128>(periods_.a) * dy;
        return numerator0 % det == 0 && numerator1 % det == 0;
    }

    int find(Vector point) const {
        for (int index = 0; index < static_cast<int>(representatives_.size()); ++index) {
            if (equivalent(point, representatives_[index])) return index;
        }
        return -1;
    }

    Matrix periods_;
    Int order_;
    std::vector<Vector> representatives_;
};

Geometry gp_make_relabelled_geometry(Matrix periods,
                                     const GPLanding::Labeler& labeler,
                                     const std::vector<Vector>& representatives) {
    Geometry geometry(periods);
    if (representatives.size() != static_cast<std::size_t>(geometry.n)) {
        throw std::logic_error("relabelled geometry has wrong representative count");
    }
    geometry.primal_edges.clear();
    geometry.matching_edges.clear();
    const std::array<Vector, 4> steps = {{{1, 0}, {0, 1}, {1, 1}, {1, -1}}};
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const Vector source = representatives[vertex];
        for (std::size_t index = 0; index < steps.size(); ++index) {
            const Vector step = steps[index];
            const Edge edge{vertex, labeler({source.x + step.x, source.y + step.y}),
                            static_cast<int>(step.x), static_cast<int>(step.y)};
            if (index < 2) geometry.primal_edges.push_back(edge);
            geometry.matching_edges.push_back(edge);
        }
    }
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = make_incident(geometry.n, geometry.matching_edges);
    return geometry;
}

std::vector<std::uint8_t> gp_transport_configuration(
        const std::vector<std::uint8_t>& source, const QuotientCoordinates& source_labels,
        const GPLanding::Labeler& target_labeler, int target_n) {
    std::vector<std::uint8_t> target(target_n, 0);
    std::vector<int> seen(target_n, 0);
    for (int label = 0; label < static_cast<int>(source.size()); ++label) {
        const int target_label = target_labeler(source_labels.representative(label));
        if (seen[target_label]) throw std::logic_error("configuration transport not bijective");
        seen[target_label] = 1;
        target[target_label] = source[label];
    }
    return target;
}

std::array<std::uint64_t, 4> gp_compare_all_configurations(
        Geometry& reference, Geometry& comparison,
        const GPLanding::Labeler& comparison_labeler, int radius) {
    if (reference.n != comparison.n || reference.n >= 63) {
        throw std::logic_error("tiny exact comparison requires equal N<63");
    }
    GPCrossClassifier ref_primal(reference, false);
    GPCrossClassifier ref_matching(reference, true);
    GPCrossClassifier cmp_primal(comparison, false);
    GPCrossClassifier cmp_matching(comparison, true);
    std::vector<GPLanding> ref_landings;
    ref_landings.emplace_back(reference, radius, true);
    std::vector<GPLanding> cmp_landings;
    cmp_landings.emplace_back(comparison, radius, true, comparison_labeler);
    std::array<std::uint64_t, 4> totals{{0, 0, 0, 0}};
    const std::uint64_t stop = 1ULL << reference.n;
    for (std::uint64_t mask = 0; mask < stop; ++mask) {
        std::vector<std::uint8_t> first(reference.n, 0);
        for (int vertex = 0; vertex < reference.n; ++vertex) {
            first[vertex] = (mask >> vertex) & 1U;
        }
        auto second = gp_transport_configuration(
            first, reference.quotient, comparison_labeler, comparison.n);
        const GPSample first_value = gp_evaluate(
            reference, ref_primal, ref_matching, ref_landings, first);
        const GPSample second_value = gp_evaluate(
            comparison, cmp_primal, cmp_matching, cmp_landings, second);
        if (!gp_sample_equal(first_value, second_value)) {
            throw std::logic_error("tiny exact relabelling oracle mismatch");
        }
        totals[0] += first_value.primal_pivotal;
        totals[1] += first_value.matching_pivotal;
        totals[2] += first_value.primal_h4[0] != 0;
        totals[3] += first_value.matching_h4[0] != 0;
    }
    return totals;
}

void gp_self_test() {
    const std::array<Vector, 8> directions = {{
        {1, 0}, {1, 1}, {0, 1}, {-1, 1},
        {-1, 0}, {-1, -1}, {0, -1}, {1, -1},
    }};
    for (int sector = 0; sector < 8; ++sector) {
        if (gp_landing_sector(directions[sector].x, directions[sector].y) != sector) {
            throw std::logic_error("axis/diagonal landing-sector oracle failed");
        }
    }
    // Exact gate 83e98fc: either direction orbit alone aliases scalar and
    // spin 4 (rank one); axis plus diagonal has response [[1,1],[1,-1]],
    // determinant -2 and rank two.
    const std::array<std::array<int, 2>, 2> orbit_response = {{{{1, 1}}, {{1, -1}}}};
    const int orbit_determinant =
        orbit_response[0][0] * orbit_response[1][1] -
        orbit_response[0][1] * orbit_response[1][0];
    if (orbit_determinant != -2) {
        throw std::logic_error("axis/diagonal scalar-spin4 rank gate failed");
    }

    // Primitive cyclic P253 equivalence.  For P=[[3,-1],[1,3]], the legacy
    // cyclic label is 3*x+y mod 10.  Exhaust every configuration after the
    // exact physical-coordinate transport.
    const Matrix primitive_periods{3, -1, 1, 3};
    Geometry primitive_hnf = make_geometry(primitive_periods);
    const Geometry frozen_r2 = make_geometry({8, -1, 1, 8});
    if (GPLanding(frozen_r2, 2, true).boundary_orbits() != 3) {
        throw std::logic_error("R2 landing does not retain axis and diagonal orbits");
    }
    GPLanding::Labeler cyclic_label = [](Vector point) {
        return static_cast<int>(positive_mod(3 * point.x + point.y, 10));
    };
    std::vector<Vector> cyclic_representatives(10);
    for (int label = 0; label < 10; ++label) cyclic_representatives[label] = {7 * label, 0};
    Geometry primitive_cyclic = gp_make_relabelled_geometry(
        primitive_periods, cyclic_label, cyclic_representatives);
    const auto primitive_totals = gp_compare_all_configurations(
        primitive_hnf, primitive_cyclic, cyclic_label, 1);

    // Nonprimitive Smith (2,4) quotient.  Compare HNF production labelling to
    // an independently enumerated direct coset labelling for all 2^8 fields.
    const Matrix nonprimitive_periods{2, -2, 2, 2};
    Geometry nonprimitive_hnf = make_geometry(nonprimitive_periods);
    const DirectLabels direct_labels(nonprimitive_periods);
    GPLanding::Labeler direct_label = [&](Vector point) {
        return direct_labels.label(point);
    };
    std::vector<Vector> direct_representatives(nonprimitive_hnf.n);
    for (int label = 0; label < nonprimitive_hnf.n; ++label) {
        direct_representatives[label] = direct_labels.representative(label);
    }
    Geometry nonprimitive_direct = gp_make_relabelled_geometry(
        nonprimitive_periods, direct_label, direct_representatives);
    const auto nonprimitive_totals = gp_compare_all_configurations(
        nonprimitive_hnf, nonprimitive_direct, direct_label, 1);
    if (nonprimitive_totals[0] + nonprimitive_totals[1] == 0) {
        throw std::logic_error("nonprimitive oracle exercised no pivotal flags");
    }

    // The same nonprimitive lattice under a unimodular period-basis change.
    const Matrix transformed_periods{2, 0, 2, 4};
    Geometry transformed = make_geometry(transformed_periods);
    GPLanding::Labeler transformed_label = [&](Vector point) {
        return transformed.quotient.label(point);
    };
    const auto transformed_totals = gp_compare_all_configurations(
        nonprimitive_hnf, transformed, transformed_label, 1);
    if (transformed_totals != nonprimitive_totals) {
        throw std::logic_error("nonprimitive period-basis oracle totals differ");
    }

    std::cout << "self-test passed: primitive N=10 configurations=1024 pivotal="
              << primitive_totals[0] << ',' << primitive_totals[1]
              << "; nonprimitive N=8 Smith=(2,4) configurations=256 pivotal="
              << nonprimitive_totals[0] << ',' << nonprimitive_totals[1]
              << " h4_nonzero=" << nonprimitive_totals[2] << ','
              << nonprimitive_totals[3]
              << "; axis+diagonal scalar/spin4 rank=2\n";
}

std::vector<Geometry> gp_geometries(const GPOptions& options) {
    std::vector<Geometry> geometries;
    geometries.reserve(options.designs.size());
    for (const auto& design : options.designs) {
        geometries.push_back(make_geometry(design.periods));
    }
    return geometries;
}

void gp_validate_landings(const GPOptions& options,
                          const std::vector<Geometry>& geometries) {
    for (const Geometry& geometry : geometries) {
        for (const int radius : options.radii) {
            const GPLanding landing(
                geometry, radius, options.cutoff == "euclidean");
            if (radius == 2 && landing.boundary_orbits() != 3) {
                throw std::invalid_argument(
                    "R2 landing must retain both axis and diagonal C4 orbits");
            }
        }
    }
}

int run_gp(int argc, char** argv) {
    const GPOptions options = gp_parse_options(argc, argv);
    if (options.self_test) {
        gp_self_test();
        return 0;
    }
    const std::vector<Geometry> geometries = gp_geometries(options);
    gp_validate_landings(options, geometries);
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
    std::vector<std::vector<GPBatchStats>> output(
        options.batches,
        std::vector<GPBatchStats>(
            options.designs.size(), GPBatchStats(options.radii.size())));
    const auto started = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::size_t design = 0; design < geometries.size(); ++design) {
            GPCrossClassifier primal_classifier(geometries[design], false);
            GPCrossClassifier matching_classifier(geometries[design], true);
            std::vector<GPLanding> landings;
            landings.reserve(options.radii.size());
            for (const int radius : options.radii) {
                landings.emplace_back(
                    geometries[design], radius, options.cutoff == "euclidean");
            }
            for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                const std::uint64_t replica = first + offset;
                auto black = gp_bernoulli_configuration(
                    geometries[design].n, options.seed, replica, options.p);
                const std::uint64_t digest = gp_field_digest(black, replica);
                output[batch][design].add(gp_evaluate(
                    geometries[design], primal_classifier, matching_classifier,
                    landings, black), digest);
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
    csv << "label,n,m00,m01,m10,m11,smith1,smith2,radius,batch,"
           "counter_first,counter_last_exclusive,samples,common_field_digest,"
           "primal_pivotal,matching_pivotal,primal_h4,matching_h4,h4_plus,h4_minus\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset + per_batch * batch;
        for (std::size_t design = 0; design < options.designs.size(); ++design) {
            const auto& common = output[batch][design];
            const Matrix& matrix = options.designs[design].periods;
            const QuotientCoordinates& quotient = geometries[design].quotient;
            for (std::size_t index = 0; index < options.radii.size(); ++index) {
                const auto& row = common.radius[index];
                csv << options.designs[design].label << ',' << geometries[design].n << ','
                    << matrix.a << ',' << matrix.b << ',' << matrix.c << ',' << matrix.d << ','
                    << quotient.smith1 << ',' << quotient.smith2 << ','
                    << options.radii[index] << ',' << batch << ',' << first << ','
                    << first + per_batch << ',' << common.samples << ','
                    << common.common_field_digest << ',' << common.primal_pivotal << ','
                    << common.matching_pivotal << ',' << row.primal_h4 << ','
                    << row.matching_h4 << ',' << row.h4_plus << ',' << row.h4_minus << '\n';
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
             << "  \"schema\": \"matching-one/general-period-multiradius-pivotal/v1\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"binary_sha256\": \"" << json_escape(options.binary_sha256) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
#ifdef _OPENMP
             << "  \"openmp\": true,\n"
#else
             << "  \"openmp\": false,\n"
#endif
             << "  \"threads_requested\": " << options.threads << ",\n"
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
        const QuotientCoordinates& quotient = geometries[index].quotient;
        metadata << "{\"label\":\"" << json_escape(options.designs[index].label)
                 << "\",\"N\":" << geometries[index].n
                 << ",\"period_matrix\":" << matrix_json(options.designs[index].periods)
                 << ",\"HNF\":[[" << quotient.h11 << ',' << quotient.h12
                 << "],[0," << quotient.h22 << "]]"
                 << ",\"smith_invariants\":[" << quotient.smith1 << ','
                 << quotient.smith2 << "]}";
    }
    metadata << "],\n"
             << "  \"period_matrix_convention\": \"row-major matrix; columns are lifted period vectors\",\n"
             << "  \"quotient_coordinates\": \"column-HNF representatives; label=rx+h11*ry\",\n"
             << "  \"rng\": \"P253 counter-derived SplitMix64 Bernoulli field\",\n"
             << "  \"graph_semantics\": \"black NN primal; complemented white NN+NNN matching\",\n"
             << "  \"root_toggle\": \"vertex zero forced absent/present for each colour\",\n"
             << "  \"landing_registry\": \"absolute lattice axis minus diagonal; no per-design sign relabeling\",\n"
             << "  \"boundary_rule\": \"euclidean: (R-1)^2 < x^2+y^2 <= R^2; chebyshev: max(|x|,|y|)=R\",\n"
             << "  \"cross_radius_coupling\": \"same configuration and pivotal flag\",\n"
             << "  \"cross_design_coupling\": \"same seed/counter HNF-label field; prefix-coupled across N\",\n"
             << "  \"common_field_digest\": \"xor of per-replica full-bitfield hashes; equality is required within equal-N pairs\",\n"
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
        return run_gp(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
