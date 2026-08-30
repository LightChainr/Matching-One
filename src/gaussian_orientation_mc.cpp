// Discovery-scale same-N Gaussian-integer orientation tomography.
//
// Primitive (a,b) represents the exact square torus with periods
// (a,b),(-b,a), N=a^2+b^2 and cyclic vertex label j=a*x+b*y (mod N).
// Two representations of the same N use identical counter-keyed occupancies.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Edge {
    int i;
    int j;
    int dx;
    int dy;
};

struct Winding {
    std::int64_t x = 0;
    std::int64_t y = 0;
};

struct Channels {
    bool direction_0 = false;
    bool direction_1 = false;
    bool either = false;
    bool both = false;
    bool cross = false;
};

const std::array<const char*, 5> kChannelNames = {
    "cross", "both", "either", "direction_0", "direction_1",
};

bool channel_value(const Channels& channels, std::size_t index) {
    switch (index) {
        case 0: return channels.cross;
        case 1: return channels.both;
        case 2: return channels.either;
        case 3: return channels.direction_0;
        case 4: return channels.direction_1;
        default: throw std::logic_error("invalid channel index");
    }
}

struct Geometry {
    int n;
    int a;
    int b;
    double theta;
    double cos4;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
};

struct PairDesign {
    int n;
    int a1;
    int b1;
    int a2;
    int b2;
};

const std::vector<PairDesign> kDesigns = {
    {65, 8, 1, 7, 4},
    {85, 9, 2, 7, 6},
    {130, 11, 3, 9, 7},
    {145, 12, 1, 9, 8},
    {170, 13, 1, 11, 7},
    {185, 13, 4, 11, 8},
    {265, 16, 3, 12, 11},
    // Keep the prospective third doubling pair in Gaussian-lineage order.
    {290, 13, 11, 17, 1},
};

int positive_mod(int value, int modulus) {
    value %= modulus;
    return value < 0 ? value + modulus : value;
}

Winding primitive(Winding value) {
    const auto divisor = std::gcd(std::llabs(value.x), std::llabs(value.y));
    if (divisor == 0) return value;
    value.x /= divisor;
    value.y /= divisor;
    if (value.x < 0 || (value.x == 0 && value.y < 0)) {
        value.x = -value.x;
        value.y = -value.y;
    }
    return value;
}

void extend_basis(std::vector<Winding>& basis, Winding value) {
    if ((value.x == 0 && value.y == 0) || basis.size() == 2) return;
    value = primitive(value);
    if (basis.empty()) {
        basis.push_back(value);
        return;
    }
    const Winding first = basis.front();
    if (first.x * value.y != first.y * value.x) basis.push_back(value);
}

class HomologyUnionFind {
  public:
    HomologyUnionFind(int n, int a, int b)
        : n_(n), a_(a), b_(b), parent_(n), size_(n), delta_x_(n), delta_y_(n),
          basis_(n) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        for (auto& basis : basis_) basis.clear();
    }

    struct FindResult {
        int root;
        std::int64_t dx;
        std::int64_t dy;
    };

    FindResult find(int x) {
        if (parent_[x] == x) return {x, 0, 0};
        const int old_parent = parent_[x];
        const FindResult up = find(old_parent);
        delta_x_[x] += up.dx;
        delta_y_[x] += up.dy;
        parent_[x] = up.root;
        return {up.root, delta_x_[x], delta_y_[x]};
    }

    Winding period_coordinates(std::int64_t dx, std::int64_t dy) const {
        // Inverse of [[a,-b],[b,a]] applied exactly to (dx,dy).
        const std::int64_t num0 = static_cast<std::int64_t>(a_) * dx +
                                  static_cast<std::int64_t>(b_) * dy;
        const std::int64_t num1 = -static_cast<std::int64_t>(b_) * dx +
                                  static_cast<std::int64_t>(a_) * dy;
        if (num0 % n_ != 0 || num1 % n_ != 0) {
            throw std::logic_error("closed displacement is not in the Gaussian period lattice");
        }
        return {num0 / n_, num1 / n_};
    }

    void add_edge(const Edge& edge) {
        FindResult fi = find(edge.i);
        FindResult fj = find(edge.j);
        std::int64_t root_dx = fi.dx + edge.dx - fj.dx;
        std::int64_t root_dy = fi.dy + edge.dy - fj.dy;
        if (fi.root == fj.root) {
            extend_basis(basis_[fi.root], period_coordinates(root_dx, root_dy));
            return;
        }
        if (size_[fi.root] < size_[fj.root]) {
            std::swap(fi, fj);
            root_dx = -root_dx;
            root_dy = -root_dy;
        }
        parent_[fj.root] = fi.root;
        delta_x_[fj.root] = root_dx;
        delta_y_[fj.root] = root_dy;
        size_[fi.root] += size_[fj.root];
        for (const Winding winding : basis_[fj.root]) extend_basis(basis_[fi.root], winding);
        basis_[fj.root].clear();
    }

    const std::vector<Winding>& component_basis(int x) {
        return basis_[find(x).root];
    }

  private:
    int n_;
    int a_;
    int b_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<std::int64_t> delta_x_;
    std::vector<std::int64_t> delta_y_;
    std::vector<std::vector<Winding>> basis_;
};

Geometry make_geometry(int a, int b) {
    if (a <= 0 || b < 0 || std::gcd(a, b) != 1) {
        throw std::invalid_argument("Gaussian representation requires a>0, b>=0, gcd(a,b)=1");
    }
    const std::int64_t n64 = static_cast<std::int64_t>(a) * a +
                             static_cast<std::int64_t>(b) * b;
    if (n64 > std::numeric_limits<int>::max()) throw std::invalid_argument("N is too large");
    const int n = static_cast<int>(n64);
    const double theta = std::atan2(static_cast<double>(b), static_cast<double>(a));
    Geometry geometry{n, a, b, theta, std::cos(4.0 * theta), {}, {}};
    geometry.primal_edges.reserve(2 * n);
    geometry.matching_edges.reserve(4 * n);
    const std::vector<std::tuple<int, int, int>> steps = {
        {a, 1, 0}, {b, 0, 1}, {a + b, 1, 1}, {a - b, 1, -1},
    };
    for (int i = 0; i < n; ++i) {
        for (std::size_t index = 0; index < steps.size(); ++index) {
            const auto [residue, dx, dy] = steps[index];
            const Edge edge{i, positive_mod(i + residue, n), dx, dy};
            if (index < 2) geometry.primal_edges.push_back(edge);
            geometry.matching_edges.push_back(edge);
        }
    }
    return geometry;
}

Channels classify(const std::vector<std::uint8_t>& active, const std::vector<Edge>& edges,
                  HomologyUnionFind& union_find) {
    union_find.reset();
    for (const Edge& edge : edges) {
        if (active[edge.i] && active[edge.j]) union_find.add_edge(edge);
    }
    Channels output;
    for (int vertex = 0; vertex < static_cast<int>(active.size()); ++vertex) {
        if (!active[vertex]) continue;
        const std::vector<Winding>& basis = union_find.component_basis(vertex);
        output.cross = output.cross || basis.size() == 2;
        for (const Winding winding : basis) {
            output.direction_0 = output.direction_0 || winding.x != 0;
            output.direction_1 = output.direction_1 || winding.y != 0;
        }
    }
    output.either = output.direction_0 || output.direction_1;
    // Deliberately configuration-level: the two directions may occur in
    // distinct rank-one components.  Cross remains a rank-two component event.
    output.both = output.direction_0 && output.direction_1;
    return output;
}

// A deliberately independent lifted-coordinate traversal used only by --self-test.
Channels classify_bfs_reference(const Geometry& geometry,
                                const std::vector<std::uint8_t>& active,
                                const std::vector<Edge>& edges) {
    std::vector<std::vector<std::tuple<int, int, int>>> adjacency(geometry.n);
    for (const Edge& edge : edges) {
        adjacency[edge.i].push_back({edge.j, edge.dx, edge.dy});
        adjacency[edge.j].push_back({edge.i, -edge.dx, -edge.dy});
    }
    std::vector<std::uint8_t> visited(geometry.n, 0);
    std::vector<std::int64_t> lift_x(geometry.n, 0), lift_y(geometry.n, 0);
    Channels output;
    for (int start = 0; start < geometry.n; ++start) {
        if (!active[start] || visited[start]) continue;
        std::queue<int> pending;
        std::vector<Winding> basis;
        visited[start] = 1;
        pending.push(start);
        while (!pending.empty()) {
            const int i = pending.front();
            pending.pop();
            for (const auto [j, dx, dy] : adjacency[i]) {
                if (!active[j]) continue;
                const std::int64_t candidate_x = lift_x[i] + dx;
                const std::int64_t candidate_y = lift_y[i] + dy;
                if (!visited[j]) {
                    visited[j] = 1;
                    lift_x[j] = candidate_x;
                    lift_y[j] = candidate_y;
                    pending.push(j);
                } else {
                    const std::int64_t cycle_x = candidate_x - lift_x[j];
                    const std::int64_t cycle_y = candidate_y - lift_y[j];
                    const std::int64_t num0 = static_cast<std::int64_t>(geometry.a) * cycle_x +
                                              static_cast<std::int64_t>(geometry.b) * cycle_y;
                    const std::int64_t num1 = -static_cast<std::int64_t>(geometry.b) * cycle_x +
                                              static_cast<std::int64_t>(geometry.a) * cycle_y;
                    if (num0 % geometry.n || num1 % geometry.n) {
                        throw std::logic_error("BFS reference found a nonperiod cycle");
                    }
                    extend_basis(basis, {num0 / geometry.n, num1 / geometry.n});
                }
            }
        }
        output.cross = output.cross || basis.size() == 2;
        for (const Winding winding : basis) {
            output.direction_0 = output.direction_0 || winding.x != 0;
            output.direction_1 = output.direction_1 || winding.y != 0;
        }
    }
    output.either = output.direction_0 || output.direction_1;
    output.both = output.direction_0 && output.direction_1;
    return output;
}

int component_count(const std::vector<std::uint8_t>& active, HomologyUnionFind& union_find) {
    int count = 0;
    for (int vertex = 0; vertex < static_cast<int>(active.size()); ++vertex) {
        if (active[vertex] && union_find.find(vertex).root == vertex) ++count;
    }
    return count;
}

struct SiteMotifs {
    int V = 0;
    int E = 0;
    int F0 = 0;
    int nnn_pos = 0;
    int nnn_neg = 0;
    int path3_x = 0;
    int path3_y = 0;
    int corners = 0;
    int right_angle = 0;
};

SiteMotifs count_site_motifs(const std::vector<std::uint8_t>& active, const Geometry& geometry) {
    SiteMotifs motifs;
    const int n = geometry.n;
    const int a = geometry.a;
    const int b = geometry.b;
    for (int i = 0; i < n; ++i) motifs.V += active[i];
    for (const Edge& edge : geometry.primal_edges) {
        if (active[edge.i] && active[edge.j]) ++motifs.E;
    }
    for (int i = 0; i < n; ++i) {
        const int ix = positive_mod(i + a, n);
        const int iy = positive_mod(i + b, n);
        const int ixy = positive_mod(i + a + b, n);
        const int ixmy = positive_mod(i + (a - b), n);
        const int ixx = positive_mod(i + a + a, n);
        const int iyy = positive_mod(i + b + b, n);
        const int s0 = active[i];
        const int s1 = active[ix];
        const int s2 = active[iy];
        const int s3 = active[ixy];
        if (s0 && s1 && s2 && s3) ++motifs.F0;
        if (s0 && active[ixy]) ++motifs.nnn_pos;
        if (s0 && active[ixmy]) ++motifs.nnn_neg;
        if (i != ix && i != ixx && ix != ixx && s0 && s1 && active[ixx]) ++motifs.path3_x;
        if (i != iy && i != iyy && iy != iyy && s0 && s2 && active[iyy]) ++motifs.path3_y;
        motifs.corners += (s0 && s1 && s2) + (s0 && s1 && s3) + (s0 && s2 && s3) + (s1 && s2 && s3);
        // The fixed-K Issue #40 oracle declares exactly the translated
        // {i,i+x,i+y} family, with multiplicity N (not all four corners).
        motifs.right_angle += s0 && s1 && s2;
    }
    return motifs;
}

constexpr int kEulerObs = 20;
double falling_ratio(int occupied, int n, int order) {
    if (order <= 0 || occupied < order || n < order) return 0.0;
    double value = 1.0;
    for (int i = 0; i < order; ++i) {
        value *= static_cast<double>(occupied - i) / static_cast<double>(n - i);
    }
    return value;
}

const std::array<const char*, kEulerObs> kEulerNames = {
    "q", "C_black", "C_white", "V", "E", "F0",
    "nnn_pos", "nnn_neg", "path3_x", "path3_y", "corners",
    "E_mc", "F0_mc", "nnn_pos_mc", "nnn_neg_mc",
    "path3_x_mc", "path3_y_mc", "corners_mc",
    "right_angle", "right_angle_mc",
};

struct EulerGeometryAccum {
    std::array<double, kEulerObs> sum{};
    std::array<std::array<double, kEulerObs>, kEulerObs> gram{};
    double wrapping_l1 = 0.0;
    double identity_l1 = 0.0;
};

struct EulerCrossAccum {
    // first_second_gram[i][j] = sum_r first_r[i] * second_r[j].
    std::array<std::array<double, kEulerObs>, kEulerObs> first_second_gram{};
};

void accumulate_euler(EulerGeometryAccum& acc, const std::array<double, kEulerObs>& values,
                      double wrapping_l1, double identity_l1) {
    for (int i = 0; i < kEulerObs; ++i) {
        acc.sum[i] += values[i];
        for (int j = 0; j < kEulerObs; ++j) acc.gram[i][j] += values[i] * values[j];
    }
    acc.wrapping_l1 += wrapping_l1;
    acc.identity_l1 += identity_l1;
}

void accumulate_euler_cross(EulerCrossAccum& acc,
                            const std::array<double, kEulerObs>& first,
                            const std::array<double, kEulerObs>& second) {
    for (int i = 0; i < kEulerObs; ++i) {
        for (int j = 0; j < kEulerObs; ++j) {
            acc.first_second_gram[i][j] += first[i] * second[j];
        }
    }
}

std::array<double, kEulerObs> pack_euler(int q, int c_black, int c_white,
                                        const SiteMotifs& motifs, int n) {
    const double k2 = falling_ratio(motifs.V, n, 2);
    const double k3 = falling_ratio(motifs.V, n, 3);
    const double k4 = falling_ratio(motifs.V, n, 4);
    return {
        static_cast<double>(q),
        static_cast<double>(c_black),
        static_cast<double>(c_white),
        static_cast<double>(motifs.V),
        static_cast<double>(motifs.E),
        static_cast<double>(motifs.F0),
        static_cast<double>(motifs.nnn_pos),
        static_cast<double>(motifs.nnn_neg),
        static_cast<double>(motifs.path3_x),
        static_cast<double>(motifs.path3_y),
        static_cast<double>(motifs.corners),
        static_cast<double>(motifs.E) - 2.0 * n * k2,
        static_cast<double>(motifs.F0) - n * k4,
        static_cast<double>(motifs.nnn_pos) - n * k2,
        static_cast<double>(motifs.nnn_neg) - n * k2,
        static_cast<double>(motifs.path3_x) - n * k3,
        static_cast<double>(motifs.path3_y) - n * k3,
        static_cast<double>(motifs.corners) - 4.0 * n * k3,
        static_cast<double>(motifs.right_angle),
        static_cast<double>(motifs.right_angle) - n * k3,
    };
}

double wrapping_spread(const Channels& primal, const Channels& matching) {
    const int q0 = static_cast<int>(channel_value(primal, 0)) -
                   static_cast<int>(channel_value(matching, 0));
    double spread = 0.0;
    for (std::size_t channel = 1; channel < kChannelNames.size(); ++channel) {
        const int q = static_cast<int>(channel_value(primal, channel)) -
                      static_cast<int>(channel_value(matching, channel));
        spread += std::abs(q - q0);
    }
    return spread;
}


std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

double counter_uniform(std::uint64_t seed, int n, std::uint64_t replica, int site) {
    // N separates designs; the two representations of one N share every U_j.
    const std::uint64_t key = seed ^ splitmix64(static_cast<std::uint64_t>(n)) ^
                              splitmix64(replica + 0xd1b54a32d192ed03ULL) ^
                              splitmix64(static_cast<std::uint64_t>(site) +
                                         0x94d049bb133111ebULL);
    return static_cast<double>(splitmix64(key) >> 11) * 0x1.0p-53;
}

void self_test() {
    for (const auto [a, b] : std::vector<std::pair<int, int>>{{2, 1}, {3, 2}}) {
        const Geometry geometry = make_geometry(a, b);
        if (positive_mod(a * a + b * b, geometry.n) != 0 ||
            positive_mod(a * (-b) + b * a, geometry.n) != 0) {
            throw std::runtime_error("cyclic label period regression failed");
        }
        HomologyUnionFind primal(geometry.n, a, b), matching(geometry.n, a, b);
        const std::uint64_t configurations = 1ULL << geometry.n;
        std::vector<std::uint8_t> active(geometry.n), complement(geometry.n);
        for (std::uint64_t mask = 0; mask < configurations; ++mask) {
            for (int site = 0; site < geometry.n; ++site) {
                active[site] = (mask >> site) & 1U;
                complement[site] = !active[site];
            }
            const Channels p = classify(active, geometry.primal_edges, primal);
            const Channels p_ref = classify_bfs_reference(geometry, active, geometry.primal_edges);
            const Channels m = classify(active, geometry.matching_edges, matching);
            const Channels m_ref = classify_bfs_reference(geometry, active, geometry.matching_edges);
            if (p.direction_0 != p_ref.direction_0 || p.direction_1 != p_ref.direction_1 ||
                p.either != p_ref.either || p.both != p_ref.both || p.cross != p_ref.cross ||
                m.direction_0 != m_ref.direction_0 || m.direction_1 != m_ref.direction_1 ||
                m.either != m_ref.either || m.both != m_ref.both || m.cross != m_ref.cross) {
                throw std::runtime_error("union-find/BFS exact configuration mismatch");
            }
            const Channels white_matching = classify(complement, geometry.matching_edges, matching);
            const int reference_difference = static_cast<int>(channel_value(p, 0)) -
                                             static_cast<int>(channel_value(white_matching, 0));
            for (std::size_t channel = 1; channel < kChannelNames.size(); ++channel) {
                const int difference = static_cast<int>(channel_value(p, channel)) -
                                       static_cast<int>(channel_value(white_matching, channel));
                if (difference != reference_difference) {
                    throw std::runtime_error("matching-channel configuration identity failed");
                }
            }
            const int c_black = component_count(active, primal);
            const int c_white = component_count(complement, matching);
            const SiteMotifs motifs = count_site_motifs(active, geometry);
            const int residual = (c_black - c_white) - (
                reference_difference + motifs.V - motifs.E + motifs.F0);
            if (residual != 0) {
                throw std::runtime_error("Euler cluster identity failed");
            }
        }
    }
    const double value = counter_uniform(17, 65, 23, 5);
    (void)counter_uniform(17, 85, 999, 7);
    if (value != counter_uniform(17, 65, 23, 5) || !(value >= 0.0 && value < 1.0)) {
        throw std::runtime_error("counter RNG regression failed");
    }
    std::cout << "self-test passed: exhaustive N=5,13 union-find/BFS and matching channels; "
                 "Euler identity; cyclic labels; counter RNG\n";
}

struct Options {
    std::uint64_t samples = 200000;
    int batches = 40;
    double p_ref = 0.592746050790;
    std::uint64_t seed = 20260828;
    int threads = 0;
    int only_n = 0;
    std::uint64_t replica_offset = 0;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
    bool euler_motifs = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N          replicas per same-N pair (default 200000)\n"
        << "  --batches B          equal batches (default 40)\n"
        << "  --p-ref P            frozen site probability (default 0.592746050790)\n"
        << "  --seed S             unsigned 64-bit seed (default 20260828)\n"
        << "  --replica-offset K   first RNG replica counter (default 0)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --n N                run only one frozen confirmation size (default all)\n"
        << "  --git-commit SHA     provenance string recorded in metadata\n"
        << "  --output-prefix PATH writes PATH.batches.csv and PATH.metadata.json\n"
        << "  --self-test           exhaustive reference checks and exit\n"
        << "  --euler-motifs        also write Euler/motif batch moments JSONL\n"
        << "  --help                show this help\n";
    std::exit(status);
}

template <typename T>
T parse_number(const std::string& text, const std::string& option) {
    std::istringstream input(text);
    T value{};
    input >> value;
    if (!input || !input.eof()) throw std::invalid_argument("invalid value for " + option);
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (++i >= argc) usage(argv[0], 2);
            return argv[i];
        };
        if (arg == "--samples") options.samples = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(next(), arg);
        else if (arg == "--p-ref") options.p_ref = parse_number<double>(next(), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(next(), arg);
        else if (arg == "--n") options.only_n = parse_number<int>(next(), arg);
        else if (arg == "--git-commit") options.git_commit = next();
        else if (arg == "--output-prefix") options.output_prefix = next();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--euler-motifs") options.euler_motifs = true;
        else if (arg == "--help") usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix is required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be positive and divisible by batches>=2");
    }
    if (!(options.p_ref > 0.0 && options.p_ref < 1.0)) {
        throw std::invalid_argument("p-ref must lie strictly between zero and one");
    }
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (options.only_n != 0 && std::none_of(kDesigns.begin(), kDesigns.end(),
            [&](const PairDesign& design) { return design.n == options.only_n; })) {
        throw std::invalid_argument(
            "--n must be one of 65, 85, 130, 145, 170, 185, 265, 290");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("replica counter range overflows uint64");
    }
    return options;
}

struct BatchCounts {
    std::uint64_t samples = 0;
    // estimator order: first primal, first matching, second primal, second matching.
    std::array<std::array<std::uint64_t, 5>, 4> sums{};
};

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        if (ch == '\\' || ch == '"') out << '\\' << ch;
        else if (ch == '\n') out << "\\n";
        else out << ch;
    }
    return out.str();
}

std::string utc_now() {
    const std::time_t value = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &value);
#else
    gmtime_r(&value, &tm);
#endif
    std::ostringstream out;
    out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

void write_euler_geometry(std::ostream& out, const char* label, int a, int b,
                         const EulerGeometryAccum& acc) {
    out << "\"" << label << "\": {\"a\": " << a << ", \"b\": " << b << ", \"sum\": [";
    for (int i = 0; i < kEulerObs; ++i) {
        if (i) out << ", ";
        out << std::setprecision(17) << acc.sum[i];
    }
    out << "], \"gram\": [";
    for (int i = 0; i < kEulerObs; ++i) {
        if (i) out << ", ";
        out << "[";
        for (int j = 0; j < kEulerObs; ++j) {
            if (j) out << ", ";
            out << std::setprecision(17) << acc.gram[i][j];
        }
        out << "]";
    }
    out << "], \"wrapping_l1\": " << std::setprecision(17) << acc.wrapping_l1
        << ", \"identity_l1\": " << acc.identity_l1 << "}";
}

std::array<double, kEulerObs> fill_euler(EulerGeometryAccum& acc, const Geometry& geometry,
                const std::vector<std::uint8_t>& black,
                const std::vector<std::uint8_t>& white,
                const Channels& primal, const Channels& matching,
                HomologyUnionFind& primal_uf, HomologyUnionFind& matching_uf) {
    const int q = static_cast<int>(primal.either) - static_cast<int>(matching.either);
    const int c_black = component_count(black, primal_uf);
    const int c_white = component_count(white, matching_uf);
    const SiteMotifs motifs = count_site_motifs(black, geometry);
    const int residual = (c_black - c_white) - (q + motifs.V - motifs.E + motifs.F0);
    const auto values = pack_euler(q, c_black, c_white, motifs, geometry.n);
    accumulate_euler(acc, values, wrapping_spread(primal, matching), std::abs(residual));
    return values;
}

void run_design(const PairDesign& design, const Options& options,
                std::ofstream& batches_file, std::ofstream* motif_file) {
    const Geometry first = make_geometry(design.a1, design.b1);
    const Geometry second = make_geometry(design.a2, design.b2);
    if (first.n != design.n || second.n != design.n) {
        throw std::logic_error("prescribed same-N design is inconsistent");
    }
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<BatchCounts> counts(options.batches);
    std::vector<EulerGeometryAccum> euler_first(options.euler_motifs ? options.batches : 0);
    std::vector<EulerGeometryAccum> euler_second(options.euler_motifs ? options.batches : 0);
    std::vector<EulerCrossAccum> euler_cross(options.euler_motifs ? options.batches : 0);

#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        BatchCounts local;
        local.samples = per_batch;
        std::vector<std::uint8_t> black(design.n), white(design.n);
        HomologyUnionFind f_primal(design.n, first.a, first.b);
        HomologyUnionFind f_matching(design.n, first.a, first.b);
        HomologyUnionFind s_primal(design.n, second.a, second.b);
        HomologyUnionFind s_matching(design.n, second.a, second.b);
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            for (int site = 0; site < design.n; ++site) {
                black[site] = counter_uniform(options.seed, design.n, replica, site) < options.p_ref;
                white[site] = !black[site];
            }
            const Channels fp = classify(black, first.primal_edges, f_primal);
            const Channels fm = classify(white, first.matching_edges, f_matching);
            const Channels sp = classify(black, second.primal_edges, s_primal);
            const Channels sm = classify(white, second.matching_edges, s_matching);
            for (std::size_t channel = 0; channel < kChannelNames.size(); ++channel) {
                local.sums[0][channel] += channel_value(fp, channel);
                local.sums[1][channel] += channel_value(fm, channel);
                local.sums[2][channel] += channel_value(sp, channel);
                local.sums[3][channel] += channel_value(sm, channel);
            }
            if (options.euler_motifs) {
                const auto first_values = fill_euler(
                    euler_first[batch], first, black, white, fp, fm, f_primal, f_matching);
                const auto second_values = fill_euler(
                    euler_second[batch], second, black, white, sp, sm, s_primal, s_matching);
                accumulate_euler_cross(euler_cross[batch], first_values, second_values);
            }
        }
        counts[batch] = local;
    }

    auto write_row = [&](int batch, std::size_t channel) {
        const BatchCounts& row = counts[batch];
        batches_file << design.n << ',' << batch << ',' << counts[batch].samples << ','
                     << std::setprecision(17) << options.p_ref << ',' << kChannelNames[channel] << ','
                     << design.a1 << ',' << design.b1 << ',' << design.a2 << ',' << design.b2
                     << ',' << row.sums[0][channel] << ',' << row.sums[1][channel]
                     << ',' << row.sums[2][channel] << ',' << row.sums[3][channel] << '\n';
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        for (std::size_t channel = 0; channel < kChannelNames.size(); ++channel) {
            write_row(batch, channel);
        }
    }
    if (motif_file != nullptr) {
        const std::uint64_t per_batch = options.samples / options.batches;
        for (int batch = 0; batch < options.batches; ++batch) {
            *motif_file << "{\"n\": " << design.n << ", \"batch\": " << batch
                        << ", \"samples\": " << per_batch
                        << ", \"p_ref\": " << std::setprecision(17) << options.p_ref
                        << ", \"names\": [";
            for (int i = 0; i < kEulerObs; ++i) {
                if (i) *motif_file << ", ";
                *motif_file << "\"" << kEulerNames[i] << "\"";
            }
            *motif_file << "], ";
            write_euler_geometry(*motif_file, "first", design.a1, design.b1, euler_first[batch]);
            *motif_file << ", ";
            write_euler_geometry(*motif_file, "second", design.a2, design.b2, euler_second[batch]);
            *motif_file << ", \"cross_gram_semantics\": \"sum first[i]*second[j] over same replicas\", "
                        << "\"cross_gram\": [";
            for (int i = 0; i < kEulerObs; ++i) {
                if (i) *motif_file << ", ";
                *motif_file << "[";
                for (int j = 0; j < kEulerObs; ++j) {
                    if (j) *motif_file << ", ";
                    *motif_file << std::setprecision(17)
                                << euler_cross[batch].first_second_gram[i][j];
                }
                *motif_file << "]";
            }
            *motif_file << "]";
            *motif_file << "}\n";
        }
    }
    std::cout << "completed N=" << design.n << " pair (" << design.a1 << ',' << design.b1
              << ")/ (" << design.a2 << ',' << design.b2 << ") samples=" << options.samples
              << '\n';
}

int run(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    if (options.self_test) {
        self_test();
        return 0;
    }
    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path batch_path = options.output_prefix.string() + ".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream batches_file(batch_path);
    if (!batches_file) throw std::runtime_error("cannot open batch output");
    batches_file << "n,batch,samples,p_ref,channel,a1,b1,a2,b2,first_primal_sum,"
                    "first_matching_sum,second_primal_sum,second_matching_sum\n";
    const std::filesystem::path motif_path = options.output_prefix.string() + ".motifs.jsonl";
    std::ofstream motif_file;
    if (options.euler_motifs) {
        motif_file.open(motif_path);
        if (!motif_file) throw std::runtime_error("cannot open motif output");
    }
    const auto started = std::chrono::steady_clock::now();
    for (const PairDesign& design : kDesigns) {
        if (options.only_n == 0 || options.only_n == design.n) {
            run_design(design, options, batches_file,
                       options.euler_motifs ? &motif_file : nullptr);
        }
    }
    batches_file.close();
    if (options.euler_motifs) motif_file.close();

    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::ostringstream command;
    for (int i = 0; i < argc; ++i) {
        if (i) command << ' ';
        command << argv[i];
    }
    std::ofstream metadata(metadata_path);
    if (!metadata) throw std::runtime_error("cannot open metadata output");
    metadata << "{\n"
             << "  \"engine\": \"same-N Gaussian site-orientation discovery\",\n"
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
             << "  \"samples_per_pair\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"p_ref\": " << std::setprecision(17) << options.p_ref << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"rng\": \"stateless SplitMix64-derived mapping (seed,N,replica,cyclic-site)\",\n"
             << "  \"coupling\": \"same cyclic labels j in Z/NZ share U_j across representations\",\n"
             << "  \"channels\": [\"cross\", \"both\", \"either\", "
                "\"direction_0\", \"direction_1\"],\n"
             << "  \"euler_motifs\": " << (options.euler_motifs ? "true" : "false") << ",\n"
             << "  \"euler_observable_names\": [\"q\", \"C_black\", \"C_white\", \"V\", \"E\", \"F0\", "
                "\"nnn_pos\", \"nnn_neg\", \"path3_x\", \"path3_y\", \"corners\", "
                "\"E_mc\", \"F0_mc\", \"nnn_pos_mc\", \"nnn_neg_mc\", "
                "\"path3_x_mc\", \"path3_y_mc\", \"corners_mc\", "
                "\"right_angle\", \"right_angle_mc\"],\n"
             << "  \"cross_geometry_joint_gram\": "
             << (options.euler_motifs ? "true" : "false") << ",\n"
             << "  \"elapsed_seconds\": " << seconds << ",\n"
             << "  \"designs\": [\n";
    bool first_row = true;
    for (const PairDesign& design : kDesigns) {
        if (options.only_n != 0 && options.only_n != design.n) continue;
        const Geometry first = make_geometry(design.a1, design.b1);
        const Geometry second = make_geometry(design.a2, design.b2);
        if (!first_row) metadata << ",\n";
        first_row = false;
        metadata << "    {\"N\": " << design.n << ", \"first\": [" << design.a1 << ','
                 << design.b1 << "], \"second\": [" << design.a2 << ',' << design.b2
                 << "], \"theta_first\": " << first.theta << ", \"theta_second\": "
                 << second.theta << ", \"cos4_first\": " << first.cos4
                 << ", \"cos4_second\": " << second.cos4 << "}";
    }
    metadata << "\n  ]\n}\n";
    std::cout << "wrote " << batch_path << "\nwrote " << metadata_path;
    if (options.euler_motifs) std::cout << "\nwrote " << motif_path;
    std::cout << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
