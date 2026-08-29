// Fixed-root, landing-marked pivotal H4 pilot for the frozen N=65 pair.
//
// The unmarked pivotal total is retained only as a Russo control.  The new
// information is the axis-minus-diagonal landing mark, which is not present
// in threshold-rank histograms.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Edge { int i, j, dx, dy; };
struct Winding { std::int64_t x = 0, y = 0; };

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

class HomologyUnionFind {
  public:
    HomologyUnionFind(int n, int a, int b)
        : n_(n), a_(a), b_(b), parent_(n), size_(n), dx_(n), dy_(n),
          rank_(n), basis_(n) { reset(); }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(dx_.begin(), dx_.end(), 0);
        std::fill(dy_.begin(), dy_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
    }

    struct Find { int root; std::int64_t dx, dy; };

    Find find(int vertex) {
        if (parent_[vertex] == vertex) return {vertex, 0, 0};
        const int old = parent_[vertex];
        const Find above = find(old);
        dx_[vertex] += above.dx;
        dy_[vertex] += above.dy;
        parent_[vertex] = above.root;
        return {above.root, dx_[vertex], dy_[vertex]};
    }

    Winding winding(std::int64_t dx, std::int64_t dy) const {
        const std::int64_t u = static_cast<std::int64_t>(a_) * dx
                             + static_cast<std::int64_t>(b_) * dy;
        const std::int64_t v = -static_cast<std::int64_t>(b_) * dx
                             + static_cast<std::int64_t>(a_) * dy;
        if (u % n_ || v % n_) throw std::logic_error("cycle outside period lattice");
        return {u / n_, v / n_};
    }

    void extend(int root, Winding value) {
        if ((value.x == 0 && value.y == 0) || rank_[root] == 2) return;
        value = primitive(value);
        if (rank_[root] == 0) {
            basis_[root][0] = value;
            rank_[root] = 1;
        } else {
            const auto first = basis_[root][0];
            if (first.x * value.y != first.y * value.x) {
                basis_[root][1] = value;
                rank_[root] = 2;
            }
        }
    }

    void add_edge(const Edge& edge) {
        Find left = find(edge.i), right = find(edge.j);
        std::int64_t root_dx = left.dx + edge.dx - right.dx;
        std::int64_t root_dy = left.dy + edge.dy - right.dy;
        if (left.root == right.root) {
            extend(left.root, winding(root_dx, root_dy));
            return;
        }
        if (size_[left.root] < size_[right.root]) {
            std::swap(left, right);
            root_dx = -root_dx;
            root_dy = -root_dy;
        }
        parent_[right.root] = left.root;
        dx_[right.root] = root_dx;
        dy_[right.root] = root_dy;
        size_[left.root] += size_[right.root];
        for (int k = 0; k < rank_[right.root]; ++k) extend(left.root, basis_[right.root][k]);
        rank_[right.root] = 0;
    }

    bool crosses(int vertex) { return rank_[find(vertex).root] == 2; }

  private:
    int n_, a_, b_;
    std::vector<int> parent_, size_;
    std::vector<std::int64_t> dx_, dy_;
    std::vector<int> rank_;
    std::vector<std::array<Winding, 2>> basis_;
};

struct Geometry {
    int n, a, b;
    std::string orientation;
    std::vector<Edge> primal, matching;

    int vertex(int x, int y) const { return positive_mod(a * x + b * y, n); }
};

Geometry make_geometry(int a, int b, const std::string& orientation) {
    if (a <= 0 || b < 0 || std::gcd(a, b) != 1) throw std::invalid_argument("bad Gaussian pair");
    Geometry out{a * a + b * b, a, b, orientation, {}, {}};
    const std::array<std::tuple<int, int, int>, 4> steps = {{
        {a, 1, 0}, {b, 0, 1}, {a + b, 1, 1}, {a - b, 1, -1},
    }};
    for (int vertex = 0; vertex < out.n; ++vertex) {
        for (int k = 0; k < 4; ++k) {
            const auto [residue, dx, dy] = steps[k];
            const Edge edge{vertex, positive_mod(vertex + residue, out.n), dx, dy};
            if (k < 2) out.primal.push_back(edge);
            out.matching.push_back(edge);
        }
    }
    return out;
}

bool crosses(const Geometry& geometry, const std::vector<std::uint8_t>& active,
             bool matching) {
    HomologyUnionFind uf(geometry.n, geometry.a, geometry.b);
    const auto& edges = matching ? geometry.matching : geometry.primal;
    for (const auto& edge : edges) {
        if (active[edge.i] && active[edge.j]) uf.add_edge(edge);
    }
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        if (active[vertex] && uf.crosses(vertex)) return true;
    }
    return false;
}

struct Landing { int axis = 0, diagonal = 0, both = 0, landed = 0, h4 = 0; };

const std::array<std::pair<int, int>, 4> kPrimalSteps = {{
    {1, 0}, {-1, 0}, {0, 1}, {0, -1},
}};
const std::array<std::pair<int, int>, 8> kMatchingSteps = {{
    {1, 0}, {-1, 0}, {0, 1}, {0, -1},
    {1, 1}, {1, -1}, {-1, 1}, {-1, -1},
}};

int sector(int x, int y, int shift) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    int base = static_cast<int>(std::floor((std::atan2(y, x) + pi / 8) / (pi / 4)));
    return positive_mod(base - shift, 8);
}

std::vector<int> components(const Geometry& geometry,
                            const std::vector<std::uint8_t>& active,
                            int radius, bool matching, bool enabled_value,
                            int registry_shift) {
    std::vector<std::pair<int, int>> points;
    std::map<std::pair<int, int>, int> index;
    std::vector<int> seen_vertex(geometry.n, 0);
    for (int y = -radius; y <= radius; ++y) {
        for (int x = -radius; x <= radius; ++x) {
            if (x == 0 && y == 0) continue;
            const int id = static_cast<int>(points.size());
            points.push_back({x, y});
            index[{x, y}] = id;
            const int vertex = geometry.vertex(x, y);
            if (seen_vertex[vertex]) throw std::logic_error("annulus is not quotient-injective");
            seen_vertex[vertex] = 1;
        }
    }
    std::vector<std::uint8_t> enabled(points.size(), 0), visited(points.size(), 0);
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        const auto [x, y] = points[i];
        enabled[i] = (static_cast<bool>(active[geometry.vertex(x, y)]) == enabled_value);
    }
    std::vector<int> masks;
    for (int start = 0; start < static_cast<int>(points.size()); ++start) {
        if (!enabled[start] || visited[start]) continue;
        std::vector<int> stack{start};
        visited[start] = 1;
        bool inner = false;
        int mask = 0;
        while (!stack.empty()) {
            const int current = stack.back(); stack.pop_back();
            const auto [x, y] = points[current];
            if ((!matching && std::abs(x) + std::abs(y) == 1)
                || (matching && std::max(std::abs(x), std::abs(y)) == 1)) inner = true;
            if (std::max(std::abs(x), std::abs(y)) == radius) {
                mask |= 1 << sector(x, y, registry_shift);
            }
            const int count = matching ? 8 : 4;
            for (int k = 0; k < count; ++k) {
                const auto step = matching ? kMatchingSteps[k] : kPrimalSteps[k];
                const auto found = index.find({x + step.first, y + step.second});
                if (found != index.end() && enabled[found->second] && !visited[found->second]) {
                    visited[found->second] = 1;
                    stack.push_back(found->second);
                }
            }
        }
        if (inner && mask) masks.push_back(mask);
    }
    return masks;
}

bool distinct_pair(const std::vector<int>& masks, int first, int second) {
    for (int i = 0; i < static_cast<int>(masks.size()); ++i) {
        for (int j = 0; j < static_cast<int>(masks.size()); ++j) {
            if (i != j && (masks[i] & (1 << first)) && (masks[j] & (1 << second))) return true;
        }
    }
    return false;
}

Landing landing_mark(const Geometry& geometry, const std::vector<std::uint8_t>& active,
                     int radius, bool open_matching, int registry_shift = 0) {
    const auto opened = components(geometry, active, radius, open_matching, true, registry_shift);
    const auto closed = components(geometry, active, radius, !open_matching, false, registry_shift);
    const bool axis = (distinct_pair(opened, 0, 4) && distinct_pair(closed, 2, 6))
                   || (distinct_pair(opened, 2, 6) && distinct_pair(closed, 0, 4));
    const bool diagonal = (distinct_pair(opened, 1, 5) && distinct_pair(closed, 3, 7))
                       || (distinct_pair(opened, 3, 7) && distinct_pair(closed, 1, 5));
    return {static_cast<int>(axis), static_cast<int>(diagonal),
            static_cast<int>(axis && diagonal), static_cast<int>(axis || diagonal),
            static_cast<int>(axis) - static_cast<int>(diagonal)};
}

struct Contribution { int pivotal = 0; Landing landing; };

Contribution pivotal(const Geometry& geometry, std::vector<std::uint8_t> active,
                     int radius, bool matching) {
    active[0] = 0;
    const bool without = crosses(geometry, active, matching);
    active[0] = 1;
    const bool with = crosses(geometry, active, matching);
    active[0] = 0;
    const int delta = static_cast<int>(with) - static_cast<int>(without);
    if (delta < 0 || delta > 1) throw std::logic_error("nonmonotone cross event");
    return {delta, delta ? landing_mark(geometry, active, radius, matching) : Landing{}};
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

std::vector<std::uint8_t> black_configuration(int n, std::uint64_t seed,
                                              std::uint64_t counter, double p) {
    std::vector<std::uint8_t> active(n, 0);
    const auto counter_key = splitmix64(counter + 0xd1b54a32d192ed03ULL);
    constexpr double inverse_53 = 1.0 / 9007199254740992.0;
    for (int vertex = 1; vertex < n; ++vertex) {
        const auto value = splitmix64(seed ^ counter_key ^ splitmix64(static_cast<std::uint64_t>(vertex)));
        const double uniform = static_cast<double>(value >> 11) * inverse_53;
        active[vertex] = uniform < p;
    }
    return active;
}

struct Sums {
    std::uint64_t samples = 0;
    std::uint64_t primal_pivotal = 0, matching_pivotal = 0;
    std::uint64_t primal_axis = 0, primal_diagonal = 0, primal_both = 0, primal_landed = 0;
    std::uint64_t matching_axis = 0, matching_diagonal = 0, matching_both = 0, matching_landed = 0;
    std::int64_t primal_h4 = 0, matching_h4 = 0;

    void add(const Contribution& primal, const Contribution& matching) {
        ++samples;
        primal_pivotal += primal.pivotal; matching_pivotal += matching.pivotal;
        primal_axis += primal.landing.axis; primal_diagonal += primal.landing.diagonal;
        primal_both += primal.landing.both; primal_landed += primal.landing.landed;
        primal_h4 += primal.landing.h4;
        matching_axis += matching.landing.axis; matching_diagonal += matching.landing.diagonal;
        matching_both += matching.landing.both; matching_landed += matching.landing.landed;
        matching_h4 += matching.landing.h4;
    }
};

struct Options {
    std::uint64_t samples = 200000, seed = 2026106201, offset = 12000000000ULL;
    int batches = 100, radius = 3, threads = 0;
    double p = 0.592746050790;
    std::string output_prefix, git_commit = "working-tree";
};

Options parse(int argc, char** argv) {
    Options out;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("missing value for " + arg);
            return argv[i];
        };
        if (arg == "--samples") out.samples = std::stoull(next());
        else if (arg == "--batches") out.batches = std::stoi(next());
        else if (arg == "--seed") out.seed = std::stoull(next());
        else if (arg == "--replica-offset") out.offset = std::stoull(next());
        else if (arg == "--radius") out.radius = std::stoi(next());
        else if (arg == "--threads") out.threads = std::stoi(next());
        else if (arg == "--p") out.p = std::stod(next());
        else if (arg == "--git-commit") out.git_commit = next();
        else if (arg == "--output-prefix") out.output_prefix = next();
        else throw std::invalid_argument("unknown argument " + arg);
    }
    if (out.output_prefix.empty() || out.samples == 0 || out.batches <= 1
        || out.samples % out.batches || out.radius != 3 || out.p <= 0 || out.p >= 1) {
        throw std::invalid_argument("require output, positive divisible samples, batches>1, radius=3, 0<p<1");
    }
    return out;
}

void write_csv(const std::string& path, const std::vector<Geometry>& geometries,
               const std::vector<std::array<Sums, 2>>& batches) {
    std::ofstream out(path);
    out << "n,a,b,orientation,batch,samples,primal_pivotal,matching_pivotal,"
           "primal_axis,primal_diagonal,primal_both,primal_landed,primal_h4,"
           "matching_axis,matching_diagonal,matching_both,matching_landed,matching_h4\n";
    for (int batch = 0; batch < static_cast<int>(batches.size()); ++batch) {
        for (int g = 0; g < 2; ++g) {
            const auto& geometry = geometries[g]; const auto& s = batches[batch][g];
            out << geometry.n << ',' << geometry.a << ',' << geometry.b << ','
                << geometry.orientation << ',' << batch << ',' << s.samples << ','
                << s.primal_pivotal << ',' << s.matching_pivotal << ','
                << s.primal_axis << ',' << s.primal_diagonal << ',' << s.primal_both << ','
                << s.primal_landed << ',' << s.primal_h4 << ','
                << s.matching_axis << ',' << s.matching_diagonal << ',' << s.matching_both << ','
                << s.matching_landed << ',' << s.matching_h4 << '\n';
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse(argc, argv);
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
        const std::vector<Geometry> geometries = {
            make_geometry(8, 1, "first"), make_geometry(7, 4, "second")
        };
        // Fail before sampling if the frozen R=3 patch aliases in either quotient.
        for (const auto& geometry : geometries) {
            std::vector<std::uint8_t> empty(geometry.n, 0);
            (void)landing_mark(geometry, empty, options.radius, false);
        }
        const auto start = std::chrono::steady_clock::now();
        std::vector<std::array<Sums, 2>> batches(options.batches);
        const std::uint64_t per_batch = options.samples / options.batches;
#pragma omp parallel for schedule(static)
        for (int batch = 0; batch < options.batches; ++batch) {
            for (std::uint64_t k = 0; k < per_batch; ++k) {
                const std::uint64_t counter = options.offset + batch * per_batch + k;
                const auto black = black_configuration(65, options.seed, counter, options.p);
                std::vector<std::uint8_t> white(65, 0);
                for (int vertex = 1; vertex < 65; ++vertex) white[vertex] = !black[vertex];
                for (int g = 0; g < 2; ++g) {
                    batches[batch][g].add(
                        pivotal(geometries[g], black, options.radius, false),
                        pivotal(geometries[g], white, options.radius, true)
                    );
                }
            }
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        const std::filesystem::path prefix(options.output_prefix);
        std::filesystem::create_directories(prefix.parent_path());
        write_csv(options.output_prefix + ".batches.csv", geometries, batches);
        std::ofstream meta(options.output_prefix + ".metadata.json");
        meta << "{\n  \"schema\": \"matching-one/marked-pivotal-h4-pilot/v1\",\n"
             << "  \"N\": 65,\n  \"representations\": [[8,1],[7,4]],\n"
             << "  \"radius\": " << options.radius << ",\n"
             << "  \"p\": \"" << std::setprecision(12) << options.p << "\",\n"
             << "  \"samples\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.offset << ",\n"
             << "  \"replica_counter_last_exclusive\": " << options.offset + options.samples << ",\n"
             << "  \"git_commit\": \"" << options.git_commit << "\",\n"
             << "  \"elapsed_seconds\": " << std::setprecision(10) << elapsed << ",\n"
             << "  \"total_pivotal_role\": \"Russo control only\",\n"
             << "  \"primary_observables\": [\"mu4\",\"a4\"]\n}\n";
        std::cout << "completed " << options.samples
                  << " N65 marked-pivotal H4 replicas in " << elapsed << " s\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
