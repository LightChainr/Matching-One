// High-throughput bidirectional threshold-rank Newman--Ziff engine for
// ordinary axis-aligned L x L square tori.
//
// This fills the production gap in Issue #47: the existing Gaussian-cyclic
// engine requires gcd(a,b)=1 and therefore cannot represent axis L>1.
//
// Output uses the same histogram/moment schema as threshold_rank_orientation_mc
// with a=L, b=0 and orientation="axis".  K conventions are identical:
//
//   K_plus  = first occupied rank with rank-2 primal homology;
//   K_minus = first k where the white matching complement has lost rank-2
//             homology, reconstructed by the reverse sweep.
//
// Build:
//   g++ -O3 -std=c++17 -fopenmp src/threshold_rank_axis_mc.cpp -o build/threshold_rank_axis_mc

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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

struct Geometry {
    int L;
    int n;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
    std::vector<std::vector<int>> primal_incident;
    std::vector<std::vector<int>> matching_incident;
};

std::vector<std::vector<int>> make_incident(int n, const std::vector<Edge>& edges) {
    std::vector<std::vector<int>> incident(n);
    for (int index = 0; index < static_cast<int>(edges.size()); ++index) {
        incident[edges[index].i].push_back(index);
        if (edges[index].j != edges[index].i) incident[edges[index].j].push_back(index);
    }
    return incident;
}

Geometry make_axis_geometry(int L) {
    if (L < 1) throw std::invalid_argument("L must be positive");
    const std::int64_t n64 = static_cast<std::int64_t>(L) * L;
    if (n64 > std::numeric_limits<int>::max()) throw std::invalid_argument("L^2 exceeds int range");
    Geometry g{L, static_cast<int>(n64), {}, {}, {}, {}};
    g.primal_edges.reserve(2 * g.n);
    g.matching_edges.reserve(4 * g.n);
    auto id = [L](int x, int y) {
        x %= L; if (x < 0) x += L;
        y %= L; if (y < 0) y += L;
        return x + L * y;
    };
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = id(x, y);
            const std::array<Edge, 4> edges = {{
                {i, id(x + 1, y), 1, 0},
                {i, id(x, y + 1), 0, 1},
                {i, id(x + 1, y + 1), 1, 1},
                {i, id(x + 1, y - 1), 1, -1},
            }};
            g.primal_edges.push_back(edges[0]);
            g.primal_edges.push_back(edges[1]);
            for (const Edge& edge : edges) g.matching_edges.push_back(edge);
        }
    }
    g.primal_incident = make_incident(g.n, g.primal_edges);
    g.matching_incident = make_incident(g.n, g.matching_edges);
    return g;
}

class AxisHomologyUnionFind {
  public:
    AxisHomologyUnionFind(int n, int L)
        : n_(n), L_(L), parent_(n), size_(n), delta_x_(n), delta_y_(n),
          rank_(n), basis_(n) { reset(); }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
    }

    struct FindResult { int root; std::int64_t dx; std::int64_t dy; };

    FindResult find(int vertex) {
        if (parent_[vertex] == vertex) return {vertex, 0, 0};
        const int old_parent = parent_[vertex];
        const FindResult above = find(old_parent);
        delta_x_[vertex] += above.dx;
        delta_y_[vertex] += above.dy;
        parent_[vertex] = above.root;
        return {above.root, delta_x_[vertex], delta_y_[vertex]};
    }

    Winding period_coordinates(std::int64_t dx, std::int64_t dy) const {
        if (dx % L_ != 0 || dy % L_ != 0) {
            throw std::logic_error("cycle displacement is outside axis period lattice");
        }
        return {dx / L_, dy / L_};
    }

    void extend(int root, Winding value) {
        if ((value.x == 0 && value.y == 0) || rank_[root] == 2) return;
        value = primitive(value);
        if (rank_[root] == 0) {
            basis_[root][0] = value;
            rank_[root] = 1;
            return;
        }
        const Winding first = basis_[root][0];
        if (first.x * value.y != first.y * value.x) {
            basis_[root][1] = value;
            rank_[root] = 2;
        }
    }

    void add_edge(const Edge& edge) {
        FindResult first = find(edge.i);
        FindResult second = find(edge.j);
        std::int64_t root_dx = first.dx + edge.dx - second.dx;
        std::int64_t root_dy = first.dy + edge.dy - second.dy;
        if (first.root == second.root) {
            extend(first.root, period_coordinates(root_dx, root_dy));
            return;
        }
        if (size_[first.root] < size_[second.root]) {
            std::swap(first, second);
            root_dx = -root_dx;
            root_dy = -root_dy;
        }
        parent_[second.root] = first.root;
        delta_x_[second.root] = root_dx;
        delta_y_[second.root] = root_dy;
        size_[first.root] += size_[second.root];
        for (std::uint8_t index = 0; index < rank_[second.root]; ++index) {
            extend(first.root, basis_[second.root][index]);
        }
        rank_[second.root] = 0;
    }

    bool component_crosses(int vertex) { return rank_[find(vertex).root] == 2; }

  private:
    int n_;
    int L_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<std::int64_t> delta_x_;
    std::vector<std::int64_t> delta_y_;
    std::vector<std::uint8_t> rank_;
    std::vector<std::array<Winding, 2>> basis_;
};

class ThresholdEngine {
  public:
    explicit ThresholdEngine(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n), uf_(geometry.n, geometry.L) {}

    int first_cross(const std::vector<int>& permutation, bool matching, bool reverse) {
        std::fill(active_.begin(), active_.end(), 0);
        uf_.reset();
        const auto& edges = matching ? geometry_.matching_edges : geometry_.primal_edges;
        const auto& incident = matching ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) uf_.add_edge(edge);
            }
            if (uf_.component_crosses(vertex)) return offset + 1;
        }
        throw std::logic_error("fully occupied axis graph did not cross wrap");
    }

    std::pair<int, int> ranks(const std::vector<int>& permutation) {
        const int k_plus = first_cross(permutation, false, false);
        const int reverse_white = first_cross(permutation, true, true);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (k_minus > k_plus) throw std::logic_error("K_minus exceeds K_plus");
        return {k_minus, k_plus};
    }

  private:
    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    AxisHomologyUnionFind uf_;
};

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

class SplitMixStream {
  public:
    explicit SplitMixStream(std::uint64_t state) : state_(state) {}
    std::uint64_t next() {
        state_ += 0x9e3779b97f4a7c15ULL;
        std::uint64_t value = state_;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }
    std::uint64_t below(std::uint64_t bound) {
        const std::uint64_t remainder =
            (std::numeric_limits<std::uint64_t>::max() % bound + 1) % bound;
        const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max() - remainder;
        while (true) {
            const std::uint64_t value = next();
            if (remainder == 0 || value <= maximum) return value % bound;
        }
    }
  private:
    std::uint64_t state_;
};

void counter_permutation(int n, std::uint64_t seed, std::uint64_t replica,
                         std::vector<int>& permutation) {
    permutation.resize(n);
    std::iota(permutation.begin(), permutation.end(), 0);
    const std::uint64_t stream_key = splitmix64(seed ^ splitmix64(replica + 0xd1b54a32d192ed03ULL));
    SplitMixStream generator(stream_key);
    for (int stop = n - 1; stop > 0; --stop) {
        const int other = static_cast<int>(generator.below(static_cast<std::uint64_t>(stop + 1)));
        std::swap(permutation[stop], permutation[other]);
    }
}

struct RankCounts {
    std::vector<std::uint64_t> minus;
    std::vector<std::uint64_t> plus;
    std::uint64_t samples = 0;
    std::uint64_t sum_minus = 0, sum_plus = 0, sum_minus2 = 0, sum_plus2 = 0;
    std::uint64_t sum_product = 0, sum_gap = 0, sum_gap2 = 0;
    explicit RankCounts(int n = 0) : minus(n + 1), plus(n + 1) {}
    void add(int k_minus, int k_plus) {
        if (!(1 <= k_minus && k_minus <= k_plus && k_plus < static_cast<int>(plus.size()))) {
            throw std::logic_error("invalid threshold rank pair");
        }
        ++minus[k_minus]; ++plus[k_plus]; ++samples;
        sum_minus += k_minus; sum_plus += k_plus;
        sum_minus2 += static_cast<std::uint64_t>(k_minus) * k_minus;
        sum_plus2 += static_cast<std::uint64_t>(k_plus) * k_plus;
        sum_product += static_cast<std::uint64_t>(k_minus) * k_plus;
        const auto gap = static_cast<std::uint64_t>(k_plus - k_minus);
        sum_gap += gap; sum_gap2 += gap * gap;
    }
};

void self_test() {
    const Geometry g = make_axis_geometry(2);
    ThresholdEngine engine(g);
    RankCounts counts(g.n);
    std::vector<int> permutation(g.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto value = engine.ranks(permutation);
        counts.add(value.first, value.second);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (counts.samples != 24 || counts.minus[2] != 16 || counts.minus[3] != 8 ||
        counts.plus[3] != 24) {
        throw std::runtime_error("axis L=2 all-permutation histogram regression failed");
    }
    if (counts.minus[1] || counts.minus[4] || counts.plus[1] || counts.plus[2] || counts.plus[4]) {
        throw std::runtime_error("axis L=2 unexpected threshold-rank mass");
    }
    std::cout << "self-test passed: axis L=2 all 24 permutations, exact K histograms\n";
}

struct Options {
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

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --L L                axis torus linear size (required)\n"
        << "  --samples N          replicas (default 1000000)\n"
        << "  --batches B          equal batches (default 100)\n"
        << "  --seed S             unsigned 64-bit seed\n"
        << "  --replica-offset K   first sample counter\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix P    writes .hist.csv/.moments.csv/.metadata.json\n"
        << "  --self-test           exact L=2 regression and exit\n";
    std::exit(status);
}

template <typename T>
T parse_number(const std::string& text, const std::string& option) {
    std::istringstream input(text);
    T value{}; input >> value;
    if (!input || !input.eof()) throw std::invalid_argument("invalid value for " + option);
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string { if (++i >= argc) usage(argv[0], 2); return argv[i]; };
        if (arg == "--L") options.L = parse_number<int>(next(), arg);
        else if (arg == "--samples") options.samples = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(next(), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(next(), arg);
        else if (arg == "--git-commit") options.git_commit = next();
        else if (arg == "--output-prefix") options.output_prefix = next();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if (options.L < 2) throw std::invalid_argument("--L must be at least 2 for production");
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

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (const char ch : value) {
        if (ch == '\\' || ch == '"') output << '\\' << ch;
        else if (ch == '\n') output << "\\n";
        else output << ch;
    }
    return output.str();
}

std::string utc_now() {
    const std::time_t value = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm time{};
#ifdef _WIN32
    gmtime_s(&time, &value);
#else
    gmtime_r(&value, &time);
#endif
    std::ostringstream output;
    output << std::put_time(&time, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

int run(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    if (options.self_test) { self_test(); return 0; }
    const Geometry geometry = make_axis_geometry(options.L);
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<RankCounts> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) output.emplace_back(geometry.n);
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        RankCounts local(geometry.n);
        ThresholdEngine engine(geometry);
        std::vector<int> permutation;
        const std::uint64_t begin = options.replica_offset + static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(geometry.n, options.seed, replica, permutation);
            const auto value = engine.ranks(permutation);
            local.add(value.first, value.second);
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
    hist << "n,a,b,orientation,batch,samples,kind,k,count\n";
    moments << "n,a,b,orientation,batch,samples,sum_kminus,sum_kplus,sum_kminus2,sum_kplus2,sum_product,sum_gap,sum_gap2\n";
    for (int batch = 0; batch < options.batches; ++batch) {
        const RankCounts& counts = output[batch];
        for (int rank = 1; rank <= geometry.n; ++rank) {
            if (counts.minus[rank]) hist << geometry.n << ',' << options.L << ",0,axis," << batch << ',' << counts.samples << ",minus," << rank << ',' << counts.minus[rank] << '\n';
            if (counts.plus[rank]) hist << geometry.n << ',' << options.L << ",0,axis," << batch << ',' << counts.samples << ",plus," << rank << ',' << counts.plus[rank] << '\n';
        }
        moments << geometry.n << ',' << options.L << ",0,axis," << batch << ',' << counts.samples << ','
                << counts.sum_minus << ',' << counts.sum_plus << ',' << counts.sum_minus2 << ','
                << counts.sum_plus2 << ',' << counts.sum_product << ',' << counts.sum_gap << ',' << counts.sum_gap2 << '\n';
    }
    hist.close(); moments.close();

    std::ostringstream command;
    for (int i = 0; i < argc; ++i) { if (i) command << ' '; command << argv[i]; }
    std::ofstream meta(meta_path);
    meta << "{\n"
         << "  \"engine\": \"axis threshold-rank Newman-Ziff\",\n"
         << "  \"generated_utc\": \"" << utc_now() << "\",\n"
         << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
         << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
         << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
         << "  \"L\": " << options.L << ",\n"
         << "  \"N\": " << geometry.n << ",\n"
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
    std::cout << "completed axis L=" << options.L << " N=" << geometry.n
              << " samples=" << options.samples << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try { return run(argc, argv); }
    catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
