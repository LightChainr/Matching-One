// High-throughput bidirectional threshold-rank Newman--Ziff engine.
//
// Each same-N orientation pair shares a counter-keyed site permutation.  For
// every orientation we retain batchwise integer K_minus/K_plus histograms and
// exact joint first/second moments.  The convention matches
// scripts/threshold_rank_nz.py:
//
//   K_plus  = first k where black B_k has a rank-2 primal component;
//   K_minus = first k where white W_k loses rank-2 matching wrapping;
//             if the reverse white sweep first crosses at r, K_minus=N-r+1.

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

struct Geometry {
    int n;
    int a;
    int b;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
    std::vector<std::vector<int>> primal_incident;
    std::vector<std::vector<int>> matching_incident;
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
    // Frozen norm-5 H4/H12 discriminator, in Gaussian-lineage order.
    {325, 17, 6, 18, 1},
    {425, 16, 13, 19, 8},
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

class HomologyUnionFind {
  public:
    HomologyUnionFind(int n, int a, int b)
        : n_(n), a_(a), b_(b), parent_(n), size_(n), delta_x_(n), delta_y_(n),
          rank_(n), basis_(n) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
    }

    struct FindResult {
        int root;
        std::int64_t dx;
        std::int64_t dy;
    };

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
        const std::int64_t num0 = static_cast<std::int64_t>(a_) * dx +
                                  static_cast<std::int64_t>(b_) * dy;
        const std::int64_t num1 = -static_cast<std::int64_t>(b_) * dx +
                                  static_cast<std::int64_t>(a_) * dy;
        if (num0 % n_ != 0 || num1 % n_ != 0) {
            throw std::logic_error("cycle displacement is outside Gaussian period lattice");
        }
        return {num0 / n_, num1 / n_};
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

    std::uint8_t component_rank(int vertex) { return rank_[find(vertex).root]; }

    bool component_direction_0(int vertex) {
        const int root = find(vertex).root;
        for (std::uint8_t index = 0; index < rank_[root]; ++index) {
            if (basis_[root][index].x != 0) return true;
        }
        return false;
    }

    bool component_direction_1(int vertex) {
        const int root = find(vertex).root;
        for (std::uint8_t index = 0; index < rank_[root]; ++index) {
            if (basis_[root][index].y != 0) return true;
        }
        return false;
    }

  private:
    int n_;
    int a_;
    int b_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<std::int64_t> delta_x_;
    std::vector<std::int64_t> delta_y_;
    std::vector<std::uint8_t> rank_;
    std::vector<std::array<Winding, 2>> basis_;
};

std::vector<std::vector<int>> make_incident(int n, const std::vector<Edge>& edges) {
    std::vector<std::vector<int>> incident(n);
    for (int index = 0; index < static_cast<int>(edges.size()); ++index) {
        const Edge& edge = edges[index];
        incident[edge.i].push_back(index);
        if (edge.j != edge.i) incident[edge.j].push_back(index);
    }
    return incident;
}

Geometry make_geometry(int a, int b) {
    if (a <= 0 || b < 0 || std::gcd(a, b) != 1) {
        throw std::invalid_argument("Gaussian representation requires a>0, b>=0, gcd=1");
    }
    const std::int64_t n64 = static_cast<std::int64_t>(a) * a +
                             static_cast<std::int64_t>(b) * b;
    if (n64 > std::numeric_limits<int>::max()) throw std::invalid_argument("N too large");
    Geometry geometry{static_cast<int>(n64), a, b, {}, {}, {}, {}};
    geometry.primal_edges.reserve(2 * geometry.n);
    geometry.matching_edges.reserve(4 * geometry.n);
    const std::array<std::tuple<int, int, int>, 4> steps = {{
        {a, 1, 0}, {b, 0, 1}, {a + b, 1, 1}, {a - b, 1, -1},
    }};
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        for (std::size_t index = 0; index < steps.size(); ++index) {
            const auto [residue, dx, dy] = steps[index];
            const Edge edge{vertex, positive_mod(vertex + residue, geometry.n), dx, dy};
            if (index < 2) geometry.primal_edges.push_back(edge);
            geometry.matching_edges.push_back(edge);
        }
    }
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = make_incident(geometry.n, geometry.matching_edges);
    return geometry;
}

class ThresholdEngine {
  public:
    explicit ThresholdEngine(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n),
          union_find_(geometry.n, geometry.a, geometry.b) {}

    int first_cross(const std::vector<int>& permutation, bool matching, bool reverse) {
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        const std::vector<Edge>& edges = matching ? geometry_.matching_edges
                                                  : geometry_.primal_edges;
        const std::vector<std::vector<int>>& incident = matching
            ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) union_find_.add_edge(edge);
            }
            if (union_find_.component_crosses(vertex)) return offset + 1;
        }
        throw std::logic_error("fully occupied graph did not cross wrap");
    }

    std::pair<int, int> ranks(const std::vector<int>& permutation) {
        const int k_plus = first_cross(permutation, false, false);
        const int reverse_white = first_cross(permutation, true, true);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (k_minus > k_plus) {
            throw std::logic_error("K_minus exceeds K_plus");
        }
        return {k_minus, k_plus};
    }

  private:
    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    HomologyUnionFind union_find_;
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
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
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
    const std::uint64_t stream_key = splitmix64(
        seed ^ splitmix64(replica + 0xd1b54a32d192ed03ULL));
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
    std::uint64_t sum_minus = 0;
    std::uint64_t sum_plus = 0;
    std::uint64_t sum_minus2 = 0;
    std::uint64_t sum_plus2 = 0;
    std::uint64_t sum_product = 0;
    std::uint64_t sum_gap = 0;
    std::uint64_t sum_gap2 = 0;

    explicit RankCounts(int n = 0) : minus(n + 1), plus(n + 1) {}

    void add(int k_minus, int k_plus) {
        if (!(1 <= k_minus && k_minus <= k_plus &&
              k_plus < static_cast<int>(plus.size()))) {
            throw std::logic_error("invalid threshold rank pair");
        }
        ++minus[k_minus];
        ++plus[k_plus];
        ++samples;
        sum_minus += k_minus;
        sum_plus += k_plus;
        sum_minus2 += static_cast<std::uint64_t>(k_minus) * k_minus;
        sum_plus2 += static_cast<std::uint64_t>(k_plus) * k_plus;
        sum_product += static_cast<std::uint64_t>(k_minus) * k_plus;
        const std::uint64_t gap = static_cast<std::uint64_t>(k_plus - k_minus);
        sum_gap += gap;
        sum_gap2 += gap * gap;
    }
};

struct PairBatch {
    RankCounts first;
    RankCounts second;
    explicit PairBatch(int n = 0) : first(n), second(n) {}
};

void self_test() {
    const Geometry geometry = make_geometry(2, 1);
    ThresholdEngine engine(geometry);
    RankCounts counts(geometry.n);
    std::vector<int> permutation(geometry.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto ranks = engine.ranks(permutation);
        counts.add(ranks.first, ranks.second);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (counts.samples != 120 || counts.minus[3] != 120 || counts.plus[4] != 120) {
        throw std::runtime_error("N=5 all-permutation rank histogram regression failed");
    }
    counter_permutation(5, 17, 0, permutation);
    const std::vector<int> expected = {4, 3, 1, 0, 2};
    if (permutation != expected || engine.ranks(permutation) != std::make_pair(3, 4)) {
        throw std::runtime_error("Python-compatible counter/permutation regression failed");
    }
    std::cout << "self-test passed: N=5 all permutations, K_minus<=K_plus, "
                 "Python-compatible counter permutation\n";
}

struct Options {
    std::uint64_t samples = 1000000;
    int batches = 100;
    std::uint64_t seed = 20260828;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    int only_n = 0;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N          replicas per orientation pair (default 1000000)\n"
        << "  --batches B          equal batches (default 100)\n"
        << "  --seed S             unsigned 64-bit seed (default 20260828)\n"
        << "  --replica-offset K   first sample counter (default 0)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --n N                only N=65,85,130,145,170,185,265,325,425 (default all)\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix PATH writes .hist.csv, .moments.csv, .metadata.json\n"
        << "  --self-test           exact tiny regression and exit\n";
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
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(next(), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(next(), arg);
        else if (arg == "--n") options.only_n = parse_number<int>(next(), arg);
        else if (arg == "--git-commit") options.git_commit = next();
        else if (arg == "--output-prefix") options.output_prefix = next();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be positive and divisible by batches>=2");
    }
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (options.only_n != 0 && std::none_of(kDesigns.begin(), kDesigns.end(),
            [&](const PairDesign& design) { return design.n == options.only_n; })) {
        throw std::invalid_argument("--n must be one of 65,85,130,145,170,185,265,325,425");
    }
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
    const std::time_t value = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
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

void run_design(const PairDesign& design, const Options& options,
                std::ofstream& histogram, std::ofstream& moments) {
    const Geometry first_geometry = make_geometry(design.a1, design.b1);
    const Geometry second_geometry = make_geometry(design.a2, design.b2);
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<PairBatch> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) output.emplace_back(design.n);

#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        PairBatch local(design.n);
        ThresholdEngine first_engine(first_geometry);
        ThresholdEngine second_engine(second_geometry);
        std::vector<int> permutation;
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(design.n, options.seed, replica, permutation);
            const auto first = first_engine.ranks(permutation);
            const auto second = second_engine.ranks(permutation);
            local.first.add(first.first, first.second);
            local.second.add(second.first, second.second);
        }
        output[batch] = std::move(local);
    }

    auto write_orientation = [&](int batch, const char* orientation, int a, int b,
                                 const RankCounts& counts) {
        for (int rank = 1; rank <= design.n; ++rank) {
            if (counts.minus[rank]) {
                histogram << design.n << ',' << a << ',' << b << ',' << orientation << ','
                          << batch << ',' << counts.samples << ",minus," << rank << ','
                          << counts.minus[rank] << '\n';
            }
            if (counts.plus[rank]) {
                histogram << design.n << ',' << a << ',' << b << ',' << orientation << ','
                          << batch << ',' << counts.samples << ",plus," << rank << ','
                          << counts.plus[rank] << '\n';
            }
        }
        moments << design.n << ',' << a << ',' << b << ',' << orientation << ',' << batch
                << ',' << counts.samples << ',' << counts.sum_minus << ',' << counts.sum_plus
                << ',' << counts.sum_minus2 << ',' << counts.sum_plus2 << ','
                << counts.sum_product << ',' << counts.sum_gap << ',' << counts.sum_gap2 << '\n';
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write_orientation(batch, "first", design.a1, design.b1, output[batch].first);
        write_orientation(batch, "second", design.a2, design.b2, output[batch].second);
    }
    std::cout << "completed N=" << design.n << " pair (" << design.a1 << ',' << design.b1
              << ")/(" << design.a2 << ',' << design.b2 << ") samples=" << options.samples
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
    const std::filesystem::path histogram_path = options.output_prefix.string() + ".hist.csv";
    const std::filesystem::path moments_path = options.output_prefix.string() + ".moments.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    std::ofstream histogram(histogram_path), moments(moments_path);
    if (!histogram || !moments) throw std::runtime_error("cannot open output files");
    histogram << "n,a,b,orientation,batch,samples,kind,k,count\n";
    moments << "n,a,b,orientation,batch,samples,sum_kminus,sum_kplus,sum_kminus2,"
               "sum_kplus2,sum_product,sum_gap,sum_gap2\n";
    const auto started = std::chrono::steady_clock::now();
    for (const PairDesign& design : kDesigns) {
        if (options.only_n == 0 || options.only_n == design.n) {
            run_design(design, options, histogram, moments);
        }
    }
    histogram.close();
    moments.close();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    std::ofstream metadata(metadata_path);
    if (!metadata) throw std::runtime_error("cannot open metadata output");
    metadata << "{\n"
             << "  \"engine\": \"same-N Gaussian threshold-rank Newman-Ziff\",\n"
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
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"rng\": \"counter-derived SplitMix64 stream plus unbiased Fisher-Yates\",\n"
             << "  \"coupling\": \"same cyclic permutation shared by same-N orientations\",\n"
             << "  \"channel\": \"rank-2 cross wrapping\",\n"
             << "  \"K_plus\": \"first black primal cross rank, 1-based\",\n"
             << "  \"K_minus\": \"first black rank after white matching cross is lost; N-r+1\",\n"
             << "  \"sparse_joint_histogram\": false,\n"
             << "  \"per_batch_joint_moments\": true,\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"designs\": [\n";
    bool first = true;
    for (const PairDesign& design : kDesigns) {
        if (options.only_n != 0 && options.only_n != design.n) continue;
        if (!first) metadata << ",\n";
        first = false;
        metadata << "    {\"N\": " << design.n << ", \"first\": ["
                 << design.a1 << ',' << design.b1 << "], \"second\": ["
                 << design.a2 << ',' << design.b2 << "]}";
    }
    metadata << "\n  ],\n"
             << "  \"histogram_csv\": \"" << json_escape(histogram_path.string()) << "\",\n"
             << "  \"moments_csv\": \"" << json_escape(moments_path.string()) << "\"\n"
             << "}\n";
    std::cout << "wrote " << histogram_path << "\nwrote " << moments_path
              << "\nwrote " << metadata_path << '\n';
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
