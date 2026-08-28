// Discovery-stage paired Monte Carlo for the square-site matching function.
//
// This deliberately uses a narrow, fixed-p common-random-number design rather
// than a full Newman-Ziff implementation.  It is intended to establish whether
// the axis/diamond Pell signal is large enough to justify a production engine.
//
// Build (Linux/macOS with an OpenMP-capable compiler):
//   g++ -O3 -std=c++17 -fopenmp src/pell_matching_mc.cpp -o pell_matching_mc
//
// Example:
//   ./pell_matching_mc --axis 17 --diamond 12 --samples 10000 --batches 20
//       --p-ref 0.59274605 --h 0.001 --seed 20260828 --threads 16
//       --output-prefix results/pell/a17_d12

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
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

struct Geometry {
    std::string name;
    int L;
    int n;
    double physical_period;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
};

class WrapUnionFind {
  public:
    explicit WrapUnionFind(int n)
        : parent_(n), size_(n), delta_x_(n), delta_y_(n), wrap_(n) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(wrap_.begin(), wrap_.end(), false);
    }

    struct FindResult {
        int root;
        int dx;
        int dy;
    };

    FindResult find(int x) {
        int node = x;
        int total_x = 0;
        int total_y = 0;
        while (parent_[node] != node) {
            total_x += delta_x_[node];
            total_y += delta_y_[node];
            node = parent_[node];
        }
        const int root = node;

        node = x;
        int remaining_x = total_x;
        int remaining_y = total_y;
        while (parent_[node] != node) {
            const int next = parent_[node];
            const int step_x = delta_x_[node];
            const int step_y = delta_y_[node];
            parent_[node] = root;
            delta_x_[node] = remaining_x;
            delta_y_[node] = remaining_y;
            remaining_x -= step_x;
            remaining_y -= step_y;
            node = next;
        }
        return {root, total_x, total_y};
    }

    void add_edge(int i, int j, int edge_dx, int edge_dy) {
        const FindResult fi = find(i);
        const FindResult fj = find(j);
        // Required position(root_j) - position(root_i).
        const int root_dx = fi.dx + edge_dx - fj.dx;
        const int root_dy = fi.dy + edge_dy - fj.dy;

        if (fi.root == fj.root) {
            if (root_dx != 0 || root_dy != 0) {
                wrap_[fi.root] = true;
            }
            return;
        }

        if (size_[fi.root] >= size_[fj.root]) {
            parent_[fj.root] = fi.root;
            delta_x_[fj.root] = root_dx;
            delta_y_[fj.root] = root_dy;
            size_[fi.root] += size_[fj.root];
            wrap_[fi.root] = wrap_[fi.root] || wrap_[fj.root];
        } else {
            parent_[fi.root] = fj.root;
            delta_x_[fi.root] = -root_dx;
            delta_y_[fi.root] = -root_dy;
            size_[fj.root] += size_[fi.root];
            wrap_[fj.root] = wrap_[fi.root] || wrap_[fj.root];
        }
    }

    bool component_wraps(int x) {
        return wrap_[find(x).root];
    }

  private:
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<int> delta_x_;
    std::vector<int> delta_y_;
    std::vector<std::uint8_t> wrap_;
};

Geometry axis_geometry(int L) {
    if (L < 2) {
        throw std::invalid_argument("axis L must be at least 2");
    }
    const std::int64_t n64 = static_cast<std::int64_t>(L) * L;
    if (n64 > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("axis geometry is too large for 32-bit site ids");
    }
    Geometry g{"axis", L, static_cast<int>(n64), static_cast<double>(L), {}, {}};
    g.primal_edges.reserve(2 * g.n);
    g.matching_edges.reserve(4 * g.n);
    auto id = [L](int x, int y) { return x + L * y; };
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = id(x, y);
            const Edge east{i, id((x + 1) % L, y), 1, 0};
            const Edge north{i, id(x, (y + 1) % L), 0, 1};
            g.primal_edges.push_back(east);
            g.primal_edges.push_back(north);
            g.matching_edges.push_back(east);
            g.matching_edges.push_back(north);
            g.matching_edges.push_back({i, id((x + 1) % L, (y + 1) % L), 1, 1});
            g.matching_edges.push_back({i, id((x + 1) % L, (y + L - 1) % L), 1, -1});
        }
    }
    return g;
}

Geometry diamond_geometry(int L) {
    if (L < 2) {
        throw std::invalid_argument("diamond L must be at least 2");
    }
    const std::int64_t n64 = 2LL * L * L;
    if (n64 > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("diamond geometry is too large for 32-bit site ids");
    }
    const int period = 2 * L;
    Geometry g{"diamond", L, static_cast<int>(n64), std::sqrt(2.0) * L, {}, {}};
    std::vector<int> ids(static_cast<std::size_t>(period) * period, -1);
    std::vector<std::pair<int, int>> coords;
    coords.reserve(g.n);
    // Preserve the reference implementation's (u outer, v inner) ordering.
    for (int u = 0; u < period; ++u) {
        for (int v = 0; v < period; ++v) {
            if (((u - v) & 1) == 0) {
                ids[static_cast<std::size_t>(u) * period + v] =
                    static_cast<int>(coords.size());
                coords.emplace_back(u, v);
            }
        }
    }
    auto mod = [period](int x) {
        x %= period;
        return x < 0 ? x + period : x;
    };
    auto id = [&](int u, int v) {
        const int result = ids[static_cast<std::size_t>(mod(u)) * period + mod(v)];
        if (result < 0) {
            throw std::logic_error("diamond edge reached an invalid parity site");
        }
        return result;
    };
    g.primal_edges.reserve(2 * g.n);
    g.matching_edges.reserve(4 * g.n);
    for (int i = 0; i < g.n; ++i) {
        const auto [u, v] = coords[i];
        const Edge east{i, id(u + 1, v - 1), 1, -1};
        const Edge north{i, id(u + 1, v + 1), 1, 1};
        g.primal_edges.push_back(east);
        g.primal_edges.push_back(north);
        g.matching_edges.push_back(east);
        g.matching_edges.push_back(north);
        g.matching_edges.push_back({i, id(u + 2, v), 2, 0});
        g.matching_edges.push_back({i, id(u, v + 2), 0, 2});
    }
    return g;
}

bool wraps(const std::vector<std::uint8_t>& active, const std::vector<Edge>& edges,
           WrapUnionFind& uf) {
    uf.reset();
    for (const Edge& edge : edges) {
        if (active[edge.i] && active[edge.j]) {
            uf.add_edge(edge.i, edge.j, edge.dx, edge.dy);
        }
    }
    for (int i = 0; i < static_cast<int>(active.size()); ++i) {
        if (active[i] && uf.component_wraps(i)) {
            return true;
        }
    }
    return false;
}

int matching_observable(const Geometry& geometry, const std::vector<double>& uniforms,
                        double p, std::vector<std::uint8_t>& active,
                        WrapUnionFind& primal_uf, WrapUnionFind& matching_uf) {
    active.resize(geometry.n);
    for (int i = 0; i < geometry.n; ++i) {
        active[i] = uniforms[i] < p;
    }
    const bool black_wrap = wraps(active, geometry.primal_edges, primal_uf);
    for (std::uint8_t& value : active) {
        value = !value;
    }
    const bool white_matching_wrap = wraps(active, geometry.matching_edges, matching_uf);
    return static_cast<int>(black_wrap) - static_cast<int>(white_matching_wrap);
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

double counter_uniform(std::uint64_t seed, std::uint64_t replica, std::uint64_t site) {
    // Stateless: results do not depend on thread count, scheduling, or batching.
    const std::uint64_t key = seed ^ splitmix64(replica + 0xd1b54a32d192ed03ULL) ^
                              splitmix64(site + 0x94d049bb133111ebULL);
    return static_cast<double>(splitmix64(key) >> 11) * 0x1.0p-53;
}

std::vector<long long> exact_power_coefficients(const Geometry& geometry) {
    if (geometry.n > 24) {
        throw std::invalid_argument("exact regression is limited to N<=24");
    }
    std::vector<long long> bernstein(geometry.n + 1, 0);
    std::vector<std::uint8_t> black(geometry.n), white(geometry.n);
    WrapUnionFind primal_uf(geometry.n), matching_uf(geometry.n);
    const std::uint64_t configurations = 1ULL << geometry.n;
    for (std::uint64_t mask = 0; mask < configurations; ++mask) {
        for (int i = 0; i < geometry.n; ++i) {
            black[i] = (mask >> i) & 1U;
            white[i] = !black[i];
        }
        const int value = static_cast<int>(wraps(black, geometry.primal_edges, primal_uf)) -
                          static_cast<int>(wraps(white, geometry.matching_edges, matching_uf));
        bernstein[__builtin_popcountll(mask)] += value;
    }
    std::vector<long long> power(geometry.n + 1, 0);
    for (int k = 0; k <= geometry.n; ++k) {
        long long choose = 1;
        for (int j = 0; j <= geometry.n - k; ++j) {
            const int degree = k + j;
            power[degree] += bernstein[k] * ((j & 1) ? -choose : choose);
            if (j < geometry.n - k) {
                choose = choose * (geometry.n - k - j) / (j + 1);
            }
        }
    }
    while (power.size() > 1 && power.back() == 0) {
        power.pop_back();
    }
    return power;
}

void require_equal(const std::string& label, const std::vector<long long>& actual,
                   const std::vector<long long>& expected) {
    if (actual != expected) {
        std::ostringstream message;
        message << label << " exact polynomial regression failed\nactual: ";
        for (long long value : actual) message << value << ' ';
        message << "\nexpected: ";
        for (long long value : expected) message << value << ' ';
        throw std::runtime_error(message.str());
    }
}

void self_test() {
    require_equal("axis L=2", exact_power_coefficients(axis_geometry(2)),
                  {-1, 0, 4, 0, -2});
    require_equal("axis L=3", exact_power_coefficients(axis_geometry(3)),
                  {-1, 0, 0, 6, 0, 0, 0, -18, 18, -4});
    require_equal("diamond L=2", exact_power_coefficients(diamond_geometry(2)),
                  {-1, 0, 0, 0, 28, -48, 24, 0, -2});

    // Counter RNG must be insensitive to call order.
    const double first = counter_uniform(17, 23, 5);
    (void)counter_uniform(17, 999, 7);
    if (first != counter_uniform(17, 23, 5) || !(first >= 0.0 && first < 1.0)) {
        throw std::runtime_error("counter RNG reproducibility regression failed");
    }
    std::cout << "self-test passed: axis L=2,3; diamond L=2; counter RNG\n";
}

struct Options {
    int axis_L = 17;
    int diamond_L = 12;
    std::uint64_t samples = 10000;
    int batches = 20;
    double p_ref = 0.59274605;
    double h = 0.001;
    std::uint64_t seed = 20260828;
    int threads = 0;
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --axis L             axis torus side (default 17)\n"
        << "  --diamond L          diamond torus parameter (default 12)\n"
        << "  --samples N          total paired replicas (default 10000)\n"
        << "  --batches B          equal independent batches (default 20)\n"
        << "  --p-ref P            center probability (default 0.59274605)\n"
        << "  --h H                common-p central-difference step (default 0.001)\n"
        << "  --seed S             unsigned 64-bit seed (default 20260828)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --output-prefix PATH writes PATH.batches.csv and PATH.metadata.json\n"
        << "  --self-test           exact tiny-torus regression and exit\n"
        << "  --help                show this help\n";
    std::exit(status);
}

template <typename T>
T parse_number(const std::string& text, const std::string& option) {
    std::istringstream input(text);
    T value;
    input >> value;
    if (!input || !input.eof()) {
        throw std::invalid_argument("invalid value for " + option + ": " + text);
    }
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("missing value for " + arg);
            return argv[i];
        };
        if (arg == "--axis") options.axis_L = parse_number<int>(value(), arg);
        else if (arg == "--diamond") options.diamond_L = parse_number<int>(value(), arg);
        else if (arg == "--samples") options.samples = parse_number<std::uint64_t>(value(), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(value(), arg);
        else if (arg == "--p-ref") options.p_ref = parse_number<double>(value(), arg);
        else if (arg == "--h") options.h = parse_number<double>(value(), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(value(), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(value(), arg);
        else if (arg == "--output-prefix") options.output_prefix = value();
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help" || arg == "-h") usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    return options;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const char c : text) {
        switch (c) {
            case '\\': escaped << "\\\\"; break;
            case '"': escaped << "\\\""; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default: escaped << c;
        }
    }
    return escaped.str();
}

struct BatchRow {
    int batch;
    std::uint64_t first_replica;
    std::uint64_t samples;
    double p;
    long long axis_sum;
    long long diamond_sum;
};

int runtime_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.self_test) {
            self_test();
            return 0;
        }
        if (options.output_prefix.empty()) {
            throw std::invalid_argument("--output-prefix is required for a simulation");
        }
        if (options.batches < 2) throw std::invalid_argument("--batches must be at least 2");
        if (options.samples < static_cast<std::uint64_t>(options.batches)) {
            throw std::invalid_argument("--samples must be at least --batches");
        }
        if (options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
            throw std::invalid_argument("--samples must be divisible by --batches for jackknife analysis");
        }
        if (!(options.h > 0.0) || !(options.p_ref - options.h > 0.0) ||
            !(options.p_ref + options.h < 1.0)) {
            throw std::invalid_argument("require 0 < p_ref-h < p_ref+h < 1");
        }
        if (options.threads < 0) throw std::invalid_argument("--threads cannot be negative");
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
#else
        if (options.threads > 1) {
            throw std::invalid_argument("binary was compiled without OpenMP support");
        }
#endif

        const Geometry axis = axis_geometry(options.axis_L);
        const Geometry diamond = diamond_geometry(options.diamond_L);
        const std::vector<double> probabilities = {
            options.p_ref - options.h, options.p_ref, options.p_ref + options.h};
        const std::uint64_t per_batch = options.samples / options.batches;
        std::vector<BatchRow> rows;
        rows.reserve(static_cast<std::size_t>(options.batches) * probabilities.size());

        const auto started = std::chrono::steady_clock::now();
        for (int batch = 0; batch < options.batches; ++batch) {
            const std::uint64_t first = static_cast<std::uint64_t>(batch) * per_batch;
            const int thread_count = runtime_threads();
            std::vector<std::vector<long long>> thread_axis(
                thread_count, std::vector<long long>(probabilities.size(), 0));
            std::vector<std::vector<long long>> thread_diamond(
                thread_count, std::vector<long long>(probabilities.size(), 0));

#pragma omp parallel
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                const int max_n = std::max(axis.n, diamond.n);
                std::vector<double> uniforms(max_n);
                std::vector<std::uint8_t> axis_active(axis.n), diamond_active(diamond.n);
                WrapUnionFind axis_primal(axis.n), axis_matching(axis.n);
                WrapUnionFind diamond_primal(diamond.n), diamond_matching(diamond.n);

#pragma omp for schedule(static)
                for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                    const std::uint64_t replica = first + offset;
                    for (int site = 0; site < max_n; ++site) {
                        uniforms[site] = counter_uniform(options.seed, replica, site);
                    }
                    for (std::size_t point = 0; point < probabilities.size(); ++point) {
                        thread_axis[tid][point] += matching_observable(
                            axis, uniforms, probabilities[point], axis_active,
                            axis_primal, axis_matching);
                        thread_diamond[tid][point] += matching_observable(
                            diamond, uniforms, probabilities[point], diamond_active,
                            diamond_primal, diamond_matching);
                    }
                }
            }

            for (std::size_t point = 0; point < probabilities.size(); ++point) {
                long long axis_sum = 0;
                long long diamond_sum = 0;
                for (int thread = 0; thread < thread_count; ++thread) {
                    axis_sum += thread_axis[thread][point];
                    diamond_sum += thread_diamond[thread][point];
                }
                rows.push_back({batch, first, per_batch, probabilities[point],
                                axis_sum, diamond_sum});
            }
            std::cerr << "completed batch " << (batch + 1) << '/' << options.batches << '\n';
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        const std::filesystem::path parent = options.output_prefix.parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        const std::filesystem::path csv_path = options.output_prefix.string() + ".batches.csv";
        const std::filesystem::path json_path = options.output_prefix.string() + ".metadata.json";

        std::ofstream csv(csv_path);
        if (!csv) throw std::runtime_error("cannot open output: " + csv_path.string());
        csv << "batch,first_replica,samples,p,axis_sum,diamond_sum,paired_difference_sum,"
               "axis_mean,diamond_mean,paired_difference_mean\n";
        csv << std::setprecision(17);
        for (const BatchRow& row : rows) {
            const long long difference = row.diamond_sum - row.axis_sum;
            csv << row.batch << ',' << row.first_replica << ',' << row.samples << ',' << row.p
                << ',' << row.axis_sum << ',' << row.diamond_sum << ',' << difference << ','
                << static_cast<double>(row.axis_sum) / row.samples << ','
                << static_cast<double>(row.diamond_sum) / row.samples << ','
                << static_cast<double>(difference) / row.samples << '\n';
        }
        csv.close();
        if (!csv) throw std::runtime_error("failed while writing: " + csv_path.string());

        std::ofstream metadata(json_path);
        if (!metadata) throw std::runtime_error("cannot open output: " + json_path.string());
        metadata << std::setprecision(17)
                 << "{\n"
                 << "  \"engine\": \"pell_matching_mc_fixed_p_v1\",\n"
                 << "  \"design\": \"three-point fixed-p common-random-number discovery scan\",\n"
                 << "  \"axis_L\": " << axis.L << ",\n"
                 << "  \"axis_sites\": " << axis.n << ",\n"
                 << "  \"axis_physical_period\": " << axis.physical_period << ",\n"
                 << "  \"diamond_L\": " << diamond.L << ",\n"
                 << "  \"diamond_sites\": " << diamond.n << ",\n"
                 << "  \"diamond_physical_period\": " << diamond.physical_period << ",\n"
                 << "  \"pell_residual\": " << (1LL * axis.n - diamond.n) << ",\n"
                 << "  \"samples\": " << options.samples << ",\n"
                 << "  \"batches\": " << options.batches << ",\n"
                 << "  \"p_ref\": " << options.p_ref << ",\n"
                 << "  \"h\": " << options.h << ",\n"
                 << "  \"seed\": " << options.seed << ",\n"
                 << "  \"rng\": \"stateless SplitMix64-derived counter mapping (seed, replica, site)\",\n"
                 << "  \"common_random_numbers\": true,\n"
                 << "  \"independent_replicas\": true,\n"
                 << "  \"effective_sample_size_per_probability\": " << options.samples << ",\n"
                 << "  \"autocorrelation_assumption\": \"none; replicas are independently counter-keyed\",\n"
                 << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
                 << "  \"threads\": " << runtime_threads() << ",\n"
                 << "  \"openmp\": "
#ifdef _OPENMP
                 << "true,\n"
#else
                 << "false,\n"
#endif
                 << "  \"elapsed_seconds\": " << elapsed << ",\n"
                 << "  \"batches_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
                 << "}\n";
        metadata.close();
        if (!metadata) throw std::runtime_error("failed while writing: " + json_path.string());

        std::cout << "wrote " << csv_path << '\n'
                  << "wrote " << json_path << '\n'
                  << "elapsed_seconds=" << std::setprecision(6) << elapsed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
