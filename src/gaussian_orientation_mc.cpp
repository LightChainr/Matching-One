// Same-N Gaussian orientation Monte Carlo (C01).
//
// Extends the PR #21 common-random-number engine and the C00 2x2 homology
// union-find.  Two primitive orientations of equal N share a cyclic Z/NZ
// occupation field.  Observables are the five wrapping channels on both the
// primal and matching graphs.
//
// Build:
//   g++ -O3 -std=c++17 -fopenmp src/gaussian_orientation_mc.cpp \
//       -o gaussian_orientation_mc
//
// Example:
//   ./gaussian_orientation_mc --rep1 8,1 --rep2 7,4 --t 1 \
//       --mode site --p 0.592746050790 --samples 200000 --batches 20 \
//       --seed 20260828 --replica-begin 0 --threads 8 \
//       --output-prefix results/server-20260828/C01/pilot_n65

#include "counter_rng.hpp"
#include "homology_union_find.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using matching::HomologyUnionFind;
using matching::PeriodMatrix;
using matching::WrappingChannels;
using matching::counter_uniform;
using matching::philox_kats_pass;

constexpr int kChannels = 5;
constexpr int kVars = 20;  // 2 orientations * (primal+matching) * 5 channels

struct Edge {
    int i;
    int j;
    int dx;
    int dy;
};

struct GaussianGeometry {
    int a = 0;
    int b = 0;
    int n = 0;
    PeriodMatrix period;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
};

int mod_n(int x, int n) {
    x %= n;
    return x < 0 ? x + n : x;
}

GaussianGeometry gaussian_geometry(int a, int b) {
    if (a <= 0 || b < 0) {
        throw std::invalid_argument("require a>0 and b>=0");
    }
    if (std::gcd(a, b) != 1) {
        throw std::invalid_argument("cyclic Gaussian representation requires gcd(a,b)=1");
    }
    const int n = a * a + b * b;
    if (n < 2) throw std::invalid_argument("N must be at least 2");
    GaussianGeometry g;
    g.a = a;
    g.b = b;
    g.n = n;
    g.period = PeriodMatrix::gaussian(a, b);
    const int ra = mod_n(a, n);
    const int rb = mod_n(b, n);
    const int rapb = mod_n(a + b, n);
    const int ramb = mod_n(a - b, n);
    g.primal_edges.reserve(static_cast<std::size_t>(2 * n));
    g.matching_edges.reserve(static_cast<std::size_t>(4 * n));
    for (int j = 0; j < n; ++j) {
        const Edge east{j, (j + ra) % n, 1, 0};
        const Edge north{j, (j + rb) % n, 0, 1};
        const Edge diag{j, (j + rapb) % n, 1, 1};
        const Edge anti{j, (j + ramb) % n, 1, -1};
        g.primal_edges.push_back(east);
        g.primal_edges.push_back(north);
        g.matching_edges.push_back(east);
        g.matching_edges.push_back(north);
        g.matching_edges.push_back(diag);
        g.matching_edges.push_back(anti);
    }
    return g;
}

WrappingChannels wrapping(HomologyUnionFind& uf, const std::vector<std::uint8_t>& active,
                          const std::vector<Edge>& edges) {
    uf.reset();
    for (const Edge& edge : edges) {
        if (active[static_cast<std::size_t>(edge.i)] &&
            active[static_cast<std::size_t>(edge.j)]) {
            uf.add_edge(edge.i, edge.j, edge.dx, edge.dy);
        }
    }
    return uf.channels();
}

void pack_pair(int* out, const WrappingChannels& primal, const WrappingChannels& matching) {
    for (int c = 0; c < kChannels; ++c) out[c] = primal.flag(c);
    for (int c = 0; c < kChannels; ++c) out[kChannels + c] = matching.flag(c);
}

std::uint32_t bond_index(int src, int dx, int dy) {
    const std::uint32_t code =
        static_cast<std::uint32_t>(dx + 2) * 5u + static_cast<std::uint32_t>(dy + 2);
    return static_cast<std::uint32_t>(src) * 32u + code;
}

WrappingChannels wrapping_bond(HomologyUnionFind& uf, const std::vector<Edge>& edges, int n,
                               int t, std::uint64_t seed, std::uint64_t replica, double p) {
    uf.reset();
    for (const Edge& edge : edges) {
        const int src = static_cast<int>((static_cast<std::int64_t>(t) * edge.i) % n);
        const std::uint32_t index = bond_index(src, edge.dx, edge.dy);
        if (counter_uniform(seed, replica, index, 1) < p) {
            uf.add_edge(edge.i, edge.j, edge.dx, edge.dy);
        }
    }
    return uf.channels();
}

void classify_site(const GaussianGeometry& geometry, HomologyUnionFind& uf,
                   const std::vector<double>& uniforms, int t, double p,
                   std::vector<std::uint8_t>& active, std::vector<std::uint8_t>& white,
                   int* packed10) {
    const int n = geometry.n;
    active.resize(static_cast<std::size_t>(n));
    white.resize(static_cast<std::size_t>(n));
    if (t == 1) {
        for (int j = 0; j < n; ++j) {
            active[static_cast<std::size_t>(j)] = uniforms[static_cast<std::size_t>(j)] < p;
        }
    } else {
        for (int j = 0; j < n; ++j) {
            const int src = static_cast<int>((static_cast<std::int64_t>(t) * j) % n);
            active[static_cast<std::size_t>(j)] =
                uniforms[static_cast<std::size_t>(src)] < p;
        }
    }
    for (int j = 0; j < n; ++j) {
        white[static_cast<std::size_t>(j)] =
            static_cast<std::uint8_t>(!active[static_cast<std::size_t>(j)]);
    }
    const WrappingChannels primal = wrapping(uf, active, geometry.primal_edges);
    const WrappingChannels matching = wrapping(uf, white, geometry.matching_edges);
    pack_pair(packed10, primal, matching);
}

void classify_bond(const GaussianGeometry& geometry, HomologyUnionFind& uf, int t,
                   std::uint64_t seed, std::uint64_t replica, double p, int* packed10) {
    const WrappingChannels primal = wrapping_bond(
        uf, geometry.primal_edges, geometry.n, t, seed, replica, p);
    const WrappingChannels matching = wrapping_bond(
        uf, geometry.matching_edges, geometry.n, t, seed, replica, p);
    pack_pair(packed10, primal, matching);
}

void require_equal_counts(const std::string& label, const std::vector<int>& actual,
                          const std::vector<int>& expected) {
    if (actual != expected) {
        std::ostringstream message;
        message << label << " exhaustive counts failed\nactual:";
        for (int value : actual) message << ' ' << value;
        message << "\nexpected:";
        for (int value : expected) message << ' ' << value;
        throw std::runtime_error(message.str());
    }
}

std::vector<int> exhaustive_primal_counts(const GaussianGeometry& geometry) {
    if (geometry.n > 16) {
        throw std::invalid_argument("exhaustive enumeration is limited to N<=16");
    }
    HomologyUnionFind uf(geometry.n, geometry.period);
    std::vector<std::uint8_t> active(static_cast<std::size_t>(geometry.n));
    std::vector<int> counts(5, 0);  // rank0, rank1, rank2, d0, d1
    const std::uint64_t n_config = 1ULL << geometry.n;
    for (std::uint64_t mask = 0; mask < n_config; ++mask) {
        for (int i = 0; i < geometry.n; ++i) {
            active[static_cast<std::size_t>(i)] = (mask >> i) & 1U;
        }
        const WrappingChannels channels = wrapping(uf, active, geometry.primal_edges);
        counts[channels.max_rank] += 1;
        counts[3] += static_cast<int>(channels.direction_0);
        counts[4] += static_cast<int>(channels.direction_1);
        if (channels.either != (channels.max_rank > 0)) {
            throw std::runtime_error("either/rank inconsistency");
        }
        if (channels.cross && channels.max_rank != 2) {
            throw std::runtime_error("cross/rank inconsistency");
        }
        if (channels.both && !(channels.direction_0 && channels.direction_1)) {
            throw std::runtime_error("both/direction inconsistency");
        }
    }
    return counts;
}

void self_test() {
    if (!philox_kats_pass()) {
        throw std::runtime_error("Philox4x32-10 KAT or counter-uniform regression failed");
    }

    const auto g21 = gaussian_geometry(2, 1);
    HomologyUnionFind uf(1, g21.period);
    const auto w1 = uf.winding(2, 1);
    const auto w2 = uf.winding(-1, 2);
    if (w1.x != 1 || w1.y != 0 || w2.x != 0 || w2.y != 1) {
        throw std::runtime_error("Gaussian (2,1) generator windings are wrong");
    }
    require_equal_counts("gaussian (2,1)", exhaustive_primal_counts(g21),
                         {16, 10, 6, 11, 11});
    require_equal_counts("gaussian (3,2)", exhaustive_primal_counts(gaussian_geometry(3, 2)),
                         {4629, 2340, 1223, 2471, 2471});

    // Occupying every site must wrap in both generators (full lattice).
    HomologyUnionFind full(g21.n, g21.period);
    std::vector<std::uint8_t> all(static_cast<std::size_t>(g21.n), 1);
    const WrappingChannels occupied = wrapping(full, all, g21.primal_edges);
    if (!occupied.cross || !occupied.either || !occupied.both) {
        throw std::runtime_error("fully occupied N=5 torus must be cross wrapping");
    }

    std::cout << "self-test passed: Philox KAT; Gaussian (2,1)/(3,2) exhaustive; windings\n";
}

struct Options {
    int a1 = 8;
    int b1 = 1;
    int a2 = 7;
    int b2 = 4;
    std::vector<int> t_values{1};
    std::string mode = "site";
    double p = 0.592746050790;
    std::uint64_t samples = 10000;
    int batches = 20;
    std::uint64_t seed = 20260828;
    std::uint64_t replica_begin = 0;
    int threads = 0;
    std::filesystem::path output_prefix;
    bool self_test = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --rep1 a,b           first Gaussian orientation (default 8,1)\n"
        << "  --rep2 a,b           second Gaussian orientation (default 7,4)\n"
        << "  --t t[,t...]         unit multipliers gcd(t,N)=1 (default 1)\n"
        << "  --mode site|bond     occupation type (default site)\n"
        << "  --p P                occupation probability (default 0.592746050790)\n"
        << "  --samples N          independent replicas (default 10000)\n"
        << "  --batches B          equal batches (default 20)\n"
        << "  --seed S             u64 Philox key (default 20260828)\n"
        << "  --replica-begin R    first replica counter (default 0)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --output-prefix PATH writes PATH.t*.batches.csv and metadata\n"
        << "  --self-test          exact tiny-torus regressions and exit\n"
        << "  --help               show this help\n";
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

std::pair<int, int> parse_pair(const std::string& text, const std::string& option) {
    const auto comma = text.find(',');
    if (comma == std::string::npos) {
        throw std::invalid_argument("expected a,b for " + option);
    }
    return {parse_number<int>(text.substr(0, comma), option),
            parse_number<int>(text.substr(comma + 1), option)};
}

std::vector<int> parse_t_list(const std::string& text) {
    std::vector<int> values;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) continue;
        values.push_back(parse_number<int>(item, "--t"));
    }
    if (values.empty()) throw std::invalid_argument("--t must contain at least one integer");
    return values;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("missing value for " + arg);
            return argv[i];
        };
        if (arg == "--rep1") {
            const auto pair = parse_pair(value(), arg);
            options.a1 = pair.first;
            options.b1 = pair.second;
        } else if (arg == "--rep2") {
            const auto pair = parse_pair(value(), arg);
            options.a2 = pair.first;
            options.b2 = pair.second;
        } else if (arg == "--t") {
            options.t_values = parse_t_list(value());
        } else if (arg == "--mode") {
            options.mode = value();
        } else if (arg == "--p") {
            options.p = parse_number<double>(value(), arg);
        } else if (arg == "--samples") {
            options.samples = parse_number<std::uint64_t>(value(), arg);
        } else if (arg == "--batches") {
            options.batches = parse_number<int>(value(), arg);
        } else if (arg == "--seed") {
            options.seed = parse_number<std::uint64_t>(value(), arg);
        } else if (arg == "--replica-begin") {
            options.replica_begin = parse_number<std::uint64_t>(value(), arg);
        } else if (arg == "--threads") {
            options.threads = parse_number<int>(value(), arg);
        } else if (arg == "--output-prefix") {
            options.output_prefix = value();
        } else if (arg == "--self-test") {
            options.self_test = true;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0], 0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
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

int runtime_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

struct Accumulator {
    std::vector<long long> sum;
    std::vector<long long> gram;
    Accumulator() : sum(kVars, 0), gram(kVars * kVars, 0) {}
    void add(const int* x) {
        for (int i = 0; i < kVars; ++i) {
            if (x[i] == 0) continue;
            sum[static_cast<std::size_t>(i)] += 1;
            for (int j = 0; j < kVars; ++j) {
                gram[static_cast<std::size_t>(i * kVars + j)] += x[j];
            }
        }
    }
    void add_from(const Accumulator& other) {
        for (std::size_t i = 0; i < sum.size(); ++i) sum[i] += other.sum[i];
        for (std::size_t i = 0; i < gram.size(); ++i) gram[i] += other.gram[i];
    }
};

const char* kVarNames[kVars] = {
    "o1_primal_cross", "o1_primal_both", "o1_primal_either", "o1_primal_direction_0",
    "o1_primal_direction_1", "o1_matching_cross", "o1_matching_both", "o1_matching_either",
    "o1_matching_direction_0", "o1_matching_direction_1", "o2_primal_cross", "o2_primal_both",
    "o2_primal_either", "o2_primal_direction_0", "o2_primal_direction_1", "o2_matching_cross",
    "o2_matching_both", "o2_matching_either", "o2_matching_direction_0", "o2_matching_direction_1",
};

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
        if (options.mode != "site" && options.mode != "bond") {
            throw std::invalid_argument("--mode must be site or bond");
        }
        if (options.batches < 2) throw std::invalid_argument("--batches must be at least 2");
        if (options.samples < static_cast<std::uint64_t>(options.batches)) {
            throw std::invalid_argument("--samples must be at least --batches");
        }
        if (options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
            throw std::invalid_argument("--samples must be divisible by --batches");
        }
        if (!(options.p > 0.0) || !(options.p < 1.0)) {
            throw std::invalid_argument("require 0 < p < 1");
        }
        if (options.threads < 0) throw std::invalid_argument("--threads cannot be negative");
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
#else
        if (options.threads > 1) {
            throw std::invalid_argument("binary was compiled without OpenMP support");
        }
#endif

        const GaussianGeometry geom1 = gaussian_geometry(options.a1, options.b1);
        const GaussianGeometry geom2 = gaussian_geometry(options.a2, options.b2);
        if (geom1.n != geom2.n) {
            throw std::invalid_argument("same-N pairing requires a1^2+b1^2 = a2^2+b2^2");
        }
        const int n = geom1.n;
        for (int t : options.t_values) {
            if (t <= 0 || std::gcd(t, n) != 1) {
                throw std::invalid_argument("each t must satisfy t>0 and gcd(t,N)=1");
            }
        }

        const std::uint64_t per_batch = options.samples / static_cast<std::uint64_t>(options.batches);
        const int n_t = static_cast<int>(options.t_values.size());
        std::vector<std::vector<Accumulator>> batch_acc(
            static_cast<std::size_t>(options.batches),
            std::vector<Accumulator>(static_cast<std::size_t>(n_t)));

        const auto started = std::chrono::steady_clock::now();
        for (int batch = 0; batch < options.batches; ++batch) {
            const std::uint64_t first =
                options.replica_begin + static_cast<std::uint64_t>(batch) * per_batch;
            const int thread_count = runtime_threads();
            std::vector<std::vector<Accumulator>> thread_acc(
                static_cast<std::size_t>(thread_count),
                std::vector<Accumulator>(static_cast<std::size_t>(n_t)));

#pragma omp parallel
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                HomologyUnionFind uf1(n, geom1.period);
                HomologyUnionFind uf2(n, geom2.period);
                std::vector<double> uniforms(static_cast<std::size_t>(n));
                std::vector<std::uint8_t> active(static_cast<std::size_t>(n));
                std::vector<std::uint8_t> white(static_cast<std::size_t>(n));
                int packed[kVars];

#pragma omp for schedule(static)
                for (std::uint64_t offset = 0; offset < per_batch; ++offset) {
                    const std::uint64_t replica = first + offset;
                    if (options.mode == "site") {
                        for (int site = 0; site < n; ++site) {
                            uniforms[static_cast<std::size_t>(site)] =
                                counter_uniform(options.seed, replica, static_cast<std::uint32_t>(site), 0);
                        }
                    }
                    int first10[10];
                    if (options.mode == "site") {
                        classify_site(geom1, uf1, uniforms, 1, options.p, active, white, first10);
                    } else {
                        classify_bond(geom1, uf1, 1, options.seed, replica, options.p, first10);
                    }
                    for (int ti = 0; ti < n_t; ++ti) {
                        const int t = options.t_values[static_cast<std::size_t>(ti)];
                        int second10[10];
                        if (options.mode == "site") {
                            classify_site(geom2, uf2, uniforms, t, options.p, active, white, second10);
                        } else {
                            classify_bond(geom2, uf2, t, options.seed, replica, options.p, second10);
                        }
                        for (int i = 0; i < 10; ++i) packed[i] = first10[i];
                        for (int i = 0; i < 10; ++i) packed[10 + i] = second10[i];
                        thread_acc[static_cast<std::size_t>(tid)][static_cast<std::size_t>(ti)].add(packed);
                    }
                }
            }

            for (int ti = 0; ti < n_t; ++ti) {
                for (int thread = 0; thread < thread_count; ++thread) {
                    batch_acc[static_cast<std::size_t>(batch)][static_cast<std::size_t>(ti)].add_from(
                        thread_acc[static_cast<std::size_t>(thread)][static_cast<std::size_t>(ti)]);
                }
            }
            std::cerr << "completed batch " << (batch + 1) << '/' << options.batches
                      << " mode=" << options.mode << " N=" << n << '\n';
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        const std::filesystem::path parent = options.output_prefix.parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);

        std::ostringstream t_joined;
        for (std::size_t i = 0; i < options.t_values.size(); ++i) {
            if (i) t_joined << ',';
            t_joined << options.t_values[i];
        }

        for (int ti = 0; ti < n_t; ++ti) {
            const int t = options.t_values[static_cast<std::size_t>(ti)];
            const std::string stem = options.output_prefix.string() + ".t" + std::to_string(t);
            const std::filesystem::path csv_path = stem + ".batches.csv";
            const std::filesystem::path json_path = stem + ".moments.json";

            std::ofstream csv(csv_path);
            if (!csv) throw std::runtime_error("cannot open output: " + csv_path.string());
            csv << "batch,first_replica,samples,t";
            for (int i = 0; i < kVars; ++i) csv << ',' << kVarNames[i];
            csv << '\n';
            for (int batch = 0; batch < options.batches; ++batch) {
                const std::uint64_t first =
                    options.replica_begin + static_cast<std::uint64_t>(batch) * per_batch;
                const Accumulator& acc = batch_acc[static_cast<std::size_t>(batch)][static_cast<std::size_t>(ti)];
                csv << batch << ',' << first << ',' << per_batch << ',' << t;
                for (int i = 0; i < kVars; ++i) csv << ',' << acc.sum[static_cast<std::size_t>(i)];
                csv << '\n';
            }
            csv.close();
            if (!csv) throw std::runtime_error("failed while writing: " + csv_path.string());

            std::ofstream json(json_path);
            if (!json) throw std::runtime_error("cannot open output: " + json_path.string());
            json << "{\n  \"t\": " << t << ",\n  \"variables\": [";
            for (int i = 0; i < kVars; ++i) {
                if (i) json << ", ";
                json << '"' << kVarNames[i] << '"';
            }
            json << "],\n  \"batches\": [\n";
            for (int batch = 0; batch < options.batches; ++batch) {
                const std::uint64_t first =
                    options.replica_begin + static_cast<std::uint64_t>(batch) * per_batch;
                const Accumulator& acc = batch_acc[static_cast<std::size_t>(batch)][static_cast<std::size_t>(ti)];
                json << "    {\"batch\": " << batch << ", \"first_replica\": " << first
                     << ", \"samples\": " << per_batch << ", \"sum\": [";
                for (int i = 0; i < kVars; ++i) {
                    if (i) json << ", ";
                    json << acc.sum[static_cast<std::size_t>(i)];
                }
                json << "], \"gram\": [";
                for (int i = 0; i < kVars * kVars; ++i) {
                    if (i) json << ", ";
                    json << acc.gram[static_cast<std::size_t>(i)];
                }
                json << "]}";
                if (batch + 1 != options.batches) json << ',';
                json << '\n';
            }
            json << "  ]\n}\n";
            json.close();
            if (!json) throw std::runtime_error("failed while writing: " + json_path.string());
            std::cout << "wrote " << csv_path << '\n' << "wrote " << json_path << '\n';
        }

        const std::filesystem::path meta_path = options.output_prefix.string() + ".metadata.json";
        std::ofstream metadata(meta_path);
        if (!metadata) throw std::runtime_error("cannot open output: " + meta_path.string());
        metadata << std::setprecision(17)
                 << "{\n"
                 << "  \"engine\": \"gaussian_orientation_mc_c01_v1\",\n"
                 << "  \"design\": \"same-N Gaussian CRN orientation scan\",\n"
                 << "  \"mode\": \"" << json_escape(options.mode) << "\",\n"
                 << "  \"N\": " << n << ",\n"
                 << "  \"rep1\": [" << geom1.a << ", " << geom1.b << "],\n"
                 << "  \"rep2\": [" << geom2.a << ", " << geom2.b << "],\n"
                 << "  \"period_matrix_1\": [[" << geom1.period.a00 << ", " << geom1.period.a01
                 << "], [" << geom1.period.a10 << ", " << geom1.period.a11 << "]],\n"
                 << "  \"period_matrix_2\": [[" << geom2.period.a00 << ", " << geom2.period.a01
                 << "], [" << geom2.period.a10 << ", " << geom2.period.a11 << "]],\n"
                 << "  \"t_values\": [" << t_joined.str() << "],\n"
                 << "  \"coupling\": \"U_j^{(1)}=U_j; U_j^{(2)}=U_{t j mod N}\",\n"
                 << "  \"p\": " << options.p << ",\n"
                 << "  \"samples\": " << options.samples << ",\n"
                 << "  \"batches\": " << options.batches << ",\n"
                 << "  \"seed\": " << options.seed << ",\n"
                 << "  \"replica_begin\": " << options.replica_begin << ",\n"
                 << "  \"replica_end\": " << (options.replica_begin + options.samples) << ",\n"
                 << "  \"rng\": \"Philox4x32-10 counter uniform; key=(seed_lo, seed_hi); "
                    "ctr=(index, replica_lo, replica_hi, stream); stream 0 site, stream 1 bond\",\n"
                 << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
                 << "  \"threads\": " << runtime_threads() << ",\n"
                 << "  \"openmp\": "
#ifdef _OPENMP
                 << "true,\n"
#else
                 << "false,\n"
#endif
                 << "  \"elapsed_seconds\": " << elapsed << "\n"
                 << "}\n";
        metadata.close();
        if (!metadata) throw std::runtime_error("failed while writing: " + meta_path.string());

        std::cout << "wrote " << meta_path << '\n'
                  << "elapsed_seconds=" << std::setprecision(6) << elapsed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
