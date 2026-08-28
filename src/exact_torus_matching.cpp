#include "torus_connectivity.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using matching::DisplacementDSU;
using matching::Graph;
using matching::Microcanonical;
using matching::Observables;
using matching::analyze_config;
using matching::euler_identity_holds;
using matching::make_graph;
using Big = boost::multiprecision::cpp_int;

struct Options {
    std::vector<int> sizes;
    int threads = 8;
    std::string outdir = "results/issue-7";
    bool dump_configs = false;
    bool check_euler = true;
};

void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0
              << " [--L N]... [--all] [--threads N] [--outdir DIR] [--dump-configs]\n";
}

bool parse_args(int argc, char** argv, Options& opt) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };
        if (a == "--L") {
            const char* v = need("--L");
            if (!v) {
                return false;
            }
            opt.sizes.push_back(std::stoi(v));
        } else if (a == "--all") {
            opt.sizes = {2, 3, 4, 5};
        } else if (a == "--threads") {
            const char* v = need("--threads");
            if (!v) {
                return false;
            }
            opt.threads = std::stoi(v);
        } else if (a == "--outdir") {
            const char* v = need("--outdir");
            if (!v) {
                return false;
            }
            opt.outdir = v;
        } else if (a == "--dump-configs") {
            opt.dump_configs = true;
        } else if (a == "--no-euler") {
            opt.check_euler = false;
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return false;
        } else {
            std::cerr << "unknown argument: " << a << "\n";
            usage(argv[0]);
            return false;
        }
    }
    if (opt.sizes.empty()) {
        std::cerr << "specify --L N or --all\n";
        return false;
    }
    return true;
}

long peak_rss_kb() {
    std::ifstream in("/proc/self/status");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("VmHWM:", 0) == 0) {
            std::istringstream iss(line.substr(6));
            long kb = 0;
            iss >> kb;
            return kb;
        }
    }
    return 0;
}

double cpu_seconds() {
    struct rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_utime.tv_sec) + static_cast<double>(u.ru_utime.tv_usec) * 1e-6 +
           static_cast<double>(u.ru_stime.tv_sec) + static_cast<double>(u.ru_stime.tv_usec) * 1e-6;
}

std::string pad_L(int L) {
    std::ostringstream oss;
    oss << 'L' << std::setw(2) << std::setfill('0') << L;
    return oss.str();
}

std::vector<Big> bernstein_to_monomial(const std::vector<uint64_t>& A, int n) {
    std::vector<Big> coeff(static_cast<std::size_t>(n + 1), 0);
    for (int k = 0; k <= n; ++k) {
        if (A[static_cast<std::size_t>(k)] == 0) {
            continue;
        }
        const Big ak = A[static_cast<std::size_t>(k)];
        Big binom = 1;
        const int nk = n - k;
        for (int j = 0; j <= nk; ++j) {
            if (j > 0) {
                binom *= (nk - j + 1);
                binom /= j;
            }
            Big term = ak * binom;
            if (j & 1) {
                term = -term;
            }
            coeff[static_cast<std::size_t>(k + j)] += term;
        }
    }
    return coeff;
}

std::vector<Big> sub_poly(const std::vector<Big>& a, const std::vector<Big>& b) {
    std::vector<Big> out = a;
    if (b.size() > out.size()) {
        out.resize(b.size(), 0);
    }
    for (std::size_t i = 0; i < b.size(); ++i) {
        out[i] -= b[i];
    }
    return out;
}

int first_mismatch(const std::vector<Big>& diff) {
    for (std::size_t i = 0; i < diff.size(); ++i) {
        if (diff[i] != 0) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int max_coeff_bit_length(const std::vector<Big>& c) {
    int bits = 0;
    for (const Big& v : c) {
        if (v == 0) {
            continue;
        }
        Big a = v < 0 ? -v : v;
        int b = 0;
        while (a > 0) {
            a >>= 1;
            ++b;
        }
        if (b > bits) {
            bits = b;
        }
    }
    return bits;
}

std::string coeffs_json(const std::vector<Big>& c) {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < c.size(); ++i) {
        if (i) {
            oss << ", ";
        }
        oss << '"' << c[i] << '"';
    }
    oss << "]";
    return oss.str();
}

void write_microcanonical_csv(const std::string& path, const Microcanonical& mc) {
    std::ofstream out(path);
    out << "k,configuration_count,sum_clusters_G,sum_clusters_Gstar_complement,"
           "wrap_H_G,wrap_V_G,wrap_E_G,wrap_B_G,"
           "wrap_H_Gstar,wrap_V_Gstar,wrap_E_Gstar,wrap_B_Gstar\n";
    for (int k = 0; k <= mc.n; ++k) {
        const std::size_t i = static_cast<std::size_t>(k);
        out << k << ',' << mc.configuration_count[i] << ',' << mc.sum_clusters_G[i] << ','
            << mc.sum_clusters_Gstar_complement[i] << ',' << mc.wrap_H_G[i] << ',' << mc.wrap_V_G[i]
            << ',' << mc.wrap_E_G[i] << ',' << mc.wrap_B_G[i] << ',' << mc.wrap_H_Gstar[i] << ','
            << mc.wrap_V_Gstar[i] << ',' << mc.wrap_E_Gstar[i] << ',' << mc.wrap_B_Gstar[i] << '\n';
    }
}

struct IdentityReport {
    std::string overall = "PASS";
    std::vector<Big> M;
    int degree = 0;
    int max_bits = 0;
    int first_mismatch_any = -1;
    std::string first_failing_class;
};

IdentityReport write_identity_json(const std::string& path, const Microcanonical& mc) {
    const int n = mc.n;
    auto N_poly = bernstein_to_monomial(mc.sum_clusters_G, n);
    auto Nhat_poly = bernstein_to_monomial(mc.sum_clusters_Gstar_complement, n);
    std::vector<Big> chi(static_cast<std::size_t>(n + 1), 0);
    chi[1] += n;
    chi[2] -= 2 * n;
    if (n >= 4) {
        chi[4] += n;
    }
    auto M = sub_poly(sub_poly(N_poly, Nhat_poly), chi);

    struct Item {
        const char* name;
        const std::vector<uint64_t>* g;
        const std::vector<uint64_t>* gs;
    };
    const Item items[] = {
        {"H", &mc.wrap_H_G, &mc.wrap_H_Gstar},
        {"V", &mc.wrap_V_G, &mc.wrap_V_Gstar},
        {"E", &mc.wrap_E_G, &mc.wrap_E_Gstar},
        {"B", &mc.wrap_B_G, &mc.wrap_B_Gstar},
    };

    IdentityReport rep;
    rep.M = M;
    rep.degree = n;
    rep.max_bits = max_coeff_bit_length(M);

    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"L\": " << mc.L << ",\n";
    oss << "  \"N\": " << n << ",\n";
    oss << "  \"degree\": " << n << ",\n";
    oss << "  \"coefficients\": " << coeffs_json(M) << ",\n";
    oss << "  \"max_coefficient_bit_length\": " << rep.max_bits << ",\n";
    oss << "  \"wrapping_convention\": {\n";
    oss << "    \"h\": \"horizontal (nonzero x winding)\",\n";
    oss << "    \"c\": \"vertical/column (nonzero y winding)\",\n";
    oss << "    \"e\": \"either (H or V)\",\n";
    oss << "    \"b\": \"both (H and V)\"\n";
    oss << "  },\n";
    oss << "  \"cluster_form\": "
           "\"M_L(p) = N_L(p) - Nhat_L(1-p) - L^2 (p - 2 p^2 + p^4)\",\n";
    oss << "  \"identities\": {\n";

    bool overall_pass = true;
    bool first_item = true;
    for (const Item& it : items) {
        auto g_poly = bernstein_to_monomial(*it.g, n);
        auto gs_poly = bernstein_to_monomial(*it.gs, n);
        auto diff_wrap = sub_poly(g_poly, gs_poly);
        auto residual = sub_poly(M, diff_wrap);
        const int mm = first_mismatch(residual);
        const char* status = (mm < 0) ? "PASS" : "FAIL";
        if (mm >= 0) {
            overall_pass = false;
            if (rep.first_mismatch_any < 0) {
                rep.first_mismatch_any = mm;
                rep.first_failing_class = it.name;
            }
        }
        if (!first_item) {
            oss << ",\n";
        }
        first_item = false;
        oss << "    \"" << it.name << "\": {\n";
        oss << "      \"identity\": \"R_G^" << it.name << "(p) - R_Gstar^" << it.name
            << "(1-p)\",\n";
        oss << "      \"status\": \"" << status << "\",\n";
        oss << "      \"first_mismatching_coefficient\": ";
        if (mm < 0) {
            oss << "null";
        } else {
            oss << mm;
        }
        oss << ",\n";
        oss << "      \"coefficients\": " << coeffs_json(diff_wrap) << ",\n";
        oss << "      \"residual\": " << coeffs_json(residual) << "\n";
        oss << "    }";
    }
    oss << "\n  },\n";
    oss << "  \"euler_mismatches\": " << mc.euler_mismatches << ",\n";
    if (mc.euler_mismatches > 0) {
        overall_pass = false;
        oss << "  \"first_euler_mask\": " << mc.first_euler_mask << ",\n";
    }
    oss << "  \"identity\": \"" << (overall_pass ? "PASS" : "FAIL") << "\",\n";
    oss << "  \"first_mismatching_coefficient\": ";
    if (rep.first_mismatch_any < 0) {
        oss << "null";
    } else {
        oss << rep.first_mismatch_any;
    }
    oss << "\n}\n";

    std::ofstream out(path);
    out << oss.str();
    rep.overall = overall_pass ? "PASS" : "FAIL";
    return rep;
}

void write_counterexample(const std::string& dir, int L, uint64_t mask, const Observables& o,
                          const Graph& g, const Graph& gstar, const IdentityReport& ident) {
    std::ostringstream path;
    path << dir << "/counterexamples/" << pad_L(L) << "_mask" << mask << ".json";
    std::ofstream out(path.str());
    const int n = L * L;
    out << "{\n";
    out << "  \"L\": " << L << ",\n";
    out << "  \"configuration_bit_pattern\": " << mask << ",\n";
    out << "  \"k\": " << o.k << ",\n";
    out << "  \"occupied_coordinates\": [";
    bool first = true;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = y * L + x;
            if (!matching::bit_set(mask, i)) {
                continue;
            }
            if (!first) {
                out << ", ";
            }
            first = false;
            out << "[" << x << ", " << y << "]";
        }
    }
    out << "],\n";
    out << "  \"matching_coordinates\": [";
    first = true;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = y * L + x;
            if (matching::bit_set(mask, i)) {
                continue;
            }
            if (!first) {
                out << ", ";
            }
            first = false;
            out << "[" << x << ", " << y << "]";
        }
    }
    out << "],\n";
    out << "  \"edge_list_G\": [";
    first = true;
    for (int i = 0; i < n; ++i) {
        if (!matching::bit_set(mask, i)) {
            continue;
        }
        for (const auto& e : g.adj[static_cast<std::size_t>(i)]) {
            if (!matching::bit_set(mask, e.to)) {
                continue;
            }
            if (!first) {
                out << ", ";
            }
            first = false;
            out << "{\"from\": " << i << ", \"to\": " << e.to << ", \"dx\": " << e.dx
                << ", \"dy\": " << e.dy << "}";
        }
    }
    out << "],\n";
    out << "  \"edge_list_Gstar\": [";
    first = true;
    for (int i = 0; i < n; ++i) {
        if (matching::bit_set(mask, i)) {
            continue;
        }
        for (const auto& e : gstar.adj[static_cast<std::size_t>(i)]) {
            if (matching::bit_set(mask, e.to)) {
                continue;
            }
            if (!first) {
                out << ", ";
            }
            first = false;
            out << "{\"from\": " << i << ", \"to\": " << e.to << ", \"dx\": " << e.dx
                << ", \"dy\": " << e.dy << "}";
        }
    }
    out << "],\n";
    out << "  \"DSU_winding_state\": {\n";
    out << "    \"clusters_G\": " << o.clusters_G << ",\n";
    out << "    \"clusters_Gstar\": " << o.clusters_Gstar << ",\n";
    out << "    \"H_G\": " << static_cast<int>(o.H_G) << ",\n";
    out << "    \"V_G\": " << static_cast<int>(o.V_G) << ",\n";
    out << "    \"E_G\": " << static_cast<int>(o.E_G) << ",\n";
    out << "    \"B_G\": " << static_cast<int>(o.B_G) << ",\n";
    out << "    \"H_Gstar\": " << static_cast<int>(o.H_Gstar) << ",\n";
    out << "    \"V_Gstar\": " << static_cast<int>(o.V_Gstar) << ",\n";
    out << "    \"E_Gstar\": " << static_cast<int>(o.E_Gstar) << ",\n";
    out << "    \"B_Gstar\": " << static_cast<int>(o.B_Gstar) << "\n";
    out << "  },\n";
    out << "  \"oracle_winding_state\": null,\n";
    out << "  \"first_mismatching_polynomial_coefficient\": ";
    if (ident.first_mismatch_any < 0) {
        out << "null";
    } else {
        out << ident.first_mismatch_any;
    }
    out << ",\n";
    out << "  \"first_failing_class\": \"" << ident.first_failing_class << "\"\n";
    out << "}\n";
    std::cerr << "wrote counterexample " << path.str() << "\n";
}

void dump_configs(const std::string& path, int L, const Graph& g, const Graph& gstar) {
    std::ofstream out(path);
    out << "mask,k,clusters_G,clusters_Gstar,H_G,V_G,E_G,B_G,H_Gstar,V_Gstar,E_Gstar,B_Gstar\n";
    const int n = L * L;
    const uint64_t total = 1ull << n;
    DisplacementDSU dsu(n);
    for (uint64_t mask = 0; mask < total; ++mask) {
        const Observables o = analyze_config(mask, g, gstar, dsu);
        out << mask << ',' << o.k << ',' << o.clusters_G << ',' << o.clusters_Gstar << ','
            << static_cast<int>(o.H_G) << ',' << static_cast<int>(o.V_G) << ','
            << static_cast<int>(o.E_G) << ',' << static_cast<int>(o.B_G) << ','
            << static_cast<int>(o.H_Gstar) << ',' << static_cast<int>(o.V_Gstar) << ','
            << static_cast<int>(o.E_Gstar) << ',' << static_cast<int>(o.B_Gstar) << '\n';
    }
}

Microcanonical enumerate_L(int L, int nthreads, bool check_euler) {
    const Graph g = make_graph(L, false);
    const Graph gstar = make_graph(L, true);
    const int n = L * L;
    const uint64_t total = 1ull << n;

#ifdef _OPENMP
    if (nthreads < 1) {
        nthreads = 1;
    }
    omp_set_num_threads(nthreads);
#else
    nthreads = 1;
#endif

    std::vector<Microcanonical> acc(static_cast<std::size_t>(nthreads), Microcanonical(L));

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
        const int nt = omp_get_num_threads();
#else
        const int tid = 0;
        const int nt = 1;
#endif
        const uint64_t begin = (total * static_cast<uint64_t>(tid)) / static_cast<uint64_t>(nt);
        const uint64_t end = (total * static_cast<uint64_t>(tid + 1)) / static_cast<uint64_t>(nt);
        DisplacementDSU dsu(n);
        Microcanonical& local = acc[static_cast<std::size_t>(tid)];
        for (uint64_t mask = begin; mask < end; ++mask) {
            const Observables o = analyze_config(mask, g, gstar, dsu);
            local.add(o);
            if (check_euler && !euler_identity_holds(L, mask, o, g)) {
                local.euler_mismatches += 1;
                if (mask < local.first_euler_mask) {
                    local.first_euler_mask = mask;
                }
            }
        }
    }

    Microcanonical sum(L);
    for (const Microcanonical& a : acc) {
        sum.absorb(a);
    }
    return sum;
}

void append_performance(const std::string& path, int L, uint64_t configs, double wall, double cpu,
                        long rss_kb) {
    const bool exists = static_cast<bool>(std::ifstream(path));
    std::ofstream out(path, std::ios::app);
    if (!exists) {
        out << "L,configurations,wall_seconds,cpu_seconds,peak_rss_kb,configs_per_sec\n";
    }
    const double cps = wall > 0.0 ? static_cast<double>(configs) / wall : 0.0;
    out << L << ',' << configs << ',' << std::setprecision(10) << wall << ',' << cpu << ',' << rss_kb
        << ',' << cps << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        return 2;
    }

#ifdef _OPENMP
    omp_set_num_threads(opt.threads);
#endif

    int rc = 0;
    for (int L : opt.sizes) {
        if (L < 2 || L > 6) {
            std::cerr << "unsupported L=" << L << " (engine allows 2..6)\n";
            rc = 1;
            continue;
        }
        const int n = L * L;
        const uint64_t total = 1ull << n;
        std::cerr << "enumerating L=" << L << " N=" << n << " configs=" << total
                  << " threads=" << opt.threads << "\n";

        const double cpu0 = cpu_seconds();
        const auto t0 = std::chrono::steady_clock::now();
        Microcanonical mc = enumerate_L(L, opt.threads, opt.check_euler);
        const auto t1 = std::chrono::steady_clock::now();
        const double cpu1 = cpu_seconds();
        const double wall = std::chrono::duration<double>(t1 - t0).count();
        const double cpu = cpu1 - cpu0;
        const long rss = peak_rss_kb();

        const std::string tag = pad_L(L);
        const std::string csv = opt.outdir + "/" + tag + "_microcanonical.csv";
        const std::string js = opt.outdir + "/" + tag + "_identity.json";
        write_microcanonical_csv(csv, mc);
        const IdentityReport ident = write_identity_json(js, mc);
        append_performance(opt.outdir + "/performance.csv", L, total, wall, cpu, rss);

        if (opt.dump_configs && L <= 4) {
            const Graph g = make_graph(L, false);
            const Graph gstar = make_graph(L, true);
            dump_configs(opt.outdir + "/" + tag + "_configs.csv", L, g, gstar);
        }

        if (ident.overall != "PASS" || mc.euler_mismatches > 0) {
            rc = 1;
            const Graph g = make_graph(L, false);
            const Graph gstar = make_graph(L, true);
            DisplacementDSU dsu(n);
            uint64_t mask = (mc.euler_mismatches > 0) ? mc.first_euler_mask : 0ull;
            if (mc.euler_mismatches == 0) {
                for (uint64_t m = 0; m < total; ++m) {
                    const Observables o = analyze_config(m, g, gstar, dsu);
                    if (!euler_identity_holds(L, m, o, g)) {
                        mask = m;
                        break;
                    }
                }
            }
            const Observables o = analyze_config(mask, g, gstar, dsu);
            write_counterexample(opt.outdir, L, mask, o, g, gstar, ident);
            std::cerr << tag << " FAIL identity=" << ident.overall
                      << " first_mismatch=" << ident.first_mismatch_any
                      << " class=" << ident.first_failing_class
                      << " euler_mismatches=" << mc.euler_mismatches << "\n";
        } else {
            std::cerr << tag << " PASS wall=" << wall << "s cpu=" << cpu
                      << "s configs/sec=" << (wall > 0 ? total / wall : 0.0) << " rss_kb=" << rss
                      << "\n";
        }
    }
    return rc;
}
