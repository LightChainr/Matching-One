#include "counter_rng.hpp"
#include "torus_connectivity.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using matching::DisplacementDSU;
using matching::Edge;
using matching::Graph;
using matching::Observables;
using matching::PhiloxStream;
using matching::analyze_config;
using matching::fisher_yates;
using matching::make_graph;
using matching::philox4x32_10;
using matching::philox4x32_10_official_kats;

constexpr uint32_t kMasterSeed = 0x4D314E5Au;  // "M1NZ"

struct CSR {
    int n = 0;
    std::vector<int> off;
    std::vector<int> to;
    std::vector<int> dx;
    std::vector<int> dy;
};

CSR make_undirected_csr(const Graph& g) {
    std::vector<std::vector<Edge>> full = g.adj;
    for (int i = 0; i < g.n; ++i) {
        for (const Edge& e : g.adj[static_cast<std::size_t>(i)]) {
            full[static_cast<std::size_t>(e.to)].push_back(Edge{i, -e.dx, -e.dy});
        }
    }
    CSR c;
    c.n = g.n;
    c.off.assign(static_cast<std::size_t>(g.n + 1), 0);
    int m = 0;
    for (int i = 0; i < g.n; ++i) {
        c.off[static_cast<std::size_t>(i)] = m;
        m += static_cast<int>(full[static_cast<std::size_t>(i)].size());
    }
    c.off[static_cast<std::size_t>(g.n)] = m;
    c.to.resize(static_cast<std::size_t>(m));
    c.dx.resize(static_cast<std::size_t>(m));
    c.dy.resize(static_cast<std::size_t>(m));
    for (int i = 0; i < g.n; ++i) {
        int z = c.off[static_cast<std::size_t>(i)];
        for (const Edge& e : full[static_cast<std::size_t>(i)]) {
            c.to[static_cast<std::size_t>(z)] = e.to;
            c.dx[static_cast<std::size_t>(z)] = e.dx;
            c.dy[static_cast<std::size_t>(z)] = e.dy;
            ++z;
        }
    }
    return c;
}

struct Engine {
    int n = 0;
    DisplacementDSU dsu;
    std::vector<uint8_t> occ;
    const CSR* csr = nullptr;
    int n_clusters = 0;
    uint8_t wrap_h = 0;
    uint8_t wrap_v = 0;

    void setup(int n_sites) {
        n = n_sites;
        dsu.reset(n);
        occ.assign(static_cast<std::size_t>(n), 0);
    }

    void begin(const CSR* c) {
        csr = c;
        dsu.reset(n);
        std::memset(occ.data(), 0, static_cast<std::size_t>(n));
        n_clusters = 0;
        wrap_h = 0;
        wrap_v = 0;
    }

    void add(int i) {
        occ[static_cast<std::size_t>(i)] = 1;
        n_clusters += 1;
        const int z0 = csr->off[static_cast<std::size_t>(i)];
        const int z1 = csr->off[static_cast<std::size_t>(i + 1)];
        for (int z = z0; z < z1; ++z) {
            const int j = csr->to[static_cast<std::size_t>(z)];
            if (!occ[static_cast<std::size_t>(j)]) {
                continue;
            }
            const auto a = dsu.find(i);
            const auto b = dsu.find(j);
            if (a.root != b.root) {
                n_clusters -= 1;
            }
            dsu.unite(i, j, csr->dx[static_cast<std::size_t>(z)], csr->dy[static_cast<std::size_t>(z)]);
            const auto r = dsu.find(i);
            wrap_h = static_cast<uint8_t>(wrap_h | dsu.wrap_h[static_cast<std::size_t>(r.root)]);
            wrap_v = static_cast<uint8_t>(wrap_v | dsu.wrap_v[static_cast<std::size_t>(r.root)]);
        }
    }
};

struct Acc {
    int n = 0;
    uint64_t replicas = 0;
    std::vector<uint64_t> cl_g, cl_gs;
    std::vector<uint64_t> h_g, v_g, e_g, b_g;
    std::vector<uint64_t> h_gs, v_gs, e_gs, b_gs;
    std::vector<uint64_t> cl_g2, cl_gs2, cl_xy;
    std::vector<uint64_t> h_xy, v_xy, e_xy, b_xy;

    explicit Acc(int n_sites = 0) { reset(n_sites); }

    void reset(int n_sites) {
        n = n_sites;
        replicas = 0;
        const std::size_t m = static_cast<std::size_t>(n + 1);
        auto z = [&](std::vector<uint64_t>& v) { v.assign(m, 0); };
        z(cl_g);
        z(cl_gs);
        z(h_g);
        z(v_g);
        z(e_g);
        z(b_g);
        z(h_gs);
        z(v_gs);
        z(e_gs);
        z(b_gs);
        z(cl_g2);
        z(cl_gs2);
        z(cl_xy);
        z(h_xy);
        z(v_xy);
        z(e_xy);
        z(b_xy);
    }

    void add_k(int k, int cg, int cgs, uint8_t hg, uint8_t vg, uint8_t hgs, uint8_t vgs) {
        const uint8_t eg = static_cast<uint8_t>(hg | vg);
        const uint8_t bg = static_cast<uint8_t>(hg & vg);
        const uint8_t egs = static_cast<uint8_t>(hgs | vgs);
        const uint8_t bgs = static_cast<uint8_t>(hgs & vgs);
        const std::size_t i = static_cast<std::size_t>(k);
        cl_g[i] += static_cast<uint64_t>(cg);
        cl_gs[i] += static_cast<uint64_t>(cgs);
        h_g[i] += hg;
        v_g[i] += vg;
        e_g[i] += eg;
        b_g[i] += bg;
        h_gs[i] += hgs;
        v_gs[i] += vgs;
        e_gs[i] += egs;
        b_gs[i] += bgs;
        cl_g2[i] += static_cast<uint64_t>(cg) * static_cast<uint64_t>(cg);
        cl_gs2[i] += static_cast<uint64_t>(cgs) * static_cast<uint64_t>(cgs);
        cl_xy[i] += static_cast<uint64_t>(cg) * static_cast<uint64_t>(cgs);
        h_xy[i] += static_cast<uint64_t>(hg) * static_cast<uint64_t>(hgs);
        v_xy[i] += static_cast<uint64_t>(vg) * static_cast<uint64_t>(vgs);
        e_xy[i] += static_cast<uint64_t>(eg) * static_cast<uint64_t>(egs);
        b_xy[i] += static_cast<uint64_t>(bg) * static_cast<uint64_t>(bgs);
    }

    void absorb(const Acc& o) {
        replicas += o.replicas;
        for (std::size_t i = 0; i < cl_g.size(); ++i) {
            cl_g[i] += o.cl_g[i];
            cl_gs[i] += o.cl_gs[i];
            h_g[i] += o.h_g[i];
            v_g[i] += o.v_g[i];
            e_g[i] += o.e_g[i];
            b_g[i] += o.b_g[i];
            h_gs[i] += o.h_gs[i];
            v_gs[i] += o.v_gs[i];
            e_gs[i] += o.e_gs[i];
            b_gs[i] += o.b_gs[i];
            cl_g2[i] += o.cl_g2[i];
            cl_gs2[i] += o.cl_gs2[i];
            cl_xy[i] += o.cl_xy[i];
            h_xy[i] += o.h_xy[i];
            v_xy[i] += o.v_xy[i];
            e_xy[i] += o.e_xy[i];
            b_xy[i] += o.b_xy[i];
        }
    }
};

struct Snap {
    std::vector<int> cl_g, cl_gs;
    std::vector<uint8_t> h_g, v_g, h_gs, v_gs;
    void resize(int n) {
        const std::size_t m = static_cast<std::size_t>(n + 1);
        cl_g.assign(m, 0);
        cl_gs.assign(m, 0);
        h_g.assign(m, 0);
        v_g.assign(m, 0);
        h_gs.assign(m, 0);
        v_gs.assign(m, 0);
    }
};

void run_pass(Engine& e, const CSR& csr, const int* sites, int n, Snap& snap, bool gstar) {
    e.begin(&csr);
    if (!gstar) {
        snap.cl_g[0] = 0;
        snap.h_g[0] = 0;
        snap.v_g[0] = 0;
    } else {
        snap.cl_gs[0] = 0;
        snap.h_gs[0] = 0;
        snap.v_gs[0] = 0;
    }
    for (int k = 1; k <= n; ++k) {
        e.add(sites[k - 1]);
        if (!gstar) {
            snap.cl_g[static_cast<std::size_t>(k)] = e.n_clusters;
            snap.h_g[static_cast<std::size_t>(k)] = e.wrap_h;
            snap.v_g[static_cast<std::size_t>(k)] = e.wrap_v;
        } else {
            snap.cl_gs[static_cast<std::size_t>(k)] = e.n_clusters;
            snap.h_gs[static_cast<std::size_t>(k)] = e.wrap_h;
            snap.v_gs[static_cast<std::size_t>(k)] = e.wrap_v;
        }
    }
}

void accumulate_matched(Acc& acc, const Snap& snap, int n) {
    for (int k = 0; k <= n; ++k) {
        const int m = n - k;
        acc.add_k(k, snap.cl_g[static_cast<std::size_t>(k)], snap.cl_gs[static_cast<std::size_t>(m)],
                  snap.h_g[static_cast<std::size_t>(k)], snap.v_g[static_cast<std::size_t>(k)],
                  snap.h_gs[static_cast<std::size_t>(m)], snap.v_gs[static_cast<std::size_t>(m)]);
    }
    acc.replicas += 1;
}

struct Options {
    std::vector<int> sizes;
    std::string mode = "shared";  // shared | independent | both
    int threads = 8;
    uint32_t seed = kMasterSeed;
    int batches = 20;
    int replicas_per_batch = -1;
    int batch_begin = 0;
    int batch_count = -1;
    std::string outdir = "results/issue-9";
    std::string exact_dir = "results/issue-7";
    bool rng_kat = false;
    bool calibrate = false;
    bool campaign = false;
    bool repro = false;
};

void usage(const char* a0) {
    std::cerr
        << "Usage: " << a0
        << " [--L N]... [--mode shared|independent|both] [--threads N] [--seed U32]\n"
        << "       [--batches N] [--replicas-per-batch N] [--batch-begin I] [--batch-count N]\n"
        << "       [--outdir DIR] [--exact-dir DIR] [--rng-kat] [--calibrate] [--campaign] [--repro]\n";
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
        } else if (a == "--mode") {
            const char* v = need("--mode");
            if (!v) {
                return false;
            }
            opt.mode = v;
        } else if (a == "--threads") {
            const char* v = need("--threads");
            if (!v) {
                return false;
            }
            opt.threads = std::stoi(v);
        } else if (a == "--seed") {
            const char* v = need("--seed");
            if (!v) {
                return false;
            }
            opt.seed = static_cast<uint32_t>(std::stoul(v, nullptr, 0));
        } else if (a == "--batches") {
            const char* v = need("--batches");
            if (!v) {
                return false;
            }
            opt.batches = std::stoi(v);
        } else if (a == "--replicas-per-batch") {
            const char* v = need("--replicas-per-batch");
            if (!v) {
                return false;
            }
            opt.replicas_per_batch = std::stoi(v);
        } else if (a == "--batch-begin") {
            const char* v = need("--batch-begin");
            if (!v) {
                return false;
            }
            opt.batch_begin = std::stoi(v);
        } else if (a == "--batch-count") {
            const char* v = need("--batch-count");
            if (!v) {
                return false;
            }
            opt.batch_count = std::stoi(v);
        } else if (a == "--outdir") {
            const char* v = need("--outdir");
            if (!v) {
                return false;
            }
            opt.outdir = v;
        } else if (a == "--exact-dir") {
            const char* v = need("--exact-dir");
            if (!v) {
                return false;
            }
            opt.exact_dir = v;
        } else if (a == "--rng-kat") {
            opt.rng_kat = true;
        } else if (a == "--calibrate") {
            opt.calibrate = true;
        } else if (a == "--campaign") {
            opt.campaign = true;
        } else if (a == "--repro") {
            opt.repro = true;
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return false;
        } else {
            std::cerr << "unknown argument: " << a << "\n";
            usage(argv[0]);
            return false;
        }
    }
    if (opt.threads < 1) {
        opt.threads = 1;
    }
    if (opt.threads > 8) {
        opt.threads = 8;
    }
    return true;
}

int default_rpb(int L) {
    switch (L) {
        case 16:
            return 5000;
        case 24:
            return 3000;
        case 32:
            return 2000;
        case 48:
            return 1000;
        case 64:
            return 750;
        case 96:
            return 300;
        case 128:
            return 200;
        case 192:
            return 80;
        case 256:
            return 50;
        default:
            return 100;
    }
}

std::string pad_L(int L) {
    std::ostringstream oss;
    oss << std::setw(3) << std::setfill('0') << L;
    return oss.str();
}

std::string hex32(uint32_t x) {
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setw(8) << std::setfill('0') << x;
    return oss.str();
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

uint64_t factorial_u64(int n) {
    uint64_t r = 1;
    for (int i = 2; i <= n; ++i) {
        r *= static_cast<uint64_t>(i);
    }
    return r;
}

struct ExactMC {
    int L = 0;
    int n = 0;
    std::vector<uint64_t> count, cl_g, cl_gs, h_g, v_g, e_g, b_g, h_gs, v_gs, e_gs, b_gs;
    bool ok = false;
};

ExactMC load_exact(const std::string& path, int L) {
    ExactMC e;
    e.L = L;
    e.n = L * L;
    const std::size_t m = static_cast<std::size_t>(e.n + 1);
    e.count.assign(m, 0);
    e.cl_g.assign(m, 0);
    e.cl_gs.assign(m, 0);
    e.h_g.assign(m, 0);
    e.v_g.assign(m, 0);
    e.e_g.assign(m, 0);
    e.b_g.assign(m, 0);
    e.h_gs.assign(m, 0);
    e.v_gs.assign(m, 0);
    e.e_gs.assign(m, 0);
    e.b_gs.assign(m, 0);
    std::ifstream in(path);
    if (!in) {
        return e;
    }
    std::string header;
    std::getline(in, header);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream iss(line);
        int k = 0;
        iss >> k;
        if (k < 0 || k > e.n) {
            continue;
        }
        const std::size_t i = static_cast<std::size_t>(k);
        iss >> e.count[i] >> e.cl_g[i] >> e.cl_gs[i] >> e.h_g[i] >> e.v_g[i] >> e.e_g[i] >> e.b_g[i] >>
            e.h_gs[i] >> e.v_gs[i] >> e.e_gs[i] >> e.b_gs[i];
    }
    e.ok = true;
    return e;
}

void write_microcanonical_csv(const std::string& path, const Acc& a) {
    std::ofstream out(path);
    out << "k,replicas,sum_clusters_G,sum_clusters_Gstar,"
           "sum_wrap_H_G,sum_wrap_V_G,sum_wrap_E_G,sum_wrap_B_G,"
           "sum_wrap_H_Gstar,sum_wrap_V_Gstar,sum_wrap_E_Gstar,sum_wrap_B_Gstar,"
           "sum_clusters_G_sq,sum_clusters_Gstar_sq,sum_clusters_G_Gstar,"
           "sum_H_G_H_Gstar,sum_V_G_V_Gstar,sum_E_G_E_Gstar,sum_B_G_B_Gstar\n";
    for (int k = 0; k <= a.n; ++k) {
        const std::size_t i = static_cast<std::size_t>(k);
        out << k << ',' << a.replicas << ',' << a.cl_g[i] << ',' << a.cl_gs[i] << ',' << a.h_g[i] << ','
            << a.v_g[i] << ',' << a.e_g[i] << ',' << a.b_g[i] << ',' << a.h_gs[i] << ',' << a.v_gs[i]
            << ',' << a.e_gs[i] << ',' << a.b_gs[i] << ',' << a.cl_g2[i] << ',' << a.cl_gs2[i] << ','
            << a.cl_xy[i] << ',' << a.h_xy[i] << ',' << a.v_xy[i] << ',' << a.e_xy[i] << ',' << a.b_xy[i]
            << '\n';
    }
}

void write_batch_bin(const std::string& path, int L, const std::vector<Acc>& batches) {
    std::ofstream out(path, std::ios::binary);
    const char magic[4] = {'N', 'Z', 'B', '1'};
    out.write(magic, 4);
    const uint32_t uL = static_cast<uint32_t>(L);
    const uint32_t nb = static_cast<uint32_t>(batches.size());
    const uint32_t nk = static_cast<uint32_t>(L * L + 1);
    out.write(reinterpret_cast<const char*>(&uL), 4);
    out.write(reinterpret_cast<const char*>(&nb), 4);
    out.write(reinterpret_cast<const char*>(&nk), 4);
    const uint64_t rpb = batches.empty() ? 0 : batches[0].replicas;
    out.write(reinterpret_cast<const char*>(&rpb), 8);
    for (const Acc& a : batches) {
        for (uint32_t k = 0; k < nk; ++k) {
            const std::size_t i = static_cast<std::size_t>(k);
            const uint64_t row[10] = {a.cl_g[i], a.cl_gs[i], a.h_g[i],  a.v_g[i],  a.e_g[i],
                                      a.b_g[i],  a.h_gs[i],  a.v_gs[i], a.e_gs[i], a.b_gs[i]};
            out.write(reinterpret_cast<const char*>(row), sizeof(row));
        }
    }
}

void append_performance(const std::string& path, int L, const std::string& mode, uint64_t replicas,
                        double wall, double cpu, long rss, int threads) {
    const bool exists = static_cast<bool>(std::ifstream(path));
    std::ofstream out(path, std::ios::app);
    if (!exists) {
        out << "L,mode,replicas,threads,wall_seconds,cpu_seconds,peak_rss_kb,site_updates,"
               "site_updates_per_sec,replicas_per_sec\n";
    }
    const uint64_t n = static_cast<uint64_t>(L) * static_cast<uint64_t>(L);
    const uint64_t updates = replicas * 2ull * n;
    const double ups = wall > 0.0 ? static_cast<double>(updates) / wall : 0.0;
    const double rps = wall > 0.0 ? static_cast<double>(replicas) / wall : 0.0;
    out << L << ',' << mode << ',' << replicas << ',' << threads << ',' << std::setprecision(12) << wall
        << ',' << cpu << ',' << rss << ',' << updates << ',' << ups << ',' << rps << '\n';
}

bool run_rng_kat(const std::string& outdir) {
    const auto kats = philox4x32_10_official_kats();
    bool all = true;
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"algorithm\": \"Philox4x32-10\",\n";
    oss << "  \"source\": \"Random123 tests/kat_vectors (DEShawResearch/random123); Salmon et al. SC11\",\n";
    oss << "  \"counter_layout\": {\n";
    oss << "    \"key0\": \"global_seed\",\n";
    oss << "    \"key1\": \"batch_id\",\n";
    oss << "    \"ctr0\": \"draw_counter\",\n";
    oss << "    \"ctr1\": \"replica_id\",\n";
    oss << "    \"ctr2\": \"stream_id (0=G/shared, 1=independent G*)\",\n";
    oss << "    \"ctr3\": \"0\"\n";
    oss << "  },\n";
    oss << "  \"kats\": [\n";
    for (std::size_t i = 0; i < kats.size(); ++i) {
        const auto got = philox4x32_10(kats[i].ctr, kats[i].key);
        bool pass = true;
        for (int w = 0; w < 4; ++w) {
            if (got[static_cast<std::size_t>(w)] != kats[i].expected[static_cast<std::size_t>(w)]) {
                pass = false;
            }
        }
        all = all && pass;
        oss << "    {\n";
        oss << "      \"ctr\": [";
        for (int w = 0; w < 4; ++w) {
            if (w) {
                oss << ", ";
            }
            oss << '"' << hex32(kats[i].ctr[static_cast<std::size_t>(w)]) << '"';
        }
        oss << "],\n      \"key\": [\"" << hex32(kats[i].key[0]) << "\", \"" << hex32(kats[i].key[1])
            << "\"],\n";
        oss << "      \"expected\": [";
        for (int w = 0; w < 4; ++w) {
            if (w) {
                oss << ", ";
            }
            oss << '"' << hex32(kats[i].expected[static_cast<std::size_t>(w)]) << '"';
        }
        oss << "],\n      \"got\": [";
        for (int w = 0; w < 4; ++w) {
            if (w) {
                oss << ", ";
            }
            oss << '"' << hex32(got[static_cast<std::size_t>(w)]) << '"';
        }
        oss << "],\n      \"pass\": " << (pass ? "true" : "false") << "\n    }";
        if (i + 1 < kats.size()) {
            oss << ",";
        }
        oss << "\n";
    }
    oss << "  ],\n";

    PhiloxStream a, b;
    a.reset(0x12345678u, 3u, 9u, 0u);
    b.reset(0x12345678u, 3u, 9u, 0u);
    bool stream_ok = true;
    for (int i = 0; i < 64; ++i) {
        if (a.next_u32() != b.next_u32()) {
            stream_ok = false;
        }
    }
    PhiloxStream c;
    c.reset(0x12345678u, 3u, 9u, 1u);
    bool distinct_stream = false;
    a.reset(0x12345678u, 3u, 9u, 0u);
    for (int i = 0; i < 16; ++i) {
        if (a.next_u32() != c.next_u32()) {
            distinct_stream = true;
        }
    }
    all = all && stream_ok && distinct_stream;
    oss << "  \"stream_determinism_same_tuple\": " << (stream_ok ? "true" : "false") << ",\n";
    oss << "  \"independent_stream_id_distinct\": " << (distinct_stream ? "true" : "false") << ",\n";
    oss << "  \"all_pass\": " << (all ? "true" : "false") << "\n";
    oss << "}\n";
    std::filesystem::create_directories(outdir);
    std::ofstream out(outdir + "/rng_validation.json");
    out << oss.str();
    std::cerr << "RNG KAT " << (all ? "PASS" : "FAIL") << " wrote " << outdir << "/rng_validation.json\n";
    return all;
}

uint64_t mask_from_sites(const int* sites, int k) {
    uint64_t m = 0;
    for (int i = 0; i < k; ++i) {
        m |= (1ull << sites[i]);
    }
    return m;
}

bool snap_matches_exact(int L, const Graph& g, const Graph& gs, DisplacementDSU& dsu, const Snap& snap,
                        const int* perm, int n, std::string& err) {
    for (int k = 0; k <= n; ++k) {
        const uint64_t mask = mask_from_sites(perm, k);
        const Observables o = analyze_config(mask, g, gs, dsu);
        const int m = n - k;
        if (o.k != k || o.clusters_G != snap.cl_g[static_cast<std::size_t>(k)] ||
            o.H_G != snap.h_g[static_cast<std::size_t>(k)] || o.V_G != snap.v_g[static_cast<std::size_t>(k)] ||
            o.clusters_Gstar != snap.cl_gs[static_cast<std::size_t>(m)] ||
            o.H_Gstar != snap.h_gs[static_cast<std::size_t>(m)] ||
            o.V_Gstar != snap.v_gs[static_cast<std::size_t>(m)]) {
            std::ostringstream oss;
            oss << "L=" << L << " k=" << k << " mask=" << mask << " exact(clG,HG,VG,clGs,HGs,VGs)=("
                << o.clusters_G << "," << int(o.H_G) << "," << int(o.V_G) << "," << o.clusters_Gstar << ","
                << int(o.H_Gstar) << "," << int(o.V_Gstar) << ") nz=("
                << snap.cl_g[static_cast<std::size_t>(k)] << "," << int(snap.h_g[static_cast<std::size_t>(k)])
                << "," << int(snap.v_g[static_cast<std::size_t>(k)]) << ","
                << snap.cl_gs[static_cast<std::size_t>(m)] << "," << int(snap.h_gs[static_cast<std::size_t>(m)])
                << "," << int(snap.v_gs[static_cast<std::size_t>(m)]) << ")";
            err = oss.str();
            return false;
        }
    }
    return true;
}

void fill_hand_perm(std::vector<int>& perm, int n, int L, int which) {
    perm.resize(static_cast<std::size_t>(n));
    if (which == 0) {
        std::iota(perm.begin(), perm.end(), 0);
    } else if (which == 1) {
        for (int i = 0; i < n; ++i) {
            perm[static_cast<std::size_t>(i)] = n - 1 - i;
        }
    } else if (which == 2) {
        int t = 0;
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                perm[static_cast<std::size_t>(t++)] = y * L + x;
            }
        }
    } else {
        int t = 0;
        for (int s = 0; s < 2 * L; ++s) {
            for (int y = 0; y < L; ++y) {
                const int x = s - y;
                if (x >= 0 && x < L) {
                    perm[static_cast<std::size_t>(t++)] = y * L + x;
                }
            }
        }
    }
}

bool run_calibrate(const Options& opt) {
    std::filesystem::create_directories(opt.outdir);
    bool all = true;
    std::ostringstream json;
    json << "{\n";
    json << "  \"exact_source\": \"" << opt.exact_dir << "\",\n";
    json << "  \"note\": \"Exact microcanonical expectations are Issue #7 totals / C(N,k). "
            "NZ shared mode uses one Fisher-Yates permutation for G (forward) and G* (reverse / "
            "complement). Reverse permutation is matching complement, not antithetic.\",\n";

    json << "  \"subset_enumeration\": {\n";
    bool first_item = true;
    for (int L : {2, 3, 4}) {
        const Graph g = make_graph(L, false);
        const Graph gs = make_graph(L, true);
        const CSR csr_g = make_undirected_csr(g);
        const CSR csr_gs = make_undirected_csr(gs);
        const int n = L * L;
        Engine eng;
        eng.setup(n);
        DisplacementDSU dsu(n);
        uint64_t mismatches = 0;
        uint64_t first_bad = ~0ull;
        const uint64_t total = 1ull << n;
        for (uint64_t mask = 0; mask < total; ++mask) {
            eng.begin(&csr_g);
            for (int i = 0; i < n; ++i) {
                if ((mask >> i) & 1ull) {
                    eng.add(i);
                }
            }
            Engine engs;
            engs.setup(n);
            engs.begin(&csr_gs);
            for (int i = 0; i < n; ++i) {
                if (((mask >> i) & 1ull) == 0ull) {
                    engs.add(i);
                }
            }
            const Observables o = analyze_config(mask, g, gs, dsu);
            if (eng.n_clusters != o.clusters_G || eng.wrap_h != o.H_G || eng.wrap_v != o.V_G ||
                engs.n_clusters != o.clusters_Gstar || engs.wrap_h != o.H_Gstar || engs.wrap_v != o.V_Gstar) {
                mismatches += 1;
                if (mask < first_bad) {
                    first_bad = mask;
                }
            }
        }
        const char* st = mismatches == 0 ? "PASS" : "FAIL";
        if (mismatches) {
            all = false;
        }
        if (!first_item) {
            json << ",\n";
        }
        first_item = false;
        json << "    \"L" << L << "\": {\"status\": \"" << st << "\", \"configs\": " << total
             << ", \"mismatches\": " << mismatches;
        if (mismatches) {
            json << ", \"first_bad_mask\": " << first_bad;
        }
        json << "}";
        std::cerr << "subset L=" << L << " " << st << " mismatches=" << mismatches << "\n";
    }
    json << "\n  },\n";

    json << "  \"hand_crafted_and_prefix\": {\n";
    first_item = true;
    for (int L : {2, 3, 4}) {
        const Graph g = make_graph(L, false);
        const Graph gs = make_graph(L, true);
        const CSR csr_g = make_undirected_csr(g);
        const CSR csr_gs = make_undirected_csr(gs);
        const int n = L * L;
        Engine eng;
        eng.setup(n);
        Snap snap;
        snap.resize(n);
        DisplacementDSU dsu(n);
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::vector<int> rev(static_cast<std::size_t>(n));
        int fails = 0;
        std::string first_err;
        auto check_perm = [&](const std::vector<int>& p) {
            run_pass(eng, csr_g, p.data(), n, snap, false);
            for (int i = 0; i < n; ++i) {
                rev[static_cast<std::size_t>(i)] = p[static_cast<std::size_t>(n - 1 - i)];
            }
            run_pass(eng, csr_gs, rev.data(), n, snap, true);
            std::string err;
            if (!snap_matches_exact(L, g, gs, dsu, snap, p.data(), n, err)) {
                fails += 1;
                if (first_err.empty()) {
                    first_err = err;
                }
            }
        };
        for (int w = 0; w < 4; ++w) {
            fill_hand_perm(perm, n, L, w);
            check_perm(perm);
        }
        const int nrand = (L == 2) ? 0 : (L == 3 ? 200 : 80);
        PhiloxStream rng;
        for (int r = 0; r < nrand; ++r) {
            rng.reset(opt.seed, 0u, static_cast<uint32_t>(r), 0u);
            fisher_yates(perm, rng);
            check_perm(perm);
        }
        if (L == 2) {
            std::iota(perm.begin(), perm.end(), 0);
            do {
                check_perm(perm);
            } while (std::next_permutation(perm.begin(), perm.end()));
        }
        const char* st = fails == 0 ? "PASS" : "FAIL";
        if (fails) {
            all = false;
        }
        if (!first_item) {
            json << ",\n";
        }
        first_item = false;
        json << "    \"L" << L << "\": {\"status\": \"" << st << "\", \"failures\": " << fails;
        if (fails) {
            json << ", \"first_error\": \"" << first_err << "\"";
        }
        json << "}";
        std::cerr << "prefix L=" << L << " " << st << " failures=" << fails << "\n";
    }
    json << "\n  },\n";

    json << "  \"exhaustive_permutations\": {\n";
    first_item = true;
    for (int L : {2, 3}) {
        const Graph g = make_graph(L, false);
        const Graph gs = make_graph(L, true);
        const CSR csr_g = make_undirected_csr(g);
        const CSR csr_gs = make_undirected_csr(gs);
        const int n = L * L;
        std::ostringstream ep;
        ep << opt.exact_dir << "/L" << std::setw(2) << std::setfill('0') << L << "_microcanonical.csv";
        const ExactMC exact = load_exact(ep.str(), L);
        Engine eng;
        eng.setup(n);
        Snap snap;
        snap.resize(n);
        Acc acc(n);
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::vector<int> rev(static_cast<std::size_t>(n));
        std::iota(perm.begin(), perm.end(), 0);
        do {
            run_pass(eng, csr_g, perm.data(), n, snap, false);
            for (int i = 0; i < n; ++i) {
                rev[static_cast<std::size_t>(i)] = perm[static_cast<std::size_t>(n - 1 - i)];
            }
            run_pass(eng, csr_gs, rev.data(), n, snap, true);
            accumulate_matched(acc, snap, n);
        } while (std::next_permutation(perm.begin(), perm.end()));
        bool pass = exact.ok;
        int first_bad_k = -1;
        std::string bad_field;
        for (int k = 0; k <= n && pass; ++k) {
            const uint64_t fac = factorial_u64(k) * factorial_u64(n - k);
            const std::size_t i = static_cast<std::size_t>(k);
            auto chk = [&](uint64_t got, uint64_t tot, const char* name) {
                if (got != tot * fac) {
                    pass = false;
                    first_bad_k = k;
                    bad_field = name;
                }
            };
            chk(acc.cl_g[i], exact.cl_g[i], "clusters_G");
            chk(acc.cl_gs[i], exact.cl_gs[i], "clusters_Gstar");
            chk(acc.h_g[i], exact.h_g[i], "H_G");
            chk(acc.v_g[i], exact.v_g[i], "V_G");
            chk(acc.e_g[i], exact.e_g[i], "E_G");
            chk(acc.b_g[i], exact.b_g[i], "B_G");
            chk(acc.h_gs[i], exact.h_gs[i], "H_Gstar");
            chk(acc.v_gs[i], exact.v_gs[i], "V_Gstar");
            chk(acc.e_gs[i], exact.e_gs[i], "E_Gstar");
            chk(acc.b_gs[i], exact.b_gs[i], "B_Gstar");
        }
        if (!pass) {
            all = false;
        }
        if (!first_item) {
            json << ",\n";
        }
        first_item = false;
        json << "    \"L" << L << "\": {\"status\": \"" << (pass ? "PASS" : "FAIL")
             << "\", \"permutations\": " << acc.replicas << ", \"loaded_exact\": "
             << (exact.ok ? "true" : "false");
        if (!pass) {
            json << ", \"first_bad_k\": " << first_bad_k << ", \"field\": \"" << bad_field << "\"";
        }
        json << "}";
        std::cerr << "exhaustive perms L=" << L << " " << (pass ? "PASS" : "FAIL")
                  << " nperms=" << acc.replicas << "\n";
    }
    json << "\n  },\n";

    json << "  \"mc_convergence\": {\n";
    {
        const int L = 4;
        const Graph g = make_graph(L, false);
        const Graph gs = make_graph(L, true);
        const CSR csr_g = make_undirected_csr(g);
        const CSR csr_gs = make_undirected_csr(gs);
        const int n = L * L;
        std::ostringstream ep;
        ep << opt.exact_dir << "/L" << std::setw(2) << std::setfill('0') << L << "_microcanonical.csv";
        const ExactMC exact = load_exact(ep.str(), L);
        const int nrep = 40000;
        Acc acc(n);
        Engine eng;
        eng.setup(n);
        Snap snap;
        snap.resize(n);
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::vector<int> rev(static_cast<std::size_t>(n));
        PhiloxStream rng;
        for (int r = 0; r < nrep; ++r) {
            rng.reset(opt.seed, 99u, static_cast<uint32_t>(r), 0u);
            fisher_yates(perm, rng);
            run_pass(eng, csr_g, perm.data(), n, snap, false);
            for (int i = 0; i < n; ++i) {
                rev[static_cast<std::size_t>(i)] = perm[static_cast<std::size_t>(n - 1 - i)];
            }
            run_pass(eng, csr_gs, rev.data(), n, snap, true);
            accumulate_matched(acc, snap, n);
        }
        double max_wrap_abs = 0.0;
        double max_cl_abs = 0.0;
        int max_wrap_k = 0;
        int max_cl_k = 0;
        std::string max_wrap_field;
        for (int k = 0; k <= n; ++k) {
            const std::size_t i = static_cast<std::size_t>(k);
            const double den = static_cast<double>(exact.count[i]);
            const double rn = static_cast<double>(acc.replicas);
            auto wrap_err = [&](uint64_t got, uint64_t tot, const char* name) {
                const double e = got / rn - tot / den;
                const double a = std::fabs(e);
                if (a > max_wrap_abs) {
                    max_wrap_abs = a;
                    max_wrap_k = k;
                    max_wrap_field = name;
                }
            };
            wrap_err(acc.h_g[i], exact.h_g[i], "H_G");
            wrap_err(acc.v_g[i], exact.v_g[i], "V_G");
            wrap_err(acc.e_g[i], exact.e_g[i], "E_G");
            wrap_err(acc.b_g[i], exact.b_g[i], "B_G");
            wrap_err(acc.h_gs[i], exact.h_gs[i], "H_Gstar");
            wrap_err(acc.v_gs[i], exact.v_gs[i], "V_Gstar");
            wrap_err(acc.e_gs[i], exact.e_gs[i], "E_Gstar");
            wrap_err(acc.b_gs[i], exact.b_gs[i], "B_Gstar");
            const double clg = acc.cl_g[i] / rn - exact.cl_g[i] / den;
            const double clgs = acc.cl_gs[i] / rn - exact.cl_gs[i] / den;
            if (std::fabs(clg) > max_cl_abs) {
                max_cl_abs = std::fabs(clg);
                max_cl_k = k;
            }
            if (std::fabs(clgs) > max_cl_abs) {
                max_cl_abs = std::fabs(clgs);
                max_cl_k = k;
            }
        }
        const double se_wrap = 1.0 / std::sqrt(static_cast<double>(nrep));
        const bool pass = exact.ok && max_wrap_abs < 5.0 * se_wrap && max_cl_abs < 0.05;
        if (!pass) {
            all = false;
        }
        json << "    \"L4\": {\"status\": \"" << (pass ? "PASS" : "FAIL") << "\", \"replicas\": " << nrep
             << ", \"max_wrap_abs_err\": " << std::setprecision(12) << max_wrap_abs
             << ", \"max_wrap_k\": " << max_wrap_k << ", \"max_wrap_field\": \"" << max_wrap_field
             << "\", \"max_cluster_abs_err\": " << max_cl_abs << ", \"max_cluster_k\": " << max_cl_k
             << ", \"bound_5_se\": " << (5.0 * se_wrap) << "}";
        std::cerr << "MC L=4 " << (pass ? "PASS" : "FAIL") << " nrep=" << nrep
                  << " max_wrap_abs=" << max_wrap_abs << " max_cl_abs=" << max_cl_abs << "\n";
    }
    json << "\n  },\n";
    json << "  \"overall\": \"" << (all ? "PASS" : "FAIL") << "\"\n";
    json << "}\n";
    std::ofstream out(opt.outdir + "/exact_calibration.json");
    out << json.str();
    std::cerr << "calibration " << (all ? "PASS" : "FAIL") << " wrote " << opt.outdir
              << "/exact_calibration.json\n";
    return all;
}

void run_one_replica(Engine& eng, Snap& snap, std::vector<int>& perm, std::vector<int>& perm2,
                     std::vector<int>& rev, const CSR& csr_g, const CSR& csr_gs, uint32_t seed,
                     uint32_t batch, uint32_t replica, bool shared, Acc& acc) {
    const int n = csr_g.n;
    PhiloxStream rng;
    rng.reset(seed, batch, replica, 0u);
    fisher_yates(perm, rng);
    run_pass(eng, csr_g, perm.data(), n, snap, false);
    if (shared) {
        for (int i = 0; i < n; ++i) {
            rev[static_cast<std::size_t>(i)] = perm[static_cast<std::size_t>(n - 1 - i)];
        }
        run_pass(eng, csr_gs, rev.data(), n, snap, true);
    } else {
        PhiloxStream rng2;
        rng2.reset(seed, batch, replica, 1u);
        fisher_yates(perm2, rng2);
        run_pass(eng, csr_gs, perm2.data(), n, snap, true);
    }
    accumulate_matched(acc, snap, n);
}

Acc run_size(int L, const std::string& mode, const Options& opt, std::vector<Acc>& batches_out,
             double& wall, double& cpu, long& rss) {
    const bool shared = (mode == "shared");
    const Graph g = make_graph(L, false);
    const Graph gs = make_graph(L, true);
    const CSR csr_g = make_undirected_csr(g);
    const CSR csr_gs = make_undirected_csr(gs);
    const int n = L * L;
    const int rpb =
        opt.replicas_per_batch > 0 ? opt.replicas_per_batch : default_rpb(L);
    int nbatches = opt.batch_count > 0 ? opt.batch_count : opt.batches;
    const int b0 = opt.batch_begin;
    if (nbatches < 1) {
        nbatches = 1;
    }
    const int nthreads = opt.threads;

#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    batches_out.assign(static_cast<std::size_t>(nbatches), Acc(n));
    Acc pooled(n);

    const double cpu0 = cpu_seconds();
    const auto t0 = std::chrono::steady_clock::now();

    for (int bi = 0; bi < nbatches; ++bi) {
        const uint32_t batch_id = static_cast<uint32_t>(b0 + bi);
        std::vector<Acc> local(static_cast<std::size_t>(nthreads), Acc(n));

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
        {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            Engine eng;
            eng.setup(n);
            Snap snap;
            snap.resize(n);
            std::vector<int> perm(static_cast<std::size_t>(n));
            std::vector<int> perm2(static_cast<std::size_t>(n));
            std::vector<int> rev(static_cast<std::size_t>(n));
            Acc& acc = local[static_cast<std::size_t>(tid)];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int r = 0; r < rpb; ++r) {
                run_one_replica(eng, snap, perm, perm2, rev, csr_g, csr_gs, opt.seed, batch_id,
                                static_cast<uint32_t>(r), shared, acc);
            }
        }

        Acc& bout = batches_out[static_cast<std::size_t>(bi)];
        for (const Acc& a : local) {
            bout.absorb(a);
        }
        pooled.absorb(bout);
        std::cerr << "  L=" << L << " mode=" << mode << " batch=" << batch_id << " replicas=" << bout.replicas
                  << "\n";
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double cpu1 = cpu_seconds();
    wall = std::chrono::duration<double>(t1 - t0).count();
    cpu = cpu1 - cpu0;
    rss = peak_rss_kb();
    return pooled;
}

int run_production(const Options& opt) {
    std::filesystem::create_directories(opt.outdir);
    std::filesystem::create_directories(opt.outdir + "/batch_bin");
    std::vector<int> sizes = opt.sizes;
    if (sizes.empty() && opt.campaign) {
        sizes = {16, 24, 32, 48, 64, 96, 128, 192, 256};
    }
    if (sizes.empty()) {
        std::cerr << "no sizes; pass --L or --campaign\n";
        return 2;
    }
    std::vector<std::string> modes;
    if (opt.mode == "both") {
        modes = {"shared", "independent"};
    } else {
        modes = {opt.mode};
    }
    int rc = 0;
    for (int L : sizes) {
        if (L < 2) {
            std::cerr << "skip L=" << L << "\n";
            continue;
        }
        for (const std::string& mode : modes) {
            if (mode != "shared" && mode != "independent") {
                std::cerr << "unknown mode " << mode << "\n";
                rc = 1;
                continue;
            }
            std::cerr << "running L=" << L << " mode=" << mode << " threads=" << opt.threads << "\n";
            std::vector<Acc> batches;
            double wall = 0, cpu = 0;
            long rss = 0;
            Acc pooled = run_size(L, mode, opt, batches, wall, cpu, rss);
            const std::string tag = pad_L(L);
            std::string csv = opt.outdir + "/microcanonical_L" + tag + ".csv";
            if (mode == "independent") {
                csv = opt.outdir + "/microcanonical_L" + tag + "_independent.csv";
            }
            if (opt.repro) {
                csv = opt.outdir + "/repro_L" + tag + "_b" + std::to_string(opt.batch_begin) + "_t" +
                      std::to_string(opt.threads) + ".csv";
            }
            write_microcanonical_csv(csv, pooled);
            const std::string bin =
                opt.outdir + "/batch_bin/L" + tag + "_" + mode + ".bin";
            if (!opt.repro) {
                write_batch_bin(bin, L, batches);
            }
            append_performance(opt.outdir + "/performance.csv", L, mode, pooled.replicas, wall, cpu, rss,
                               opt.threads);
            std::cerr << "done L=" << L << " mode=" << mode << " replicas=" << pooled.replicas
                      << " wall=" << wall << "s site_updates/s="
                      << (wall > 0 ? (static_cast<double>(pooled.replicas) * 2.0 * L * L / wall) : 0.0)
                      << " file=" << csv << "\n";
        }
    }
    return rc;
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
    if (opt.rng_kat || opt.campaign) {
        if (!run_rng_kat(opt.outdir)) {
            rc = 1;
        }
    }
    if (opt.calibrate || opt.campaign) {
        if (!run_calibrate(opt)) {
            rc = 1;
        }
    }
    if (opt.campaign || opt.repro || !opt.sizes.empty()) {
        const int prc = run_production(opt);
        if (prc) {
            rc = prc;
        }
    }
    if (!opt.rng_kat && !opt.calibrate && !opt.campaign && !opt.repro && opt.sizes.empty()) {
        usage(argv[0]);
        return 2;
    }
    return rc;
}
