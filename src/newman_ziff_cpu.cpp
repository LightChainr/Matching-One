// C05 / P33: threshold-rank Newman-Ziff on the shared C00 homology DSU.
//
// Reuses issue-9 Philox Fisher-Yates (counter_rng.hpp) and the C00/PR21
// HomologyUnionFind (exact adj(P)/det(P) windings).  One DSU, not a fork.
//
// Off-by-one convention (frozen by exact tiny tests):
//   Add sites in increasing U order (a uniform permutation).
//   K_plus  = min k in {0,...,N+1} such that the primal CROSS-wraps after k
//             occupied black sites. Empty never cross-wraps (K_plus>=1);
//             never wraps => K_plus=N+1.
//   K_minus = min k such that the white matching complement does NOT
//             CROSS-wrap. If m* is the first reverse-permutation occupation
//             at which matching CROSS-wraps, K_minus = N - m* + 1.
//             White never wraps => K_minus=0; always wraps => K_minus=N+1.
//   For 1<=K<=N the transition occurs when the K-th site is occupied:
//             T | K=k  ~ Beta(k, N+1-k)  (k-th uniform order statistic).
//   Cross channel requires K_minus <= K_plus configuration-by-configuration.
//
// Matching function from the rank pair:
//   M(p) = P(K_plus <= m) - P(K_minus > m),  m ~ Binomial(N,p).
// Production records CROSS ranks only. Either is an exact-test diagnostic,
// not a second replication of M.
//
// Build:
//   g++ -O3 -std=c++17 -fopenmp src/newman_ziff_cpu.cpp -o build/newman_ziff_cpu

#include "counter_rng.hpp"
#include "homology_union_find.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using matching::HomologyUnionFind;
using matching::PeriodMatrix;
using matching::PhiloxStream;
using matching::fisher_yates;
using matching::philox4x32_10;
using matching::philox4x32_10_official_kats;

constexpr uint32_t kMasterSeed = 0xC0504E5Au;  // C05 NZ

struct Edge {
    int i = 0;
    int j = 0;
    int dx = 0;
    int dy = 0;
};

struct CSR {
    int n = 0;
    std::vector<int> off;
    std::vector<int> to;
    std::vector<int> dx;
    std::vector<int> dy;
};

CSR make_csr(int n, const std::vector<Edge>& directed) {
    std::vector<std::vector<std::tuple<int, int, int>>> adj(static_cast<std::size_t>(n));
    for (const Edge& e : directed) {
        adj[static_cast<std::size_t>(e.i)].push_back({e.j, e.dx, e.dy});
        adj[static_cast<std::size_t>(e.j)].push_back({e.i, -e.dx, -e.dy});
    }
    CSR c;
    c.n = n;
    c.off.assign(static_cast<std::size_t>(n + 1), 0);
    int m = 0;
    for (int i = 0; i < n; ++i) {
        c.off[static_cast<std::size_t>(i)] = m;
        m += static_cast<int>(adj[static_cast<std::size_t>(i)].size());
    }
    c.off[static_cast<std::size_t>(n)] = m;
    c.to.resize(static_cast<std::size_t>(m));
    c.dx.resize(static_cast<std::size_t>(m));
    c.dy.resize(static_cast<std::size_t>(m));
    for (int i = 0; i < n; ++i) {
        int z = c.off[static_cast<std::size_t>(i)];
        for (const auto [j, dx, dy] : adj[static_cast<std::size_t>(i)]) {
            c.to[static_cast<std::size_t>(z)] = j;
            c.dx[static_cast<std::size_t>(z)] = dx;
            c.dy[static_cast<std::size_t>(z)] = dy;
            ++z;
        }
    }
    return c;
}

struct Geometry {
    std::string name;
    int n = 0;
    int a = 0;
    int b = 0;
    int L = 0;
    PeriodMatrix period;
    CSR primal;
    CSR matching;
    double theta = 0.0;
    double cos4 = 0.0;
};

int positive_mod(int value, int modulus) {
    value %= modulus;
    return value < 0 ? value + modulus : value;
}

Geometry axis_geometry(int L) {
    if (L < 2) throw std::invalid_argument("axis L must be >= 2");
    Geometry g;
    g.name = "axis_L" + std::to_string(L);
    g.L = L;
    g.a = L;
    g.b = 0;
    g.n = L * L;
    g.period = PeriodMatrix::diagonal(L, L);
    g.theta = 0.0;
    g.cos4 = 1.0;
    std::vector<Edge> primal;
    std::vector<Edge> matching;
    primal.reserve(static_cast<std::size_t>(2 * g.n));
    matching.reserve(static_cast<std::size_t>(4 * g.n));
    const int offs[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const int i = y * L + x;
            for (int k = 0; k < 4; ++k) {
                const int dx = offs[k][0];
                const int dy = offs[k][1];
                const int nx = positive_mod(x + dx, L);
                const int ny = positive_mod(y + dy, L);
                const int j = ny * L + nx;
                const Edge e{i, j, dx, dy};
                if (k < 2) primal.push_back(e);
                matching.push_back(e);
            }
        }
    }
    g.primal = make_csr(g.n, primal);
    g.matching = make_csr(g.n, matching);
    return g;
}

Geometry gaussian_geometry(int a, int b) {
    if (a <= 0 || b < 0 || std::gcd(a, b) != 1) {
        throw std::invalid_argument("Gaussian representation requires a>0, b>=0, gcd=1");
    }
    Geometry g;
    g.a = a;
    g.b = b;
    g.n = a * a + b * b;
    g.L = 0;
    g.period = PeriodMatrix::gaussian(a, b);
    g.theta = std::atan2(static_cast<double>(b), static_cast<double>(a));
    g.cos4 = std::cos(4.0 * g.theta);
    std::ostringstream nm;
    nm << "g_" << a << "_" << b;
    g.name = nm.str();
    std::vector<Edge> primal;
    std::vector<Edge> matching;
    primal.reserve(static_cast<std::size_t>(2 * g.n));
    matching.reserve(static_cast<std::size_t>(4 * g.n));
    const int ra = positive_mod(a, g.n);
    const int rb = positive_mod(b, g.n);
    const int rapb = positive_mod(a + b, g.n);
    const int ramb = positive_mod(a - b, g.n);
    for (int j = 0; j < g.n; ++j) {
        const Edge east{j, (j + ra) % g.n, 1, 0};
        const Edge north{j, (j + rb) % g.n, 0, 1};
        const Edge diag{j, (j + rapb) % g.n, 1, 1};
        const Edge anti{j, (j + ramb) % g.n, 1, -1};
        primal.push_back(east);
        primal.push_back(north);
        matching.push_back(east);
        matching.push_back(north);
        matching.push_back(diag);
        matching.push_back(anti);
    }
    g.primal = make_csr(g.n, primal);
    g.matching = make_csr(g.n, matching);
    return g;
}

struct Engine {
    HomologyUnionFind uf;
    std::vector<uint8_t> occ;
    const CSR* csr = nullptr;

    explicit Engine(int n, PeriodMatrix period) : uf(n, period), occ(static_cast<std::size_t>(n), 0) {}

    void begin(const CSR* c) {
        csr = c;
        uf.reset();
        std::memset(occ.data(), 0, occ.size());
    }

    void add(int i) {
        occ[static_cast<std::size_t>(i)] = 1;
        const int z0 = csr->off[static_cast<std::size_t>(i)];
        const int z1 = csr->off[static_cast<std::size_t>(i + 1)];
        for (int z = z0; z < z1; ++z) {
            const int j = csr->to[static_cast<std::size_t>(z)];
            if (!occ[static_cast<std::size_t>(j)]) continue;
            uf.add_edge(i, j, csr->dx[static_cast<std::size_t>(z)], csr->dy[static_cast<std::size_t>(z)]);
        }
    }

    bool cross() const { return uf.channels().cross; }
};

// Returns (K_minus, K_plus) for the CROSS channel under the frozen convention.
std::pair<int, int> threshold_ranks(Engine& eng, const Geometry& g, const int* perm) {
    const int n = g.n;
    eng.begin(&g.primal);
    int k_plus = n + 1;
    for (int k = 1; k <= n; ++k) {
        eng.add(perm[k - 1]);
        if (eng.cross()) {
            k_plus = k;
            break;
        }
    }
    eng.begin(&g.matching);
    int m_star = n + 1;
    for (int m = 1; m <= n; ++m) {
        eng.add(perm[n - m]);
        if (eng.cross()) {
            m_star = m;
            break;
        }
    }
    const int k_minus = n - m_star + 1;
    return {k_minus, k_plus};
}

using JointKey = std::uint64_t;
JointKey joint_key(int km, int kp) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(km)) << 32) |
           static_cast<std::uint32_t>(kp);
}
std::pair<int, int> split_key(JointKey key) {
    return {static_cast<int>(key >> 32), static_cast<int>(key & 0xffffffffu)};
}

struct Acc {
    int n = 0;
    std::uint64_t replicas = 0;
    std::uint64_t kminus_le_kplus = 0;
    std::uint64_t kminus_gt_kplus = 0;
    std::vector<std::uint64_t> km_hist;
    std::vector<std::uint64_t> kp_hist;
    std::map<JointKey, std::uint64_t> joint;
    // first/second joint moments (integer sums over replicas)
    std::int64_t sum_km = 0, sum_kp = 0, sum_km2 = 0, sum_kp2 = 0, sum_kmkp = 0;
    std::int64_t sum_gap = 0, sum_gap2 = 0;

    explicit Acc(int n_sites = 0) { reset(n_sites); }

    void reset(int n_sites) {
        n = n_sites;
        replicas = 0;
        kminus_le_kplus = 0;
        kminus_gt_kplus = 0;
        const std::size_t m = static_cast<std::size_t>(n + 2);  // 0..N+1
        km_hist.assign(m, 0);
        kp_hist.assign(m, 0);
        joint.clear();
        sum_km = sum_kp = sum_km2 = sum_kp2 = sum_kmkp = 0;
        sum_gap = sum_gap2 = 0;
    }

    void add_ranks(int km, int kp) {
        replicas += 1;
        if (km < 0 || km > n + 1 || kp < 0 || kp > n + 1) {
            throw std::runtime_error("rank out of range");
        }
        km_hist[static_cast<std::size_t>(km)] += 1;
        kp_hist[static_cast<std::size_t>(kp)] += 1;
        joint[joint_key(km, kp)] += 1;
        const std::int64_t km64 = km;
        const std::int64_t kp64 = kp;
        sum_km += km64;
        sum_kp += kp64;
        sum_km2 += km64 * km64;
        sum_kp2 += kp64 * kp64;
        sum_kmkp += km64 * kp64;
        const std::int64_t gap = kp64 - km64;
        sum_gap += gap;
        sum_gap2 += gap * gap;
        if (km <= kp) kminus_le_kplus += 1;
        else kminus_gt_kplus += 1;
    }

    void absorb(const Acc& o) {
        replicas += o.replicas;
        kminus_le_kplus += o.kminus_le_kplus;
        kminus_gt_kplus += o.kminus_gt_kplus;
        for (std::size_t i = 0; i < km_hist.size(); ++i) {
            km_hist[i] += o.km_hist[i];
            kp_hist[i] += o.kp_hist[i];
        }
        for (const auto& kv : o.joint) joint[kv.first] += kv.second;
        sum_km += o.sum_km;
        sum_kp += o.sum_kp;
        sum_km2 += o.sum_km2;
        sum_kp2 += o.sum_kp2;
        sum_kmkp += o.sum_kmkp;
        sum_gap += o.sum_gap;
        sum_gap2 += o.sum_gap2;
    }
};

struct CoupleAcc {
    std::uint64_t replicas = 0;
    std::int64_t sum_km1 = 0, sum_km2 = 0, sum_kp1 = 0, sum_kp2 = 0;
    std::int64_t sum_km1km2 = 0, sum_kp1kp2 = 0, sum_km1kp2 = 0, sum_kp1km2 = 0;
    std::int64_t sum_gap1gap2 = 0;

    void add(int km1, int kp1, int km2, int kp2) {
        replicas += 1;
        const std::int64_t a = km1, b = kp1, c = km2, d = kp2;
        sum_km1 += a;
        sum_kp1 += b;
        sum_km2 += c;
        sum_kp2 += d;
        sum_km1km2 += a * c;
        sum_kp1kp2 += b * d;
        sum_km1kp2 += a * d;
        sum_kp1km2 += b * c;
        sum_gap1gap2 += (b - a) * (d - c);
    }

    void absorb(const CoupleAcc& o) {
        replicas += o.replicas;
        sum_km1 += o.sum_km1;
        sum_km2 += o.sum_km2;
        sum_kp1 += o.sum_kp1;
        sum_kp2 += o.sum_kp2;
        sum_km1km2 += o.sum_km1km2;
        sum_kp1kp2 += o.sum_kp1kp2;
        sum_km1kp2 += o.sum_km1kp2;
        sum_kp1km2 += o.sum_kp1km2;
        sum_gap1gap2 += o.sum_gap1gap2;
    }
};

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

bool rng_kat_pass(std::ostream& log) {
    const auto kats = philox4x32_10_official_kats();
    bool all = true;
    for (const auto& kat : kats) {
        const auto got = philox4x32_10(kat.ctr, kat.key);
        bool pass = got == kat.expected;
        all = all && pass;
        log << "  kat " << (pass ? "PASS" : "FAIL") << " expected=" << hex32(kat.expected[0])
            << " got=" << hex32(got[0]) << "\n";
    }
    PhiloxStream a, b;
    a.reset(0x12345678u, 3u, 9u, 0u);
    b.reset(0x12345678u, 3u, 9u, 0u);
    bool stream_ok = true;
    for (int i = 0; i < 64; ++i) {
        if (a.next_u32() != b.next_u32()) stream_ok = false;
    }
    all = all && stream_ok;
    log << "  stream_determinism " << (stream_ok ? "PASS" : "FAIL") << "\n";
    return all;
}

// Microcanonical wrapping indicator from a k-subset prefix of perm.
bool subset_cross(Engine& eng, const CSR& csr, const int* sites, int k) {
    eng.begin(&csr);
    for (int i = 0; i < k; ++i) eng.add(sites[i]);
    return eng.cross();
}

std::vector<double> bernstein_from_hists(const Acc& acc) {
    // a_k / C(N,k) wait: a_k is the INTEGER sum of D over k-subsets.
    // From permutations: each k-subset appears k!(N-k)! times.
    // sum_perms D_at_k = a_k * k! * (N-k)!
    // D_at_k = 1{K+ <= k} - 1{K- > k}
    // We return the mean D at each k, i.e. E[D|occupation k] = a_k / C(N,k).
    const int n = acc.n;
    std::vector<double> mean_d(static_cast<std::size_t>(n + 1), 0.0);
    if (acc.replicas == 0) return mean_d;
    std::uint64_t cdf_plus = 0;
    std::uint64_t sf_minus = acc.replicas;  // P(K_minus > -1) = 1
    // At occupation k: P(K+ <= k) and P(K- > k)
    for (int k = 0; k <= n; ++k) {
        cdf_plus += acc.kp_hist[static_cast<std::size_t>(k)];
        sf_minus -= acc.km_hist[static_cast<std::size_t>(k)];
        // after subtracting hist[k], sf_minus = count(K_minus > k)
        const double p_plus = static_cast<double>(cdf_plus) / static_cast<double>(acc.replicas);
        const double p_white = static_cast<double>(sf_minus) / static_cast<double>(acc.replicas);
        mean_d[static_cast<std::size_t>(k)] = p_plus - p_white;
    }
    return mean_d;
}

double binomial_convolve(const std::vector<double>& qk, int n, double p) {
    if (p <= 0.0) return qk.front();
    if (p >= 1.0) return qk.back();
    const double ratio_up = p / (1.0 - p);
    int mode = static_cast<int>((n + 1) * p);
    if (mode < 0) mode = 0;
    if (mode > n) mode = n;
    double acc = 1.0;
    double s = qk[static_cast<std::size_t>(mode)];
    double tot = 1.0;
    double w = 1.0;
    for (int k = mode + 1; k <= n; ++k) {
        w *= static_cast<double>(n - k + 1) / static_cast<double>(k) * ratio_up;
        if (w == 0.0) break;
        tot += w;
        s += w * qk[static_cast<std::size_t>(k)];
    }
    w = 1.0;
    const double ratio_dn = (1.0 - p) / p;
    for (int k = mode - 1; k >= 0; --k) {
        w *= static_cast<double>(k + 1) / static_cast<double>(n - k) * ratio_dn;
        if (w == 0.0) break;
        tot += w;
        s += w * qk[static_cast<std::size_t>(k)];
    }
    return s / tot;
}

std::vector<int> exact_subset_mean_d_cross(const Geometry& g) {
    // Return integer a_k = sum_{|C|=k} D_cross(C).  N must be small.
    const int n = g.n;
    if (n > 16) throw std::invalid_argument("subset enumeration limited to N<=16");
    std::vector<int> a(static_cast<std::size_t>(n + 1), 0);
    Engine eng(n, g.period);
    std::vector<int> sites(static_cast<std::size_t>(n));
    const std::uint64_t total = 1ull << n;
    for (std::uint64_t mask = 0; mask < total; ++mask) {
        int k = 0;
        for (int i = 0; i < n; ++i) {
            if ((mask >> i) & 1ull) sites[static_cast<std::size_t>(k++)] = i;
        }
        const bool black = subset_cross(eng, g.primal, sites.data(), k);
        int w = 0;
        for (int i = 0; i < n; ++i) {
            if (((mask >> i) & 1ull) == 0ull) sites[static_cast<std::size_t>(w++)] = i;
        }
        const bool white = subset_cross(eng, g.matching, sites.data(), n - k);
        a[static_cast<std::size_t>(k)] += static_cast<int>(black) - static_cast<int>(white);
    }
    return a;
}

double eval_power(const std::vector<int>& coeff, double p) {
    double r = 0.0;
    for (int i = static_cast<int>(coeff.size()) - 1; i >= 0; --i) {
        r = r * p + static_cast<double>(coeff[static_cast<std::size_t>(i)]);
    }
    return r;
}

bool run_exact_tests(const std::string& outdir, std::string& json_out) {
    std::filesystem::create_directories(outdir);
    std::ostringstream j;
    bool all = true;
    auto flag = [&](bool p) {
        if (!p) all = false;
        return p ? "PASS" : "FAIL";
    };

    j << "{\n";
    j << "  \"off_by_one_convention\": {\n";
    j << "    \"K_plus\": \"smallest black occupation k at which primal CROSS wrapping is true; "
         "K_plus=N+1 if never\",\n";
    j << "    \"K_minus\": \"N-m_star+1 where m_star is first matching occupation (reverse perm) "
         "with CROSS wrapping; K_minus=0 if white never wraps\",\n";
    j << "    \"beta\": \"T|K=k ~ Beta(k, N+1-k) for k=1..N; K=0 => T=0; K=N+1 => T not in [0,1)\",\n";
    j << "    \"matching_function\": \"M(p)=P(K_plus<=m)-P(K_minus>m), m~Binomial(N,p)\",\n";
    j << "    \"channel\": \"cross; either is diagnostic only, not a second M replication\"\n";
    j << "  },\n";

    std::ostringstream kat_log;
    const bool kat = rng_kat_pass(kat_log);
    j << "  \"rng_kat\": {\"status\": \"" << flag(kat) << "\"},\n";
    std::cerr << "RNG KAT " << (kat ? "PASS" : "FAIL") << "\n" << kat_log.str();

    j << "  \"exhaustive\": {\n";
    bool first_item = true;
    struct Case {
        Geometry g;
        const char* label;
        std::vector<int> power;  // published either-wrap matching polynomial; identity => same as cross
        bool have_power;
    };
    std::vector<Case> cases;
    {
        Case c{axis_geometry(2), "axis_L2", {-1, 0, 4, 0, -2}, true};
        cases.push_back(std::move(c));
    }
    {
        Case c{axis_geometry(3), "axis_L3", {-1, 0, 0, 6, 0, 0, 0, -18, 18, -4}, true};
        cases.push_back(std::move(c));
    }
    {
        Case c{gaussian_geometry(2, 1), "gaussian_2_1", {}, false};
        cases.push_back(std::move(c));
    }

    for (Case& cs : cases) {
        const Geometry& g = cs.g;
        const int n = g.n;
        Acc acc(n);
        Engine eng(n, g.period);
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::iota(perm.begin(), perm.end(), 0);
        do {
            const auto [km, kp] = threshold_ranks(eng, g, perm.data());
            acc.add_ranks(km, kp);
        } while (std::next_permutation(perm.begin(), perm.end()));

        const std::vector<int> a_subset = exact_subset_mean_d_cross(g);
        // Convert permutation means to integer a_k via a_k = mean_d[k] * C(N,k)
        // but we compare permutation-implied a_k with subset a_k exactly:
        // sum_perms D(k) = a_k * k! * (N-k)!
        // replicas = N!, so mean_d[k] * N! = a_k * k! * (N-k)!
        // a_k = mean_d[k] * C(N,k)
        const std::vector<double> mean_d = bernstein_from_hists(acc);
        bool poly_ok = true;
        std::vector<int> a_from_perm(static_cast<std::size_t>(n + 1), 0);
        for (int k = 0; k <= n; ++k) {
            // C(n,k) * mean_d should be an integer equal to a_subset[k]
            double comb = 1.0;
            for (int i = 0; i < k; ++i) {
                comb *= static_cast<double>(n - i) / static_cast<double>(i + 1);
            }
            const double ak = mean_d[static_cast<std::size_t>(k)] * comb;
            const int aki = static_cast<int>(std::llround(ak));
            a_from_perm[static_cast<std::size_t>(k)] = aki;
            if (std::fabs(ak - static_cast<double>(a_subset[static_cast<std::size_t>(k)])) > 1e-6) {
                poly_ok = false;
            }
        }

        bool ineq = acc.kminus_gt_kplus == 0;
        bool power_ok = true;
        double max_m_err = 0.0;
        if (cs.have_power) {
            for (double p : {0.1, 0.3, 0.5, 0.541196100146197, 0.586511455112676, 0.7, 0.9}) {
                const double from_hist = binomial_convolve(mean_d, n, p);
                const double from_power = eval_power(cs.power, p);
                max_m_err = std::max(max_m_err, std::fabs(from_hist - from_power));
            }
            power_ok = max_m_err < 1e-10;
        }

        // Numerical derivative at p=0.5 vs finite difference of reconstructed M
        auto M = [&](double p) { return binomial_convolve(mean_d, n, p); };
        const double h = 1e-5;
        const double m1 = (M(0.5 + h) - M(0.5 - h)) / (2.0 * h);
        double m1_ref = 0.0;
        if (cs.have_power) {
            m1_ref = (eval_power(cs.power, 0.5 + h) - eval_power(cs.power, 0.5 - h)) / (2.0 * h);
        }
        const bool deriv_ok = !cs.have_power || std::fabs(m1 - m1_ref) < 1e-6;

        const bool pass = poly_ok && ineq && power_ok && deriv_ok;
        flag(pass);
        if (!first_item) j << ",\n";
        first_item = false;
        j << "    \"" << cs.label << "\": {\"status\": \"" << (pass ? "PASS" : "FAIL")
          << "\", \"N\": " << n << ", \"permutations\": " << acc.replicas
          << ", \"Kminus_le_Kplus\": " << acc.kminus_le_kplus
          << ", \"Kminus_gt_Kplus\": " << acc.kminus_gt_kplus
          << ", \"reconstruct_subset_ak\": " << (poly_ok ? "true" : "false")
          << ", \"power_max_abs_err\": " << std::setprecision(16) << max_m_err
          << ", \"derivative_p05_ok\": " << (deriv_ok ? "true" : "false")
          << ", \"a_k_cross\": [";
        for (int k = 0; k <= n; ++k) {
            if (k) j << ", ";
            j << a_subset[static_cast<std::size_t>(k)];
        }
        j << "]}";
        std::cerr << "exact " << cs.label << " " << (pass ? "PASS" : "FAIL")
                  << " nperms=" << acc.replicas << " ineq_fail=" << acc.kminus_gt_kplus
                  << " max_M_err=" << max_m_err << "\n";
    }
    j << "\n  },\n";

    // Axis L=4: 2^16 subset a_k vs Monte Carlo reconstruction (not 16! perms).
    {
        const Geometry g = axis_geometry(4);
        const int n = g.n;
        const std::vector<int> a_subset = exact_subset_mean_d_cross(g);
        std::vector<double> mean_d(static_cast<std::size_t>(n + 1), 0.0);
        std::vector<double> comb(static_cast<std::size_t>(n + 1), 1.0);
        for (int k = 1; k <= n; ++k) {
            comb[static_cast<std::size_t>(k)] =
                comb[static_cast<std::size_t>(k - 1)] * static_cast<double>(n - k + 1) /
                static_cast<double>(k);
        }
        for (int k = 0; k <= n; ++k) {
            mean_d[static_cast<std::size_t>(k)] =
                static_cast<double>(a_subset[static_cast<std::size_t>(k)]) /
                comb[static_cast<std::size_t>(k)];
        }
        const std::vector<int> power = {-1, 0, 0, 0, 8, 0, 32, -64, 172, -704, 1104, -608, -56, 128, 16, -32, 6};
        double max_err = 0.0;
        for (double p : {0.2, 0.4, 0.5, 0.590672112331028, 0.7}) {
            max_err = std::max(max_err, std::fabs(binomial_convolve(mean_d, n, p) - eval_power(power, p)));
        }
        Acc acc(n);
        Engine eng(n, g.period);
        std::vector<int> perm(static_cast<std::size_t>(n));
        PhiloxStream rng;
        const int nrep = 20000;
        for (int r = 0; r < nrep; ++r) {
            rng.reset(kMasterSeed, 99u, static_cast<uint32_t>(r), 0u);
            fisher_yates(perm, rng);
            const auto [km, kp] = threshold_ranks(eng, g, perm.data());
            acc.add_ranks(km, kp);
        }
        const std::vector<double> mc_d = bernstein_from_hists(acc);
        double max_mc = 0.0;
        for (int k = 0; k <= n; ++k) {
            max_mc = std::max(max_mc, std::fabs(mc_d[static_cast<std::size_t>(k)] -
                                                mean_d[static_cast<std::size_t>(k)]));
        }
        const bool ineq = acc.kminus_gt_kplus == 0;
        const bool pass = max_err < 1e-10 && ineq && max_mc < 0.05;
        flag(pass);
        j << "  \"axis_L4\": {\"status\": \"" << (pass ? "PASS" : "FAIL")
          << "\", \"subset_vs_published_power_max_abs\": " << std::setprecision(16) << max_err
          << ", \"mc_replicas\": " << nrep << ", \"max_meanD_abs_err\": " << max_mc
          << ", \"Kminus_gt_Kplus\": " << acc.kminus_gt_kplus << "},\n";
        std::cerr << "exact axis_L4 " << (pass ? "PASS" : "FAIL") << " power_err=" << max_err
                  << " mc_meanD_err=" << max_mc << " ineq_fail=" << acc.kminus_gt_kplus << "\n";
    }

    j << "  \"overall\": \"" << (all ? "PASS" : "FAIL") << "\"\n";
    j << "}\n";
    json_out = j.str();
    std::ofstream out(outdir + "/exact_tests.json");
    out << json_out;
    std::cerr << "exact tests " << (all ? "PASS" : "FAIL") << " wrote " << outdir
              << "/exact_tests.json\n";
    return all;
}

void write_hist_csv(const std::string& path, const std::vector<std::uint64_t>& hist, const char* name) {
    std::ofstream out(path);
    out << "k," << name << "\n";
    for (std::size_t k = 0; k < hist.size(); ++k) out << k << ',' << hist[k] << '\n';
}

void write_joint_csv(const std::string& path, const Acc& acc) {
    std::ofstream out(path);
    out << "kminus,kplus,count\n";
    for (const auto& kv : acc.joint) {
        const auto [km, kp] = split_key(kv.first);
        out << km << ',' << kp << ',' << kv.second << '\n';
    }
}

void write_moments_row(std::ofstream& out, const std::string& geom, int batch, const Acc& a) {
    out << geom << ',' << batch << ',' << a.replicas << ',' << a.n << ',' << a.sum_km << ','
        << a.sum_kp << ',' << a.sum_km2 << ',' << a.sum_kp2 << ',' << a.sum_kmkp << ','
        << a.sum_gap << ',' << a.sum_gap2 << ',' << a.kminus_le_kplus << ',' << a.kminus_gt_kplus
        << '\n';
}

struct Options {
    bool exact = false;
    std::vector<int> axis_L;
    std::vector<std::pair<int, int>> pair;  // two (a,b)
    std::uint64_t samples = 1000000;
    int batches = 40;
    int threads = 8;
    uint32_t seed = kMasterSeed;
    int batch_begin = 0;
    std::string outdir = "results/server-20260828/C05";
};

void usage(const char* a0) {
    std::cerr << "Usage: " << a0
              << " [--exact-tests] [--axis L]... [--pair a1,b1 a2,b2]\n"
              << "       [--samples N] [--batches B] [--threads T] [--seed U32]\n"
              << "       [--batch-begin I] [--outdir DIR]\n";
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
        if (a == "--exact-tests") {
            opt.exact = true;
        } else if (a == "--axis") {
            const char* v = need("--axis");
            if (!v) return false;
            opt.axis_L.push_back(std::stoi(v));
        } else if (a == "--pair") {
            if (i + 2 >= argc) {
                std::cerr << "--pair needs two a,b tokens\n";
                return false;
            }
            auto parse_ab = [](const std::string& s) {
                const auto c = s.find(',');
                if (c == std::string::npos) throw std::invalid_argument("expected a,b");
                return std::pair<int, int>{std::stoi(s.substr(0, c)), std::stoi(s.substr(c + 1))};
            };
            opt.pair.push_back(parse_ab(argv[++i]));
            opt.pair.push_back(parse_ab(argv[++i]));
        } else if (a == "--samples") {
            const char* v = need("--samples");
            if (!v) return false;
            opt.samples = std::stoull(v);
        } else if (a == "--batches") {
            const char* v = need("--batches");
            if (!v) return false;
            opt.batches = std::stoi(v);
        } else if (a == "--threads") {
            const char* v = need("--threads");
            if (!v) return false;
            opt.threads = std::stoi(v);
        } else if (a == "--seed") {
            const char* v = need("--seed");
            if (!v) return false;
            opt.seed = static_cast<uint32_t>(std::stoul(v, nullptr, 0));
        } else if (a == "--batch-begin") {
            const char* v = need("--batch-begin");
            if (!v) return false;
            opt.batch_begin = std::stoi(v);
        } else if (a == "--outdir") {
            const char* v = need("--outdir");
            if (!v) return false;
            opt.outdir = v;
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return false;
        } else {
            std::cerr << "unknown argument: " << a << "\n";
            usage(argv[0]);
            return false;
        }
    }
    if (opt.threads < 1) opt.threads = 1;
    if (opt.threads > 8) opt.threads = 8;
    if (opt.batches < 2) opt.batches = 2;
    return true;
}

void write_geom_meta(const std::string& path, const Geometry& g, const Options& opt,
                     std::uint64_t samples, double wall, double cpu, long rss) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "{\n";
    out << "  \"name\": \"" << g.name << "\",\n";
    out << "  \"N\": " << g.n << ",\n";
    out << "  \"a\": " << g.a << ",\n";
    out << "  \"b\": " << g.b << ",\n";
    out << "  \"L\": " << g.L << ",\n";
    out << "  \"theta\": " << g.theta << ",\n";
    out << "  \"cos4\": " << g.cos4 << ",\n";
    out << "  \"period_matrix\": [[" << g.period.a00 << ", " << g.period.a01 << "], ["
        << g.period.a10 << ", " << g.period.a11 << "]],\n";
    out << "  \"det\": " << g.period.det() << ",\n";
    out << "  \"channel\": \"cross\",\n";
    out << "  \"rng_algorithm\": \"Philox4x32-10\",\n";
    out << "  \"rng_seed\": " << opt.seed << ",\n";
    out << "  \"samples\": " << samples << ",\n";
    out << "  \"batches\": " << opt.batches << ",\n";
    out << "  \"batch_begin\": " << opt.batch_begin << ",\n";
    out << "  \"threads\": " << opt.threads << ",\n";
    out << "  \"wall_seconds\": " << wall << ",\n";
    out << "  \"cpu_seconds\": " << cpu << ",\n";
    out << "  \"peak_rss_kb\": " << rss << "\n";
    out << "}\n";
}

int run_single(const Geometry& g, const Options& opt) {
    const int n = g.n;
    if (opt.samples % static_cast<std::uint64_t>(opt.batches) != 0) {
        std::cerr << "samples must be divisible by batches\n";
        return 2;
    }
    const int rpb = static_cast<int>(opt.samples / static_cast<std::uint64_t>(opt.batches));
    const int nthreads = opt.threads;
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    const std::string gdir = opt.outdir + "/" + g.name;
    std::filesystem::create_directories(gdir);
    std::ofstream moments(gdir + "/batch_moments.csv");
    moments << "geom,batch,replicas,n,sum_km,sum_kp,sum_km2,sum_kp2,sum_kmkp,sum_gap,sum_gap2,"
               "kminus_le_kplus,kminus_gt_kplus\n";
    std::ofstream bkm(gdir + "/batch_kminus_hist.csv");
    std::ofstream bkp(gdir + "/batch_kplus_hist.csv");
    bkm << "batch,k,count\n";
    bkp << "batch,k,count\n";

    Acc pooled(n);
    const double cpu0 = cpu_seconds();
    const auto t0 = std::chrono::steady_clock::now();

    for (int bi = 0; bi < opt.batches; ++bi) {
        const uint32_t batch_id = static_cast<uint32_t>(opt.batch_begin + bi);
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
            Engine eng(n, g.period);
            std::vector<int> perm(static_cast<std::size_t>(n));
            Acc& acc = local[static_cast<std::size_t>(tid)];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int r = 0; r < rpb; ++r) {
                PhiloxStream rng;
                rng.reset(opt.seed, batch_id, static_cast<uint32_t>(r), 0u);
                fisher_yates(perm, rng);
                const auto [km, kp] = threshold_ranks(eng, g, perm.data());
                acc.add_ranks(km, kp);
            }
        }
        Acc bout(n);
        for (const Acc& a : local) bout.absorb(a);
        pooled.absorb(bout);
        write_moments_row(moments, g.name, static_cast<int>(batch_id), bout);
        for (int k = 0; k <= n + 1; ++k) {
            bkm << batch_id << ',' << k << ',' << bout.km_hist[static_cast<std::size_t>(k)] << '\n';
            bkp << batch_id << ',' << k << ',' << bout.kp_hist[static_cast<std::size_t>(k)] << '\n';
        }
        std::cerr << "  " << g.name << " batch=" << batch_id << " replicas=" << bout.replicas
                  << " ineq_fail=" << bout.kminus_gt_kplus << "\n";
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();
    const double cpu = cpu_seconds() - cpu0;
    const long rss = peak_rss_kb();
    write_hist_csv(gdir + "/kminus_hist.csv", pooled.km_hist, "count");
    write_hist_csv(gdir + "/kplus_hist.csv", pooled.kp_hist, "count");
    write_joint_csv(gdir + "/joint_hist.csv", pooled);
    write_geom_meta(gdir + "/run_meta.json", g, opt, pooled.replicas, wall, cpu, rss);
    std::cerr << "done " << g.name << " replicas=" << pooled.replicas << " wall=" << wall
              << "s ineq_fail=" << pooled.kminus_gt_kplus << "\n";
    return pooled.kminus_gt_kplus ? 1 : 0;
}

int run_pair(const Geometry& g1, const Geometry& g2, const Options& opt) {
    if (g1.n != g2.n) throw std::invalid_argument("same-N pair required");
    const int n = g1.n;
    if (opt.samples % static_cast<std::uint64_t>(opt.batches) != 0) {
        std::cerr << "samples must be divisible by batches\n";
        return 2;
    }
    const int rpb = static_cast<int>(opt.samples / static_cast<std::uint64_t>(opt.batches));
    const int nthreads = opt.threads;
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
    const std::string tag = "n" + std::to_string(n);
    const std::string d1 = opt.outdir + "/" + g1.name;
    const std::string d2 = opt.outdir + "/" + g2.name;
    const std::string dc = opt.outdir + "/coupling_" + tag;
    std::filesystem::create_directories(d1);
    std::filesystem::create_directories(d2);
    std::filesystem::create_directories(dc);

    auto open_moments = [](const std::string& dir) {
        std::ofstream out(dir + "/batch_moments.csv");
        out << "geom,batch,replicas,n,sum_km,sum_kp,sum_km2,sum_kp2,sum_kmkp,sum_gap,sum_gap2,"
               "kminus_le_kplus,kminus_gt_kplus\n";
        return out;
    };
    std::ofstream m1 = open_moments(d1);
    std::ofstream m2 = open_moments(d2);
    std::ofstream couple(dc + "/batch_coupling.csv");
    couple << "batch,replicas,sum_km1,sum_km2,sum_kp1,sum_kp2,sum_km1km2,sum_kp1kp2,sum_km1kp2,"
              "sum_kp1km2,sum_gap1gap2\n";
    std::ofstream bkm1(d1 + "/batch_kminus_hist.csv");
    std::ofstream bkp1(d1 + "/batch_kplus_hist.csv");
    std::ofstream bkm2(d2 + "/batch_kminus_hist.csv");
    std::ofstream bkp2(d2 + "/batch_kplus_hist.csv");
    bkm1 << "batch,k,count\n";
    bkp1 << "batch,k,count\n";
    bkm2 << "batch,k,count\n";
    bkp2 << "batch,k,count\n";

    Acc pooled1(n), pooled2(n);
    CoupleAcc pooled_c;
    const double cpu0 = cpu_seconds();
    const auto t0 = std::chrono::steady_clock::now();

    for (int bi = 0; bi < opt.batches; ++bi) {
        const uint32_t batch_id = static_cast<uint32_t>(opt.batch_begin + bi);
        std::vector<Acc> loc1(static_cast<std::size_t>(nthreads), Acc(n));
        std::vector<Acc> loc2(static_cast<std::size_t>(nthreads), Acc(n));
        std::vector<CoupleAcc> locc(static_cast<std::size_t>(nthreads));
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
        {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            Engine e1(n, g1.period);
            Engine e2(n, g2.period);
            std::vector<int> perm(static_cast<std::size_t>(n));
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int r = 0; r < rpb; ++r) {
                PhiloxStream rng;
                rng.reset(opt.seed, batch_id, static_cast<uint32_t>(r), 0u);
                fisher_yates(perm, rng);
                const auto r1 = threshold_ranks(e1, g1, perm.data());
                const auto r2 = threshold_ranks(e2, g2, perm.data());
                loc1[static_cast<std::size_t>(tid)].add_ranks(r1.first, r1.second);
                loc2[static_cast<std::size_t>(tid)].add_ranks(r2.first, r2.second);
                locc[static_cast<std::size_t>(tid)].add(r1.first, r1.second, r2.first, r2.second);
            }
        }
        Acc b1(n), b2(n);
        CoupleAcc bc;
        for (int t = 0; t < nthreads; ++t) {
            b1.absorb(loc1[static_cast<std::size_t>(t)]);
            b2.absorb(loc2[static_cast<std::size_t>(t)]);
            bc.absorb(locc[static_cast<std::size_t>(t)]);
        }
        pooled1.absorb(b1);
        pooled2.absorb(b2);
        pooled_c.absorb(bc);
        write_moments_row(m1, g1.name, static_cast<int>(batch_id), b1);
        write_moments_row(m2, g2.name, static_cast<int>(batch_id), b2);
        couple << batch_id << ',' << bc.replicas << ',' << bc.sum_km1 << ',' << bc.sum_km2 << ','
               << bc.sum_kp1 << ',' << bc.sum_kp2 << ',' << bc.sum_km1km2 << ',' << bc.sum_kp1kp2
               << ',' << bc.sum_km1kp2 << ',' << bc.sum_kp1km2 << ',' << bc.sum_gap1gap2 << '\n';
        for (int k = 0; k <= n + 1; ++k) {
            bkm1 << batch_id << ',' << k << ',' << b1.km_hist[static_cast<std::size_t>(k)] << '\n';
            bkp1 << batch_id << ',' << k << ',' << b1.kp_hist[static_cast<std::size_t>(k)] << '\n';
            bkm2 << batch_id << ',' << k << ',' << b2.km_hist[static_cast<std::size_t>(k)] << '\n';
            bkp2 << batch_id << ',' << k << ',' << b2.kp_hist[static_cast<std::size_t>(k)] << '\n';
        }
        std::cerr << "  pair N=" << n << " batch=" << batch_id << " replicas=" << b1.replicas
                  << " ineq_fail=" << (b1.kminus_gt_kplus + b2.kminus_gt_kplus) << "\n";
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double wall = std::chrono::duration<double>(t1 - t0).count();
    const double cpu = cpu_seconds() - cpu0;
    const long rss = peak_rss_kb();
    write_hist_csv(d1 + "/kminus_hist.csv", pooled1.km_hist, "count");
    write_hist_csv(d1 + "/kplus_hist.csv", pooled1.kp_hist, "count");
    write_joint_csv(d1 + "/joint_hist.csv", pooled1);
    write_hist_csv(d2 + "/kminus_hist.csv", pooled2.km_hist, "count");
    write_hist_csv(d2 + "/kplus_hist.csv", pooled2.kp_hist, "count");
    write_joint_csv(d2 + "/joint_hist.csv", pooled2);
    write_geom_meta(d1 + "/run_meta.json", g1, opt, pooled1.replicas, wall, cpu, rss);
    write_geom_meta(d2 + "/run_meta.json", g2, opt, pooled2.replicas, wall, cpu, rss);
    std::ofstream cm(dc + "/coupling_moments.json");
    cm << std::setprecision(17);
    cm << "{\n  \"N\": " << n << ",\n  \"replicas\": " << pooled_c.replicas
       << ",\n  \"coupling\": \"CRN same_U_j: one Fisher-Yates permutation of cyclic labels\",\n"
       << "  \"sum_km1\": " << pooled_c.sum_km1 << ",\n  \"sum_km2\": " << pooled_c.sum_km2
       << ",\n  \"sum_kp1\": " << pooled_c.sum_kp1 << ",\n  \"sum_kp2\": " << pooled_c.sum_kp2
       << ",\n  \"sum_km1km2\": " << pooled_c.sum_km1km2 << ",\n  \"sum_kp1kp2\": " << pooled_c.sum_kp1kp2
       << ",\n  \"sum_km1kp2\": " << pooled_c.sum_km1kp2 << ",\n  \"sum_kp1km2\": " << pooled_c.sum_kp1km2
       << ",\n  \"sum_gap1gap2\": " << pooled_c.sum_gap1gap2 << ",\n  \"wall_seconds\": " << wall
       << "\n}\n";
    std::cerr << "done pair N=" << n << " " << g1.name << "/" << g2.name
              << " replicas=" << pooled1.replicas << " wall=" << wall << "s\n";
    return (pooled1.kminus_gt_kplus || pooled2.kminus_gt_kplus) ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) return 2;
#ifdef _OPENMP
    omp_set_num_threads(opt.threads);
#endif
    std::filesystem::create_directories(opt.outdir);
    int rc = 0;
    if (opt.exact) {
        std::string json;
        if (!run_exact_tests(opt.outdir + "/exact", json)) rc = 1;
        if (!opt.axis_L.empty() || !opt.pair.empty()) {
            if (rc) {
                std::cerr << "exact tests FAILED; not starting production histograms\n";
                return rc;
            }
        } else {
            return rc;
        }
    }
    if (!opt.axis_L.empty() && rc == 0) {
        for (int L : opt.axis_L) {
            const int prc = run_single(axis_geometry(L), opt);
            if (prc) rc = prc;
        }
    }
    if (opt.pair.size() == 2 && rc == 0) {
        const Geometry g1 = gaussian_geometry(opt.pair[0].first, opt.pair[0].second);
        const Geometry g2 = gaussian_geometry(opt.pair[1].first, opt.pair[1].second);
        const int prc = run_pair(g1, g2, opt);
        if (prc) rc = prc;
    } else if (opt.pair.size() != 0 && opt.pair.size() != 2) {
        std::cerr << "--pair requires exactly two orientations\n";
        return 2;
    }
    if (!opt.exact && opt.axis_L.empty() && opt.pair.empty()) {
        usage(argv[0]);
        return 2;
    }
    return rc;
}
