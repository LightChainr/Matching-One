// P11 engine: exact / high-precision transfer-matrix DP for square-lattice
// site percolation spanning counts A_{n,m}(k)  (Mertens 2022 definitions,
// arXiv:2109.12102; spanning = one occupied 4-neighbour cluster touching the
// virtual all-occupied row 0 and row m of the free-boundary n x m grid).
//
// Sweep design: the state map is advanced ONE CELL at a time (extend(sigma,c)
// semantics, paper Sec. 3.5), so per-row work is O(n * S) not O(2^n * S).
// State = unsigned __int128, 5 bits per cell:  0 empty, 1 top-connected,
// >=2 non-top cluster label (canonical by first occurrence; partial states
// are canonicalized at every occupied-cell step because they are map keys).
//
// Drop rules (exactly two):
//   * the top cluster loses its last cell with no new-row presence -> drop
//     (no future cluster can reach the virtual row 0 any more);
//   * a non-top cluster losing all cells simply vanishes (a finished
//     non-spanning cluster) -- the configuration stays valid.
//
// Modes:
//   crt   exact A(k) modulo primes (Chinese remainder outside); per-row
//         telemetry on stderr; optional map dump/load for checkpoints.
//   f128  __float128 evaluation of R(p), R'(p); Newton root finding for
//         p_med / p_cell with final bracket verification.
//
// All arithmetic in crt mode is exact modular integer arithmetic. In f128
// mode all rounding is IEEE quad precision; the achieved digit count is
// audited against the crt/exact route at widths where both run.

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#ifdef HAVE_QUADMATH
#include <quadmath.h>
#endif

typedef unsigned __int128 u128;
#if defined(HAVE_QUADMATH)
typedef __float128 f128;
#elif defined(HAVE_LDBL128)
// aarch64 Linux: long double is IEEE binary128 (113-bit mantissa)
typedef long double f128;
#else
typedef double f128;        // placeholder only; f128/res modes refuse to run
#endif

static const int BITS = 5, CELLMAX = 32;
static const uint64_t EMPTY = 0, TOP = 1;

static inline uint64_t cell_of(u128 s, int i) {
    return (uint64_t)((s >> (BITS * i)) & (CELLMAX - 1));
}
static inline void set_cell(u128& s, int i, uint64_t v) {
    s &= ~(((u128)(uint64_t)(CELLMAX - 1)) << (BITS * i));
    s |= ((u128)v) << (BITS * i);
}

struct StateHash {
    size_t operator()(u128 s) const {
        uint64_t lo = (uint64_t)s, hi = (uint64_t)(s >> 64);
        lo ^= lo >> 33; lo *= 0xff51afd7ed558ccdULL;
        hi ^= hi >> 33; hi *= 0xc4ceb9fe1a85ec53ULL;
        uint64_t x = lo ^ hi;
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33;
        return (size_t)x;
    }
};

static int canonical_relabel(u128& s, int n) {
    int mapv[CELLMAX];
    for (int i = 0; i < CELLMAX; ++i) mapv[i] = -1;
    int next = 0;
    u128 out = 0;
    for (int i = 0; i < n; ++i) {
        uint64_t c = cell_of(s, i);
        int v;
        if (c == EMPTY) v = 0;
        else if (c == TOP) v = 1;
        else {
            if (mapv[c] < 0) mapv[c] = 2 + (next++);
            v = mapv[c];
        }
        out |= ((u128)v) << (BITS * i);
    }
    s = out;
    return 2 + next;
}

static inline bool occurs(u128 s, int a, int b, uint64_t L) {
    for (int i = a; i < b; ++i)
        if (cell_of(s, i) == L) return true;
    return false;
}

// ---------------------------------------------------------------------------
// value policies
// ---------------------------------------------------------------------------

struct CrtVal {
    std::vector<uint64_t> c;   // coefficients mod p, compact
    void add(const CrtVal& o, uint64_t mod) {
        if (o.c.size() > c.size()) c.resize(o.c.size(), 0);
        for (size_t i = 0; i < o.c.size(); ++i) {
            c[i] += o.c[i];
            if (c[i] >= mod) c[i] -= mod;
        }
    }
    void shift1() { c.insert(c.begin(), 0); }
    void trim() { while (!c.empty() && c.back() == 0) c.pop_back(); }
};

struct F128Val {
    f128 v = 0, dv = 0;
    void add(const F128Val& o, uint64_t) { v += o.v; dv += o.dv; }
    void shift1() {}
    void trim() {}
};

// ---------------------------------------------------------------------------

template <class Val>
struct Sweeper {
    int n = 0, m = 0;
    uint64_t pmod = 0;    // prime (crt) or 0 (f128)
    f128 p = 0, qp = 0;   // f128 branch weights
    long long extends_cnt = 0, insert_cnt = 0;

    typedef std::unordered_map<u128, Val, StateHash> Map;

    // advance every partial state across cell c
    void cell_step(Map& cur, Map& next, int c) {
        next.clear();
        next.reserve(cur.size() * 2 + 16);
        for (auto& kv : cur) {
            const u128 s = kv.first;
            const Val& src = kv.second;
            const uint64_t oldc = cell_of(s, c);
            const uint64_t left = (c > 0) ? cell_of(s, c - 1) : EMPTY;
            ++extends_cnt;

            // ---- branch b=0 evaluated first (it takes src's coefficient
            //      vector by move in crt mode); b=1 then reuses tmp0's data.
            Val tmp0;
            bool b0_insert = false;
            {
                bool drop = false;
                if (oldc != EMPTY && !occurs(s, c + 1, n, oldc) &&
                    !occurs(s, 0, c, oldc)) {
                    if (oldc == TOP) drop = true;
                    // non-top: the finished cluster just vanishes
                }
                if (!drop) {
                    u128 s2 = s;
                    set_cell(s2, c, EMPTY);
                    // emptying a cell can remove a label's last occurrence and
                    // leave a gap in the label numbering; re-canonicalize so
                    // equivalent states share one map key
                    canonical_relabel(s2, n);
                    if constexpr (std::is_same_v<Val, F128Val>) {
                        tmp0.v = src.v * qp;
                        tmp0.dv = src.dv * qp - src.v;
                    } else {
                        tmp0.c = std::move(const_cast<Val&>(src).c);
                    }
                    next[s2].add(tmp0, pmod);
                    b0_insert = true;
                    ++insert_cnt;
                }
            }
            // ---- branch b=1 : new cell c occupied
            {
                u128 s2 = s;
                uint64_t lab;
                if (oldc == EMPTY && left == EMPTY) {
                    int mx = 1;
                    for (int i = 0; i < n; ++i)
                        mx = std::max(mx, (int)cell_of(s, i));
                    lab = (uint64_t)(mx + 1);
                } else if (oldc != EMPTY && left != EMPTY) {
                    if (oldc == TOP || left == TOP) {
                        lab = TOP;
                        uint64_t other = (oldc == TOP) ? left : oldc;
                        if (other != TOP)
                            for (int i = 0; i < n; ++i)
                                if (cell_of(s2, i) == other)
                                    set_cell(s2, i, TOP);
                    } else {
                        lab = left;
                        // relabel oldc across the WHOLE row: the class may
                        // also own new-row cells merged earlier in this sweep
                        for (int i = 0; i < n; ++i)
                            if (cell_of(s2, i) == oldc) set_cell(s2, i, lab);
                    }
                } else {
                    lab = (oldc != EMPTY) ? oldc : left;
                }
                set_cell(s2, c, lab);
                canonical_relabel(s2, n);
                Val tmp;
                if constexpr (std::is_same_v<Val, F128Val>) {
                    tmp.v = src.v * p;
                    tmp.dv = src.dv * p + src.v;  // d/dp (v_old * p)
                } else {
                    tmp.c = b0_insert ? tmp0.c
                                      : std::move(const_cast<Val&>(src).c);
                    tmp.shift1();
                }
                tmp.trim();
                next[s2].add(tmp, pmod);
                ++insert_cnt;
            }
        }
    }

    void row_step(Map& cur, Map& scratch) {
        for (int c = 0; c < n; ++c) {
            cell_step(cur, scratch, c);
            std::swap(cur, scratch);
        }
    }
};

// ---------------------------------------------------------------------------
// initial row-1 map: 2^n - 1 configurations, occupied sites top-connected
// ---------------------------------------------------------------------------

template <class Val>
static void init_row1(typename Sweeper<Val>::Map& map0, int n,
                      uint64_t mod, f128 p, f128 q) {
    for (uint64_t mask = 1; mask < (1ull << n); ++mask) {
        u128 s = 0;
        int k = 0;
        for (int j = 0; j < n; ++j)
            if ((mask >> j) & 1) { set_cell(s, j, TOP); ++k; }
        Val v;
        if constexpr (std::is_same_v<Val, CrtVal>) {
            v.c.assign(k + 1, 0);
            v.c[k] = 1 % mod;
        } else {
            f128 base = 1;
            for (int i = 0; i < k; ++i) base *= p;
            for (int i = 0; i < n - k; ++i) base *= q;
            f128 d = 0;
            if (k > 0) d += ((f128)k / p) * base;
            if (n - k > 0) d -= ((f128)(n - k) / q) * base;
            v.v = base; v.dv = d;
        }
        map0[s] = std::move(v);
    }
}

// ---------------------------------------------------------------------------
// checkpoint dump/load (crt): deterministic byte layout (sorted keys), so the
// file's SHA-256 (computed outside) is a deterministic state-enumeration hash
// ---------------------------------------------------------------------------

static void map_dump(const Sweeper<CrtVal>::Map& mp, const char* path,
                     int n, int row, uint64_t mod) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror("fopen"); exit(1); }
    uint32_t hdr[4] = {(uint32_t)n, (uint32_t)row, (uint32_t)mod,
                       (uint32_t)mp.size()};
    fwrite(hdr, sizeof(uint32_t), 4, f);
    std::vector<u128> keys;
    keys.reserve(mp.size());
    for (auto& kv : mp) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    for (u128 k : keys) {
        fwrite(&k, sizeof(u128), 1, f);
        const std::vector<uint64_t>& c = mp.at(k).c;
        uint32_t len = (uint32_t)c.size();
        fwrite(&len, sizeof(uint32_t), 1, f);
        if (len) fwrite(c.data(), sizeof(uint64_t), len, f);
    }
    fclose(f);
}

static void map_load(Sweeper<CrtVal>::Map& mp, const char* path,
                     int& n, int& row, uint64_t& mod) {
    FILE* f = fopen(path, "rb");
    if (!f) { perror("fopen"); exit(1); }
    uint32_t hdr[4];
    if (fread(hdr, sizeof(uint32_t), 4, f) != 4) { fprintf(stderr, "bad hdr\n"); exit(1); }
    n = hdr[0]; row = hdr[1]; mod = hdr[2];
    mp.clear();
    for (uint32_t i = 0; i < hdr[3]; ++i) {
        u128 k;
        if (fread(&k, sizeof(u128), 1, f) != 1) { fprintf(stderr, "bad key\n"); exit(1); }
        uint32_t len;
        if (fread(&len, sizeof(uint32_t), 1, f) != 1) { fprintf(stderr, "bad len\n"); exit(1); }
        CrtVal v; v.c.resize(len);
        if (len && fread(v.c.data(), sizeof(uint64_t), len, f) != len) {
            fprintf(stderr, "bad coeff\n"); exit(1);
        }
        mp[k] = std::move(v);
    }
    fclose(f);
}

static double now_s() {
    static auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// crt mode
// ---------------------------------------------------------------------------

int mode_crt(int n, int m, const char* map_in, const char* map_out,
             int dump_row) {
    // verified primes (sympy.isprime), pairwise coprime
    uint64_t primes[] = {1073741789ULL, 1073741783ULL, 1073741741ULL,
                         1073741723ULL, 1073741719ULL, 1073741717ULL,
                         1073741689ULL, 1073741671ULL, 1073741663ULL,
                         1073741651ULL};
    const int np = 10;
    int need = (n * m + 29) / 30;   // prod(p_i) > 2^(nm) >= max A(k)
    if (map_in) need = 1;  // restart: single pass with the checkpoint modulus
    if (need > np) { fprintf(stderr, "not enough primes\n"); return 1; }
    printf("{\"mode\":\"crt\",\"n\":%d,\"m\":%d,\"primes\":[", n, m);
    for (int i = 0; i < need; ++i)
        printf("%s%llu", i ? "," : "", (unsigned long long)primes[i]);
    printf("],\"passes\":[\n");
    for (int ip = 0; ip < need; ++ip) {
        uint64_t mod = primes[ip];
        Sweeper<CrtVal> sw; sw.n = n; sw.m = m; sw.pmod = mod;
        Sweeper<CrtVal>::Map cur, scratch;
        int start_row = 1;
        if (map_in) map_load(cur, map_in, n, start_row, mod);
        else init_row1<CrtVal>(cur, n, mod, 0, 0);
        double t0 = now_s();
        sw.extends_cnt = sw.insert_cnt = 0;
        for (int r = start_row + 1; r <= m; ++r) {
            sw.row_step(cur, scratch);
            fprintf(stderr,
                    "{\"pass\":%d,\"row\":%d,\"states\":%zu,\"t\":%.2f,"
                    "\"extends\":%lld,\"inserts\":%lld}\n",
                    ip, r, cur.size(), now_s() - t0, sw.extends_cnt,
                    sw.insert_cnt);
            if (map_out && r == dump_row) {
                std::string path = std::string(map_out) + ".p" +
                                   std::to_string(ip);
                map_dump(cur, path.c_str(), n, r, mod);
            }
        }
        std::vector<uint64_t> A(n * m + 1, 0);
        for (auto& kv : cur) {
            const std::vector<uint64_t>& c = kv.second.c;
            for (size_t k = 0; k < c.size(); ++k) {
                A[k] += c[k];
                if (A[k] >= mod) A[k] -= mod;
            }
        }
        printf("{\"prime\":%llu,\"A\":[", (unsigned long long)mod);
        for (int k = 0; k <= n * m; ++k)
            printf("%s%llu", k ? "," : "", (unsigned long long)A[k]);
        printf("],\"peak_states\":%zu,\"extends\":%lld,\"inserts\":%lld,"
               "\"t\":%.2f}%s\n",
               cur.size(), sw.extends_cnt, sw.insert_cnt, now_s() - t0,
               ip + 1 < need ? "," : "");
        if (map_out && !dump_row) {
            std::string path = std::string(map_out) + ".p" + std::to_string(ip);
            map_dump(cur, path.c_str(), n, m, mod);
        }
    }
    printf("]}\n");
    return 0;
}

// ---------------------------------------------------------------------------
// f128 mode: Newton root finding with bracket verification
// ---------------------------------------------------------------------------

static std::string f128_str(f128 x, int digits) {
    (void)digits;
    char buf[512];
#if defined(HAVE_QUADMATH)
    quadmath_snprintf(buf, sizeof buf, "%.45Qg", x);
#elif defined(HAVE_LDBL128)
    snprintf(buf, sizeof buf, "%.45Lg", x);
#else
    snprintf(buf, sizeof buf, "%.17g", (double)x);
#endif
    return buf;
}

static f128 f128_from_str(const char* s) {
#if defined(HAVE_QUADMATH)
    return strtoflt128(s, nullptr);
#else
    return strtold(s, nullptr);
#endif
}

int mode_f128(int n, int kind, f128 p0, int iters, f128 bracket_d) {
#ifndef HAVE_QUADMATH
    if (sizeof(f128) < 16) {
        fprintf(stderr, "f128 mode needs binary128 long double or quadmath\n");
        return 1;
    }
#endif
    auto eval = [&](f128 p, f128& f, f128& df) {
        f = 0; df = 0;
        for (int w = 0; w < (kind == 1 ? 2 : 1); ++w) {
            int nn = (w == 0) ? n : n - 1;
            if (nn < 1) continue;
            Sweeper<F128Val> sw; sw.n = nn; sw.m = nn; sw.pmod = 0;
            sw.p = p; sw.qp = 1 - p;
            Sweeper<F128Val>::Map cur, scratch;
            init_row1<F128Val>(cur, nn, 0, p, 1 - p);
            for (int r = 2; r <= nn; ++r) sw.row_step(cur, scratch);
            f128 R = 0, dR = 0;
            for (auto& kv : cur) { R += kv.second.v; dR += kv.second.dv; }
            if (w == 0) { f += R; df += dR; }
            else { f -= R; df -= dR; }
        }
        if (kind == 0) f -= (f128)0.5;
    };
    f128 p = p0, f = 0, df = 0;
    for (int it = 0; it < iters; ++it) {
        eval(p, f, df);
        f128 step = f / df;
        p -= step;
        fprintf(stderr, "{\"it\":%d,\"p\":%s,\"f\":%s,\"step\":%s}\n", it,
                f128_str(p, 45).c_str(), f128_str(f, 12).c_str(),
                f128_str(step, 12).c_str());
        if (fabsl((long double)step) < 1e-36L) break;
    }
    f128 fl, fr;
    eval(p - bracket_d, fl, df);
    eval(p + bracket_d, fr, df);
    bool ok = (fl < 0 && fr > 0) || (fl > 0 && fr < 0);
    printf("{\"mode\":\"f128\",\"kind\":\"%s\",\"n\":%d,"
           "\"p_root\":%s,\"f_left\":%s,\"f_right\":%s,"
           "\"bracket_d\":%s,\"bracket_ok\":%s}\n",
           kind == 0 ? "med" : "cell", n,
           f128_str(p, 45).c_str(), f128_str(fl, 12).c_str(),
           f128_str(fr, 12).c_str(), f128_str(bracket_d, 6).c_str(),
           ok ? "true" : "false");
    return ok ? 0 : 2;
}

// ---------------------------------------------------------------------------
// resource telemetry mode: one evaluation at fixed p with per-row stats
// ---------------------------------------------------------------------------

int mode_res(int n, f128 p) {
    Sweeper<F128Val> sw; sw.n = n; sw.m = n; sw.pmod = 0;
    sw.p = p; sw.qp = 1 - p;
    Sweeper<F128Val>::Map cur, scratch;
    init_row1<F128Val>(cur, n, 0, p, 1 - p);
    printf("{\"mode\":\"res\",\"n\":%d,\"rows\":[\n", n);
    fprintf(stderr, "{\"row\":1,\"states\":%zu,\"t\":%.3f}\n", cur.size(),
            now_s());
    for (int r = 2; r <= n; ++r) {
        double t0 = now_s();
        sw.row_step(cur, scratch);
        printf("{\"row\":%d,\"states\":%zu,\"row_t\":%.3f,"
               "\"extends\":%lld,\"inserts\":%lld}%s\n",
               r, cur.size(), now_s() - t0, sw.extends_cnt, sw.insert_cnt,
               r < n ? "," : "");
    }
    printf("]}\n");
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: <crt|f128|res> ...\n"); return 1; }
    std::string mode = argv[1];
    std::vector<std::string> a(argv + 2, argv + argc);
    auto get = [&](const char* k, std::string def = "") -> std::string {
        for (size_t i = 0; i + 1 < a.size(); ++i)
            if (a[i] == k) return a[i + 1];
        return def;
    };
    int n = atoi(get("--n", "0").c_str());
    int m = atoi(get("--m", "0").c_str()); if (!m) m = n;
    if (mode == "crt") {
        std::string mi = get("--map-in", ""), mo = get("--map-out", "");
        int dump_row = atoi(get("--dump-row", "0").c_str());
        return mode_crt(n, m, mi.empty() ? nullptr : mi.c_str(),
                        mo.empty() ? nullptr : mo.c_str(), dump_row);
    } else if (mode == "f128" || mode == "res") {
#if !defined(HAVE_QUADMATH) && !defined(HAVE_LDBL128)
        fprintf(stderr, "f128/res modes need binary128 (build with "
                        "-DHAVE_LDBL128 on aarch64 or quadmath)\n");
        return 1;
#else
        if (sizeof(f128) < 16) {
            fprintf(stderr, "long double is not binary128 here; refusing\n");
            return 1;
        }
        if (mode == "f128") {
            int kind = get("--kind", "med") == "cell" ? 1 : 0;
            f128 p0 = f128_from_str(get("--p0", "0.5927").c_str());
            int iters = atoi(get("--iters", "12").c_str());
            f128 bd = f128_from_str(get("--bracket", "1e-30").c_str());
            return mode_f128(n, kind, p0, iters, bd);
        }
        f128 p = f128_from_str(get("--p", "0.592746050792").c_str());
        return mode_res(n, p);
#endif
    }
    fprintf(stderr, "unknown mode\n");
    return 1;
}
