// Independent screening census for C(d, H) against one method interval.
//
// Written from the P2 protocol (manuscript section 3.2, Theorem 2) re-derived by
// hand: with S the scale, w_k = round(S * m^k) and c_0 = S, the integer
//
//     T(a) = sum_{k=0}^{d} c_k a_k
//
// satisfies |T(a) - S P_a(m)| <= H * rho, and if P_a has a root in [l, u] then
// |P_a(m)| <= D (u-l)/2 with D the global derivative bound. Hence every
// root-carrying a obeys |T(a)| <= B = ceil(S D (u-l)/2 + H rho), and screening on
// that inequality loses nothing.
//
// The coefficients c_k, the bound B and the scale S are all supplied by the
// driver, which computes them from the interval alone. This program never reads
// the repository's own screen.
//
// Two search paths, deliberately different code:
//
//   brute : nested enumeration of every tuple, with the LAST variable solved for
//           exactly rather than looped. Deterministic order, no sorting.
//   mitm  : split the coefficient vector into a left and a right block, sort the
//           (small) right block once, binary-search it from every left tuple.
//
// Both emit every a with |T(a)| <= B. `--shard k n` restricts the outermost
// variable so independent processes can cover disjoint pieces.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;
typedef long long ll;

static int D = 4, H = 100;
static ll C[8];          // C[0] = S, C[k] = w_k
static ll B = 0;
static int shardK = 0, shardN = 1;

static inline ll fdiv(ll a, ll b) {              // floor, b > 0
    ll q = a / b, r = a % b;
    if (r != 0 && r < 0) --q;
    return q;
}
static inline ll cdivp(ll a, ll b) {             // ceil, b > 0
    ll q = a / b, r = a % b;
    if (r != 0 && r > 0) ++q;
    return q;
}
static inline ll iabs(ll x) { return x < 0 ? -x : x; }

struct Rec { int a[8]; ll t; };
static vector<Rec> hits;
static ll scanned = 0;

static inline void emit(const int *a, ll t) {
    Rec r; memset(r.a, 0, sizeof r.a);
    for (int k = 0; k <= D; ++k) r.a[k] = a[k];
    r.t = t;
    hits.push_back(r);
}

static void brute_d4() {
    int a[8];
    for (a[0] = -H; a[0] <= H; ++a[0]) {
        if (((a[0] + H) % shardN) != shardK) continue;
        for (int a1 = -H; a1 <= H; ++a1) {
            ll r1 = C[1] * (ll)a1;
            for (int a2 = -H; a2 <= H; ++a2) {
                ll r2 = r1 + C[2] * (ll)a2;
                for (int a3 = -H; a3 <= H; ++a3) {
                    ll r3 = r2 + C[3] * (ll)a3;
                    for (int a4 = 1; a4 <= H; ++a4) {
                        ll t = C[0] * (ll)a[0] + r3 + C[4] * (ll)a4;
                        ++scanned;
                        if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = a3; a[4] = a4; emit(a, t); }
                    }
                }
            }
        }
    }
}

// Same enumeration, but the innermost variable is solved instead of looped.
// C[4] > 0, so |r3 + C[4] a4| <= B  <=>  a4 in [ceil((-B-r3)/C4), floor((B-r3)/C4)].
static void brute_d4_solve() {
    int a[8];
    for (a[0] = -H; a[0] <= H; ++a[0]) {
        if (((a[0] + H) % shardN) != shardK) continue;
        ll base0 = C[0] * (ll)a[0];
        for (int a1 = -H; a1 <= H; ++a1) {
            ll r1 = C[1] * (ll)a1;
            for (int a2 = -H; a2 <= H; ++a2) {
                ll r2 = r1 + C[2] * (ll)a2;
                for (int a3 = -H; a3 <= H; ++a3) {
                    ll r3 = r2 + C[3] * (ll)a3;
                    ll base = base0 + r3;
                    ll lo = cdivp(-B - base, C[4]);
                    ll hi = fdiv(B - base, C[4]);
                    if (lo < 1) lo = 1;
                    if (hi > H) hi = H;
                    for (ll a4 = lo; a4 <= hi; ++a4) {
                        ll t = base + C[4] * a4;
                        ++scanned;
                        if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = a3; a[4] = (int)a4; emit(a, t); }
                    }
                }
            }
        }
    }
}

// Degree d: a_d is the leading coefficient in [1, H]; a_0..a_{d-1} are in [-H, H].
// The leading term is solved exactly from the prefix instead of looped, so the
// enumeration touches every tuple in C(d,H) but wastes no inner iterations.
static void brute_solve() {
    int a[8];
    ll C0 = C[0];
    for (a[0] = -H; a[0] <= H; ++a[0]) {
        if (((a[0] + H) % shardN) != shardK) continue;
        ll b0 = C0 * (ll)a[0];
        if (D == 1) {
            ll lo = cdivp(-B - b0, C[1]), hi = fdiv(B - b0, C[1]);
            if (lo < 1) lo = 1; if (hi > H) hi = H;
            for (ll x = lo; x <= hi; ++x) { ll t = b0 + C[1] * x; ++scanned; if (iabs(t) <= B) { a[1] = (int)x; emit(a, t); } }
            continue;
        }
        for (int a1 = -H; a1 <= H; ++a1) {
            ll b1 = b0 + C[1] * (ll)a1;
            if (D == 2) {
                ll lo = cdivp(-B - b1, C[2]), hi = fdiv(B - b1, C[2]);
                if (lo < 1) lo = 1; if (hi > H) hi = H;
                for (ll x = lo; x <= hi; ++x) { ll t = b1 + C[2] * x; ++scanned; if (iabs(t) <= B) { a[1] = a1; a[2] = (int)x; emit(a, t); } }
                continue;
            }
            for (int a2 = -H; a2 <= H; ++a2) {
                ll b2 = b1 + C[2] * (ll)a2;
                if (D == 3) {
                    ll lo = cdivp(-B - b2, C[3]), hi = fdiv(B - b2, C[3]);
                    if (lo < 1) lo = 1; if (hi > H) hi = H;
                    for (ll x = lo; x <= hi; ++x) { ll t = b2 + C[3] * x; ++scanned; if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = (int)x; emit(a, t); } }
                    continue;
                }
                for (int a3 = -H; a3 <= H; ++a3) {
                    ll b3 = b2 + C[3] * (ll)a3;
                    ll lo = cdivp(-B - b3, C[4]), hi = fdiv(B - b3, C[4]);
                    if (lo < 1) lo = 1; if (hi > H) hi = H;
                    for (ll x = lo; x <= hi; ++x) { ll t = b3 + C[4] * x; ++scanned; if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = a3; a[4] = (int)x; emit(a, t); } }
                }
            }
        }
    }
}

// Degree 4, leading variable looped rather than solved: used only as a slow
// cross-check that the exact solve loses nothing.
static void brute_d4_loop() {
    int a[8];
    for (a[0] = -H; a[0] <= H; ++a[0]) {
        if (((a[0] + H) % shardN) != shardK) continue;
        ll b0 = C[0] * (ll)a[0];
        for (int a1 = -H; a1 <= H; ++a1) {
            ll b1 = b0 + C[1] * (ll)a1;
            for (int a2 = -H; a2 <= H; ++a2) {
                ll b2 = b1 + C[2] * (ll)a2;
                for (int a3 = -H; a3 <= H; ++a3) {
                    ll b3 = b2 + C[3] * (ll)a3;
                    for (int a4 = 1; a4 <= H; ++a4) {
                        ll t = b3 + C[4] * (ll)a4;
                        ++scanned;
                        if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = a3; a[4] = a4; emit(a, t); }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// mitm for d = 4: left block {0,1,2}, right block {3,4}.
// ---------------------------------------------------------------------------
struct REnt { ll val; int a3, a4; };
static bool rless(const REnt &x, const REnt &y) { return x.val < y.val; }

static void mitm_d4() {
    vector<REnt> R;
    R.reserve((size_t)(2 * H + 1) * H);
    for (int a3 = -H; a3 <= H; ++a3)
        for (int a4 = 1; a4 <= H; ++a4)
            R.push_back({C[3] * (ll)a3 + C[4] * (ll)a4, a3, a4});
    sort(R.begin(), R.end(), rless);
    vector<ll> keys(R.size());
    for (size_t i = 0; i < R.size(); ++i) keys[i] = R[i].val;

    int a[8];
    for (a[0] = -H; a[0] <= H; ++a[0]) {
        if (((a[0] + H) % shardN) != shardK) continue;
        ll b0 = C[0] * (ll)a[0];
        for (int a1 = -H; a1 <= H; ++a1) {
            ll b1 = b0 + C[1] * (ll)a1;
            for (int a2 = -H; a2 <= H; ++a2) {
                ll L = b1 + C[2] * (ll)a2;
                ll want = -L;
                size_t lo = lower_bound(keys.begin(), keys.end(), want - B) - keys.begin();
                size_t hi = upper_bound(keys.begin(), keys.end(), want + B) - keys.begin();
                for (size_t i = lo; i < hi; ++i) {
                    ll t = L + R[i].val;
                    ++scanned;
                    if (iabs(t) <= B) { a[1] = a1; a[2] = a2; a[3] = R[i].a3; a[4] = R[i].a4; emit(a, t); }
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    const char *mode = "brute";
    const char *outpath = nullptr;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--mode") && i + 1 < argc) mode = argv[++i];
        else if (!strcmp(argv[i], "--shard") && i + 2 < argc) { shardK = atoi(argv[++i]); shardN = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--d") && i + 1 < argc) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--H") && i + 1 < argc) H = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--B") && i + 1 < argc) B = atoll(argv[++i]);
        else if (!strcmp(argv[i], "--c") && i + 1 < argc) {
            char *p = argv[++i];
            int k = 0;
            for (char *tok = strtok(p, ","); tok && k <= D; tok = strtok(nullptr, ",")) C[k++] = atoll(tok);
        } else outpath = argv[i];
    }
    if (shardN < 1) shardN = 1;

    auto t0 = chrono::steady_clock::now();
    if (!strcmp(mode, "brute_loop")) { brute_d4_loop(); }
    else if (!strcmp(mode, "mitm")) {
        if (D != 4) { fprintf(stderr, "mitm implemented for d=4 only\n"); return 2; }
        mitm_d4();
    } else {
        if (D < 1 || D > 4) { fprintf(stderr, "unsupported degree %d\n", D); return 2; }
        brute_solve();
    }
    double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

    FILE *f = outpath ? fopen(outpath, "w") : stdout;
    if (!f) { fprintf(stderr, "cannot open %s\n", outpath); return 2; }
    fprintf(f, "# degree %d height %d mode %s shard %d/%d B %lld hits %zu scanned %lld seconds %.3f\n",
            D, H, mode, shardK, shardN, B, hits.size(), scanned, secs);
    for (const Rec &r : hits) {
        for (int k = 0; k <= D; ++k) fprintf(f, "%d ", r.a[k]);
        fprintf(f, "%lld\n", r.t);
    }
    if (outpath) fclose(f);
    fprintf(stderr, "degree=%d H=%d mode=%s shard=%d/%d B=%lld hits=%zu scanned=%lld t=%.3fs\n",
            D, H, mode, shardK, shardN, B, hits.size(), scanned, secs);
    return 0;
}
