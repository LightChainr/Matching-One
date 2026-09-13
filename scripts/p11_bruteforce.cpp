// Path B (fast): brute-force enumeration of all 2^(nm) occupancy patterns on
// the free-boundary n x m grid; union-find spanning check, counts by k.
// Independent of any transfer matrix. Usage: ./p11_bruteforce n [m]
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

int main(int argc, char** argv) {
    int n = atoi(argv[1]);
    int m = argc > 2 ? atoi(argv[2]) : n;
    int nm = n * m;
    vector<unsigned long long> A(nm + 1, 0);
    vector<int> par(nm + 1), stk;
    for (unsigned long long mask = 0; mask < (1ull << nm); ++mask) {
        int k = __builtin_popcountll(mask & ((1ull << nm) - 1));
        if (nm > 63) { fprintf(stderr, "nm too large\n"); return 1; }
        // union-find with path halving over occupied sites only
        for (int i = 0; i < nm; ++i) par[i] = i;
        auto find = [&](int x) {
            while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
            return x;
        };
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j) {
                int v = i * n + j;
                if (!((mask >> v) & 1ull)) continue;
                if (j + 1 < n && ((mask >> (v + 1)) & 1ull)) {
                    int a = find(v), b = find(v + 1);
                    if (a != b) par[a] = b;
                }
                if (i + 1 < m && ((mask >> (v + n)) & 1ull)) {
                    int a = find(v), b = find(v + n);
                    if (a != b) par[a] = b;
                }
            }
        bool spans = false;
        for (int j = 0; j < n && !spans; ++j) {
            if (!((mask >> j) & 1ull)) continue;            // top row
            for (int j2 = (m - 1) * n; j2 < m * n; ++j2) {  // bottom row
                if (!((mask >> j2) & 1ull)) continue;
                if (find(j) == find(j2)) { spans = true; break; }
            }
        }
        if (spans) ++A[k];
    }
    printf("{\"path\":\"brute_force_cpp\",\"n\":%d,\"m\":%d,\"A\":[", n, m);
    for (int k = 0; k <= nm; ++k)
        printf("%s%llu", k ? "," : "", A[k]);
    printf("]}\n");
    return 0;
}
