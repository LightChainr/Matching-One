// Paired rank births in one site-permutation filtration. C++17, no libraries.
// Lifted union-find convention adapted from src/threshold_rank_axis_mc.cpp.
// Global image rank is the span of ALL cycle windings, not Boolean wraps and
// not the maximum of component ranks. Coordinates are integer torus periods.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct Neighbor { int vertex, dx, dy; };

class BirthEngine {
    int L, N, degree;
    std::vector<std::array<Neighbor,6>> neighbor;
    std::vector<int> parent, size;
    std::vector<std::int64_t> dx, dy;
    std::vector<unsigned char> active;
    int rank = 0;
    std::int64_t basis_x = 0, basis_y = 0;
    struct Root { int vertex; std::int64_t x, y; };
    Root find(int v) {
        if (parent[v] == v) return {v, 0, 0};
        Root r = find(parent[v]);
        dx[v] += r.x; dy[v] += r.y; parent[v] = r.vertex;
        return {r.vertex, dx[v], dy[v]};
    }
    void edge(int u, Neighbor e) {
        Root a = find(u), b = find(e.vertex);
        std::int64_t x = a.x + e.dx - b.x;
        std::int64_t y = a.y + e.dy - b.y;
        if (a.vertex == b.vertex) {
            if (x % L || y % L) throw std::logic_error("nonperiodic cycle");
            x /= L; y /= L;
            if (!(x || y) || rank == 2) return;
            if (rank == 0) { basis_x=x; basis_y=y; rank=1; }
            else if (basis_x*y != basis_y*x) rank=2;
            return;
        }
        if (size[a.vertex] < size[b.vertex]) {
            std::swap(a,b); x=-x; y=-y;
        }
        parent[b.vertex]=a.vertex; dx[b.vertex]=x; dy[b.vertex]=y;
        size[a.vertex]+=size[b.vertex];
    }
public:
    BirthEngine(int length, bool triangular): L(length), N(L*L),
        degree(triangular?6:4), neighbor(N), parent(N), size(N), dx(N),
        dy(N), active(N) {
        if (L < 3 || L > 1024) throw std::invalid_argument("3 <= L <= 1024");
        const int sx[6]={1,-1,0,0,1,-1}, sy[6]={0,0,1,-1,1,-1};
        for (int y=0;y<L;++y) for (int x=0;x<L;++x)
            for (int j=0;j<degree;++j)
                neighbor[x+L*y][j]={((x+sx[j]+L)%L)+L*((y+sy[j]+L)%L),sx[j],sy[j]};
    }
    std::pair<int,int> births(const std::vector<int>& order) {
        std::fill(active.begin(),active.end(),0);
        rank=0; int j1=0;
        for (int k=0;k<N;++k) {
            int v=order[k]; active[v]=1; parent[v]=v; size[v]=1; dx[v]=dy[v]=0;
            for (int j=0;j<degree;++j) {
                Neighbor e=neighbor[v][j];
                if (active[e.vertex]) edge(v,e);
            }
            // Inspect only after every incident edge of the new site is added:
            // a genuine 0->2 insertion has J1=J2, regardless of edge order.
            if (rank >= 1 && j1 == 0) j1=k+1;
            if (rank == 2) return {j1,k+1};
        }
        throw std::logic_error("full torus did not reach rank two");
    }
};

// Unbiased bounded integer, unlike rng()%bound. Explicit Fisher-Yates makes
// reproduction independent of the implementation of std::shuffle.
std::uint64_t bounded(std::mt19937_64& rng, std::uint64_t bound) {
    const std::uint64_t threshold=(-bound)%bound;
    while (true) { auto r=rng(); if (r>=threshold) return r%bound; }
}

int main(int argc,char** argv) {
    try {
        if (argc != 7) {
            std::cerr << "engine square|triangular L samples seed output.json sample|exact\n";
            return 2;
        }
        std::string lattice=argv[1], mode=argv[6];
        if (lattice!="square" && lattice!="triangular") throw std::invalid_argument("lattice");
        int L=std::stoi(argv[2]), N=L*L;
        std::uint64_t samples=std::stoull(argv[3]), seed=std::stoull(argv[4]);
        if (mode!="sample" && mode!="exact") throw std::invalid_argument("mode");
        if (mode=="exact" && L!=3) throw std::invalid_argument("exact limited to L=3");
        if (samples==0 && mode=="sample") throw std::invalid_argument("samples must be positive");
        BirthEngine engine(L,lattice=="triangular");
        std::vector<int> order(N); std::iota(order.begin(),order.end(),0);
        std::mt19937_64 rng(seed);
        std::map<std::pair<int,int>,std::uint64_t> histogram;
        std::uint64_t count=0, atom=0, sum_d=0, sum_d2=0;
        auto start=std::chrono::steady_clock::now();
        do {
            if (mode=="sample") {
                std::iota(order.begin(),order.end(),0);
                for (int i=N-1;i>0;--i) std::swap(order[i],order[bounded(rng,i+1)]);
            }
            auto birth=engine.births(order); ++histogram[birth];
            std::uint64_t d=birth.second-birth.first;
            ++count; atom+=(d==0); sum_d+=d; sum_d2+=d*d;
            if (mode=="sample" && count==samples) break;
        } while (mode=="sample" || std::next_permutation(order.begin(),order.end()));
        if (mode=="exact") {
            const std::uint64_t numerator=lattice=="square"?3:2;
            const std::uint64_t d2_num=lattice=="square"?43:81;
            const std::uint64_t d2_den=lattice=="square"?14:28;
            if (count!=362880 || atom*35!=count*numerator || sum_d*2!=3*count ||
                sum_d2*d2_den!=d2_num*count)
                throw std::logic_error("L=3 exact paired-birth control failed");
        }
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream out(argv[5]); if (!out) throw std::runtime_error("cannot open output");
        out.precision(17);
        out << "{\"schema\":\"matching-one.joint-birth-batch.v1\",\"lattice\":\"" << lattice
            << "\",\"L\":" << L << ",\"N\":" << N << ",\"mode\":\"" << mode
            << "\",\"seed\":\"" << seed << "\",\"samples\":" << count
            << ",\"runtime_seconds\":" << seconds << ",\"direct_atom_count\":" << atom
            << ",\"sum_D\":" << sum_d << ",\"sum_D2\":" << sum_d2
            << ",\"histogram_columns\":[\"J1\",\"J2\",\"count\"],\"histogram\":[";
        bool comma=false;
        for (const auto& h:histogram) {
            if (comma) out << ','; comma=true;
            out << '[' << h.first.first << ',' << h.first.second << ',' << h.second << ']';
        }
        out << "]}\n";
        if (!out) throw std::runtime_error("output write failed");
        std::cout << lattice << " L=" << L << " " << mode << " samples=" << count
                  << " seconds=" << seconds << " atom=" << atom << " mean_D="
                  << static_cast<double>(sum_d)/count << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
