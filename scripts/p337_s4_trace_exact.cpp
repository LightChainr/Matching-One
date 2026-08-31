// Fixed S4 seam insertion: missing per-component first-deck mod2/mod3 data.
// The geometry, white Alexander components, site order, edge/face counting
// and binary traversal are copied from p337_closed_source_finite_exact.cpp.
// Only black rollback carries an additional Z/6 potential and cycle flags.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

class RollbackComponents {
    std::array<int, 25> parent{}, size{};
    std::vector<std::pair<int, int>> changed;
public:
    int components = 0;
    void activate(int v) { parent[v] = v; size[v] = 1; ++components; }
    int root(int v) const { while (parent[v] != v) v = parent[v]; return v; }
    std::size_t mark() const { return changed.size(); }
    void join(int a, int b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a, b);
        changed.emplace_back(b, size[a]); parent[b] = a; size[a] += size[b];
        --components;
    }
    void undo(std::size_t mark) {
        while (changed.size() > mark) {
            const auto [b, old_size] = changed.back(); changed.pop_back();
            size[parent[b]] = old_size; parent[b] = b; ++components;
        }
        --components;
    }
};

class RollbackSeamComponents {
    std::array<int, 25> parent{}, size{}, potential{}, bad2{}, bad3{};
    struct Change { int a, b, size_a, bad2_a, bad3_a, total2, total3; };
    std::vector<Change> changed;
    static int mod6(int x) { x %= 6; return x < 0 ? x+6 : x; }
    std::pair<int,int> root_shift(int v) const {
        int s=0;
        while (parent[v]!=v) { s+=potential[v]; v=parent[v]; }
        return {v,mod6(s)};
    }
public:
    int components=0, total_bad2=0, num_bad3=0;
    void activate(int v) {
        parent[v]=v; size[v]=1; potential[v]=bad2[v]=bad3[v]=0;
        ++components;
    }
    std::size_t mark() const { return changed.size(); }
    void join(int v, int u, int gain) {
        auto [a, sa]=root_shift(v);
        auto [b, sb]=root_shift(u);
        // phi(u)-phi(v)=gain, so phi(root_u)-phi(root_v)=gain+sa-sb.
        int delta=mod6(gain+sa-sb);
        if (a==b) {
            const int next2=bad2[a] || delta%2, next3=bad3[a] || delta%3;
            if (next2==bad2[a] && next3==bad3[a]) return;
            changed.push_back({a,-1,size[a],bad2[a],bad3[a],total_bad2,num_bad3});
            total_bad2+=next2-bad2[a]; num_bad3+=next3-bad3[a];
            bad2[a]=next2; bad3[a]=next3;
            return;
        }
        if (size[a]<size[b]) { std::swap(a,b); delta=mod6(-delta); }
        changed.push_back({a,b,size[a],bad2[a],bad3[a],total_bad2,num_bad3});
        total_bad2-=bad2[a]+bad2[b]; num_bad3-=bad3[a]+bad3[b];
        bad2[a]=bad2[a] || bad2[b]; bad3[a]=bad3[a] || bad3[b];
        total_bad2+=bad2[a]; num_bad3+=bad3[a];
        parent[b]=a; potential[b]=delta; size[a]+=size[b]; --components;
    }
    void undo(std::size_t mark) {
        while (changed.size()>mark) {
            const auto c=changed.back(); changed.pop_back();
            total_bad2=c.total2; num_bad3=c.total3;
            bad2[c.a]=c.bad2_a; bad3[c.a]=c.bad3_a;
            if (c.b>=0) {
                size[c.a]=c.size_a; parent[c.b]=c.b; potential[c.b]=0;
                ++components;
            }
        }
        --components; // activated vertex has no remaining joins/cycles
    }
};

class Enumerator {
    static constexpr int n=25;
    // Per k, key=(g,q,any bad2 component,number bad3 components).
    std::array<std::map<std::array<int,4>,std::uint64_t>,n+1> histogram;
    std::array<int,n> colour{};
    struct BlackEdge { int vertex, gain; };
    std::array<std::vector<BlackEdge>,n> black_previous;
    std::array<std::vector<int>,n> white_previous;
    std::array<std::vector<std::array<int,4>>,n> closing_faces;
    RollbackSeamComponents black;
    RollbackComponents white;
    int a,b;
    static int mod(int x) { x%=n; return x<0 ? x+n : x; }
    int key(int x,int y) const { return n*mod(a*x+b*y)+mod(-b*x+a*y); }

    void visit(int v,int k,int edges,int faces) {
        if (v==n) {
            const int q=black.components-white.components-(k-edges+faces);
            if (q < -1 || q > 1) throw std::logic_error("invalid digital rank value");
            const int c=black.components+white.components;
            const int bv=2*n-4*k+edges, sstar=c+faces+bv;
            const int g=2*n+1-k-sstar;
            ++histogram[k][{g,q,int(black.total_bad2>0),black.num_bad3}];
            return;
        }
        colour[v]=0;
        auto mark=white.mark(); white.activate(v);
        for (int u:white_previous[v]) if (colour[u]==0) white.join(v,u);
        visit(v+1,k,edges,faces);
        white.undo(mark);

        colour[v]=1;
        mark=black.mark(); black.activate(v);
        int extra_edges=0,extra_faces=0;
        for (const auto& edge:black_previous[v]) if (colour[edge.vertex]==1) {
            black.join(v,edge.vertex,edge.gain); ++extra_edges;
        }
        for (const auto& face:closing_faces[v]) {
            bool full=true;
            for (int u:face) if (colour[u]!=1) full=false;
            extra_faces+=full;
        }
        visit(v+1,k+1,edges+extra_edges,faces+extra_faces);
        black.undo(mark); colour[v]=-1;
    }
public:
    Enumerator(int aa,int bb):a(aa),b(bb) {
        if (!((a==5 && b==0)||(a==4 && b==3)))
            throw std::invalid_argument("fixed scientific pair is (5,0),(4,3)");
        colour.fill(-1);
        std::array<int,n*n> index; index.fill(-1);
        std::vector<std::pair<int,int>> reps;
        reps.emplace_back(0,0); index[key(0,0)]=0;
        for (std::size_t i=0;i<reps.size();++i) {
            auto [x,y]=reps[i];
            for (const auto& step:std::array<std::pair<int,int>,2>{{{1,0},{0,1}}}) {
                int xx=x+step.first,yy=y+step.second,h=key(xx,yy);
                if (index[h]<0) { index[h]=static_cast<int>(reps.size()); reps.emplace_back(xx,yy); }
            }
        }
        if (reps.size()!=n) throw std::logic_error("wrong quotient area");
        for (int v=0;v<n;++v) {
            auto [x,y]=reps[v];
            for (int dx=-1;dx<=1;++dx) for (int dy=-1;dy<=1;++dy) {
                if (dx==0 && dy==0) continue;
                int u=index[key(x+dx,y+dy)];
                if (u<v) {
                    white_previous[v].push_back(u);
                    if (dx==0 || dy==0) {
                        const int deck_numerator=a*(x+dx-reps[u].first)+b*(y+dy-reps[u].second);
                        if (deck_numerator%n) throw std::logic_error("nonintegral first deck gain");
                        black_previous[v].push_back({u,deck_numerator/n});
                    }
                }
            }
            std::array<int,4> face{v,index[key(x+1,y)],index[key(x,y+1)],index[key(x+1,y+1)]};
            closing_faces[*std::max_element(face.begin(),face.end())].push_back(face);
        }
    }
    void run(const char* path) {
        std::ifstream exists(path);
        if (exists.good()) throw std::runtime_error("output already exists");
        visit(0,0,0,0);
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot create output");
        out << "k,g,q,bad2,n_bad3,count\n";
        for (int k=0;k<=n;++k) for (const auto& [key,count]:histogram[k]) {
            out << k;
            for (int value:key) out << ',' << value;
            out << ',' << count << '\n';
        }
        out.close(); if (!out) throw std::runtime_error("output write failed");
    }
};

int main(int argc,char** argv) {
    try {
        if (argc!=4) throw std::invalid_argument("usage: enumerate a b output.csv");
        const auto start=std::chrono::steady_clock::now();
        Enumerator(std::stoi(argv[1]),std::stoi(argv[2])).run(argv[3]);
        const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "{\"configurations\":33554432,\"elapsed_seconds\":" << sec << "}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
