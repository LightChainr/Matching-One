// Prepared Q1 derivative of the forced regular completion Kreg=K2bar+K0bar.
// DO NOT RUN before the final character/topology gate, contract and root GO.
// Keep the old N25 quotient, NN/matching DSUs, Alexander q and binary traversal.
// At vacant origin, occupied-neighbor edge-nodes use black NN components;
// each vacant-neighbor incident edge-node is its own distinct singleton.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
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

struct Totals {
    std::uint64_t count=0;
    std::int64_t q=0, e=0, activation4=0, q_activation4=0, e_activation4=0;
};

class Enumerator {
    static constexpr int n=25;
    std::array<Totals,n+1> histogram{};
    std::array<int,n> colour{};
    std::array<int,4> ports{}; // N,S,E,W around the representative origin0
    std::array<std::vector<int>,n> black_previous,white_previous;
    std::array<std::vector<std::array<int,4>>,n> closing_faces;
    RollbackComponents black,white;
    int a,b;
    static int mod(int x) { x%=n; return x<0 ? x+n : x; }
    int key(int x,int y) const { return n*mod(a*x+b*y)+mod(-b*x+a*y); }

    int four_times_activation() const {
        if (colour[0]!=0) return 0;
        std::array<int,4> outside{};
        for (int i=0;i<4;++i) {
            // With origin vacant, a vacant neighbor leaves this incident
            // edge-node isolated. Never query its inactive black DSU root.
            outside[i]=colour[ports[i]]==1 ? black.root(ports[i]) : n+i;
        }
        int equal_pairs=0;
        for (int i=0;i<4;++i) for (int j=i+1;j<4;++j)
            equal_pairs+=outside[i]==outside[j];
        if (equal_pairs==0) return 4; // four distinct outside components
        if (equal_pairs==1) { // exactly 2+1+1
            const bool opposite=outside[0]==outside[1] || outside[2]==outside[3];
            return opposite ? 4 : 2;
        }
        if (equal_pairs==2) { // exactly 2+2; two opposite or adjacent pairs
            const bool opposite=outside[0]==outside[1] && outside[2]==outside[3];
            return opposite ? -2 : -1;
        }
        return 0; // 3+1 (three equal pairs) or all4 (six equal pairs)
    }

    void visit(int v,int k,int edges,int faces) {
        if (v==n) {
            const int q=black.components-white.components-(k-edges+faces);
            if (q < -1 || q > 1) throw std::logic_error("invalid digital rank value");
            const int e=q*q, a4=four_times_activation();
            auto& row=histogram[k];
            ++row.count; row.q+=q; row.e+=e; row.activation4+=a4;
            row.q_activation4+=q*a4; row.e_activation4+=e*a4;
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
        for (int u:black_previous[v]) if (colour[u]==1) {
            black.join(v,u); ++extra_edges;
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
        ports={index[key(0,1)],index[key(0,-1)],index[key(1,0)],index[key(-1,0)]};
        for (int v=0;v<n;++v) {
            auto [x,y]=reps[v];
            for (int dx=-1;dx<=1;++dx) for (int dy=-1;dy<=1;++dy) {
                if (dx==0 && dy==0) continue;
                int u=index[key(x+dx,y+dy)];
                if (u<v) {
                    white_previous[v].push_back(u);
                    if (dx==0 || dy==0) black_previous[v].push_back(u);
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
        out << "k,count,sum_q,sum_e,sum_a4,sum_qa4,sum_ea4\n";
        for (int k=0;k<=n;++k) {
            const auto& row=histogram[k];
            out << k << ',' << row.count << ',' << row.q << ',' << row.e << ','
                << row.activation4 << ',' << row.q_activation4 << ',' << row.e_activation4 << '\n';
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
