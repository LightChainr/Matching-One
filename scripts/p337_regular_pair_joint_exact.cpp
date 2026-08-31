// Prepared exact joint source for d_logQ d_epsilon^2 U of canonical Kreg.
// DO NOT RUN before the root's frozen joint-source contract and explicit GO.
// Preserve the original N25 quotient, black/white rollback, q/E and traversal.
// Only origin-to-other-vacant-site signed Bell8 kernel cross moments are new.
// Adjacent ports sharing one physical vacant edge use the same edge-node ID.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

using Kernel = std::vector<std::int16_t>;

Kernel read_kernel(const char* path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read exact kernel TSV");
    Kernel kernel(1U<<24,0); // Sparse TSV omissions are exact zeros.
    std::string line,field;
    int key_column=-1,g_column=-1;
    while (std::getline(in,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty() || line.front()=='#') continue;
        std::stringstream stream(line);
        std::vector<std::string> fields;
        while (std::getline(stream,field,'\t')) fields.push_back(field);
        if (key_column<0) {
            for (std::size_t i=0;i<fields.size();++i) {
                if (fields[i]=="key" || fields[i]=="packed_key") key_column=static_cast<int>(i);
                if (fields[i]=="g16") g_column=static_cast<int>(i);
            }
            if (key_column<0 || g_column<0) throw std::runtime_error("kernel needs key,g16 TSV columns");
            continue;
        }
        if (fields.size()<=static_cast<std::size_t>(std::max(key_column,g_column)))
            throw std::runtime_error("incomplete kernel TSV row");
        const auto packed=std::stoull(fields[key_column]);
        const auto g16=std::stoll(fields[g_column]);
        if (packed>=kernel.size() || g16<std::numeric_limits<std::int16_t>::min() ||
            g16>std::numeric_limits<std::int16_t>::max())
            throw std::runtime_error("kernel entry out of fixed packed-key/int16 range");
        kernel[packed]=static_cast<std::int16_t>(g16);
    }
    if (key_column<0) throw std::runtime_error("empty kernel TSV");
    return kernel;
}

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
    std::int64_t q=0, e=0;
    // Ordered as total, adjacent, nonadjacent; all source sums stay int64.
    std::array<std::int64_t,3> b16{},qb16{},eb16{};
};

class Enumerator {
    static constexpr int n=25;
    std::array<Totals,n+1> histogram{};
    std::array<int,n> colour{};
    std::array<std::array<int,4>,n> ports{},port_edges{}; // N,E,S,W
    std::array<bool,n> adjacent_to_origin{};
    const Kernel& kernel;
    std::array<std::vector<int>,n> black_previous,white_previous;
    std::array<std::vector<std::array<int,4>>,n> closing_faces;
    RollbackComponents black,white;
    int a,b;
    static int mod(int x) { x%=n; return x<0 ? x+n : x; }
    int key(int x,int y) const { return n*mod(a*x+b*y)+mod(-b*x+a*y); }

    std::array<std::int64_t,2> sum_origin_pair_kernels() const {
        std::array<std::int64_t,2> sums{}; // adjacent, nonadjacent
        if (colour[0]!=0) return sums;
        std::array<int,n> roots{};
        for (int v=0;v<n;++v) if (colour[v]==1) roots[v]=black.root(v);
        std::array<int,4> origin_outside{};
        for (int i=0;i<4;++i) {
            const int u=ports[0][i];
            origin_outside[i]=colour[u]==1 ? roots[u] : n+port_edges[0][i];
        }
        for (int y=1;y<n;++y) if (colour[y]==0) {
            std::array<int,8> outside{},labels{};
            std::uint32_t packed=0;
            int next_label=0;
            for (int i=0;i<8;++i) {
                if (i<4) outside[i]=origin_outside[i];
                else {
                    const int d=i-4,u=ports[y][d];
                    outside[i]=colour[u]==1 ? roots[u] : n+port_edges[y][d];
                }
                int previous=0;
                while (previous<i && outside[previous]!=outside[i]) ++previous;
                labels[i]=previous<i ? labels[previous] : next_label++;
                packed|=static_cast<std::uint32_t>(labels[i])<<(3*i);
            }
            sums[adjacent_to_origin[y] ? 0 : 1]+=kernel[packed];
        }
        return sums;
    }

    void visit(int v,int k,int edges,int faces) {
        if (v==n) {
            const int q=black.components-white.components-(k-edges+faces);
            if (q < -1 || q > 1) throw std::logic_error("invalid digital rank value");
            const int e=q*q;
            const auto pair_sums=sum_origin_pair_kernels();
            const std::array<std::int64_t,3> b16{pair_sums[0]+pair_sums[1],pair_sums[0],pair_sums[1]};
            auto& row=histogram[k];
            ++row.count; row.q+=q; row.e+=e;
            for (int group=0;group<3;++group) {
                row.b16[group]+=b16[group];
                row.qb16[group]+=q*b16[group];
                row.eb16[group]+=e*b16[group];
            }
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
    Enumerator(int aa,int bb,const Kernel& lookup):kernel(lookup),a(aa),b(bb) {
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
        // Unique undirected physical edges are the N/E edges out of each site.
        // Opposite S/W ports refer to that same edge at its other endpoint.
        for (int v=0;v<n;++v) {
            auto [x,y]=reps[v];
            ports[v]={index[key(x,y+1)],index[key(x+1,y)],index[key(x,y-1)],index[key(x-1,y)]};
            port_edges[v]={2*v,2*v+1,2*ports[v][2],2*ports[v][3]+1};
        }
        for (int u:ports[0]) adjacent_to_origin[u]=true;
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
        out << "k,count,sum_q,sum_e,sum_b16,sum_qb16,sum_eb16,"
               "sum_b16_adj,sum_qb16_adj,sum_eb16_adj,"
               "sum_b16_far,sum_qb16_far,sum_eb16_far\n";
        for (int k=0;k<=n;++k) {
            const auto& row=histogram[k];
            out << k << ',' << row.count << ',' << row.q << ',' << row.e;
            for (int group=0;group<3;++group)
                out << ',' << row.b16[group] << ',' << row.qb16[group] << ',' << row.eb16[group];
            out << '\n';
        }
        out.close(); if (!out) throw std::runtime_error("output write failed");
    }
};

int main(int argc,char** argv) {
    try {
        if (argc!=5) throw std::invalid_argument("usage: enumerate a b kernel.tsv output.csv");
        const auto start=std::chrono::steady_clock::now();
        const auto kernel=read_kernel(argv[3]);
        Enumerator(std::stoi(argv[1]),std::stoi(argv[2]),kernel).run(argv[4]);
        const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "{\"configurations\":33554432,\"elapsed_seconds\":" << sec << "}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

