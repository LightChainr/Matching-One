// Completion geometry on independent uniform-permutation filtrations.
// Lifted union-find convention retained from ../birth-gap-20260929/engine.cpp.
// The old engine and archives are unchanged. C++17, standard library only.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct Neighbor { int vertex, dx, dy; };
struct Observation { int j1=0, j2=0, dx_b=0, dy_b=0, nu_b=-1, tau=0, nu_tau=-1; };

class Geometry {
    int L, N, degree;
    std::vector<std::array<Neighbor,6>> neighbor;
    std::vector<int> parent, size;
    std::vector<std::int64_t> dx, dy;
    std::vector<unsigned char> active;
    int rank=0;
    std::int64_t bx=0, by=0;
    struct Root { int vertex; std::int64_t x,y; };
    Root find(int v) {
        if (parent[v]==v) return {v,0,0};
        Root r=find(parent[v]);
        dx[v]+=r.x; dy[v]+=r.y; parent[v]=r.vertex;
        return {r.vertex,dx[v],dy[v]};
    }
    void edge(int v, Neighbor e) {
        Root a=find(v), b=find(e.vertex);
        std::int64_t x=a.x+e.dx-b.x, y=a.y+e.dy-b.y;
        if (a.vertex==b.vertex) {
            if (x%L || y%L) throw std::logic_error("nonperiodic cycle");
            x/=L; y/=L;
            if (!(x||y) || rank==2) return;
            if (!rank) {bx=x;by=y;rank=1;}
            else if (bx*y!=by*x) rank=2;
        } else {
            if (size[a.vertex]<size[b.vertex]) {std::swap(a,b);x=-x;y=-y;}
            parent[b.vertex]=a.vertex;dx[b.vertex]=x;dy[b.vertex]=y;
            size[a.vertex]+=size[b.vertex];
        }
    }
public:
    Geometry(int length,bool triangular):L(length),N(L*L),degree(triangular?6:4),
        neighbor(N),parent(N),size(N),dx(N),dy(N),active(N) {
        if (L<3 || L>1024) throw std::invalid_argument("L outside [3,1024]");
        const int sx[]={1,-1,0,0,1,-1},sy[]={0,0,1,-1,1,-1};
        for(int y=0;y<L;++y)for(int x=0;x<L;++x)for(int j=0;j<degree;++j)
            neighbor[x+L*y][j]={((x+sx[j]+L)%L)+L*((y+sy[j]+L)%L),sx[j],sy[j]};
    }
    void reset(){std::fill(active.begin(),active.end(),0);rank=0;bx=by=0;}
    int insert(int v) {
        if(active[v])throw std::logic_error("double insertion");
        active[v]=1;parent[v]=v;size[v]=1;dx[v]=dy[v]=0;
        for(int j=0;j<degree;++j)if(active[neighbor[v][j].vertex])edge(v,neighbor[v][j]);
        return rank;
    }
    std::pair<int,int> direction()const{
        if(rank!=1)return {0,0};
        auto g=std::gcd(bx,by);auto x=bx/g,y=by/g;
        if(x<0 || (x==0 && y<0)){x=-x;y=-y;}
        return {static_cast<int>(x),static_cast<int>(y)};
    }
    int completion_count() {
        if(rank!=1)return -1;
        int result=0;
        // A new cycle through v arises precisely from two neighbours in
        // the same old component. Cycles between different old components
        // cannot arise from a single star. Existing windings are all along b.
        for(int v=0;v<N;++v)if(!active[v]){
            std::array<Root,6> roots; int used=0; bool completes=false;
            for(int j=0;j<degree && !completes;++j){
                Neighbor e=neighbor[v][j];if(!active[e.vertex])continue;
                Root r=find(e.vertex);r.x=e.dx-r.x;r.y=e.dy-r.y;
                for(int h=0;h<used;++h)if(roots[h].vertex==r.vertex){
                    auto x=r.x-roots[h].x,y=r.y-roots[h].y;
                    if(x%L || y%L)throw std::logic_error("hypothetical nonperiodic cycle");
                    if(bx*y!=by*x){completes=true;break;}
                }
                roots[used++]=r;
            }
            result+=completes;
        }
        return result;
    }
    Observation observe(const std::vector<int>&order,int b,int tau){
        reset();Observation o;o.tau=tau;
        for(int k=1;k<=N;++k){
            insert(order[k-1]);
            if(rank>=1 && !o.j1)o.j1=k;
            if(rank==2){o.j2=k;return o;}
            if(k==b){
                auto d=direction();o.dx_b=d.first;o.dy_b=d.second;
                o.nu_b=completion_count();
            }
            if(k==tau)o.nu_tau=(tau==b?o.nu_b:completion_count());
        }
        throw std::logic_error("full torus not rank two");
    }
    int control(std::ostream&out){
        if(L!=3)throw std::logic_error("control L3 only");
        out<<"mask,rank,dx,dy,nu2\n";int controls=0;
        for(int mask=0;mask<(1<<N);++mask){
            reset();for(int v=0;v<N;++v)if(mask>>v&1)insert(v);
            auto d=direction();int nu=completion_count();
            if(rank==1){
                int brute=0;
                for(int v=0;v<N;++v)if(!(mask>>v&1)){
                    Geometry copy=*this;brute+=(copy.insert(v)==2);++controls;
                }
                if(nu!=brute)throw std::logic_error("completion probe mismatch");
            }
            out<<mask<<','<<rank<<','<<d.first<<','<<d.second<<','<<nu<<'\n';
        }
        return controls;
    }
};

std::uint64_t bounded(std::mt19937_64&rng,std::uint64_t bound){
    std::uint64_t threshold=(-bound)%bound;
    while(true){auto r=rng();if(r>=threshold)return r%bound;}
}

int main(int argc,char**argv){
    try{
        if(argc<3)throw std::invalid_argument("mode lattice ...");
        std::string mode=argv[1],lattice=argv[2];
        if(lattice!="square" && lattice!="triangular")throw std::invalid_argument("lattice");
        if(mode=="control"){
            if(argc!=4)throw std::invalid_argument("control square|triangular output.csv");
            Geometry g(3,lattice=="triangular");std::ofstream out(argv[3]);
            if(!out)throw std::runtime_error("cannot open control output");
            int checks=g.control(out);out.close();
            if(!out)throw std::runtime_error("control write failed");
            std::cout<<"512 subsets; "<<checks<<" single-site completion probes agree\n";return 0;
        }
        if(mode!="sample" || argc!=10)
            throw std::invalid_argument("sample square|triangular L samples permutation_seed time_seed b c output.csv");
        int L=std::stoi(argv[3]),N=L*L,b=std::stoi(argv[7]),c=std::stoi(argv[8]);
        auto samples=std::stoull(argv[4]),seed=std::stoull(argv[5]),time_seed=std::stoull(argv[6]);
        if(!samples || b<1 || b>=c || c>N)throw std::invalid_argument("sample/time bounds");
        Geometry g(L,lattice=="triangular");std::vector<int>order(N);
        std::mt19937_64 rng(seed),time_rng(time_seed);
        std::ofstream out(argv[9]);if(!out)throw std::runtime_error("cannot open output");
        out<<"J1,J2,dx_b,dy_b,nu_b,tau,nu_tau\n";
        auto start=std::chrono::steady_clock::now();
        for(std::uint64_t i=0;i<samples;++i){
            std::iota(order.begin(),order.end(),0);
            for(int v=N-1;v>0;--v)std::swap(order[v],order[bounded(rng,v+1)]);
            int tau=b+bounded(time_rng,c-b);
            auto o=g.observe(order,b,tau);
            out<<o.j1<<','<<o.j2<<','<<o.dx_b<<','<<o.dy_b<<','<<o.nu_b<<','<<o.tau<<','<<o.nu_tau<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("output write failed");
        double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"samples="<<samples<<" seconds="<<sec<<'\n';return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
