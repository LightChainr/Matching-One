// Conditional mean of the existing safe-insertion target, with no safe probes.
// Source permutations and lifted union-find are reused without modification.
#define main original_completion_engine_main
#include "../completion-hazard-production-20260929/engine.cpp"
#undef main
#include <unordered_map>
#include <unordered_set>

using EdgeKey = std::uint64_t;
EdgeKey pair_key(int a,int b){if(a>b)std::swap(a,b);return (EdgeKey(a)<<32)|unsigned(b);}

struct Contact {
    std::array<int,6> component{};
    std::array<std::int64_t,6> height{};
    int n=0;
};

struct PairResult {
    int nu=0;
    std::uint64_t edges=0, candidates=0;
};

class ProjectedPairs {
    int L,N,z;
    std::vector<std::array<Neighbor,6>> neighbor;
    std::vector<unsigned char> occupied,safe;
    std::vector<int> component,stack;
    std::vector<std::int64_t> potential;
    std::vector<Contact> contact;
public:
    std::unordered_set<EdgeKey> edges;
    ProjectedPairs(int length,bool triangular):L(length),N(L*L),z(triangular?6:4),
        neighbor(N),occupied(N),safe(N),component(N),potential(N),contact(N){
        const int sx[]={1,-1,0,0,1,-1},sy[]={0,0,1,-1,1,-1};
        for(int y=0;y<L;++y)for(int x=0;x<L;++x)for(int j=0;j<z;++j)
            neighbor[x+L*y][j]={((x+sx[j]+L)%L)+L*((y+sy[j]+L)%L),sx[j],sy[j]};
        stack.reserve(N);
    }
    PairResult count(const std::vector<int>& order,int k,std::pair<int,int> d){
        auto ell=[d](int x,int y){return std::int64_t(d.first)*y-std::int64_t(d.second)*x;};
        std::fill(occupied.begin(),occupied.end(),0);
        std::fill(safe.begin(),safe.end(),0);
        std::fill(component.begin(),component.end(),-1);
        for(int i=0;i<k;++i)occupied[order[i]]=1;
        // Rank one guarantees a consistent transverse potential even when
        // more than one occupied component carries the persistent direction.
        for(int root=0;root<N;++root)if(occupied[root] && component[root]<0){
            component[root]=root;potential[root]=0;stack.clear();stack.push_back(root);
            while(!stack.empty()){
                int u=stack.back();stack.pop_back();
                for(int j=0;j<z;++j){auto e=neighbor[u][j];if(!occupied[e.vertex])continue;
                    auto h=potential[u]+ell(e.dx,e.dy);
                    if(component[e.vertex]<0){component[e.vertex]=root;potential[e.vertex]=h;stack.push_back(e.vertex);}
                    else if(component[e.vertex]!=root || potential[e.vertex]!=h)
                        throw std::logic_error("inconsistent rank-one transverse potential");
                }
            }
        }
        PairResult out;
        for(int v=0;v<N;++v)if(!occupied[v]){
            auto& row=contact[v];row.n=0;bool completes=false;
            for(int j=0;j<z && !completes;++j){
                auto e=neighbor[v][j];if(!occupied[e.vertex])continue;
                int c=component[e.vertex];auto a=potential[e.vertex]-ell(e.dx,e.dy);
                int h=0;while(h<row.n && row.component[h]!=c)++h;
                if(h<row.n){if(row.height[h]!=a)completes=true;}
                else {row.component[row.n]=c;row.height[row.n++]=a;}
            }
            if(completes)++out.nu;
            else safe[v]=1;
        }
        // Each root pair has groups of vertices with equal relative height.
        // Only cross-group pairs are enumerated; every one is a true edge.
        using HeightGroups=std::unordered_map<std::int64_t,std::vector<int>>;
        std::unordered_map<EdgeKey,HeightGroups> buckets;
        for(int v=0;v<N;++v)if(safe[v]){
            const auto& row=contact[v];
            for(int i=0;i<row.n;++i)for(int j=i+1;j<row.n;++j){
                int c=row.component[i],cc=row.component[j];
                auto h=row.height[i]-row.height[j];
                if(c>cc){std::swap(c,cc);h=-h;}
                buckets[pair_key(c,cc)][h].push_back(v);
            }
        }
        edges.clear();
        for(const auto& bucket:buckets){
            const auto& groups=bucket.second;
            for(auto a=groups.begin();a!=groups.end();++a){
                auto b=a;++b;
                for(;b!=groups.end();++b)
                    for(int v:a->second)for(int w:b->second){
                        edges.insert(pair_key(v,w));++out.candidates;
                    }
            }
        }
        // A direct lattice edge provides one more relative-height constraint.
        // One shared old component is enough for a conflict with that edge.
        for(int v=0;v<N;++v)if(safe[v]){
            const auto& a=contact[v];
            for(int j=0;j<z;++j){
                auto e=neighbor[v][j];int w=e.vertex;
                if(w<=v || !safe[w])continue;
                const auto& b=contact[w];bool mismatch=false;
                for(int p=0;p<a.n && !mismatch;++p)for(int q=0;q<b.n;++q)
                    if(a.component[p]==b.component[q] && a.height[p]-b.height[q]+ell(e.dx,e.dy)!=0){mismatch=true;break;}
                if(mismatch)edges.insert(pair_key(v,w));
            }
        }
        out.edges=edges.size();return out;
    }
};

int main(int argc,char**argv){
    try{
        if(argc==2 && std::string(argv[1])=="control"){
            const std::array<int,4> masks{{4403,4405,4371,4402}};
            const std::array<int,4> expected{{2,3,4,6}};
            for(int i=0;i<4;++i){
                Geometry g(4,i>=2);ProjectedPairs counter(4,i>=2);std::vector<int> order;
                for(int v=0;v<16;++v)if(masks[i]>>v&1){order.push_back(v);g.insert(v);}
                int k=order.size();auto r=counter.count(order,k,g.direction());
                std::unordered_set<EdgeKey> brute;
                for(int v=0;v<16;++v)if(!(masks[i]>>v&1)){
                    Geometry one=g;if(one.insert(v)==2)continue;
                    for(int w=v+1;w<16;++w)if(!(masks[i]>>w&1)){
                        Geometry alone=g;if(alone.insert(w)==2)continue;
                        Geometry two=one;if(two.insert(w)==2)brute.insert(pair_key(v,w));
                    }
                }
                if(r.nu!=0 || r.edges!=unsigned(expected[i]) || counter.edges!=brute)
                    throw std::logic_error("physical pair witness mismatch");
                std::cout<<"mask="<<masks[i]<<" nu="<<r.nu<<" e="<<r.edges<<" nonlocal_candidates="<<r.candidates<<'\n';
            }
            return 0;
        }
        if(argc!=6)throw std::invalid_argument("L prefixes existing_permutation_seed b output.csv");
        int L=std::stoi(argv[1]),N=L*L,b=std::stoi(argv[4]);
        auto samples=std::stoull(argv[2]),seed=std::stoull(argv[3]);
        if(L<3 || b<1 || b>N-2 || !samples)throw std::invalid_argument("bounds");
        Geometry g(L,false);ProjectedPairs counter(L,false);std::vector<int>order(N);
        std::mt19937_64 rng(seed);std::ofstream out(argv[5]);
        if(!out)throw std::runtime_error("cannot open output");
        out<<"J1,rank_b,dx_b,dy_b,nu_b,synergy_edges,nonlocal_candidates\n";
        auto start=std::chrono::steady_clock::now();
        for(std::uint64_t i=0;i<samples;++i){
            std::iota(order.begin(),order.end(),0);
            for(int v=N-1;v>0;--v)std::swap(order[v],order[bounded(rng,v+1)]);
            g.reset();int j1=0,rank=0;
            for(int k=1;k<=b;++k){rank=g.insert(order[k-1]);if(rank && !j1)j1=k;if(rank==2)break;}
            if(rank!=1){out<<j1<<','<<rank<<",0,0,-1,-1,-1\n";continue;}
            auto d=g.direction();auto r=counter.count(order,b,d);
            out<<j1<<",1,"<<d.first<<','<<d.second<<','<<r.nu<<','<<r.edges<<','<<r.candidates<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("write failed");
        std::cout<<"prefixes="<<samples<<" seconds="
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
