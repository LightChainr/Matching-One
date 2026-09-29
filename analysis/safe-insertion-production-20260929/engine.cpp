// Independent-prefix safe-insertion probes. Reuse unchanged tested geometry.
#define main original_completion_engine_main
#include "../completion-hazard-production-20260929/engine.cpp"
#undef main

struct SafeObservation {
    int j1=0, rank_b=0, dx=0, dy=0, nu=-1;
    std::array<int,4> y{{-1,-1,-1,-1}};
};

SafeObservation probe(Geometry& g,const std::vector<int>& order,int b,
                      std::mt19937_64& probe_rng) {
    g.reset(); SafeObservation o;
    int N=static_cast<int>(order.size());
    for(int k=1;k<=b;++k){
        o.rank_b=g.insert(order[k-1]);
        if(o.rank_b>=1 && !o.j1)o.j1=k;
        if(o.rank_b==2)return o; // rank cannot subsequently decrease
    }
    if(o.rank_b!=1)return o;
    auto d=g.direction();o.dx=d.first;o.dy=d.second;
    o.nu=g.completion_count();
    if(o.nu==N-b)return o; // no safe insertion; conditional probe undefined
    for(int j=0;j<4;++j){
        // An independent probe stream chooses uniformly among vacancies;
        // reject completions, yielding uniform S(A). Copies share no edits.
        while(true){
            int v=order[b+bounded(probe_rng,N-b)];
            Geometry copy=g;
            if(copy.insert(v)==2)continue;
            if(copy.direction()!=d)throw std::logic_error("safe direction changed");
            o.y[j]=copy.completion_count()-o.nu;
            if(o.y[j]<0)throw std::logic_error("completion count decreased");
            break;
        }
    }
    return o;
}

int main(int argc,char**argv){
    try{
        if(argc==2 && std::string(argv[1])=="control"){
            // One known physical pair, not a repeat of the L4 census.
            const std::array<int,2> masks{{4403,4405}};
            const std::array<std::array<int,3>,2> expected{{{{6,4,0}},{{6,2,2}}}};
            for(int i=0;i<2;++i){
                Geometry g(4,false);int rank=0;
                for(int v=0;v<16;++v)if(masks[i]>>v&1)rank=g.insert(v);
                if(rank!=1 || g.direction()!=std::make_pair(0,1) || g.completion_count()!=0)
                    throw std::logic_error("physical witness prefix mismatch");
                std::array<int,3> counts{{0,0,0}};
                for(int v=0;v<16;++v)if(!(masks[i]>>v&1)){
                    Geometry copy=g;
                    if(copy.insert(v)!=1)throw std::logic_error("witness unsafe");
                    int nu=copy.completion_count();
                    if(nu<0 || nu>2)throw std::logic_error("witness degree bounds");
                    ++counts[nu];
                }
                if(counts!=expected[i])throw std::logic_error("witness successor counts");
                std::cout<<"mask="<<masks[i]<<" next_nu_counts="<<counts[0]<<','<<counts[1]<<','<<counts[2]<<'\n';
            }
            return 0;
        }
        if(argc!=7)throw std::invalid_argument("L samples permutation_seed probe_seed b output.csv");
        int L=std::stoi(argv[1]),N=L*L,b=std::stoi(argv[5]);
        auto samples=std::stoull(argv[2]),seed=std::stoull(argv[3]),probe_seed=std::stoull(argv[4]);
        if(!samples || b<1 || b>N-2)throw std::invalid_argument("sample/count bounds");
        Geometry g(L,false);std::vector<int>order(N);
        std::mt19937_64 rng(seed),probe_rng(probe_seed);
        std::ofstream out(argv[6]);if(!out)throw std::runtime_error("cannot open output");
        out<<"J1,rank_b,dx_b,dy_b,nu_b,Y0,Y1,Y2,Y3\n";
        auto start=std::chrono::steady_clock::now();
        for(std::uint64_t i=0;i<samples;++i){
            std::iota(order.begin(),order.end(),0);
            for(int v=N-1;v>0;--v)std::swap(order[v],order[bounded(rng,v+1)]);
            auto o=probe(g,order,b,probe_rng);
            out<<o.j1<<','<<o.rank_b<<','<<o.dx<<','<<o.dy<<','<<o.nu;
            for(auto y:o.y)out<<','<<y;
            out<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("write failed");
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"samples="<<samples<<" seconds="<<seconds<<'\n';return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
