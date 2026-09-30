// One source kick, original uniform continuation, existing permutation seeds.
// Reuses the exact projected-pair and lifted topology implementations.
#define MATCHING_ONE_PAIR_LIBRARY
#include "../exact-completion-pairs-20260929/engine.cpp"

int main(int argc,char**argv){
    try {
        if(argc!=7)throw std::invalid_argument("L prefixes existing_seed b lag output.csv");
        const int L=std::stoi(argv[1]),N=L*L,b=std::stoi(argv[4]),h=std::stoi(argv[5]);
        const auto samples=std::stoull(argv[2]),seed=std::stoull(argv[3]);
        if(!samples || b<1 || h<2 || b+h>N)throw std::invalid_argument("bounds");
        Geometry g(L,false);ProjectedPairs counter(L,false);
        std::vector<int>order(N),degree(N);std::mt19937_64 rng(seed);
        std::ofstream out(argv[6]);if(!out)throw std::runtime_error("output");
        out<<"J1,rank_b,dx_b,dy_b,nu_b,synergy_edges,sum_degree_squared,survival_h,window_degree_sum\n";
        const auto start=std::chrono::steady_clock::now();
        for(std::uint64_t row=0;row<samples;++row){
            std::iota(order.begin(),order.end(),0);
            for(int v=N-1;v>0;--v)std::swap(order[v],order[bounded(rng,v+1)]);
            g.reset();int j1=0,rank=0;
            for(int k=1;k<=b;++k){rank=g.insert(order[k-1]);if(rank && !j1)j1=k;if(rank==2)break;}
            if(rank!=1){out<<j1<<','<<rank<<",0,0,-1,-1,-1,-1,-1\n";continue;}
            const auto d=g.direction();const auto r=counter.count(order,b,d);
            std::fill(degree.begin(),degree.end(),0);
            for(auto edge:counter.edges){++degree[int(edge>>32)];++degree[int(edge&0xffffffffULL)];}
            std::uint64_t squared=0,window_sum=0;
            for(int value:degree)squared+=std::uint64_t(value)*value;
            for(int k=b;k<b+h;++k){
                const int v=order[k];window_sum+=degree[v];
                if(rank==1)rank=g.insert(v);
            }
            out<<j1<<",1,"<<d.first<<','<<d.second<<','<<r.nu<<','<<r.edges<<','
               <<squared<<','<<(rank==1)<<','<<window_sum<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("write failed");
        std::cout<<"prefixes="<<samples<<" lag="<<h<<" seconds="
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
