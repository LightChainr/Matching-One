// First-birth prefix replay, followed by a uniform-safe reference chain.
// First-proposal survival flags exactly couple it to the killed uniform chain.
#define MATCHING_ONE_PAIR_LIBRARY
#include "../exact-completion-pairs-20260929/engine.cpp"

int main(int argc,char**argv){
    try{
        if(argc!=7)throw std::invalid_argument("L samples permutation_seed continuation_seed b output.csv");
        int L=std::stoi(argv[1]),N=L*L,b=std::stoi(argv[5]);
        auto samples=std::stoull(argv[2]),seed=std::stoull(argv[3]),continuation=std::stoull(argv[4]);
        if(L<3 || b<1 || b>N-2 || !samples)throw std::invalid_argument("bounds");
        Geometry g(L,false);ProjectedPairs counter(L,false);
        std::vector<int>order(N),reference,vacancies;reference.reserve(N);vacancies.reserve(N);
        std::mt19937_64 rng(seed),future(continuation);
        std::ofstream out(argv[6]);if(!out)throw std::runtime_error("cannot open output");
        out<<"J1,entry_rank,dx_entry,dy_entry,reference_rank_b,natural_alive,nu_b,synergy_edges,first_rejection_count,proposal_count\n";
        auto start=std::chrono::steady_clock::now();
        for(std::uint64_t sample=0;sample<samples;++sample){
            std::iota(order.begin(),order.end(),0);
            for(int v=N-1;v>0;--v)std::swap(order[v],order[bounded(rng,v+1)]);
            g.reset();int j1=0,rank=0;
            for(int k=1;k<=b;++k){rank=g.insert(order[k-1]);if(rank){j1=k;break;}}
            if(rank!=1){out<<j1<<','<<rank<<",0,0,"<<rank<<",0,-1,-1,0,0\n";continue;}
            auto direction=g.direction();int alive=1,first_rejection=0;
            std::uint64_t proposals=0;
            reference.assign(order.begin(),order.begin()+j1);
            vacancies.assign(order.begin()+j1,order.end());
            for(int k=j1;k<b;++k){
                const int m=vacancies.size();bool accepted=false;
                // A fresh random ordering, without replacement. Rejected
                // vacancies remain vacant and eligible on the next count.
                for(int r=m-1;r>=0;--r){
                    std::swap(vacancies[r],vacancies[bounded(future,r+1)]);
                    const int v=vacancies[r];++proposals;
                    const bool completes=g.completes_if_inserted(v);
                    if(r==m-1 && completes && alive){alive=0;first_rejection=k+1;}
                    if(completes)continue;
                    if(g.insert(v)!=1)throw std::logic_error("unsafe accepted site");
                    reference.push_back(v);
                    std::swap(vacancies[r],vacancies.back());vacancies.pop_back();
                    accepted=true;break;
                }
                if(!accepted){rank=-1;alive=0;break;} // exact cemetery, no retry cap
            }
            if(rank==-1){
                out<<j1<<",1,"<<direction.first<<','<<direction.second<<",-1,0,-1,-1,"
                   <<first_rejection<<','<<proposals<<'\n';continue;
            }
            if(int(reference.size())!=b || g.direction()!=direction)throw std::logic_error("reference endpoint");
            auto pairs=counter.count(reference,b,direction);
            out<<j1<<",1,"<<direction.first<<','<<direction.second<<",1,"<<alive<<','<<pairs.nu<<','<<pairs.edges<<','
               <<first_rejection<<','<<proposals<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("write failed");
        std::cout<<"prefixes="<<samples<<" seconds="
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
