// One frozen positive-policy contrast experiment; never run before F0.
#define P334_CONTACT_LIBRARY_ONLY
#include "vendor/p334_next_label_contact_coordinates.cpp"

namespace prospective {
constexpr int baseline_draws=32, contrast_draws=64;
constexpr std::uint64_t fork_seed=202608311638334ULL;
struct Label { int v,safe,e[2],c[2],rank[2],h[2]; };
struct Row { std::uint64_t index,counter; int batch,k0,r[2]; double feature[2][4],weight[2]; };

int checkpoint_rank(Checkpoint& z) {
    std::vector<Vector> lines;
    for(int v=0;v<z.geometry.n;++v) if(z.active[v]) {
        auto m=z.uf.component_mark(v);
        if(m.rank==2)return 2;
        if(m.rank==1)lines.push_back(m.line);
    }
    return image_rank(lines);
}
int rank_after_r0(Checkpoint& z,int v) {
    std::map<int,Vector> anchor; std::vector<Vector> windings;
    for(int ei:z.geometry.primal_incident[v]) {
        const auto& e=z.geometry.primal_edges[ei]; int other=e.i==v?e.j:e.i;
        if(!z.active[other])continue;
        auto f=z.uf.find(other); Vector step=e.i==v?Vector{e.dx,e.dy}:Vector{-e.dx,-e.dy};
        Vector a{step.x-f.dx,step.y-f.dy}; auto inserted=anchor.emplace(f.root,a);
        if(!inserted.second)windings.push_back(z.geometry.quotient.winding(a.x-inserted.first->second.x,a.y-inserted.first->second.y));
    }
    return image_rank(windings);
}
std::uint64_t key(int n,std::uint64_t index,int kind,int draw) {
    if(index>=(1ULL<<20)||kind<0||kind>15||draw<0||draw>=256)throw std::runtime_error("frozen RNG address overflow");
    // Old8 and old-prefix New64 used addresses below bit41. Same bijective
    // mapping/fork seed; bit63 separates this experiment's starting-key domain.
    std::uint64_t address=(1ULL<<63)|(static_cast<std::uint64_t>(n)<<48)|(index<<20)|(static_cast<std::uint64_t>(kind)<<16)|(static_cast<std::uint64_t>(draw)<<2);
    return splitmix64(fork_seed^splitmix64(address));
}
std::vector<std::string> split(const std::string& line) {
    std::stringstream in(line); std::vector<std::string> out;std::string x;
    while(std::getline(in,x,','))out.push_back(x); return out;
}
std::string stem(int n,int shard) {
    std::ostringstream out;out<<"N"<<n<<".shard"<<std::setfill('0')<<std::setw(3)<<shard;return out.str();
}
std::vector<Label> census(Checkpoint& a,Checkpoint& b,const std::vector<int>& vacant,Row& row) {
    int d=vacant.size();std::vector<Label> ls;std::map<int,std::array<long long,3>> sums;
    for(int v:vacant) {
        int r0=rank_after_r0(a,v),r1=rank_after_r0(b,v);
        auto x=a.read(v,0,r0),y=b.read(v,0,r1);
        if(x.e>4||y.e>4||x.cycles<0||x.cycles>3||y.cycles<0||y.cycles>3)throw std::runtime_error("degree/loop bound drift");
        Label l{v,(r0==0&&r1==0),{x.e,y.e},{x.c,y.c},{r0,r1},{0,0}};
        if(l.safe) {auto& s=sums[5*x.e+y.e];++s[0];s[1]+=x.cycles;s[2]+=y.cycles;}
        ls.push_back(l);
    }
    for(auto& l:ls)if(l.safe) {
        const auto& s=sums.at(5*l.e[0]+l.e[1]);
        for(int o=0;o<2;++o) {
            l.h[o]=static_cast<int>(s[0]*(l.e[o]-l.c[o])-s[o+1]);
            if(std::abs(l.h[o])>3*d)throw std::runtime_error("H bound drift");
            row.feature[o][0]+=1.0/d;
            row.feature[o][1]+=static_cast<double>(l.h[o])*l.h[o]/(static_cast<double>(d)*d*d);
            row.feature[o][2]+=static_cast<double>(l.e[o])/d;
            row.feature[o][3]+=static_cast<double>(l.e[o]-l.c[o])/d;
            row.weight[o]+=static_cast<double>(std::max(l.h[o],0))/(static_cast<double>(d)*d);
        }
    }
    for(const auto& group:sums)for(int o=0;o<2;++o) {
        long long check=0;for(const auto& l:ls)if(l.safe&&5*l.e[0]+l.e[1]==group.first)check+=l.h[o];
        if(check)throw std::runtime_error("class score not exactly centered");
    }
    return ls;
}
void features(int n,int shard,int per_shard,std::uint64_t counter_begin,std::uint64_t seed,const std::filesystem::path& out,const std::string& freeze) {
    auto first_geometry=make_geometry({n,n==325?57:132,0,1}),second_geometry=make_geometry({n,n==325?18:268,0,1});
    Checkpoint first(first_geometry),second(second_geometry);
    std::filesystem::create_directories(out);
    auto base=stem(n,shard);
    auto ppath=out/(base+".prefix.csv.gz"),cpath=out/(base+".census.csv.gz");
    if(std::filesystem::exists(ppath)||std::filesystem::exists(cpath))throw std::runtime_error("refuse to overwrite prefix/census");
    gzFile p=gzopen(ppath.c_str(),"wb1"),c=gzopen(cpath.c_str(),"wb1");if(!p||!c)throw std::runtime_error("cannot open blind files");
    gzprintf(p,"index,batch,counter,k0,rank_first,rank_second,first_mass,first_energy,first_degree,first_loop,second_mass,second_energy,second_degree,second_loop,W_first,W_second\n");
    gzprintf(c,"index,label,safe,e_first,c_first,rank_first,e_second,c_second,rank_second,hnum_first,hnum_second\n");
    std::vector<int> perm,vacant;std::array<int,9> cells{};std::uint64_t sites=0;
    for(int j=0;j<per_shard;++j) {
        std::uint64_t index=static_cast<std::uint64_t>(shard)*per_shard+j;
        Row row{};row.index=index;row.counter=counter_begin+index;row.batch=shard;row.k0=n==325?193:252;
        counter_permutation(n,seed,row.counter,perm);first.load(perm,row.k0);second.load(perm,row.k0);
        row.r[0]=checkpoint_rank(first);row.r[1]=checkpoint_rank(second);++cells[3*row.r[0]+row.r[1]];
        if(row.r[0]==0&&row.r[1]==0) {
            vacant.assign(perm.begin()+row.k0,perm.end());std::sort(vacant.begin(),vacant.end());
            auto labels=census(first,second,vacant,row);
            for(const auto& l:labels) {
                if(gzprintf(c,"%llu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",(unsigned long long)index,l.v,l.safe,l.e[0],l.c[0],l.rank[0],l.e[1],l.c[1],l.rank[1],l.h[0],l.h[1])<=0)throw std::runtime_error("census write error");
                ++sites;
            }
        }
        if(gzprintf(p,"%llu,%d,%llu,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(unsigned long long)row.index,row.batch,(unsigned long long)row.counter,row.k0,row.r[0],row.r[1],row.feature[0][0],row.feature[0][1],row.feature[0][2],row.feature[0][3],row.feature[1][0],row.feature[1][1],row.feature[1][2],row.feature[1][3],row.weight[0],row.weight[1])<=0)throw std::runtime_error("prefix write error");
    }
    if(gzclose(p)!=Z_OK||gzclose(c)!=Z_OK)throw std::runtime_error("blind gzip close error");
    std::ofstream meta(out/(base+".features.json"));
    meta<<"{\"freeze_commit\":\""<<freeze<<"\",\"N\":"<<n<<",\"shard\":"<<shard<<",\"all_prefixes\":"<<per_shard<<",\"cell_counts\":[";
    for(int i=0;i<9;++i)meta<<(i?",":"")<<cells[i];
    meta<<"],\"census_label_rows\":"<<sites<<",\"tail_paths\":0,\"future_rank_engine_calls\":0}\n";
}
void shuffle(std::vector<int>& order,std::uint64_t seed) {
    SplitMixStream rng(seed);for(int stop=static_cast<int>(order.size())-1;stop>0;--stop)std::swap(order[stop],order[rng.below(stop+1)]);
}
int choose(const std::vector<Label>& labels,int source,bool positive,int cls,SplitMixStream& rng) {
    std::uint64_t total=0;
    for(const auto& l:labels)if(l.safe&&(cls<0||5*l.e[0]+l.e[1]==cls))total+=std::max(positive?l.h[source]:-l.h[source],0);
    if(total==0)throw std::runtime_error("empty contrast mass");auto pick=rng.below(total);
    for(std::size_t i=0;i<labels.size();++i) {const auto& l=labels[i];if(!l.safe||(cls>=0&&5*l.e[0]+l.e[1]!=cls))continue;
        auto w=static_cast<std::uint64_t>(std::max(positive?l.h[source]:-l.h[source],0));if(pick<w)return i;pick-=w;}
    throw std::runtime_error("weighted label draw failed");
}
void tails(int n,int shard,int per_shard,std::uint64_t counter_begin,std::uint64_t seed,const std::filesystem::path& out,const std::string& freeze,const std::string& prediction_hash) {
    auto base=stem(n,shard);auto target=out/(base+".tails.csv.gz");
    if(std::filesystem::exists(target))throw std::runtime_error("refuse to overwrite tails");
    std::ifstream lock(out/(base+".prediction.sha256"));std::string saved;lock>>saved;
    if(saved!=prediction_hash||saved.size()!=64)throw std::runtime_error("prediction must be saved before any tail");
    auto first_geometry=make_geometry({n,n==325?57:132,0,1}),second_geometry=make_geometry({n,n==325?18:268,0,1});
    ThresholdEngine first(first_geometry),second(second_geometry);
    std::map<std::uint64_t,std::vector<Label>> labels;
    gzFile input=gzopen((out/(base+".census.csv.gz")).c_str(),"rb");if(!input)throw std::runtime_error("missing census");
    char buffer[4096];gzgets(input,buffer,sizeof(buffer));
    while(gzgets(input,buffer,sizeof(buffer))) {auto f=split(buffer);if(f.size()!=11)throw std::runtime_error("census schema");
        Label l{std::stoi(f[1]),std::stoi(f[2]),{std::stoi(f[3]),std::stoi(f[6])},{std::stoi(f[4]),std::stoi(f[7])},{std::stoi(f[5]),std::stoi(f[8])},{std::stoi(f[9]),std::stoi(f[10])}};
        labels[std::stoull(f[0])].push_back(l);
    }gzclose(input);
    gzFile output=gzopen(target.c_str(),"wb1");if(!output)throw std::runtime_error("cannot open tails");
    gzprintf(output,"index,batch,counter,group,source,draw,arm,label,e_first,e_second,first_k1,first_k2,second_k1,second_k2\n");
    std::vector<int> old,vacant,order,permutation;std::uint64_t paths=0;
    const int k0=n==325?193:252;
    for(const auto& item:labels) {
        auto index=item.first;const auto& ls=item.second;
        if(index/static_cast<std::uint64_t>(per_shard)!=static_cast<std::uint64_t>(shard)||ls.size()!=static_cast<std::size_t>(n-k0))throw std::runtime_error("incomplete shard label census");
        auto counter=counter_begin+index;counter_permutation(n,seed,counter,old);
        vacant.assign(old.begin()+k0,old.end());std::sort(vacant.begin(),vacant.end());
        for(std::size_t i=0;i<ls.size();++i)if(ls[i].v!=vacant[i])throw std::runtime_error("frozen prefix/census label mismatch");
        auto emit=[&](char group,int source,int draw,int arm,int li,const std::vector<int>& remaining) {
            permutation.assign(old.begin(),old.begin()+k0);permutation.push_back(ls[li].v);
            for(int v:remaining)if(v!=ls[li].v)permutation.push_back(v);
            auto f=first.ranks(permutation),s=second.ranks(permutation);
            if(group=='I'&&(f.first<=k0+1||s.first<=k0+1||!ls[li].safe))throw std::runtime_error("intervention changed immediate rank");
            if(gzprintf(output,"%llu,%d,%llu,%c,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",(unsigned long long)index,shard,(unsigned long long)counter,group,source,draw,arm,ls[li].v,ls[li].e[0],ls[li].e[1],f.first,f.second,s.first,s.second)<=0)throw std::runtime_error("tail write error");
            ++paths;
        };
        for(int q=0;q<baseline_draws;++q) {
            order=vacant;shuffle(order,key(n,index,1,q));int label=order.front();int li=std::lower_bound(vacant.begin(),vacant.end(),label)-vacant.begin();emit('B',-1,q,0,li,order);
        }
        for(int source=0;source<2;++source) {
            long long total=0;for(const auto& l:ls)total+=std::max(l.h[source],0);if(!total)continue;
            for(int q=0;q<contrast_draws;++q) {
                SplitMixStream rng(key(n,index,4+source,q));int plus=choose(ls,source,true,-1,rng);
                int cls=5*ls[plus].e[0]+ls[plus].e[1];int minus=choose(ls,source,false,cls,rng);
                if(ls[plus].e[0]!=ls[minus].e[0]||ls[plus].e[1]!=ls[minus].e[1])throw std::runtime_error("paired coarse class drift");
                order=vacant;shuffle(order,key(n,index,8+source,q));emit('I',source,q,0,plus,order);emit('I',source,q,1,minus,order);
            }
        }
    }
    if(gzclose(output)!=Z_OK)throw std::runtime_error("tail gzip close error");
    std::ofstream meta(out/(base+".tails.json"));meta<<"{\"freeze_commit\":\""<<freeze<<"\",\"prediction_sha256\":\""<<prediction_hash<<"\",\"N\":"<<n<<",\"shard\":"<<shard<<",\"cell00_prefixes\":"<<labels.size()<<",\"baseline_draws\":32,\"contrasts_per_source\":64,\"tail_paths\":"<<paths<<",\"sampling\":\"Rao-Blackwellized same-class positive-policy coupling\"}\n";
}
}
#ifndef P334_PROSPECTIVE_LIBRARY_ONLY
int main(int argc,char** argv) {
    try {
        if(argc!=9&&argc!=10)throw std::runtime_error("usage: prospective features|tails N shard per_shard counter_begin prefix_seed out freeze_commit [prediction_sha256]");
        std::string mode=argv[1],freeze=argv[8];int n=std::stoi(argv[2]),shard=std::stoi(argv[3]),count=std::stoi(argv[4]);
        if((n!=325&&n!=425)||shard<0||shard>=60||count!=5000||freeze.size()!=40)throw std::runtime_error("outside frozen production domain");
        if(mode=="features"&&argc==9)prospective::features(n,shard,count,std::stoull(argv[5]),std::stoull(argv[6]),argv[7],freeze);
        else if(mode=="tails"&&argc==10)prospective::tails(n,shard,count,std::stoull(argv[5]),std::stoull(argv[6]),argv[7],freeze,argv[9]);
        else throw std::runtime_error("invalid production phase");
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
#endif
