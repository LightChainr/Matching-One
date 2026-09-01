// Prospective lag1 production. DO NOT RUN before root's actual freeze approval.
// Backend geometry/UF are verbatim bfab0330. Fresh tagged N/seed/counter domain.
#ifndef MATCHING_NORM4_BACKEND
#error "Provide packaged immutable backend"
#endif
#define main unused_archived_production_main
#include MATCHING_NORM4_BACKEND
#undef main

namespace {
using Int=std::int64_t;
// q,E,s,qs,Es,count01,count02,count12,s_previous01,s_previous02,s_previous12
using Cell=std::array<Int,11>;
class NewObserver {
 public:
  explicit NewObserver(const Geometry&g):g_(g),active_(g.n),
#ifdef MATCHING_NORM4_INTEGER
    uf_(g.quotient),
#else
    uf_(g.n,g.a,g.b),
#endif
    black_(g.n+1),white_(g.n+1){}
  void observe(const std::vector<int>&order,std::vector<Cell>&cells){
    const int kp=sweep(order,false,false,black_);
    const int km=g_.n-sweep(order,true,true,white_)+1;
    if(!(1<=km && km<=kp && kp<=g_.n))throw std::logic_error("threshold ordering");
    for(int k=0;k<=g_.n;++k){
      const Int q=-1+(k>=km)+(k>=kp),e=q*q,s=black_[k]+white_[g_.n-k];
      const std::array<Int,5> v{{q,e,s,q*s,e*s}};
      for(int j=0;j<5;++j)cells[k][j]+=v[j];
    }
    auto mark=[&](int k,int type){++cells[k][5+type];cells[k][8+type]+=black_[k-1]+white_[g_.n-k+1];};
    if(km==kp)mark(km,1);else{mark(km,0);mark(kp,2);}
  }
 private:
  int sweep(const std::vector<int>&o,bool matching,bool reverse,std::vector<int>&components){
    std::fill(active_.begin(),active_.end(),0);uf_.reset();components[0]=0;
    int count=0,cross=0;
    const auto&edges=matching?g_.matching_edges:g_.primal_edges;
    const auto&incident=matching?g_.matching_incident:g_.primal_incident;
    for(int offset=0;offset<g_.n;++offset){
      const int v=o[reverse?g_.n-1-offset:offset];active_[v]=1;++count;
      for(int ei:incident[v]){const Edge&e=edges[ei];if(active_[e.i] && active_[e.j]){
        if(uf_.find(e.i).root!=uf_.find(e.j).root)--count;uf_.add_edge(e);}}
      components[offset+1]=count;if(cross==0 && uf_.component_crosses(v))cross=offset+1;
    }
    if(cross==0 || count!=1)throw std::logic_error("invalid full graph");return cross;
  }
  const Geometry&g_;std::vector<std::uint8_t>active_;HomologyUnionFind uf_;
  std::vector<int>black_,white_;
};
constexpr std::uint64_t master_seed=0x5031353450524f31ULL;
constexpr std::uint64_t first_counter=30000000000ULL;
constexpr int batches=200;
std::uint64_t domain_seed(int n){return splitmix64(master_seed^splitmix64(static_cast<std::uint64_t>(n)^0x4e4f524d345631ULL));}
}
#ifndef P154_OBSERVER_ONLY
int main(int argc,char**argv){
 try{
  if(argc!=7)throw std::invalid_argument("N batch_begin batch_end output.csv threads freeze_commit; invoke only through authorized driver");
  const int n=std::stoi(argv[1]),begin=std::stoi(argv[2]),end=std::stoi(argv[3]),threads=std::stoi(argv[5]);
  const std::string freeze=argv[6];
  if(freeze.size()!=40 || freeze.find_first_not_of("0123456789abcdef")!=std::string::npos)throw std::invalid_argument("actual freeze commit required");
  if(!(0<=begin && begin<end && end<=batches && threads>=1 && threads<=14))throw std::invalid_argument("batch/worker bounds");
#ifdef MATCHING_NORM4_INTEGER
  if(n!=340)throw std::invalid_argument("integer backend only N340");
#else
  if(n!=85)throw std::invalid_argument("primitive backend only N85");
#endif
  const std::uint64_t total=n==85?5000000ULL:160000000ULL,per_batch=total/batches,seed=domain_seed(n);
  if(std::filesystem::exists(argv[4]))throw std::runtime_error("refusing overwrite");
  auto found=std::find_if(kDesigns.begin(),kDesigns.end(),[n](const PairDesign&d){return d.n==n;});
  if(found==kDesigns.end())throw std::logic_error("design absent");const PairDesign&d=*found;
#ifdef MATCHING_NORM4_INTEGER
  const Geometry first=make_geometry(d.first),second=make_geometry(d.second);
#else
  const Geometry first=make_geometry(d.a1,d.b1),second=make_geometry(d.a2,d.b2);
#endif
  std::vector<std::array<std::vector<Cell>,2>>output(end-begin);
  #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
  for(int b=begin;b<end;++b){
    NewObserver f(first),s(second);std::vector<int>order;auto&row=output[b-begin];row[0].resize(n+1);row[1].resize(n+1);
    const auto start=first_counter+per_batch*static_cast<std::uint64_t>(b);
    for(std::uint64_t c=start;c<start+per_batch;++c){counter_permutation(n,seed,c,order);f.observe(order,row[0]);s.observe(order,row[1]);}
  }
  std::ofstream out(argv[4]);if(!out)throw std::runtime_error("output open");
  out<<"n,orientation,batch,k,samples,sum_q,sum_e,sum_s,sum_qs,sum_es,event_count01,event_count02,event_count12,sum_s_previous01,sum_s_previous02,sum_s_previous12\n";
  for(int b=begin;b<end;++b)for(int g=0;g<2;++g)for(int k=0;k<=n;++k){
    out<<n<<','<<(g==0?"first":"second")<<','<<b<<','<<k<<','<<per_batch;
    for(Int x:output[b-begin][g][k])out<<','<<x;out<<'\n';
  }
  out.close();if(!out)throw std::runtime_error("output write");
  std::cout<<"domain=P154_PROSPECTIVE_TRANSMISSION_V1_20260831 N="<<n<<" N_seed="<<seed
           <<" counters=["<<first_counter+per_batch*begin<<','<<first_counter+per_batch*end
           <<") samples="<<per_batch*(end-begin)<<" freeze_commit="<<freeze<<" threads="<<threads<<'\n';
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
#endif
