// One new joint mark on each OLD norm4 permutation. No new RNG domain.
#ifndef MATCHING_NORM4_BACKEND
#error "Supply packaged immutable bfab0330 backend"
#endif
#define main unused_archived_main
#include MATCHING_NORM4_BACKEND
#undef main

namespace {
using I = std::int64_t;
// count, q_late, E_late, s_early, q_late*s_early, E_late*s_early,
// s_early^2, s_now, q_late*s_now, E_late*s_now, conditional on early rank.
using Cell = std::array<I,10>;
using Profile = std::vector<std::array<Cell,3>>;
class Replay {
 public:
  explicit Replay(const Geometry& g):g_(g),active_(g.n),
#ifdef MATCHING_NORM4_INTEGER
    uf_(g.quotient),
#else
    uf_(g.n,g.a,g.b),
#endif
    black_(g.n+1),white_(g.n+1) {}
  void observe(const std::vector<int>& order, Profile& cells, int lag) {
    const int kp=sweep(order,false,false,black_);
    const int km=g_.n-sweep(order,true,true,white_)+1;
    if (!(1<=km && km<=kp && kp<=g_.n)) throw std::logic_error("invalid old thresholds");
    for(int k=0;k<=g_.n;++k) {
      const int l=std::max(0,k-lag),r=(l>=km)+(l>=kp);
      const I q=-1+(k>=km)+(k>=kp),e=q*q;
      const I s=black_[l]+white_[g_.n-l],now=black_[k]+white_[g_.n-k];
      const Cell v{{1,q,e,s,q*s,e*s,s*s,now,q*now,e*now}};
      for(std::size_t j=0;j<v.size();++j)cells[k][r][j]+=v[j];
    }
  }
 private:
  int sweep(const std::vector<int>& order,bool matching,bool reverse,std::vector<int>& components) {
    std::fill(active_.begin(),active_.end(),0);uf_.reset();components[0]=0;
    int count=0,cross=0;
    const auto& edges=matching?g_.matching_edges:g_.primal_edges;
    const auto& incident=matching?g_.matching_incident:g_.primal_incident;
    for(int offset=0;offset<g_.n;++offset) {
      const int v=order[reverse?g_.n-1-offset:offset];active_[v]=1;++count;
      for(int ei:incident[v]) {
        const Edge& e=edges[ei];
        if(active_[e.i] && active_[e.j]) {
          if(uf_.find(e.i).root!=uf_.find(e.j).root)--count;
          uf_.add_edge(e);
        }
      }
      components[offset+1]=count;
      if(cross==0 && uf_.component_crosses(v))cross=offset+1;
    }
    if(cross==0 || count!=1)throw std::logic_error("invalid full graph");
    return cross;
  }
  const Geometry& g_;std::vector<std::uint8_t> active_;HomologyUnionFind uf_;
  std::vector<int> black_,white_;
};
}

int main(int argc,char**argv) {
 try {
  if(argc!=4)throw std::invalid_argument("usage: replay N output.csv threads");
  const int n=std::stoi(argv[1]),threads=std::stoi(argv[3]);
  if(threads<1 || threads>16)throw std::invalid_argument("threads must be1..16");
#ifdef MATCHING_NORM4_INTEGER
  if(n!=260 && n!=340)throw std::invalid_argument("integer backend N");
#else
  if(n!=65 && n!=85 && n!=130 && n!=170)throw std::invalid_argument("primitive backend N");
#endif
  const auto it=std::find_if(kDesigns.begin(),kDesigns.end(),[n](const PairDesign&d){return d.n==n;});
  if(it==kDesigns.end())throw std::logic_error("missing old design");
  const PairDesign& d=*it;
#ifdef MATCHING_NORM4_INTEGER
  const Geometry first=make_geometry(d.first),second=make_geometry(d.second);
#else
  const Geometry first=make_geometry(d.a1,d.b1),second=make_geometry(d.a2,d.b2);
#endif
  int lag=0;while(lag*lag<n)++lag;
  const bool endpoint=n==260 || n==340;
  const std::uint64_t seed=n==260?2026105401ULL:n==340?2026105402ULL:2026104501ULL;
  const std::uint64_t start=endpoint?8200000000ULL:5100000000ULL;
  if(std::filesystem::exists(argv[2]))throw std::runtime_error("refusing overwrite");
  std::vector<std::array<Profile,2>> batches(100);
  #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
  for(int b=0;b<100;++b) {
    Replay r0(first),r1(second);std::vector<int> order;
    batches[b][0].resize(n+1);batches[b][1].resize(n+1);
    auto interval=[&](std::uint64_t begin,std::uint64_t count) {
      for(std::uint64_t counter=begin;counter<begin+count;++counter) {
        counter_permutation(n,seed,counter,order);
        r0.observe(order,batches[b][0],lag);r1.observe(order,batches[b][1],lag);
      }
    };
    interval(start+1000ULL*b,1000);
    if(endpoint)interval(start+100000ULL+9000ULL*b,9000);
  }
  std::ofstream out(argv[2]);if(!out)throw std::runtime_error("output open failed");
  out<<"n,orientation,batch,k,source_k,lag,early_rank,count,sum_q,sum_e,sum_s_early,sum_qs_early,sum_es_early,sum_s2_early,sum_s_now,sum_qs_now,sum_es_now\n";
  for(int b=0;b<100;++b)for(int g=0;g<2;++g)for(int k=0;k<=n;++k)for(int r=0;r<3;++r) {
    out<<n<<','<<(g==0?"first":"second")<<','<<b<<','<<k<<','<<std::max(0,k-lag)<<','<<lag<<','<<r;
    for(I x:batches[b][g][k][r])out<<','<<x;
    out<<'\n';
  }
  out.close();if(!out)throw std::runtime_error("output write failed");
  std::cout<<"N="<<n<<" lag="<<lag<<" old_permutations="<<(endpoint?1000000:100000)
           <<" paired_geometries=2 batches=100 threads="<<threads<<" new_counters=0\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
