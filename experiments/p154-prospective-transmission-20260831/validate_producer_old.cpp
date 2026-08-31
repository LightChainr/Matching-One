// Deterministic replay of ONLY original batch0: no prospective seed is invoked.
// Compiled separately, never used by the prospective shard driver.
#define P154_OBSERVER_ONLY 1
#include "producer.cpp"
int main(int argc,char**argv){
 try{
  if(argc!=2)throw std::invalid_argument("old-audit output.csv");
#ifdef MATCHING_NORM4_INTEGER
  const int n=340;const std::uint64_t old_seed=2026105402ULL,start=8200000000ULL;
#else
  const int n=85;const std::uint64_t old_seed=2026104501ULL,start=5100000000ULL;
#endif
  const auto d=*std::find_if(kDesigns.begin(),kDesigns.end(),[n](const PairDesign&d){return d.n==n;});
#ifdef MATCHING_NORM4_INTEGER
  Geometry first=make_geometry(d.first),second=make_geometry(d.second);
#else
  Geometry first=make_geometry(d.a1,d.b1),second=make_geometry(d.a2,d.b2);
#endif
  NewObserver a(first),b(second);std::array<std::vector<Cell>,2> v;v[0].resize(n+1);v[1].resize(n+1);std::vector<int>order;
  auto old_interval=[&](std::uint64_t offset,std::uint64_t count){for(std::uint64_t c=start+offset;c<start+offset+count;++c){
    counter_permutation(n,old_seed,c,order);a.observe(order,v[0]);b.observe(order,v[1]);}};
  old_interval(0,1000);if(n==340)old_interval(100000,9000);
  std::ofstream out(argv[1]);out<<"n,g,k,q,e,s,qs,es,count01,count02,count12,s_previous01,s_previous02,s_previous12\n";
  for(int g=0;g<2;++g)for(int k=0;k<=n;++k){out<<n<<','<<g<<','<<k;for(auto x:v[g][k])out<<','<<x;out<<'\n';}
  out.close();if(!out)throw std::runtime_error("old audit write failed");
  std::cout<<"OLD_BATCH_ONLY N="<<n<<" old_samples="<<(n==340?10000:1000)<<" prospective_samples=0 prospective_N_seed_not_used="<<domain_seed(n)<<'\n';
  return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
