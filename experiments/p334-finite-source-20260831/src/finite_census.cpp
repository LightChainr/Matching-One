// Exact finite vacant-label census on the original checkpoints, never a tail.
#define P334_CONTACT_LIBRARY_ONLY
#include "p334_next_label_contact_coordinates.cpp"

namespace {
Vector essential_line(Checkpoint& checkpoint, int old_rank) {
    if (old_rank != 1) return {0,0};
    for (int v = 0; v < checkpoint.geometry.n; ++v) {
        if (!checkpoint.active[v]) continue;
        const auto mark = checkpoint.uf.component_mark(v);
        if (mark.rank == 1) return mark.line;
    }
    throw std::runtime_error("R1 checkpoint has no existing essential component");
}

int after_rank(Checkpoint& checkpoint, int vertex, int old_rank, Vector line) {
    if (old_rank == 2) return 2;
    std::map<int,Vector> anchors;
    Vector basis = line;
    int rank = old_rank;
    for (const int edge_index: checkpoint.geometry.primal_incident[vertex]) {
        const auto& edge = checkpoint.geometry.primal_edges[edge_index];
        const int other = edge.i == vertex ? edge.j : edge.i;
        if (!checkpoint.active[other]) continue;
        const auto found = checkpoint.uf.find(other);
        const Vector step = edge.i == vertex ? Vector{edge.dx,edge.dy} : Vector{-edge.dx,-edge.dy};
        const Vector alpha{step.x-found.dx,step.y-found.dy};
        const auto inserted = anchors.emplace(found.root,alpha);
        if (inserted.second) continue;
        const auto anchor = inserted.first->second;
        const auto w = checkpoint.geometry.quotient.winding(alpha.x-anchor.x,alpha.y-anchor.y);
        if (w.x == 0 && w.y == 0) continue;
        if (rank == 0) { rank=1; basis=w; }
        else if (static_cast<__int128>(basis.x)*w.y != static_cast<__int128>(basis.y)*w.x)
            return 2;
    }
    return rank;
}
}

int main(int argc, char** argv) {
 try {
  if(argc!=3) throw std::runtime_error("usage: finite_census N output_dir");
  int n=std::stoi(argv[1]);
  if(n!=325&&n!=425) throw std::runtime_error("N325/425 only");
  std::filesystem::path out(argv[2]);
  if(std::filesystem::exists(out)) throw std::runtime_error("refuse overwrite");
  std::filesystem::create_directories(out);
  auto prefixes=original_prefixes(n);
  auto geom0=make_geometry({n,n==325?57:132,0,1});
  auto geom1=make_geometry({n,n==325?18:268,0,1});
  Checkpoint first(geom0), second(geom1);
  uint64_t seed=n==325?20260831430325ULL:20260831430425ULL;
  gzFile f=gzopen((out/"census.csv.gz").c_str(),"wb1");
  if(!f) throw std::runtime_error("output failed");
  gzprintf(f,"N,batch,counter,k0,first_rank,second_rank,first_e,second_e,L_first,L_second,count\n");
  std::vector<int> perm; uint64_t positions=0,joint=0,rows=0;
  auto start=std::chrono::steady_clock::now();
  for(const auto& item:prefixes){
   auto counter=item.first;auto p=item.second;
   int batch=p[0],k0=p[1],r0=p[2],r1=p[3];
   counter_permutation(n,seed,counter,perm);
   first.load(perm,k0);second.load(perm,k0);
   auto l0=essential_line(first,r0),l1=essential_line(second,r1);
   std::map<std::array<int,4>,int> counts;
   for(int k=k0;k<n;k++){
    int v=perm[k];positions++;
    int a=after_rank(first,v,r0,l0),b=after_rank(second,v,r1,l1);
    if(a!=r0||b!=r1)continue;
    auto m0=first.read(v,r0,a),m1=second.read(v,r1,b);
    std::array<int,4> key{m0.e,m1.e,r0==0?m0.cycles:0,r1==0?m1.cycles:0};
    counts[key]++;joint++;
   }
   for(const auto& c:counts){
    auto x=c.first;
    gzprintf(f,"%d,%d,%llu,%d,%d,%d,%d,%d,%d,%d,%d\n",n,batch,(unsigned long long)counter,k0,r0,r1,x[0],x[1],x[2],x[3],c.second);rows++;
   }
  }
  if(gzclose(f)!=Z_OK)throw std::runtime_error("gzip failed");
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::ofstream meta(out/"metadata.json");
  meta<<std::setprecision(15)<<"{\"N\":"<<n<<",\"prefixes\":"<<prefixes.size()<<",\"positions\":"<<positions<<",\"joint_safe\":"<<joint<<",\"census_rows\":"<<rows<<",\"elapsed_seconds\":"<<seconds<<",\"new_samples\":0,\"tail_replays\":0}\n";
  std::cout<<"N"<<n<<" finite-class census "<<rows<<" rows in "<<seconds<<" sec\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
