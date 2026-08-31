// Reobserve only the two existing Phase-E 20k counter blocks. No new samples.
// The included backend has the same Git blob as the original 0578105 runner.
#define main threshold_rank_integer_period_hidden_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main
#include <set>

namespace {
constexpr double p_ref = 0.59274605079;
int replay_count(int n, std::uint64_t seed, std::uint64_t replica) {
    SplitMixStream stream(splitmix64(seed ^ splitmix64(replica + 0x8cb92ba72f3d8dd7ULL)));
    int k = 0;
    for (int i = 0; i < n; ++i)
        k += static_cast<double>(stream.next() >> 11) / 9007199254740992.0 < p_ref;
    return k;
}
struct Sums {
    std::uint64_t samples=0, k1=0, k2=0, i0=0, i1=0, i2=0;
    std::uint64_t k=0, kk=0, edges=0, i0k=0, i2k=0, i0kk=0, i2kk=0, i0edges=0, i2edges=0;
    void add(int a, int b, int count, int occupied_edges) {
        const std::uint64_t z0=count<a, z2=count>=b, factorial=static_cast<std::uint64_t>(count)*(count-1);
        ++samples; k1+=a; k2+=b; i0+=z0; i2+=z2; i1+=1-z0-z2;
        k+=count; kk+=factorial; edges+=occupied_edges;
        i0k+=z0*count; i2k+=z2*count; i0kk+=z0*factorial; i2kk+=z2*factorial;
        i0edges+=z0*occupied_edges; i2edges+=z2*occupied_edges;
    }
};
void write(std::ostream& out, int n, int a, int b, const char* direction, int batch, const Sums& s) {
    out<<n<<','<<a<<','<<b<<','<<direction<<','<<batch<<','<<s.samples<<','<<s.k1<<','<<s.k2
       <<','<<s.i0<<','<<s.i1<<','<<s.i2<<','<<s.k<<','<<s.kk<<','<<s.edges
       <<','<<s.i0k<<','<<s.i2k<<','<<s.i0kk<<','<<s.i2kk<<','<<s.i0edges<<','<<s.i2edges<<'\n';
}
int count_edges(const Geometry& g, const std::vector<int>& permutation, int k) {
    std::vector<std::uint8_t> occupied(g.n,0);
    for(int i=0;i<k;++i) occupied[permutation[i]]=1;
    int total=0;
    for(const Edge& edge:g.primal_edges) total+=occupied[edge.i] && occupied[edge.j];
    return total;
}
void require_simple_nn(const Geometry& g) {
    std::set<std::pair<int,int>> seen;
    for(const Edge& e:g.primal_edges) {
        if(e.i==e.j || !seen.insert(std::minmax(e.i,e.j)).second)
            throw std::runtime_error("fixed-K formula requires a simple NN graph");
    }
    if(seen.size()!=2U*g.n) throw std::runtime_error("expected 2N undirected NN edges");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=3) throw std::runtime_error("usage: replay 65|130 output.csv");
        const int n=std::stoi(argv[1]);
        if(n!=65 && n!=130) throw std::runtime_error("only archived N65/N130 blocks are supported");
        const bool small=n==65;
        const Matrix m1=small?Matrix{8,-1,1,8}:Matrix{11,-3,3,11};
        const Matrix m2=small?Matrix{7,-4,4,7}:Matrix{9,-7,7,9};
        const Geometry g1=make_geometry(m1),g2=make_geometry(m2);
        require_simple_nn(g1); require_simple_nn(g2);
        const std::uint64_t seed=small?202615465ULL:2026154130ULL;
        const std::uint64_t offset=small?15466000000ULL:15466200000ULL;
        if(std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite a replay artifact");
        std::ofstream out(argv[2]);
        if(!out) throw std::runtime_error("cannot create output");
        out<<"n,a,b,orientation,batch,samples,sum_k1,sum_k2,sum_i0,sum_i1,sum_i2,"
              "sum_k,sum_kk,sum_edges,sum_i0k,sum_i2k,sum_i0kk,sum_i2kk,sum_i0edges,sum_i2edges\n";
        const auto start=std::chrono::steady_clock::now();
        ThresholdEngine e1(g1),e2(g2);
        std::vector<int> permutation;
        for(int batch=0;batch<100;++batch) {
            Sums first,second;
            for(int j=0;j<200;++j) {
                const std::uint64_t replica=offset+200ULL*batch+j;
                counter_permutation(n,seed,replica,permutation);
                const int k=replay_count(n,seed,replica);
                const auto r1=e1.ranks(permutation),r2=e2.ranks(permutation);
                first.add(r1.first,r1.second,k,count_edges(g1,permutation,k));
                second.add(r2.first,r2.second,k,count_edges(g2,permutation,k));
            }
            write(out,n,small?8:11,small?1:3,"first",batch,first);
            write(out,n,small?7:9,small?4:7,"second",batch,second);
        }
        std::cout<<"N="<<n<<" archived_samples=20000 batches=100 elapsed_seconds="
                 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
}
