// Reobserve the old P40 configurations by K. No new counter or p point.
// These joint sufficient statistics support analytic Bernoulli reweighting.
#ifndef MATCHING_P40_BACKEND
#error "Supply the immutable P40 backend through the replay driver"
#endif
#define main p40_archived_unused_thermal_main
#include MATCHING_P40_BACKEND
#undef main

namespace {
using ThermalInt=std::int64_t;
using ThermalCell=std::array<ThermalInt,6>;
constexpr double thermal_p=0.59274605079;
constexpr std::uint64_t thermal_seed=4020260830ULL;
constexpr std::uint64_t thermal_offset=40000000ULL;
constexpr int thermal_batches=100, thermal_batch_samples=10000;

void thermal_observe(std::vector<ThermalCell>& cells,int k,const Geometry& g,
                     const std::vector<std::uint8_t>& black,
                     const std::vector<std::uint8_t>& white,
                     HomologyUnionFind& black_uf,HomologyUnionFind& white_uf) {
    const Channels p=classify(black,g.primal_edges,black_uf);
    const Channels w=classify(white,g.matching_edges,white_uf);
    const ThermalInt q=static_cast<int>(p.either)-static_cast<int>(w.either);
    const ThermalInt e=q*q;
    const ThermalInt s=component_count(black,black_uf)+component_count(white,white_uf);
    const ThermalCell values{1,q,e,s,q*s,e*s};
    for (std::size_t i=0;i<values.size();++i) cells[k][i]+=values[i];
}

void thermal_write(std::ostream& out,const Geometry& g,const char* direction,
                   int batch,const std::vector<ThermalCell>& cells) {
    for (int k=0;k<=g.n;++k) {
        out<<g.n<<','<<g.a<<','<<g.b<<','<<direction<<','<<batch<<','<<k;
        for (const ThermalInt value:cells[k]) out<<','<<value;
        out<<'\n';
    }
}
} // namespace

int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::runtime_error("usage: p40-source-thermal-replay 65|85 output.csv");
        const int n=std::stoi(argv[1]);
        if (n!=65 && n!=85) throw std::runtime_error("only frozen P40 N65/N85 are supported");
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite a replay artifact");
        const Geometry first=n==65?make_geometry(8,1):make_geometry(9,2);
        const Geometry second=n==65?make_geometry(7,4):make_geometry(7,6);
        HomologyUnionFind fb(n,first.a,first.b),fw(n,first.a,first.b);
        HomologyUnionFind sb(n,second.a,second.b),sw(n,second.a,second.b);
        std::vector<std::uint8_t> black(n),white(n);
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create output");
        out<<"n,a,b,orientation,batch,k,samples,sum_q,sum_e,sum_s,sum_qs,sum_es\n";
        for (int batch=0;batch<thermal_batches;++batch) {
            std::vector<ThermalCell> f(n+1),s(n+1);
            const std::uint64_t begin=thermal_offset+static_cast<std::uint64_t>(batch)*thermal_batch_samples;
            for (std::uint64_t counter=begin;counter<begin+thermal_batch_samples;++counter) {
                int k=0;
                for (int site=0;site<n;++site) {
                    black[site]=counter_uniform(thermal_seed,n,counter,site)<thermal_p;
                    white[site]=!black[site];
                    k+=black[site];
                }
                thermal_observe(f,k,first,black,white,fb,fw);
                thermal_observe(s,k,second,black,white,sb,sw);
            }
            thermal_write(out,first,"first",batch,f);
            thermal_write(out,second,"second",batch,s);
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing output");
        std::cout<<"N="<<n<<": observed old 1000000 counters in 100 paired K-stratified batches; new samples=0\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
