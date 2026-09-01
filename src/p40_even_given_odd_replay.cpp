// Reobserve the frozen P40 million counters, using the exact archived backend.
// Only missing E_top/source/control products are new; no new configurations.
// The driver supplies MATCHING_P40_BACKEND after checking its immutable blob.
#ifndef MATCHING_P40_BACKEND
#error "Compile with the archived P40 backend supplied by the replay driver"
#endif
#define main p40_archived_engine_unused_main
#include MATCHING_P40_BACKEND
#undef main

namespace {
using EvenInt = std::int64_t;
constexpr double even_p = 0.59274605079;
constexpr std::uint64_t even_seed = 4020260830ULL;
constexpr std::uint64_t even_offset = 40000000ULL;
constexpr int even_batches = 100;
constexpr int even_batch_samples = 10000;
constexpr std::array<const char*, 12> even_names = {
    "q", "e", "s", "k", "kk", "edges", "chi",
    "e_s", "e_k", "e_kk", "e_edges", "e_chi"
};

struct EvenSums {
    std::array<EvenInt, even_names.size()> sum{};
    void add(const Geometry& geometry, const std::vector<std::uint8_t>& black,
             const std::vector<std::uint8_t>& white,
             HomologyUnionFind& black_uf, HomologyUnionFind& white_uf) {
        const Channels p = classify(black, geometry.primal_edges, black_uf);
        const Channels w = classify(white, geometry.matching_edges, white_uf);
        const EvenInt q = static_cast<int>(p.either) - static_cast<int>(w.either);
        const EvenInt e = q*q;
        const EvenInt cb = component_count(black, black_uf);
        const EvenInt cw = component_count(white, white_uf);
        const EvenInt s = cb+cw;
        const SiteMotifs motifs = count_site_motifs(black, geometry);
        const EvenInt k = motifs.V;
        const EvenInt kk = k*(k-1);
        const EvenInt edges = motifs.E;
        const EvenInt chi = k-edges+motifs.F0;
        const std::array<EvenInt, even_names.size()> values = {
            q,e,s,k,kk,edges,chi,e*s,e*k,e*kk,e*edges,e*chi
        };
        for (std::size_t j=0; j<values.size(); ++j) sum[j] += values[j];
    }
};

void even_write(std::ostream& out, const Geometry& g, const char* direction,
                int batch, const EvenSums& sums) {
    out << g.n << ',' << g.a << ',' << g.b << ',' << direction << ','
        << batch << ',' << even_batch_samples;
    for (EvenInt value:sums.sum) out << ',' << value;
    out << '\n';
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc!=3) throw std::runtime_error("usage: p40-even-replay 65|85 output.csv");
        const int n=std::stoi(argv[1]);
        if (n!=65 && n!=85) throw std::runtime_error("only frozen P40 N65/N85 are supported");
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite a replay file");
        const Geometry first = n==65 ? make_geometry(8,1) : make_geometry(9,2);
        const Geometry second = n==65 ? make_geometry(7,4) : make_geometry(7,6);
        HomologyUnionFind fb(n,first.a,first.b), fw(n,first.a,first.b);
        HomologyUnionFind sb(n,second.a,second.b), sw(n,second.a,second.b);
        std::vector<std::uint8_t> black(n),white(n);
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create replay output");
        out << "n,a,b,orientation,batch,samples";
        for (const char* name:even_names) out << ",sum_" << name;
        out << '\n';
        for (int batch=0; batch<even_batches; ++batch) {
            EvenSums f,s;
            const std::uint64_t begin=even_offset+static_cast<std::uint64_t>(batch)*even_batch_samples;
            for (std::uint64_t counter=begin; counter<begin+even_batch_samples; ++counter) {
                for (int site=0; site<n; ++site) {
                    black[site]=counter_uniform(even_seed,n,counter,site)<even_p;
                    white[site]=!black[site];
                }
                f.add(first,black,white,fb,fw);
                s.add(second,black,white,sb,sw);
            }
            even_write(out,first,"first",batch,f);
            even_write(out,second,"second",batch,s);
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing replay output");
        std::cout << "Reobserved N=" << n << ": 100 aligned batches, 1000000 old counters, zero new samples\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
