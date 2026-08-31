// Reconstruct only the immutable checkpoint and already sampled next labels.
// No call to ranks(), a tail RNG, a child DP, or suffix replay occurs here.
#define main p334_unused_threshold_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main
#include <zlib.h>

namespace {
struct ContactMark {
    int e = 0, c = 0, cycles = 0;
    int contact_rank = -1, null_cycle_rank = -1, architecture = -1;
};

int image_rank(const std::vector<Vector>& generators) {
    Vector basis{0, 0};
    int rank = 0;
    for (const auto w: generators) {
        if (w.x == 0 && w.y == 0) continue;
        if (rank == 0) { basis = w; rank = 1; }
        else if (static_cast<__int128>(basis.x)*w.y != static_cast<__int128>(basis.y)*w.x)
            return 2;
    }
    return rank;
}

struct Checkpoint {
    const Geometry& geometry;
    HomologyUnionFind uf;
    std::vector<std::uint8_t> active;

    explicit Checkpoint(const Geometry& g): geometry(g), uf(g.quotient), active(g.n) {}

    void load(const std::vector<int>& permutation, int k0) {
        uf.reset();
        std::fill(active.begin(), active.end(), 0);
        for (int k = 0; k < k0; ++k) {
            const int vertex = permutation[k];
            active[vertex] = 1;
            for (const int edge_index: geometry.primal_incident[vertex]) {
                const auto& edge = geometry.primal_edges[edge_index];
                if (active[edge.i] && active[edge.j]) uf.add_edge(edge);
            }
        }
    }

    ContactMark read(int vertex, int old_rank, int recorded_after) {
        if (active[vertex]) throw std::runtime_error("sampled label is occupied at checkpoint");
        std::map<int, std::vector<Vector>> contacts;
        ContactMark result;
        for (const int edge_index: geometry.primal_incident[vertex]) {
            const auto& edge = geometry.primal_edges[edge_index];
            const int other = edge.i == vertex ? edge.j : edge.i;
            if (!active[other]) continue;
            ++result.e;
            const Vector step = edge.i == vertex ? Vector{edge.dx, edge.dy} : Vector{-edge.dx, -edge.dy};
            const auto found = uf.find(other);
            // alpha=delta-p_C(u); root gauge is never shared between components.
            contacts[found.root].push_back({step.x-found.dx, step.y-found.dy});
        }
        result.c = static_cast<int>(contacts.size());
        result.cycles = result.e-result.c;
        if (old_rank != 0) return result;
        std::vector<Vector> all;
        bool one_component_rank2 = false;
        for (const auto& item: contacts) {
            const auto anchor = item.second.front();
            std::vector<Vector> local;
            for (std::size_t j = 1; j < item.second.size(); ++j) {
                const auto alpha = item.second[j];
                const auto w = geometry.quotient.winding(alpha.x-anchor.x, alpha.y-anchor.y);
                local.push_back(w);
                all.push_back(w);
            }
            one_component_rank2 = one_component_rank2 || image_rank(local) == 2;
        }
        result.contact_rank = image_rank(all);
        result.null_cycle_rank = result.cycles-result.contact_rank;
        result.architecture = result.contact_rank == 2 ? (one_component_rank2 ? 1 : 2) : 0;
        if (result.contact_rank != recorded_after || result.null_cycle_rank < 0)
            throw std::runtime_error("contact theorem/prefix-source mismatch");
        if (result.architecture == 2 && (result.e != 4 || result.c != 2))
            throw std::runtime_error("independent two-component birth is not 2+2");
        return result;
    }
};

std::vector<std::uint64_t> parse_csv(const std::string& line) {
    std::stringstream stream(line);
    std::string cell;
    std::vector<std::uint64_t> values;
    while (std::getline(stream, cell, ',')) values.push_back(std::stoull(cell));
    return values;
}

std::map<std::uint64_t, std::array<int,4>> original_prefixes(int n) {
    std::ifstream input("results/p334-full-birth-archive/N"+std::to_string(n)+".csv");
    std::string line;
    std::getline(input, line);
    if (line != "N,batch,counter,k0,first_k1,first_k2,first_rank,second_k1,second_k2,second_rank")
        throw std::runtime_error("missing original full-prefix archive");
    std::map<std::uint64_t, std::array<int,4>> result;
    while (std::getline(input,line)) {
        const auto x = parse_csv(line);
        if (x.size() != 10 || x[0] != static_cast<unsigned>(n)) throw std::runtime_error("bad original prefix");
        result.emplace(x[2], std::array<int,4>{static_cast<int>(x[1]),static_cast<int>(x[3]),static_cast<int>(x[6]),static_cast<int>(x[9])});
    }
    if (result.size() != 20000) throw std::runtime_error("not the complete original20k population");
    return result;
}

struct Counts {
    std::uint64_t r0 = 0, safe = 0, rank1 = 0, single = 0, two_two = 0, safe_cycles = 0;
    std::array<std::uint64_t,4> safe_loop_merger{{0,0,0,0}};
    void add(const ContactMark& m, int old_rank) {
        if (old_rank != 0) return;
        ++r0;
        if (m.contact_rank == 0) {
            ++safe;
            safe_cycles += m.cycles;
            ++safe_loop_merger[(m.cycles > 0 ? 2 : 0)+(m.c > 1 ? 1 : 0)];
        } else if (m.contact_rank == 1) ++rank1;
        else if (m.architecture == 1) ++single;
        else ++two_two;
    }
};
}

#ifndef P334_CONTACT_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: contact_coordinates N output_directory code_commit");
        const int n = std::stoi(argv[1]);
        if (n != 325 && n != 425) throw std::runtime_error("only the original two sizes");
        const std::filesystem::path out(argv[2]);
        if (std::filesystem::exists(out)) throw std::runtime_error("refuse to overwrite contact output");
        const auto prefixes = original_prefixes(n);
        const auto first_geometry = make_geometry({n, n == 325 ? 57 : 132, 0, 1});
        const auto second_geometry = make_geometry({n, n == 325 ? 18 : 268, 0, 1});
        Checkpoint first(first_geometry), second(second_geometry);
        const std::uint64_t seed = n == 325 ? 20260831430325ULL : 20260831430425ULL;
        std::filesystem::create_directories(out);
        std::vector<int> permutation;
        std::array<Counts,2> counts;
        std::uint64_t output_rows = 0, reconstructed_prefixes = 0, unique_label_queries = 0;
        const auto start = std::chrono::steady_clock::now();
        for (int batch = 0; batch < 20; ++batch) {
            std::ostringstream name;
            name << "N" << n << ".batch" << std::setw(2) << std::setfill('0') << batch << ".csv.gz";
            const auto source = std::filesystem::path("results/p334-nested-next-label-forks")/("N"+std::to_string(n))/name.str();
            gzFile input = gzopen(source.c_str(), "rb"), output = gzopen((out/name.str()).c_str(), "wb1");
            if (!input || !output) throw std::runtime_error("cannot open immutable input/contact output");
            char buffer[4096];
            if (!gzgets(input,buffer,sizeof(buffer))) throw std::runtime_error("missing fork header");
            gzprintf(output,"N,batch,counter,quartet,group,next_label,first_oldrank,first_rank_after,first_e,first_c,first_new_cycles,first_contact_rank,first_r0_null_cycle_rank,first_r0_rank2_arch,second_oldrank,second_rank_after,second_e,second_c,second_new_cycles,second_contact_rank,second_r0_null_cycle_rank,second_r0_rank2_arch\n");
            std::uint64_t prior = std::numeric_limits<std::uint64_t>::max();
            std::map<int,std::array<ContactMark,2>> cache;
            while (gzgets(input,buffer,sizeof(buffer))) {
                const auto x = parse_csv(buffer);
                if (x.size() != 16) throw std::runtime_error("bad frozen fork schema");
                if (x[8] != 0) continue;  // no duplicate row for the second suffix.
                const auto counter = x[2];
                const int k0 = static_cast<int>(x[3]), old0 = static_cast<int>(x[4]), old1 = static_cast<int>(x[5]);
                if (prefixes.at(counter) != std::array<int,4>{batch,k0,old0,old1})
                    throw std::runtime_error("fork/original prefix keys differ");
                if (counter != prior) {
                    counter_permutation(n,seed,counter,permutation); // exact old prefix only.
                    first.load(permutation,k0); second.load(permutation,k0);
                    prior = counter; cache.clear(); ++reconstructed_prefixes;
                }
                const int vertex = static_cast<int>(x[9]);
                auto it = cache.find(vertex);
                if (it == cache.end()) {
                    it = cache.emplace(vertex,std::array<ContactMark,2>{first.read(vertex,old0,x[10]),second.read(vertex,old1,x[11])}).first;
                    ++unique_label_queries;
                }
                const auto& a = it->second[0]; const auto& b = it->second[1];
                counts[0].add(a,old0); counts[1].add(b,old1);
                gzprintf(output,"%d,%d,%llu,%llu,%llu,%d,%d,%llu,%d,%d,%d,%d,%d,%d,%d,%llu,%d,%d,%d,%d,%d,%d\n",
                    n,batch,static_cast<unsigned long long>(counter),static_cast<unsigned long long>(x[6]),static_cast<unsigned long long>(x[7]),vertex,
                    old0,static_cast<unsigned long long>(x[10]),a.e,a.c,a.cycles,a.contact_rank,a.null_cycle_rank,a.architecture,
                    old1,static_cast<unsigned long long>(x[11]),b.e,b.c,b.cycles,b.contact_rank,b.null_cycle_rank,b.architecture);
                ++output_rows;
            }
            if (gzclose(input) != Z_OK || gzclose(output) != Z_OK) throw std::runtime_error("gzip failed");
        }
        if (reconstructed_prefixes != 20000 || output_rows != 320000) throw std::runtime_error("incomplete original domain");
        std::ofstream meta(out/"metadata.json");
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        meta << std::setprecision(15) << "{\n\"N\":" << n << ",\n\"code_commit\":\"" << argv[3]
             << "\",\n\"fork_source_commit\":\"e32a85939279b8574278024d647b56d2d1485247\",\n"
             << "\"prefix_source_commit\":\"9c495ab13e65f2bc93dc0849ee3b73f88724c4b1\",\n"
             << "\"original_seed\":" << seed << ",\n\"reconstructed_prefixes\":" << reconstructed_prefixes
             << ",\n\"sampled_label_rows\":" << output_rows << ",\n\"unique_prefix_label_queries\":" << unique_label_queries
             << ",\n\"new_samples\":0,\n\"tail_replays\":0,\n\"DP_calls\":0,\n\"elapsed_seconds\":" << seconds << ",\n\"R0_counts\":[\n";
        for (int o = 0; o < 2; ++o) {
            const auto& c = counts[o];
            meta << "{\"orientation\":\"" << (o == 0 ? "first" : "second") << "\",\"R0_label_draws\":" << c.r0
                 << ",\"safe_R0_draws\":" << c.safe << ",\"rank1_draws\":" << c.rank1
                 << ",\"rank2_single_component\":" << c.single << ",\"rank2_two_components_2plus2\":" << c.two_two
                 << ",\"safe_contractible_cycle_rank_sum\":" << c.safe_cycles << ",\"safe_loop_merger_00_01_10_11\":["
                 << c.safe_loop_merger[0] << ',' << c.safe_loop_merger[1] << ',' << c.safe_loop_merger[2] << ',' << c.safe_loop_merger[3] << "]}" << (o == 0 ? ",\n" : "\n");
        }
        meta << "]}\n";
        meta.close();
        std::cout << "N" << n << " COMPLETE " << output_rows << " rows, " << unique_label_queries << " unique prefix-labels in " << seconds << " seconds\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
#endif
