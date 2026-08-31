// Add lag-one source/event marks to the SAME archived norm-4 permutations.
// The driver includes bfab0330's original RNG, quotient labels and homology UF.
// Only three transitions are marked. Source at terminal K is the BULK cluster
// count at K-1, not the cluster count after the event and not s/N.
#ifndef MATCHING_NORM4_BACKEND
#error "Supply the immutable norm-4 production backend through the replay driver"
#endif
#define main norm4_lagged_archived_unused_production_main
#include MATCHING_NORM4_BACKEND
#undef main

namespace {
constexpr int lagged_batches = 100;
using LaggedInt = std::int64_t;
// Counts 01,02,12 followed by corresponding sums of s_previous.
using LaggedCell = std::array<LaggedInt, 6>;

class LaggedFiltrationReplay {
  public:
    explicit LaggedFiltrationReplay(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n),
#ifdef MATCHING_NORM4_INTEGER
          uf_(geometry.quotient),
#else
          uf_(geometry.n, geometry.a, geometry.b),
#endif
          black_components_(geometry.n + 1), white_components_(geometry.n + 1) {}

    void observe(const std::vector<int>& permutation, std::vector<LaggedCell>& cells) {
        const int k_plus = sweep(permutation, false, false, black_components_);
        const int reverse_white = sweep(permutation, true, true, white_components_);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (!(1 <= k_minus && k_minus <= k_plus && k_plus <= geometry_.n)) {
            throw std::logic_error("old-engine activation thresholds are inconsistent");
        }
        const auto mark = [&](int k, std::size_t event) {
            // White reverse prefix N-(K-1) is exactly the complement of the
            // black forward prefix K-1 on THIS permutation.
            const LaggedInt s_previous = black_components_[k - 1]
                + white_components_[geometry_.n - k + 1];
            ++cells[k][event];
            cells[k][event + 3] += s_previous;
        };
        if (k_minus == k_plus) {
            mark(k_minus, 1);  // direct 0->2: Delta q=2, Delta E=0
        } else {
            mark(k_minus, 0);  // 0->1: Delta q=1, Delta E=-1
            mark(k_plus, 2);   // 1->2: Delta q=1, Delta E=+1
        }
    }

  private:
    // The old source replay's complete forward/reverse activation sweeps.
    // Continue after rank two so both source component arrays remain complete.
    int sweep(const std::vector<int>& permutation, bool matching, bool reverse,
              std::vector<int>& components) {
        std::fill(active_.begin(), active_.end(), 0);
        uf_.reset();
        components[0] = 0;
        int occupied_components = 0;
        int first_cross = 0;
        const auto& edges = matching ? geometry_.matching_edges : geometry_.primal_edges;
        const auto& incident = matching ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            ++occupied_components;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) {
                    if (uf_.find(edge.i).root != uf_.find(edge.j).root) --occupied_components;
                    uf_.add_edge(edge);
                }
            }
            components[offset + 1] = occupied_components;
            if (first_cross == 0 && uf_.component_crosses(vertex)) first_cross = offset + 1;
        }
        if (first_cross == 0 || occupied_components != 1) {
            throw std::logic_error("fully occupied archived graph is not connected and rank two");
        }
        return first_cross;
    }

    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    HomologyUnionFind uf_;
    std::vector<int> black_components_;
    std::vector<int> white_components_;
};

void write_lagged(std::ostream& out, int n, int a, int b, const char* orientation,
                  int batch, std::uint64_t samples, const std::vector<LaggedCell>& cells) {
    for (int k = 0; k <= n; ++k) {
        out << n << ',' << a << ',' << b << ',' << orientation << ',' << batch
            << ',' << k << ',' << samples;
        for (const LaggedInt value : cells[k]) out << ',' << value;
        out << '\n';
    }
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("usage: norm4-lagged-source-replay N output.csv");
        const int n = std::stoi(argv[1]);
#ifdef MATCHING_NORM4_INTEGER
        if (n != 260 && n != 340) throw std::invalid_argument("integer-period variant supports N260/N340 only");
#else
        if (n != 65 && n != 85 && n != 130 && n != 170) {
            throw std::invalid_argument("primitive variant supports N65/N85/N130/N170 only");
        }
#endif
        const auto found = std::find_if(kDesigns.begin(), kDesigns.end(),
            [n](const PairDesign& design) { return design.n == n; });
        if (found == kDesigns.end()) throw std::logic_error("frozen backend lacks requested design");
        const PairDesign& design = *found;
#ifdef MATCHING_NORM4_INTEGER
        // N260/N340 must retain original HNF quotient labels. a,b below are
        // only output lineage labels; they do not construct these geometries.
        const Geometry first = make_geometry(design.first);
        const Geometry second = make_geometry(design.second);
#else
        const Geometry first = make_geometry(design.a1, design.b1);
        const Geometry second = make_geometry(design.a2, design.b2);
#endif
        const bool endpoint = n == 260 || n == 340;
        const std::uint64_t seed = n == 260 ? 2026105401ULL :
                                   n == 340 ? 2026105402ULL : 2026104501ULL;
        const std::uint64_t counter_begin = endpoint ? 8200000000ULL : 5100000000ULL;
        const std::uint64_t samples_per_batch = endpoint ? 10000ULL : 1000ULL;
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite lagged marks");
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create lagged replay output");
        out << "n,a,b,orientation,batch,k,samples,event_count01,event_count02,event_count12,"
               "sum_s_previous01,sum_s_previous02,sum_s_previous12\n";
        LaggedFiltrationReplay first_replay(first), second_replay(second);
        std::vector<int> permutation;
        for (int batch = 0; batch < lagged_batches; ++batch) {
            std::vector<LaggedCell> first_cells(n + 1), second_cells(n + 1);
            std::uint64_t observed = 0;
            const auto reobserve_interval = [&](std::uint64_t begin, std::uint64_t count) {
                for (std::uint64_t counter = begin; counter < begin + count; ++counter) {
                    counter_permutation(n, seed, counter, permutation);
                    first_replay.observe(permutation, first_cells);
                    second_replay.observe(permutation, second_cells);
                    ++observed;
                }
            };
            reobserve_interval(counter_begin + 1000ULL * batch, 1000);
            if (endpoint) {
                // Same nested union as the saved endpoint source/line arrays:
                // old1000 of batch b plus incremental9000 of batch b.
                reobserve_interval(counter_begin + 100000ULL + 9000ULL * batch, 9000);
            }
            if (observed != samples_per_batch) throw std::logic_error("frozen batch count differs");
            write_lagged(out, n, design.a1, design.b1, "first", batch, observed, first_cells);
            write_lagged(out, n, design.a2, design.b2, "second", batch, observed, second_cells);
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing lagged replay output");
        std::cout << "N=" << n << ": lag=1 event/source marks on " << samples_per_batch * lagged_batches
                  << " old permutations; 100 aligned batches; seed=" << seed
                  << "; counters=[" << counter_begin << ',' << counter_begin + samples_per_batch * lagged_batches
                  << "); endpoint nested union=" << (endpoint ? "true" : "false")
                  << "; s_previous=CB[K-1]+CW[N-K+1] bulk; new samples=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
