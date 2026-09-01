// Reobserve OLD norm-4 production permutations at every K. The default remains
// the original first 100,000; an explicit endpoint increment cannot repeat them.
// The driver includes the immutable bfab0330 backend, including its original
// quotient labels, homology arithmetic and counter-keyed Fisher--Yates RNG.
// No Bernoulli draw, probability point, or new permutation counter is added.
#ifndef MATCHING_NORM4_BACKEND
#error "Supply an immutable norm-4 production backend through the replay driver"
#endif
#define main norm4_archived_unused_production_main
#include MATCHING_NORM4_BACKEND
#undef main

namespace {
using ReplayInt = std::int64_t;
using ReplayCell = std::array<ReplayInt, 6>;
constexpr int replay_batches = 100;

// This is the old activation loop continued after the first rank-2 crossing.
// Inactive UF singleton entries never count as occupied components.  A new
// active vertex creates one component; an edge joining two roots removes one.
class FullFiltrationReplay {
  public:
    explicit FullFiltrationReplay(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n),
#ifdef MATCHING_NORM4_INTEGER
          uf_(geometry.quotient),
#else
          uf_(geometry.n, geometry.a, geometry.b),
#endif
          black_components_(geometry.n + 1), white_components_(geometry.n + 1) {}

    void observe(const std::vector<int>& permutation, std::vector<ReplayCell>& cells) {
        const int k_plus = sweep(permutation, false, false, black_components_);
        const int reverse_white = sweep(permutation, true, true, white_components_);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (!(1 <= k_minus && k_minus <= k_plus && k_plus <= geometry_.n)) {
            throw std::logic_error("old-engine activation thresholds are inconsistent");
        }
        for (int k = 0; k <= geometry_.n; ++k) {
            // Digital Alexander / the old engine's exact two-activation object.
            const ReplayInt q = -1 + (k >= k_minus) + (k >= k_plus);
            const ReplayInt e = q * q;
            const ReplayInt s = black_components_[k] + white_components_[geometry_.n - k];
            const ReplayCell values{1, q, e, s, q * s, e * s};
            for (std::size_t field = 0; field < values.size(); ++field) {
                cells[k][field] += values[field];
            }
        }
    }

  private:
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
            throw std::logic_error("fully occupied production graph is not connected and rank two");
        }
        return first_cross;
    }

    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    HomologyUnionFind uf_;
    std::vector<int> black_components_;
    std::vector<int> white_components_;
};

void write_filtration(std::ostream& out, int n, int a, int b, const char* direction,
                      int batch, const std::vector<ReplayCell>& cells) {
    for (int k = 0; k <= n; ++k) {
        out << n << ',' << a << ',' << b << ',' << direction << ',' << batch << ',' << k;
        for (const ReplayInt value : cells[k]) out << ',' << value;
        out << '\n';
    }
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 5) {
            throw std::invalid_argument("usage: norm4-source-thermal-replay N output.csv [old_offset samples]");
        }
        const int n = std::stoi(argv[1]);
        const std::uint64_t offset = argc == 5 ? std::stoull(argv[3]) : 0;
        const std::uint64_t replay_samples = argc == 5 ? std::stoull(argv[4]) : 100000;
        const std::uint64_t archive_size = n == 260 || n == 340 ? 1000000000ULL : 1900000000ULL;
        if (replay_samples == 0 || replay_samples % replay_batches != 0 ||
            offset > archive_size || replay_samples > archive_size - offset ||
            (argc == 5 && (offset < 100000 || (n != 260 && n != 340)))) {
            throw std::invalid_argument("increment must partition an unmarked old endpoint interval into 100 equal batches");
        }
        const std::uint64_t replay_batch_samples = replay_samples / replay_batches;
#ifdef MATCHING_NORM4_INTEGER
        if (n != 260 && n != 340) throw std::invalid_argument("integer-period variant supports N260/N340 only");
#else
        if (n != 65 && n != 85 && n != 130 && n != 170) {
            throw std::invalid_argument("primitive variant supports N65/N85/N130/N170 only");
        }
#endif
        const auto found = std::find_if(kDesigns.begin(), kDesigns.end(),
            [n](const PairDesign& design) { return design.n == n; });
        if (found == kDesigns.end()) throw std::logic_error("frozen backend lacks the requested design");
        const PairDesign& design = *found;
#ifdef MATCHING_NORM4_INTEGER
        const Geometry first = make_geometry(design.first);
        const Geometry second = make_geometry(design.second);
#else
        const Geometry first = make_geometry(design.a1, design.b1);
        const Geometry second = make_geometry(design.a2, design.b2);
#endif
        const std::uint64_t seed = n == 260 ? 2026105401ULL :
                                   n == 340 ? 2026105402ULL : 2026104501ULL;
        const std::uint64_t counter_begin = (n == 260 || n == 340 ? 8200000000ULL : 5100000000ULL) + offset;
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite a replay artifact");
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create replay output");
        out << "n,a,b,orientation,batch,k,samples,sum_q,sum_e,sum_s,sum_qs,sum_es\n";
        FullFiltrationReplay first_replay(first), second_replay(second);
        std::vector<int> permutation;
        for (int batch = 0; batch < replay_batches; ++batch) {
            std::vector<ReplayCell> first_cells(n + 1), second_cells(n + 1);
            const std::uint64_t begin = counter_begin + static_cast<std::uint64_t>(batch) * replay_batch_samples;
            for (std::uint64_t counter = begin; counter < begin + replay_batch_samples; ++counter) {
                counter_permutation(n, seed, counter, permutation);
                first_replay.observe(permutation, first_cells);
                second_replay.observe(permutation, second_cells);
            }
            write_filtration(out, n, design.a1, design.b1, "first", batch, first_cells);
            write_filtration(out, n, design.a2, design.b2, "second", batch, second_cells);
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing replay output");
        std::cout << "N=" << n << ": reobserved " << replay_samples
                  << " old permutations, 100 aligned batches, every K=0..N; seed=" << seed
                  << "; counters=[" << counter_begin << ',' << counter_begin + replay_samples
                  << "); new samples=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
