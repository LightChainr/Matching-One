// New conditional-line marks on the SAME archived norm-4 permutations.
// The driver supplies the immutable bfab0330 production backend, including RNG,
// quotient labels, physical edge lifts and homology UF. No new counter is used.
#ifndef MATCHING_NORM4_BACKEND
#error "Supply the immutable norm-4 production backend through the replay driver"
#endif
#define main norm4_line_archived_unused_production_main
#include MATCHING_NORM4_BACKEND
#undef main

namespace {
constexpr int line_replay_batches = 100;

struct LineReplayCell {
    std::uint64_t samples = 0;
    std::int64_t rank1 = 0;
    std::int64_t rank1_s = 0;
    long double line4_re = 0;
    long double line4_im = 0;
    long double line4_s_re = 0;
    long double line4_s_im = 0;
};

// All cycles are expressed in the physical square-lattice x/y frame. Global
// ambient span is needed: this is not the rank of a selected UF component.
// Rank-one cycles are collinear, so the first nonzero lift determines their
// unoriented line. Sign and integer multiples cancel in the fourth character.
struct PhysicalCycleLine {
    int rank = 0;
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::array<long double, 2> character{{0, 0}};

    void add(std::int64_t dx, std::int64_t dy) {
        if (rank == 2 || (dx == 0 && dy == 0)) return;
        if (rank == 1) {
            if (static_cast<__int128>(x) * dy != static_cast<__int128>(y) * dx) {
                rank = 2;
                character = {{0, 0}};
            }
            return;
        }
        rank = 1;
        x = dx;
        y = dy;
        const long double a = static_cast<long double>(x);
        const long double b = static_cast<long double>(y);
        const long double xx = a * a;
        const long double yy = b * b;
        const long double denominator = (xx + yy) * (xx + yy);
        character = {{(xx * xx - 6 * xx * yy + yy * yy) / denominator,
                      4 * a * b * (xx - yy) / denominator}};
    }
};

class FullLineFiltrationReplay {
  public:
    explicit FullLineFiltrationReplay(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n),
#ifdef MATCHING_NORM4_INTEGER
          uf_(geometry.quotient),
#else
          uf_(geometry.n, geometry.a, geometry.b),
#endif
          black_components_(geometry.n + 1), white_components_(geometry.n + 1),
          black_rank_(geometry.n + 1), black_line_(geometry.n + 1) {}

    void observe(const std::vector<int>& permutation, std::vector<LineReplayCell>& cells) {
        const int k_plus = sweep(permutation, false, false, black_components_);
        const int reverse_white = sweep(permutation, true, true, white_components_);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (!(1 <= k_minus && k_minus <= k_plus && k_plus <= geometry_.n)) {
            throw std::logic_error("old-engine activation thresholds are inconsistent");
        }
        for (int k = 0; k <= geometry_.n; ++k) {
            // In-run semantic identity, using the two sweeps already required
            // for s. This does not generate a second permutation or test run.
            const int expected_rank = (k >= k_minus) + (k >= k_plus);
            if (black_rank_[k] != expected_rank) {
                throw std::logic_error("physical ambient rank disagrees with archived thresholds at N=" +
                    std::to_string(geometry_.n) + ", K=" + std::to_string(k));
            }
            LineReplayCell& cell = cells[k];
            ++cell.samples;
            if (black_rank_[k] != 1) continue;
            const std::int64_t s = black_components_[k] + white_components_[geometry_.n - k];
            ++cell.rank1;
            cell.rank1_s += s;
            cell.line4_re += black_line_[k][0];
            cell.line4_im += black_line_[k][1];
            cell.line4_s_re += static_cast<long double>(s) * black_line_[k][0];
            cell.line4_s_im += static_cast<long double>(s) * black_line_[k][1];
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
        PhysicalCycleLine ambient;
        if (!matching) {
            black_rank_[0] = 0;
            black_line_[0] = {{0, 0}};
        }
        const auto& edges = matching ? geometry_.matching_edges : geometry_.primal_edges;
        const auto& incident = matching ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            ++occupied_components;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (!active_[edge.i] || !active_[edge.j]) continue;
                const auto first = uf_.find(edge.i);
                const auto second = uf_.find(edge.j);
                if (first.root != second.root) {
                    --occupied_components;
                } else if (!matching) {
                    // UF offsets and edge.dx/dy are PHYSICAL lifts. Do not
                    // reinterpret the backend's period-basis winding as x/y.
                    ambient.add(first.dx + edge.dx - second.dx,
                                first.dy + edge.dy - second.dy);
                }
                uf_.add_edge(edge);
            }
            const int k = offset + 1;
            components[k] = occupied_components;
            if (!matching) {
                black_rank_[k] = static_cast<std::uint8_t>(ambient.rank);
                black_line_[k] = ambient.character;
            }
            if (first_cross == 0 && uf_.component_crosses(vertex)) first_cross = k;
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
    std::vector<std::uint8_t> black_rank_;
    std::vector<std::array<long double, 2>> black_line_;
};

void write_line_filtration(std::ostream& out, int n, int a, int b, const char* direction,
                           int batch, const std::vector<LineReplayCell>& cells) {
    for (int k = 0; k <= n; ++k) {
        const LineReplayCell& cell = cells[k];
        out << n << ',' << a << ',' << b << ',' << direction << ',' << batch << ',' << k
            << ',' << cell.samples << ',' << cell.rank1 << ',' << cell.rank1_s
            << ',' << cell.line4_re << ',' << cell.line4_im
            << ',' << cell.line4_s_re << ',' << cell.line4_s_im << '\n';
    }
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("usage: norm4-source-line-replay N output.csv");
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
        const std::uint64_t samples = endpoint ? 1000000ULL : 100000ULL;
        if (std::filesystem::exists(argv[2])) throw std::runtime_error("refusing to overwrite a line replay artifact");
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create line replay output");
        out << std::setprecision(21);
        out << "n,a,b,orientation,batch,k,samples,sum_rank1,sum_rank1_s,"
               "sum_line4_re,sum_line4_im,sum_line4_s_re,sum_line4_s_im\n";
        FullLineFiltrationReplay first_replay(first), second_replay(second);
        std::vector<int> permutation;
        for (int batch = 0; batch < line_replay_batches; ++batch) {
            std::vector<LineReplayCell> first_cells(n + 1), second_cells(n + 1);
            auto reobserve_interval = [&](std::uint64_t begin, std::uint64_t count) {
                for (std::uint64_t counter = begin; counter < begin + count; ++counter) {
                    counter_permutation(n, seed, counter, permutation);
                    first_replay.observe(permutation, first_cells);
                    second_replay.observe(permutation, second_cells);
                }
            };
            reobserve_interval(counter_begin + 1000ULL * batch, 1000);
            if (endpoint) {
                // Preserve the prior endpoint's exact nested analysis batches:
                // old first100k batch b UNION added900k batch b, not new 10k chunks.
                reobserve_interval(counter_begin + 100000ULL + 9000ULL * batch, 9000);
            }
            write_line_filtration(out, n, design.a1, design.b1, "first", batch, first_cells);
            write_line_filtration(out, n, design.a2, design.b2, "second", batch, second_cells);
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing line replay output");
        std::cout << "N=" << n << ": new physical rank1-line marks on " << samples
                  << " old permutations; 100 aligned batches; every K=0..N; seed=" << seed
                  << "; counters=[" << counter_begin << ',' << counter_begin + samples
                  << "); endpoint nested union=" << (endpoint ? "true" : "false")
                  << "; s=CB+CW bulk; line4=physical (dx+i*dy)^4/|v|^4; new samples=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
