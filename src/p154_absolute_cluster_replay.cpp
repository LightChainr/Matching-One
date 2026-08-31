// Deterministic reobservation of the existing P154 mixed-plane 20k blocks.
// No new counters, samples, probabilities, geometries, or random streams.
//
// Historical metadata claims freeze commit
// 0578105d92d3822cb48f5c421bd23ff339295cc6; the resolvable freeze commit is
// 05781051b76001f2b18560d7b0914f2481412584. Preserve that distinction instead
// of rewriting old metadata. The shared backend Git blob is the source anchor:
// 22058703c12b168e844088277c9b61d64b9c1d2c
// (SHA256 7df0d9362b31111eab8fc73ab4032eea37458fde0ff2d720bded8e7b530fa94a).
//
// First moments and the upper triangle of all second moments are retained for
// (q,e,k,kk,edges,s,chi,cb,cw). Here s=cb+cw is unnormalised, q=I2-I0,
// e=I2+I0, and chi is computed independently as K - NN edges + black faces.
// This intentionally redundant vector must not be blindly covariance-inverted.
// Normalisations, centred responses, controls, and paired delete-one estimates
// belong to the downstream analysis, not this integer sufficient-statistic file.

#define main threshold_rank_integer_period_hidden_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main

#include <set>

namespace {

constexpr double absolute_p_ref = 0.59274605079;
constexpr std::size_t absolute_dimension = 9;
constexpr std::array<const char*, absolute_dimension> absolute_keys = {
    "q", "e", "k", "kk", "edges", "s", "chi", "cb", "cw"
};

int absolute_replay_count(int n, std::uint64_t seed, std::uint64_t replica) {
    // Same independent count stream as p154_local_singlet_pilot.cpp and the
    // fixed-K interaction replay. Its count is Bin(N,p), not a fixed K.
    SplitMixStream stream(splitmix64(seed ^ splitmix64(
        replica + 0x8cb92ba72f3d8dd7ULL)));
    int k = 0;
    for (int i = 0; i < n; ++i) {
        k += static_cast<double>(stream.next() >> 11) / 9007199254740992.0
             < absolute_p_ref;
    }
    return k;
}

void absolute_require_simple_nn(const Geometry& geometry) {
    std::set<std::pair<int, int>> seen;
    for (const Edge& edge : geometry.primal_edges) {
        if (edge.i == edge.j ||
            !seen.insert(std::minmax(edge.i, edge.j)).second) {
            throw std::runtime_error("archived NN graph must be simple");
        }
    }
    if (seen.size() != 2U * geometry.n) {
        throw std::runtime_error("expected 2N undirected NN edges");
    }
}

struct AbsoluteCounts {
    Int black = 0;
    Int white = 0;
    Int edges = 0;
    Int faces = 0;
    Int chi = 0;
};

AbsoluteCounts absolute_counts(const Geometry& geometry,
                               const std::vector<int>& permutation, int k) {
    std::vector<std::uint8_t> occupied(geometry.n, 0);
    for (int i = 0; i < k; ++i) occupied[permutation[i]] = 1;

    HomologyUnionFind black(geometry.quotient), white(geometry.quotient);
    AbsoluteCounts result;
    for (const Edge& edge : geometry.primal_edges) {
        if (occupied[edge.i] && occupied[edge.j]) {
            ++result.edges;
            black.add_edge(edge);
        }
    }
    // The white matching graph is NN plus both diagonal directions. Diagonals
    // are used for white connectivity, not as cells in the black Euler count.
    for (const Edge& edge : geometry.matching_edges) {
        if (!occupied[edge.i] && !occupied[edge.j]) white.add_edge(edge);
    }

    std::vector<std::uint8_t> black_seen(geometry.n, 0);
    std::vector<std::uint8_t> white_seen(geometry.n, 0);
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        // Inactive UF singletons are not clusters of the opposite colour.
        if (occupied[vertex]) {
            const int root = black.find(vertex).root;
            if (!black_seen[root]) {
                black_seen[root] = 1;
                ++result.black;
            }
        } else {
            const int root = white.find(vertex).root;
            if (!white_seen[root]) {
                white_seen[root] = 1;
                ++result.white;
            }
        }

        // Each quotient vertex indexes one translated fundamental square.
        const Vector z = geometry.quotient.representative(vertex);
        const int right = geometry.quotient.label({z.x + 1, z.y});
        const int up = geometry.quotient.label({z.x, z.y + 1});
        const int diagonal = geometry.quotient.label({z.x + 1, z.y + 1});
        result.faces += occupied[vertex] && occupied[right] &&
                        occupied[up] && occupied[diagonal];
    }
    result.chi = static_cast<Int>(k) - result.edges + result.faces;
    return result;
}

struct AbsoluteSums {
    Int samples = 0, k1 = 0, k2 = 0, i0 = 0, i1 = 0, i2 = 0;
    Int k = 0, kk = 0, edges = 0;
    Int i0k = 0, i2k = 0, i0kk = 0, i2kk = 0;
    Int i0edges = 0, i2edges = 0;
    Int faces = 0, i0s = 0, i2s = 0;
    std::array<Int, absolute_dimension> first{};
    std::array<std::array<Int, absolute_dimension>, absolute_dimension> second{};

    void add(int first_rank, int second_rank, int count,
             const AbsoluteCounts& observation, int n,
             std::uint64_t replica, const char* orientation) {
        const Int z0 = count < first_rank;
        const Int z2 = count >= second_rank;
        const Int q = z2 - z0;
        const Int e = z2 + z0;
        const Int factorial = static_cast<Int>(count) * (count - 1);
        const Int s = observation.black + observation.white;

        // Independent cell count versus cluster count/topological identity.
        // Never substitute the right-hand side as the definition of chi.
        if (observation.chi != observation.black - observation.white - q) {
            std::ostringstream message;
            message << "Euler/cluster identity failed: N=" << n
                    << " replica=" << replica << " orientation=" << orientation
                    << " cell_chi=" << observation.chi
                    << " black=" << observation.black
                    << " white=" << observation.white << " q=" << q;
            throw std::runtime_error(message.str());
        }

        ++samples;
        k1 += first_rank;
        k2 += second_rank;
        i0 += z0;
        i1 += 1 - z0 - z2;
        i2 += z2;
        k += count;
        kk += factorial;
        edges += observation.edges;
        i0k += z0 * count;
        i2k += z2 * count;
        i0kk += z0 * factorial;
        i2kk += z2 * factorial;
        i0edges += z0 * observation.edges;
        i2edges += z2 * observation.edges;
        faces += observation.faces;
        i0s += z0 * s;
        i2s += z2 * s;

        const std::array<Int, absolute_dimension> values = {
            q, e, static_cast<Int>(count), factorial, observation.edges,
            s, observation.chi, observation.black, observation.white
        };
        // All products fit signed int64 on the fixed N<=130, 200-sample batch.
        for (std::size_t i = 0; i < absolute_dimension; ++i) {
            first[i] += values[i];
            for (std::size_t j = i; j < absolute_dimension; ++j) {
                second[i][j] += values[i] * values[j];
            }
        }
    }
};

void absolute_header(std::ostream& out) {
    // Preserve every old fixed-K replay field, in its original order.
    out << "n,a,b,orientation,batch,samples,sum_k1,sum_k2,sum_i0,sum_i1,sum_i2,"
           "sum_k,sum_kk,sum_edges,sum_i0k,sum_i2k,sum_i0kk,sum_i2kk,"
           "sum_i0edges,sum_i2edges";
    for (std::size_t i : {0U, 1U, 5U, 6U, 7U, 8U}) {
        out << ",sum_" << absolute_keys[i];
    }
    out << ",sum_faces,sum_i0s,sum_i2s";
    for (std::size_t i = 0; i < absolute_dimension; ++i) {
        for (std::size_t j = i; j < absolute_dimension; ++j) {
            out << ",sum_" << absolute_keys[i] << '_' << absolute_keys[j];
        }
    }
    out << '\n';
}

void absolute_write(std::ostream& out, int n, int a, int b,
                    const char* orientation, int batch, const AbsoluteSums& s) {
    out << n << ',' << a << ',' << b << ',' << orientation << ',' << batch << ','
        << s.samples << ',' << s.k1 << ',' << s.k2 << ',' << s.i0 << ',' << s.i1
        << ',' << s.i2 << ',' << s.k << ',' << s.kk << ',' << s.edges << ','
        << s.i0k << ',' << s.i2k << ',' << s.i0kk << ',' << s.i2kk << ','
        << s.i0edges << ',' << s.i2edges;
    for (std::size_t i : {0U, 1U, 5U, 6U, 7U, 8U}) out << ',' << s.first[i];
    out << ',' << s.faces << ',' << s.i0s << ',' << s.i2s;
    for (std::size_t i = 0; i < absolute_dimension; ++i) {
        for (std::size_t j = i; j < absolute_dimension; ++j) {
            out << ',' << s.second[i][j];
        }
    }
    out << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error("usage: absolute-cluster-replay 65|130 output.csv");
        }
        const int n = std::stoi(argv[1]);
        if (n != 65 && n != 130) {
            throw std::runtime_error("only archived N65/N130 blocks are supported");
        }
        const bool small = n == 65;
        const Matrix first_matrix = small ? Matrix{8, -1, 1, 8}
                                          : Matrix{11, -3, 3, 11};
        const Matrix second_matrix = small ? Matrix{7, -4, 4, 7}
                                           : Matrix{9, -7, 7, 9};
        const Geometry first_geometry = make_geometry(first_matrix);
        const Geometry second_geometry = make_geometry(second_matrix);
        absolute_require_simple_nn(first_geometry);
        absolute_require_simple_nn(second_geometry);
        const std::uint64_t seed = small ? 202615465ULL : 2026154130ULL;
        const std::uint64_t offset = small ? 15466000000ULL : 15466200000ULL;
        if (std::filesystem::exists(argv[2])) {
            throw std::runtime_error("refusing to overwrite a replay artifact");
        }
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot create output");
        absolute_header(out);

        const auto start = std::chrono::steady_clock::now();
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        std::vector<int> permutation;
        for (int batch = 0; batch < 100; ++batch) {
            AbsoluteSums first, second;
            for (int j = 0; j < 200; ++j) {
                const std::uint64_t replica = offset + 200ULL * batch + j;
                counter_permutation(n, seed, replica, permutation);
                const int k = absolute_replay_count(n, seed, replica);
                const auto first_ranks = first_engine.ranks(permutation);
                const auto second_ranks = second_engine.ranks(permutation);
                first.add(first_ranks.first, first_ranks.second, k,
                          absolute_counts(first_geometry, permutation, k),
                          n, replica, "first");
                second.add(second_ranks.first, second_ranks.second, k,
                           absolute_counts(second_geometry, permutation, k),
                           n, replica, "second");
            }
            absolute_write(out, n, small ? 8 : 11, small ? 1 : 3, "first", batch, first);
            absolute_write(out, n, small ? 7 : 9, small ? 4 : 7, "second", batch, second);
        }
        out.close();
        if (!out) throw std::runtime_error("failed to finish replay output");
        std::cout << "N=" << n << " archived_samples=20000 batches=100"
                  << " paired_geometry_configurations=40000"
                  << " new_random_samples=0 cell_identity_failures=0 elapsed_seconds="
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start).count()
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
