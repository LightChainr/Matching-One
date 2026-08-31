// Fresh-prefix, finite +/- source-normal intervention. No archive is read.
// The original integer-period UF and K1/K2 engine are included unchanged.
#define main p334_unused_threshold_main
#include "threshold_rank_integer_period_mc.cpp"
#undef main
#include <zlib.h>

namespace {
using I128 = __int128;
using U128 = unsigned __int128;
constexpr int observable_count = 6, response_count = 24;
constexpr long double p_reference = 0.59274605079L;
constexpr std::uint64_t prefix_tag = 0xf3340001ULL;
constexpr std::uint64_t treatment_tag = 0xf3340002ULL;

struct NormalInterventionOptions {
    int n = 0, batch = -1, prefixes = 0, reps = 0;
    std::uint64_t seed = 0;
    std::string out, commit;
    bool smoke = false;
};

NormalInterventionOptions parse_normal_options(int argc, char** argv) {
    NormalInterventionOptions o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--smoke") { o.smoke = true; continue; }
        if (i+1 == argc) throw std::invalid_argument("missing value for "+key);
        const std::string value = argv[++i];
        if (key == "--n") o.n = std::stoi(value);
        else if (key == "--batch") o.batch = std::stoi(value);
        else if (key == "--prefixes") o.prefixes = std::stoi(value);
        else if (key == "--reps") o.reps = std::stoi(value);
        else if (key == "--seed") o.seed = std::stoull(value);
        else if (key == "--output") o.out = value;
        else if (key == "--code-commit") o.commit = value;
        else throw std::invalid_argument("unknown option "+key);
    }
    if ((o.n != 325 && o.n != 425) || o.batch < 0 || o.batch >= 20 ||
        o.prefixes <= 0 || o.prefixes > 25000 || o.reps <= 0 || o.reps > 8 ||
        !o.seed || o.out.empty() || o.commit.size() != 40 ||
        o.commit.find_first_not_of("0123456789abcdef") != std::string::npos)
        throw std::invalid_argument("usage: --n {325,425} --batch {0..19} --prefixes {1..25000} --reps {1..8} --seed UINT64 --output NEW_DIRECTORY --code-commit FULL_SHA [--smoke]");
    if (o.smoke && (o.prefixes > 64 || o.seed == 202608311920334ULL))
        throw std::invalid_argument("technical smoke is capped at64 prefixes and excludes the frozen production seed");
    return o;
}

struct Contact { int e = 0, c = 0, rank = 0; };

void span_add(Vector w, Vector& basis, int& rank) {
    if (rank == 2 || (w.x == 0 && w.y == 0)) return;
    if (rank == 0) { basis = w; rank = 1; }
    else if (static_cast<I128>(basis.x)*w.y != static_cast<I128>(basis.y)*w.x) rank = 2;
}

struct FreshCheckpoint {
    const Geometry& geometry;
    HomologyUnionFind uf;
    std::vector<std::uint8_t> occupied;
    explicit FreshCheckpoint(const Geometry& g): geometry(g), uf(g.quotient), occupied(g.n) {}

    int load(const std::vector<int>& permutation, int k0) {
        uf.reset();
        std::fill(occupied.begin(), occupied.end(), 0);
        for (int k = 0; k < k0; ++k) {
            const int vertex = permutation[k];
            occupied[vertex] = 1;
            for (int ei: geometry.primal_incident[vertex]) {
                const auto& edge = geometry.primal_edges[ei];
                if (occupied[edge.i] && occupied[edge.j]) uf.add_edge(edge);
            }
        }
        Vector basis{};
        int rank = 0;
        for (int k = 0; k < k0 && rank < 2; ++k) {
            const auto mark = uf.component_mark(permutation[k]);
            if (mark.rank == 2) rank = 2;
            else if (mark.rank == 1) span_add(mark.line, basis, rank);
        }
        return rank;
    }

    // Only called for old rank0. Differences use anchors within the SAME UF root.
    Contact contact(int vertex) {
        std::array<int,4> roots{};
        std::array<Vector,4> anchors{};
        Vector winding_basis{};
        Contact r;
        for (int ei: geometry.primal_incident[vertex]) {
            const auto& edge = geometry.primal_edges[ei];
            const int other = edge.i == vertex ? edge.j : edge.i;
            if (!occupied[other]) continue;
            ++r.e;
            const Vector step = edge.i == vertex ? Vector{edge.dx,edge.dy} : Vector{-edge.dx,-edge.dy};
            const auto found = uf.find(other);
            const Vector alpha{step.x-found.dx,step.y-found.dy};
            int j = 0;
            while (j < r.c && roots[j] != found.root) ++j;
            if (j == r.c) {
                if (r.c == 4) throw std::logic_error("non-square contact degree");
                roots[j] = found.root; anchors[j] = alpha; ++r.c;
            } else {
                span_add(geometry.quotient.winding(alpha.x-anchors[j].x,alpha.y-anchors[j].y),
                         winding_basis,r.rank);
            }
        }
        return r;
    }
};

struct ClassMoments {
    I128 count = 0;
    std::array<I128,2> sum{};
    std::array<I128,2> square{};
};

struct NormalLaw {
    int rank = 0;
    std::array<std::vector<I128>,2> numerator;
    std::array<I128,2> maximum{};
    I128 denominator = 1;
};

NormalLaw source_normal_law(FreshCheckpoint& first, FreshCheckpoint& second,
                           const std::vector<int>& vacant) {
    const int d = static_cast<int>(vacant.size());
    std::array<ClassMoments,25> classes{};
    std::vector<int> cls(d,-1);
    std::vector<std::array<int,2>> loops(d);
    for (int u = 0; u < d; ++u) {
        const auto a = first.contact(vacant[u]), b = second.contact(vacant[u]);
        if (a.rank || b.rank) continue;
        cls[u] = 5*a.e+b.e;
        loops[u] = {a.e-a.c,b.e-b.c};
        auto& c = classes[cls[u]];
        ++c.count;
        for (int o = 0; o < 2; ++o) {
            c.sum[o] += loops[u][o];
            c.square[o] += loops[u][o]*loops[u][o];
        }
    }
    std::vector<std::array<I128,2>> v(d), t(d);
    I128 gff = 0, gfs = 0, gss = 0;
    std::array<std::array<I128,2>,2> m{};
    for (int u = 0; u < d; ++u) {
        if (cls[u] < 0) continue;
        const auto& c = classes[cls[u]];
        for (int o = 0; o < 2; ++o) {
            v[u][o] = c.count*loops[u][o]-c.sum[o];
            t[u][o] = v[u][o]*v[u][o]-(c.count*c.square[o]-c.sum[o]*c.sum[o]);
        }
        gff += v[u][0]*v[u][0];
        gfs += v[u][0]*v[u][1];
        gss += v[u][1]*v[u][1];
        for (int o = 0; o < 2; ++o) for (int k = 0; k < 2; ++k)
            m[o][k] += t[u][o]*v[u][k];
    }
    const I128 det = gff*gss-gfs*gfs, trace = gff+gss;
    if (det < 0) throw std::logic_error("negative exact source Gram determinant");
    NormalLaw result;
    result.rank = det > 0 ? 2 : (trace > 0 ? 1 : 0);
    // rank1 Moore-Penrose: Gnum^+ = Gnum / trace(Gnum)^2.
    const I128 divisor = result.rank == 2 ? det : (result.rank == 1 ? trace*trace : 1);
    result.denominator = static_cast<I128>(d)*d*divisor;
    for (int o = 0; o < 2; ++o) {
        I128 af = 0, as = 0;
        if (result.rank == 2) {
            af = m[o][0]*gss-m[o][1]*gfs;
            as = m[o][1]*gff-m[o][0]*gfs;
        } else if (result.rank == 1) {
            af = m[o][0]*gff+m[o][1]*gfs;
            as = m[o][0]*gfs+m[o][1]*gss;
        }
        result.numerator[o].resize(d);
        std::array<I128,25> class_sum{};
        std::array<I128,2> orthogonal{};
        for (int u = 0; u < d; ++u) {
            const I128 z = divisor*t[u][o]-af*v[u][0]-as*v[u][1];
            result.numerator[o][u] = z;
            result.maximum[o] = std::max(result.maximum[o], z < 0 ? -z : z);
            if (cls[u] >= 0) class_sum[cls[u]] += z;
            for (int k = 0; k < 2; ++k) orthogonal[k] += z*v[u][k];
        }
        if (orthogonal[0] || orthogonal[1] ||
            std::any_of(class_sum.begin(),class_sum.end(),[](I128 z){return z != 0;}))
            throw std::logic_error("exact source projection identities failed");
    }
    return result;
}

U128 uniform_below128(SplitMixStream& rng, U128 bound) {
    if (bound == 0) throw std::logic_error("zero rational sampling denominator");
    if (bound == 1) return 0;
    int bits = 0;
    for (U128 x = bound-1; x; x >>= 1) ++bits;
    const U128 mask = bits == 128 ? ~U128(0) : (U128(1)<<bits)-1;
    while (true) {
        const std::uint64_t high = rng.next(), low = rng.next();
        const U128 value = ((U128(high)<<64)|low) & mask;
        if (value < bound) return value;
    }
}

std::array<int,2> coupled_labels(const NormalLaw& law, int axis,
                                const std::vector<int>& vacant, SplitMixStream& rng) {
    const U128 bound = U128(law.maximum[axis])*vacant.size();
    const U128 u = uniform_below128(rng,bound);
    std::array<U128,2> cumulative{};
    std::array<int,2> selected{-1,-1};
    for (std::size_t j = 0; j < vacant.size(); ++j) {
        for (int arm = 0; arm < 2; ++arm) {
            const I128 weight = law.maximum[axis]+(arm == 0 ? 1 : -1)*law.numerator[axis][j];
            cumulative[arm] += U128(weight);
            if (selected[arm] < 0 && u < cumulative[arm]) selected[arm] = vacant[j];
        }
    }
    if (selected[0] < 0 || selected[1] < 0 || cumulative[0] != bound || cumulative[1] != bound)
        throw std::logic_error("exact rational label sampler failed");
    return selected;
}

std::uint64_t fresh_domain(std::uint64_t seed, int n, std::uint64_t tag) {
    return splitmix64(seed ^ splitmix64((std::uint64_t(n)<<32)|tag));
}

std::uint64_t treatment_key(std::uint64_t domain, std::uint64_t counter,
                            int rep, int axis, int slot) {
    // Fixed8 slots per prefix even in smoke; changing --reps never aliases a later prefix.
    const std::uint64_t id = (((counter*8+rep)*2+axis)*2+slot);
    return splitmix64(domain ^ splitmix64(id));
}

void shared_remainder(const std::vector<int>& vacant, SplitMixStream& rng,
                      std::vector<int>& priorities) {
    priorities = vacant;
    for (int j = static_cast<int>(priorities.size())-1; j > 0; --j)
        std::swap(priorities[j],priorities[rng.below(j+1)]);
}

void arm_permutation(const std::vector<int>& prefix, int k0, int selected,
                     const std::vector<int>& priorities, std::vector<int>& out) {
    out.assign(prefix.begin(),prefix.begin()+k0);
    out.push_back(selected);
    for (int u: priorities) if (u != selected) out.push_back(u);
}

std::vector<long double> canonical_tails(int n) {
    std::vector<long double> mass(n+1), tail(n+2,0);
    long double total = 0;
    for (int k = 0; k <= n; ++k) {
        mass[k] = std::exp(std::lgamma(static_cast<long double>(n+1))-
            std::lgamma(static_cast<long double>(k+1))-
            std::lgamma(static_cast<long double>(n-k+1))+
            k*std::log(p_reference)+(n-k)*std::log1p(-p_reference));
        total += mass[k];
    }
    for (int k = n; k >= 0; --k) tail[k] = tail[k+1]+mass[k]/total;
    return tail;
}

std::array<long double,observable_count> readout(std::pair<int,int> births,
                                              int n, const std::vector<long double>& tail) {
    const long double k1 = births.first, k2 = births.second;
    return {tail[births.first]+tail[births.second]-1,
            1-tail[births.first]+tail[births.second],
            (k1+k2)/(2*(n+1)),(k2-k1)/(n+1),k1,k2};
}

std::array<std::string,response_count> response_labels() {
    std::array<std::string,response_count> labels;
    const std::array<std::string,2> role{"first","second"};
    const std::array<std::string,observable_count> obs{"A_ref","E_ref","C","W","K1","K2"};
    for (int s = 0; s < 2; ++s) for (int o = 0; o < 2; ++o) for (int f = 0; f < observable_count; ++f)
        labels[(s*2+o)*observable_count+f] = "source_"+role[s]+"__receiver_"+role[o]+"__"+obs[f];
    return labels;
}

template <typename Range> void json_numbers(std::ostream& out, const Range& values) {
    out << '[';
    bool comma = false;
    for (const auto value: values) { if (comma) out << ','; out << value; comma = true; }
    out << ']';
}
}

int main(int argc, char** argv) {
    try {
        const auto options = parse_normal_options(argc,argv);
        const int n = options.n, k0 = n == 325 ? 193 : 252;
        const std::filesystem::path out(options.out);
        if (std::filesystem::exists(out)) throw std::runtime_error("refuse to overwrite output directory");
        std::filesystem::create_directories(out);
        const auto first_geometry = make_geometry({n,n == 325 ? 57 : 132,0,1});
        const auto second_geometry = make_geometry({n,n == 325 ? 18 : 268,0,1});
        FreshCheckpoint first(first_geometry), second(second_geometry);
        ThresholdEngine first_engine(first_geometry), second_engine(second_geometry);
        const auto tail = canonical_tails(n);
        const auto labels = response_labels();
        const std::uint64_t prefix_seed = fresh_domain(options.seed,n,prefix_tag);
        const std::uint64_t treatment_seed = fresh_domain(options.seed,n,treatment_tag);
        std::vector<int> permutation, vacant, priorities, arm;
        std::array<long double,response_count> sum{};
        std::array<std::array<long double,response_count>,response_count> outer{};
        std::array<std::uint64_t,9> cells{};
        std::array<std::uint64_t,3> gram_ranks{};
        std::array<std::uint64_t,2> zero_b{}, equal_labels{}, active_pairs{};
        std::uint64_t completed_tail_permutations = 0;
        gzFile csv = gzopen((out/"prefixes.csv.gz").c_str(),"wb1");
        if (!csv) throw std::runtime_error("cannot create prefix output");
        std::ostringstream header;
        header << "N,batch,counter,k0,first_rank,second_rank,gram_rank,B_first,B_second";
        for (const auto& label: labels) header << ',' << label;
        header << '\n';
        gzputs(csv,header.str().c_str());
        const auto start = std::chrono::steady_clock::now();
        for (int offset = 0; offset < options.prefixes; ++offset) {
            const std::uint64_t counter = std::uint64_t(options.batch)*options.prefixes+offset;
            counter_permutation(n,prefix_seed,counter,permutation);
            const int first_rank = first.load(permutation,k0), second_rank = second.load(permutation,k0);
            ++cells[3*first_rank+second_rank];
            std::array<long double,response_count> value{};
            std::array<long double,2> amplitude{};
            int gram_rank = -1;
            if (first_rank == 0 && second_rank == 0) {
                vacant.assign(permutation.begin()+k0,permutation.end());
                std::sort(vacant.begin(),vacant.end());
                const auto law = source_normal_law(first,second,vacant);
                gram_rank = law.rank; ++gram_ranks[gram_rank];
                for (int source = 0; source < 2; ++source) {
                    if (!law.maximum[source]) { ++zero_b[source]; continue; }
                    amplitude[source] = static_cast<long double>(law.maximum[source])/
                                        static_cast<long double>(law.denominator);
                    for (int rep = 0; rep < options.reps; ++rep) {
                        ++active_pairs[source];
                        SplitMixStream choose(treatment_key(treatment_seed,counter,rep,source,0));
                        const auto selected = coupled_labels(law,source,vacant,choose);
                        if (selected[0] == selected[1]) { ++equal_labels[source]; continue; }
                        SplitMixStream remaining(treatment_key(treatment_seed,counter,rep,source,1));
                        shared_remainder(vacant,remaining,priorities);
                        std::array<std::array<long double,observable_count>,2> plus;
                        for (int sign = 0; sign < 2; ++sign) {
                            arm_permutation(permutation,k0,selected[sign],priorities,arm);
                            const std::array<std::pair<int,int>,2> births{first_engine.ranks(arm),second_engine.ranks(arm)};
                            ++completed_tail_permutations;
                            for (int receiver = 0; receiver < 2; ++receiver) {
                                const auto f = readout(births[receiver],n,tail);
                                if (sign == 0) plus[receiver] = f;
                                else for (int j = 0; j < observable_count; ++j)
                                    value[(source*2+receiver)*observable_count+j] +=
                                        amplitude[source]*(plus[receiver][j]-f[j])/(2*options.reps);
                            }
                        }
                    }
                }
            }
            for (int j = 0; j < response_count; ++j) {
                sum[j] += value[j];
                for (int k = 0; k < response_count; ++k) outer[j][k] += value[j]*value[k];
            }
            std::ostringstream row;
            row << std::setprecision(17) << n << ',' << options.batch << ',' << counter << ',' << k0 << ','
                << first_rank << ',' << second_rank << ',' << gram_rank << ',' << amplitude[0] << ',' << amplitude[1];
            for (auto v: value) row << ',' << v;
            row << '\n';
            if (gzputs(csv,row.str().c_str()) <= 0) throw std::runtime_error("prefix output write failed");
        }
        if (gzclose(csv) != Z_OK) throw std::runtime_error("prefix gzip close failed");
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream meta(out/"batch.json");
        meta << std::setprecision(19) << "{\n\"code_commit\":\"" << options.commit
             << "\",\n\"backend_blob\":\"fcd4dc0fbfeb498713ee99efda990187cd762b34\",\n"
             << "\"protocol_commit\":\"43079652\",\n\"N\":" << n << ",\n\"k0\":" << k0
             << ",\n\"periods\":[[" << n << ',' << (n == 325 ? 57 : 132) << ",0,1],["
             << n << ',' << (n == 325 ? 18 : 268) << ",0,1]],\n\"p_ref\":" << p_reference
             << ",\n\"batch\":" << options.batch << ",\n\"full_prefix_denominator\":" << options.prefixes
             << ",\n\"counter_start\":" << std::uint64_t(options.batch)*options.prefixes
             << ",\n\"counter_end_exclusive\":" << std::uint64_t(options.batch+1)*options.prefixes
             << ",\n\"paired_reps_per_own_axis\":" << options.reps << ",\n\"master_seed\":" << options.seed
             << ",\n\"prefix_domain_seed\":" << prefix_seed << ",\n\"treatment_domain_seed\":" << treatment_seed
             << ",\n\"technical_smoke\":" << (options.smoke ? "true" : "false")
             << ",\n\"original00_count\":" << cells[0] << ",\n\"rankcell_counts_00_to_22\":";
        json_numbers(meta,cells);
        meta << ",\n\"source_gram_rank_counts_0_1_2_in00\":"; json_numbers(meta,gram_ranks);
        meta << ",\n\"zero_B_own_axes\":"; json_numbers(meta,zero_b);
        meta << ",\n\"nonzero_B_paired_draws\":"; json_numbers(meta,active_pairs);
        meta << ",\n\"identical_label_pairs_skipped\":"; json_numbers(meta,equal_labels);
        meta << ",\n\"completed_tail_permutations_two_receivers_each\":" << completed_tail_permutations
             << ",\n\"labels\":[";
        for (int j = 0; j < response_count; ++j) meta << (j ? ",\"" : "\"") << labels[j] << '"';
        meta << "],\n\"sum_prefix_responses\":"; json_numbers(meta,sum);
        auto mean = sum;
        for (auto& v: mean) v /= options.prefixes;
        meta << ",\n\"mean_full_prefix\":"; json_numbers(meta,mean);
        meta << ",\n\"sum_outer_prefix_responses\":[";
        for (int j = 0; j < response_count; ++j) { if (j) meta << ','; json_numbers(meta,outer[j]); }
        meta << "],\n\"elapsed_seconds\":" << elapsed
             << ",\n\"compiler\":\"" << __VERSION__ << "\",\n"
             << "\"semantics\":\"Fresh iid uniform ordered prefixes; non00 exactly zero, retained in denominator. Rank-aware exact source projection (MP at rank1, zero at rank0). Rational q+/- common inverseCDF label coupling; shared uniform remaining-label priority order; both receivers share each arm. Response B*(Fplus-Fminus)/2; no Taylor approximation or old-source read.\"\n}\n";
        meta.close();
        if (!meta) throw std::runtime_error("metadata write failed");
        std::ofstream complete(out/"COMPLETED");
        complete << options.commit << '\n'; complete.close();
        if (!complete) throw std::runtime_error("completion marker write failed");
        // Operational receipt only: scientific response values are never printed.
        std::cout << "COMPLETE N=" << n << " batch=" << options.batch << " fresh_prefixes=" << options.prefixes
                  << " active00=" << cells[0] << " gram_ranks=" << gram_ranks[0] << ',' << gram_ranks[1] << ',' << gram_ranks[2]
                  << " paired_tail_permutations=" << completed_tail_permutations << " elapsed_seconds=" << elapsed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
