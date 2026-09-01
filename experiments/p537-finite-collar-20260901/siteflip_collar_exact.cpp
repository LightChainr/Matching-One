// Exact radius-one finite-collar site-flip sufficient statistics for Issue #537.
//
// Fix x=0 vacant and z=East(x).  The common N25 collar is the injective
// Chebyshev ball B_inf(z,1).  In the alternating mask forced by x=West(z),
// its labelled arms are B_N,W_E,B_S,W_W.  The four diagonal occupations are
// the complete inner landing word.  Global same-colour connectivity is kept
// separately as the outer attachment bits J_B,J_W.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
constexpr int N = 25;
constexpr int PHYSICAL_EDGES = 2 * N;
constexpr std::uint64_t BACKGROUNDS = 1ULL << 23;
constexpr std::uint32_t KEY_SPACE = 1U << 24;

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

std::int64_t integer(const std::string& text) {
    std::size_t end = 0;
    const auto value = std::stoll(text, &end);
    if (end != text.size()) throw std::invalid_argument("invalid integer: " + text);
    return value;
}

struct Kernel {
    std::vector<std::int32_t> g16 = std::vector<std::int32_t>(KEY_SPACE, 0);
    std::size_t rows = 0;
};

Kernel read_kernel(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read kernel: " + path);
    Kernel kernel;
    std::vector<std::uint32_t> seen;
    std::string line;
    int key_column = -1, g_column = -1;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tsv(line);
        if (key_column < 0) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "key" || fields[i] == "packed_key") key_column = int(i);
                if (fields[i] == "g16") g_column = int(i);
            }
            if (key_column < 0 || g_column < 0)
                throw std::runtime_error("kernel needs key (or packed_key) and g16");
            continue;
        }
        const auto wide_key = integer(fields.at(key_column));
        if (wide_key < 0 || wide_key >= KEY_SPACE) throw std::runtime_error("kernel key out of range");
        const auto key = std::uint32_t(wide_key);
        if (std::find(seen.begin(), seen.end(), key) != seen.end())
            throw std::runtime_error("duplicate kernel key");
        const auto value = integer(fields.at(g_column));
        if (value < std::numeric_limits<std::int32_t>::min() ||
            value > std::numeric_limits<std::int32_t>::max())
            throw std::runtime_error("kernel value out of int32 range");
        kernel.g16[key] = std::int32_t(value);
        seen.push_back(key);
    }
    kernel.rows = seen.size();
    if (!kernel.rows) throw std::runtime_error("empty kernel");
    return kernel;
}

struct DSU {
    std::array<int,N> parent{}, size{};
    DSU() {
        for (int i = 0; i < N; ++i) { parent[i] = i; size[i] = 1; }
    }
    int root(int v) {
        while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
        return v;
    }
    void join(int a, int b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a,b);
        parent[b] = a; size[a] += size[b];
    }
};

struct State {
    int q = 0;
    std::array<int,N> occupied_root{};
    std::array<int,N> vacant_root{};
};

struct RowKey {
    std::int8_t dx = 0, dy = 0;
    std::uint8_t k_minus = 0, rank0 = 0, rank1 = 0;
    std::uint8_t arm_mask = 0, alternating = 0, corner_mask = 0;
    std::uint8_t outer_occupied_join = 0, outer_vacant_join = 0;
    std::uint8_t local_source_contact_mask = 0, source_absent = 0;
    std::uint32_t bell0 = 0, bell1 = 0;
    bool operator==(const RowKey& x) const {
        return std::tie(dx,dy,k_minus,rank0,rank1,arm_mask,alternating,corner_mask,
                        outer_occupied_join,outer_vacant_join,local_source_contact_mask,
                        source_absent,bell0,bell1) ==
               std::tie(x.dx,x.dy,x.k_minus,x.rank0,x.rank1,x.arm_mask,x.alternating,x.corner_mask,
                        x.outer_occupied_join,x.outer_vacant_join,x.local_source_contact_mask,
                        x.source_absent,x.bell0,x.bell1);
    }
};

struct RowHash {
    std::size_t operator()(const RowKey& x) const {
        std::uint64_t h = 1469598103934665603ULL;
        auto mix = [&](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
        mix(std::uint8_t(x.dx)); mix(std::uint8_t(x.dy)); mix(x.k_minus);
        mix(x.rank0); mix(x.rank1); mix(x.arm_mask); mix(x.alternating); mix(x.corner_mask);
        mix(x.outer_occupied_join); mix(x.outer_vacant_join); mix(x.local_source_contact_mask);
        mix(x.source_absent); mix(x.bell0); mix(x.bell1);
        return std::size_t(h);
    }
};

struct Totals {
    std::uint64_t count = 0;
    std::int64_t sum_q0 = 0, sum_e0 = 0, sum_a16_0 = 0, sum_q0_a16_0 = 0, sum_e0_a16_0 = 0;
    std::int64_t sum_q1 = 0, sum_e1 = 0, sum_a16_1 = 0, sum_q1_a16_1 = 0, sum_e1_a16_1 = 0;
};

class Enumerator {
    int a, b;
    const Kernel& kernel;
    std::array<std::array<int,4>,N> neighbors{}, edge_id{}; // N,E,S,W
    std::array<std::array<int,4>,N> diagonal{};              // NE,NW,SW,SE
    std::array<std::pair<int,int>,N> coordinate{};
    std::array<bool,N> occupied{};
    int z = -1;
    std::array<int,4> z_cardinal{};
    std::array<int,4> z_corner{}; // NE,SE,SW,NW
    std::vector<int> free_sites;
    std::unordered_map<RowKey,Totals,RowHash> rows;

    static int mod(int x) { x %= N; return x < 0 ? x + N : x; }
    int quotient_key(int x, int y) const { return N*mod(a*x+b*y) + mod(-b*x+a*y); }

    State evaluate() const {
        DSU black, white;
        int k = 0, edges = 0, faces = 0;
        for (int v = 0; v < N; ++v) if (occupied[v]) ++k;
        for (int v = 0; v < N; ++v) {
            if (occupied[v]) {
                for (int direction : {0,1}) if (occupied[neighbors[v][direction]]) {
                    ++edges; black.join(v,neighbors[v][direction]);
                }
            } else {
                for (int direction : {0,1}) if (!occupied[neighbors[v][direction]])
                    white.join(v,neighbors[v][direction]);
                for (int direction : {0,1}) if (!occupied[diagonal[v][direction]])
                    white.join(v,diagonal[v][direction]);
            }
            if (occupied[v] && occupied[neighbors[v][1]] && occupied[neighbors[v][0]] &&
                occupied[diagonal[v][0]]) ++faces;
        }
        std::array<bool,N> black_seen{}, white_seen{};
        int black_components = 0, white_components = 0;
        State state;
        state.occupied_root.fill(-1);
        state.vacant_root.fill(-1);
        for (int v = 0; v < N; ++v) {
            if (occupied[v]) {
                const int root = black.root(v);
                state.occupied_root[v] = root;
                if (!black_seen[root]) { black_seen[root] = true; ++black_components; }
            } else {
                const int root = white.root(v);
                state.vacant_root[v] = root;
                if (!white_seen[root]) { white_seen[root] = true; ++white_components; }
            }
        }
        state.q = black_components - white_components - (k - edges + faces);
        if (state.q < -1 || state.q > 1) throw std::logic_error("digital rank outside {-1,0,1}");
        return state;
    }

    int outside_id(int center, int direction, const State& state) const {
        const int u = neighbors[center][direction];
        return occupied[u] ? state.occupied_root[u] : N + edge_id[center][direction];
    }

    std::uint32_t bell_key(int y, const State& state) const {
        std::array<int,N+PHYSICAL_EDGES> label;
        label.fill(-1);
        int next = 0;
        std::uint32_t key = 0;
        for (int side = 0; side < 2; ++side) {
            const int center = side ? y : 0;
            for (int direction = 0; direction < 4; ++direction) {
                const int id = outside_id(center,direction,state);
                if (label[id] < 0) label[id] = next++;
                key |= std::uint32_t(label[id]) << (3*(4*side+direction));
            }
        }
        return key;
    }

    int collar_corner_mask() const {
        int result = 0;
        for (int i = 0; i < 4; ++i) result |= int(occupied[z_corner[i]]) << i;
        return result;
    }

    // Return bit 0 for B_N and bit 1 for B_S.  Under arm_mask=5 every
    // occupied collar vertex belongs to exactly one of these inner arms.
    int local_black_arm(int v) const {
        if (!occupied[v]) return 0;
        if (v == z_cardinal[0] || v == z_corner[0] || v == z_corner[3]) return 1;
        if (v == z_cardinal[2] || v == z_corner[1] || v == z_corner[2]) return 2;
        return 0;
    }

    int local_source_contact_mask(int y, bool absent) const {
        if (absent) return 0;
        int result = 0;
        for (int center : {0,y}) for (int direction = 0; direction < 4; ++direction)
            result |= local_black_arm(neighbors[center][direction]);
        return result;
    }

    void consume(std::uint64_t mask) {
        occupied.fill(false);
        int k_minus = 0;
        for (std::size_t i = 0; i < free_sites.size(); ++i) {
            occupied[free_sites[i]] = (mask >> i) & 1ULL;
            k_minus += occupied[free_sites[i]];
        }
        // x=0 and z are off in the frozen background.
        const State state0 = evaluate();
        int arm_mask = 0;
        for (int direction = 0; direction < 4; ++direction)
            arm_mask |= int(occupied[z_cardinal[direction]]) << direction;
        const bool alternating = arm_mask == 5 || arm_mask == 10;
        const int corner_mask = collar_corner_mask();

        occupied[z] = true;
        const State state1 = evaluate();
        occupied[z] = false;
        const int e0 = state0.q*state0.q, e1 = state1.q*state1.q;

        int outer_occupied_join = 0, outer_vacant_join = 0;
        if (alternating) {
            std::array<int,2> black_direction{}, white_direction{};
            int bi = 0, wi = 0;
            for (int direction = 0; direction < 4; ++direction) {
                if (occupied[z_cardinal[direction]]) black_direction[bi++] = direction;
                else white_direction[wi++] = direction;
            }
            if (bi != 2 || wi != 2) throw std::logic_error("alternating collar arm count");
            outer_occupied_join =
                state0.occupied_root[z_cardinal[black_direction[0]]] ==
                state0.occupied_root[z_cardinal[black_direction[1]]];
            // With z occupied, state1's matching graph is exactly the off-z
            // vacant graph used to attach the two separator arms.
            outer_vacant_join =
                state1.vacant_root[z_cardinal[white_direction[0]]] ==
                state1.vacant_root[z_cardinal[white_direction[1]]];
            const int predicted = outer_occupied_join + outer_vacant_join - 1;
            if (state1.q-state0.q != predicted)
                throw std::logic_error("finite-collar identity failure");
            if (!outer_occupied_join && !outer_vacant_join)
                throw std::logic_error("forbidden double-distinct outer attachment");
        }

        constexpr std::uint32_t SOURCE_ABSENT_SENTINEL = KEY_SPACE;
        for (int y = 1; y < N; ++y) {
            if (y == z) continue;
            const bool absent = occupied[y];
            std::uint32_t bell0 = SOURCE_ABSENT_SENTINEL, bell1 = SOURCE_ABSENT_SENTINEL;
            std::int64_t a16_0 = 0, a16_1 = 0;
            if (!absent) {
                bell0 = bell_key(y,state0);
                occupied[z] = true;
                bell1 = bell_key(y,state1);
                occupied[z] = false;
                a16_0 = std::int64_t(kernel.g16[bell0]);
                a16_1 = std::int64_t(kernel.g16[bell1]);
            }
            const RowKey key{
                std::int8_t(coordinate[y].first),std::int8_t(coordinate[y].second),
                std::uint8_t(k_minus),std::uint8_t(state0.q+1),std::uint8_t(state1.q+1),
                std::uint8_t(arm_mask),std::uint8_t(alternating),std::uint8_t(corner_mask),
                std::uint8_t(outer_occupied_join),std::uint8_t(outer_vacant_join),
                std::uint8_t(local_source_contact_mask(y,absent)),std::uint8_t(absent),bell0,bell1
            };
            auto& total = rows[key];
            ++total.count;
            total.sum_q0 += state0.q; total.sum_e0 += e0;
            total.sum_a16_0 += a16_0; total.sum_q0_a16_0 += state0.q*a16_0;
            total.sum_e0_a16_0 += e0*a16_0;
            total.sum_q1 += state1.q; total.sum_e1 += e1;
            total.sum_a16_1 += a16_1; total.sum_q1_a16_1 += state1.q*a16_1;
            total.sum_e1_a16_1 += e1*a16_1;
        }
    }

public:
    Enumerator(int aa, int bb, const Kernel& table) : a(aa), b(bb), kernel(table) {
        if (!((a == 5 && b == 0) || (a == 4 && b == 3)))
            throw std::invalid_argument("fixed N25 geometries are (5,0),(4,3)");
        for (auto& x : edge_id) x.fill(-1);
        std::array<int,N*N> index; index.fill(-1);
        std::vector<std::pair<int,int>> representatives{{0,0}};
        index[quotient_key(0,0)] = 0;
        for (std::size_t i = 0; i < representatives.size(); ++i) {
            const auto [x,y] = representatives[i];
            for (const auto [dx,dy] : std::array<std::pair<int,int>,2>{{{1,0},{0,1}}}) {
                const int key = quotient_key(x+dx,y+dy);
                if (index[key] < 0) {
                    index[key] = int(representatives.size());
                    representatives.emplace_back(x+dx,y+dy);
                }
            }
        }
        if (representatives.size() != N) throw std::logic_error("wrong quotient area");
        for (int v = 0; v < N; ++v) {
            const auto [x,y] = representatives[v];
            neighbors[v] = {index[quotient_key(x,y+1)],index[quotient_key(x+1,y)],
                            index[quotient_key(x,y-1)],index[quotient_key(x-1,y)]};
            diagonal[v] = {index[quotient_key(x+1,y+1)],index[quotient_key(x-1,y+1)],
                           index[quotient_key(x-1,y-1)],index[quotient_key(x+1,y-1)]};
        }
        int next_edge = 0;
        for (int v = 0; v < N; ++v) for (int direction : {0,1}) {
            const int u = neighbors[v][direction], reverse = direction+2;
            if (edge_id[v][direction] >= 0 || edge_id[u][reverse] >= 0)
                throw std::logic_error("bad reciprocal edge assignment");
            edge_id[v][direction] = edge_id[u][reverse] = next_edge++;
        }
        if (next_edge != PHYSICAL_EDGES) throw std::logic_error("wrong edge count");

        z = neighbors[0][1];
        z_cardinal = neighbors[z];
        z_corner = {diagonal[z][0],diagonal[z][3],diagonal[z][2],diagonal[z][1]};
        std::array<bool,N> collar_seen{};
        collar_seen[z] = true;
        for (int v : z_cardinal) collar_seen[v] = true;
        for (int v : z_corner) collar_seen[v] = true;
        if (std::count(collar_seen.begin(),collar_seen.end(),true) != 9)
            throw std::logic_error("radius-one collar is not injective");
        for (int v = 0; v < N; ++v) if (v != 0 && v != z) free_sites.push_back(v);
        if (free_sites.size() != 23) throw std::logic_error("wrong free-site count");

        std::array<bool,N> found{};
        std::array<std::tuple<int,int,int>,N> best{};
        for (int x = -12; x <= 12; ++x) for (int y = -12; y <= 12; ++y) {
            const int v = index[quotient_key(x,y)];
            const auto candidate = std::make_tuple(x*x+y*y,x,y);
            if (!found[v] || candidate < best[v]) {
                found[v] = true; best[v] = candidate; coordinate[v] = {x,y};
            }
        }
        for (bool ok : found) if (!ok) throw std::logic_error("missing quotient coordinate");
        rows.reserve(1 << 20);
    }

    void run(const std::string& output) {
        if (std::ifstream(output).good()) throw std::runtime_error("output already exists");
        for (std::uint64_t mask = 0; mask < BACKGROUNDS; ++mask) consume(mask);
        std::vector<std::pair<RowKey,Totals>> ordered(rows.begin(),rows.end());
        std::sort(ordered.begin(),ordered.end(),[](const auto& lhs, const auto& rhs) {
            const auto& x = lhs.first; const auto& y = rhs.first;
            return std::tie(x.k_minus,x.rank0,x.rank1,x.arm_mask,x.corner_mask,
                            x.outer_occupied_join,x.outer_vacant_join,x.local_source_contact_mask,
                            x.source_absent,x.dx,x.dy,x.bell0,x.bell1) <
                   std::tie(y.k_minus,y.rank0,y.rank1,y.arm_mask,y.corner_mask,
                            y.outer_occupied_join,y.outer_vacant_join,y.local_source_contact_mask,
                            y.source_absent,y.dx,y.dy,y.bell0,y.bell1);
        });
        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot create output");
        out << "k_minus,rank0,rank1,arm_mask,alternating_four_arm,collar_corner_mask,"
               "outer_occupied_join,outer_vacant_join,local_source_contact_mask,"
               "source_absent,y_dx,y_dy,bell0,bell1,count,"
               "sum_q0,sum_E0,sum_a16_0,sum_q0_a16_0,sum_E0_a16_0,"
               "sum_q1,sum_E1,sum_a16_1,sum_q1_a16_1,sum_E1_a16_1\n";
        std::uint64_t total_count = 0;
        for (const auto& [key,total] : ordered) {
            total_count += total.count;
            out << int(key.k_minus) << ',' << int(key.rank0) << ',' << int(key.rank1) << ','
                << int(key.arm_mask) << ',' << int(key.alternating) << ',' << int(key.corner_mask) << ','
                << int(key.outer_occupied_join) << ',' << int(key.outer_vacant_join) << ','
                << int(key.local_source_contact_mask) << ',' << int(key.source_absent) << ','
                << int(key.dx) << ',' << int(key.dy) << ',' << key.bell0 << ',' << key.bell1 << ','
                << total.count << ',' << total.sum_q0 << ',' << total.sum_e0 << ','
                << total.sum_a16_0 << ',' << total.sum_q0_a16_0 << ',' << total.sum_e0_a16_0 << ','
                << total.sum_q1 << ',' << total.sum_e1 << ',' << total.sum_a16_1 << ','
                << total.sum_q1_a16_1 << ',' << total.sum_e1_a16_1 << '\n';
        }
        if (!out) throw std::runtime_error("output write failed");
        std::cerr << "rows=" << ordered.size() << " pair_fibres=" << total_count << '\n';
    }
};
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::invalid_argument("usage: siteflip_collar_exact a b kernel.tsv output.csv");
        const int a = int(integer(argv[1])), b = int(integer(argv[2]));
        const auto start = std::chrono::steady_clock::now();
        const Kernel kernel = read_kernel(argv[3]);
        Enumerator(a,b,kernel).run(argv[4]);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"completed\",\"schema\":\"matching-one/p537-finite-collar-producer/v1\","
                  << "\"N\":25,\"a\":" << a << ",\"b\":" << b
                  << ",\"fixed_x\":0,\"fixed_z_direction\":\"E\",\"backgrounds\":" << BACKGROUNDS
                  << ",\"kernel_rows\":" << kernel.rows << ",\"collar\":\"B_inf(z,1)\","
                  << "\"a_scale_denominator\":16,\"fixed_x_pair_scale_denominator\":25,"
                  << "\"c4_fixed_z_orbit_multiplier\":4,\"source_absent_sentinel\":" << KEY_SPACE
                  << ",\"elapsed_seconds\":" << elapsed << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
