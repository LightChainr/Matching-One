// Exact first-witness producer for the Issue #537 one-defect gate.
//
// Fix x=0 vacant and z=East(x), enumerate axis-N25 backgrounds in increasing
// bit-mask order, and stop at the first physical z flip that changes q+1,
// the x/y Bell partition, and the actual g16 source value.  Requiring the
// last condition keeps the falsifier independent of how the beta*B Schur
// counterterm is allocated to topological cells.  The certificate also records the
// joint x4+y4+z4 carrier partition and separate outer black/white partitions;
// independent Bell canonicalization is therefore checkable from the joint
// key instead of being trusted as an unrelated label.
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
    std::array<int,N> occupied_degree{};
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
        state.occupied_degree.fill(0);
        for (int v = 0; v < N; ++v) {
            if (occupied[v]) {
                const int root = black.root(v);
                state.occupied_root[v] = root;
                if (!black_seen[root]) { black_seen[root] = true; ++black_components; }
                for (int direction = 0; direction < 4; ++direction)
                    state.occupied_degree[v] += occupied[neighbors[v][direction]];
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

    std::uint64_t joint_key(int y, const State& state) const {
        std::array<int,N+PHYSICAL_EDGES> label;
        label.fill(-1);
        int next = 0;
        std::uint64_t key = 0;
        const std::array<int,3> centers{{0,y,z}};
        for (int side = 0; side < 3; ++side) {
            for (int direction = 0; direction < 4; ++direction) {
                const int id = outside_id(centers[side],direction,state);
                if (label[id] < 0) label[id] = next++;
                if (label[id] >= 16) throw std::logic_error("joint key needs more than four bits");
                key |= std::uint64_t(label[id]) << (4*(4*side+direction));
            }
        }
        return key;
    }

    static std::array<int,12> joint_labels(std::uint64_t key) {
        std::array<int,12> labels{};
        for (int i = 0; i < 12; ++i) labels[i] = int((key >> (4*i)) & 15U);
        return labels;
    }

    static std::uint32_t bell_from_joint(std::uint64_t key) {
        std::array<int,16> canonical;
        canonical.fill(-1);
        int next = 0;
        std::uint32_t bell = 0;
        for (int i = 0; i < 8; ++i) {
            const int label = int((key >> (4*i)) & 15U);
            if (canonical[label] < 0) canonical[label] = next++;
            bell |= std::uint32_t(canonical[label]) << (3*i);
        }
        return bell;
    }

    static int terminal_incidence(std::uint64_t key) {
        const auto labels = joint_labels(key);
        std::array<bool,16> source{}, thermal{};
        for (int i = 0; i < 8; ++i) source[labels[i]] = true;
        for (int i = 8; i < 12; ++i) thermal[labels[i]] = true;
        int count = 0;
        for (int label = 0; label < 16; ++label) count += source[label] && thermal[label];
        return count;
    }

    std::string source_component(int y) const {
        const auto [dx,dy] = coordinate[y];
        if (std::abs(dx)+std::abs(dy) == 1) return "nn_other";
        if ((std::abs(dx) == 2 && dy == 0) || (dx == 0 && std::abs(dy) == 2)) return "axial2";
        if (std::abs(dx) == 1 && std::abs(dy) == 1) return "diag1";
        return "unscored";
    }

    static void write_array(std::ostream& out, const std::array<int,12>& values) {
        out << '[';
        for (int i = 0; i < 12; ++i) {
            if (i) out << ',';
            out << values[i];
        }
        out << ']';
    }

    static void write_array(std::ostream& out, const std::array<int,4>& values) {
        out << '[';
        for (int i = 0; i < 4; ++i) {
            if (i) out << ',';
            out << values[i];
        }
        out << ']';
    }

    std::array<int,4> black_partition(const State& state) const {
        std::array<int,4> result{{-1,-1,-1,-1}};
        std::array<int,N> label;
        label.fill(-1);
        int next = 0;
        for (int direction = 0; direction < 4; ++direction) {
            const int v = z_cardinal[direction];
            if (!occupied[v]) continue;
            const int root = state.occupied_root[v];
            if (label[root] < 0) label[root] = next++;
            result[direction] = label[root];
        }
        return result;
    }

    std::array<int,4> white_cut_partition(const State& state1) const {
        std::array<int,4> result{{-1,-1,-1,-1}};
        std::array<int,N> label;
        label.fill(-1);
        int next = 0;
        for (int direction = 0; direction < 4; ++direction) {
            const int v = z_cardinal[direction];
            if (occupied[v]) continue;
            const int root = state1.vacant_root[v];
            if (label[root] < 0) label[root] = next++;
            result[direction] = label[root];
        }
        return result;
    }

    int source_port_occupied_mask(int y) const {
        int result = 0;
        const std::array<int,2> centers{{0,y}};
        for (int side = 0; side < 2; ++side)
            for (int direction = 0; direction < 4; ++direction)
                result |= int(occupied[neighbors[centers[side]][direction]]) << (4*side+direction);
        return result;
    }

    bool source_touches_black_landing_component(int y, const State& state) const {
        std::array<bool,N> landing_root{};
        for (int direction = 0; direction < 4; ++direction) {
            const int v = z_cardinal[direction];
            if (occupied[v]) landing_root[state.occupied_root[v]] = true;
        }
        for (int center : {0,y}) for (int direction = 0; direction < 4; ++direction) {
            const int v = neighbors[center][direction];
            if (occupied[v] && landing_root[state.occupied_root[v]]) return true;
        }
        return false;
    }

    bool off_port_degree_branch(const State& state) const {
        std::array<bool,N> landing_root{}, is_port{};
        for (int direction = 0; direction < 4; ++direction) {
            const int v = z_cardinal[direction];
            is_port[v] = true;
            if (occupied[v]) landing_root[state.occupied_root[v]] = true;
        }
        for (int v = 0; v < N; ++v)
            if (occupied[v] && !is_port[v] && landing_root[state.occupied_root[v]] &&
                state.occupied_degree[v] >= 3) return true;
        return false;
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

    void write_witness(const std::string& output, std::uint64_t mask, int y, int k_minus,
                       const State& state0, const State& state1, int arm_mask,
                       int corner_mask, int outer_occupied_join, int outer_vacant_join,
                       std::uint32_t bell0, std::uint32_t bell1,
                       std::int64_t a16_0, std::int64_t a16_1) {
        if (std::ifstream(output).good()) throw std::runtime_error("output already exists");

        const int local_contact = local_source_contact_mask(y,false);
        const int source_ports0 = source_port_occupied_mask(y);
        const bool global_contact = source_touches_black_landing_component(y,state0);
        const bool degree_branch = off_port_degree_branch(state0);
        const auto black0 = black_partition(state0);
        const std::uint64_t joint0 = joint_key(y,state0);

        occupied[z] = true;
        const int source_ports1 = source_port_occupied_mask(y);
        const auto black1 = black_partition(state1);
        const auto white_cut = white_cut_partition(state1);
        const std::uint64_t joint1 = joint_key(y,state1);
        occupied[z] = false;

        if (bell_from_joint(joint0) != bell0 || bell_from_joint(joint1) != bell1)
            throw std::logic_error("joint source restriction does not reproduce Bell keys");

        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot create witness output");
        out << "{\n"
            << "  \"schema\": \"matching-one/p537-one-defect-witness/v1\",\n"
            << "  \"status\": \"first_observable_diagonal_edge_in_axis_scan\",\n"
            << "  \"N\": 25,\n"
            << "  \"geometry\": {\"id\": \"axis\", \"a\": " << a << ", \"b\": " << b << "},\n"
            << "  \"scan_contract\": {\"background_order\": \"ascending 23-bit mask\", "
               "\"source_order\": \"vertex id 1..24 excluding z\", "
               "\"filter\": \"source present; rank, Bell, and g16 all change\", "
               "\"stopped_at_first_match\": true},\n"
            << "  \"transition_id\": \"axis-N25:x0:zE:y(" << coordinate[y].first << ','
            << coordinate[y].second << "):eta" << mask << ":0to1\",\n"
            << "  \"background_mask\": " << mask << ",\n"
            << "  \"k_minus\": " << k_minus << ",\n"
            << "  \"vertices\": {\"x\": 0, \"z\": " << z << ", \"y\": " << y
            << ", \"y_coordinate\": [" << coordinate[y].first << ',' << coordinate[y].second << "]},\n"
            << "  \"free_site_order\": [";
        for (std::size_t i = 0; i < free_sites.size(); ++i) {
            if (i) out << ',';
            out << free_sites[i];
        }
        out << "],\n  \"occupied_background_vertices\": [";
        bool first = true;
        for (int v = 0; v < N; ++v) if (occupied[v]) {
            if (!first) out << ',';
            out << v;
            first = false;
        }
        out << "],\n"
            << "  \"collar\": {\"arm_order\": [\"N\",\"E\",\"S\",\"W\"], "
            << "\"arm_mask\": " << arm_mask << ", \"alternating_four_arm\": "
            << ((arm_mask == 5 || arm_mask == 10) ? "true" : "false")
            << ", \"outer_join_identity_applicable\": "
            << ((arm_mask == 5 || arm_mask == 10) ? "true" : "false")
            << ", \"corner_mask_NE_SE_SW_NW\": " << corner_mask
            << ", \"outer_black_join_J_B\": " << outer_occupied_join
            << ", \"outer_white_join_J_W\": " << outer_vacant_join << "},\n"
            << "  \"outer_C\": {\"meaning\": \"global first-occurrence carrier partition on x4+y4+z4 ports\",\n"
            << "    \"before_key\": " << joint0 << ", \"before_labels\": ";
        write_array(out,joint_labels(joint0));
        out << ",\n    \"after_key\": " << joint1 << ", \"after_labels\": ";
        write_array(out,joint_labels(joint1));
        out << ",\n    \"terminal_incidence_before\": " << terminal_incidence(joint0)
            << ", \"terminal_incidence_after\": " << terminal_incidence(joint1) << "},\n"
            << "  \"outer_B\": {\"meaning\": \"global NN-black partition of the four z-cardinal ports; -1 means vacant\", "
               "\"before\": ";
        write_array(out,black0);
        out << ", \"after\": ";
        write_array(out,black1);
        out << "},\n"
            << "  \"outer_W\": {\"meaning\": \"global matching-white partition in the off-z cut; -1 means occupied\", "
               "\"cut_partition\": ";
        write_array(out,white_cut);
        out << "},\n"
            << "  \"source_ports\": {\"occupied_mask_before\": " << source_ports0
            << ", \"occupied_mask_after\": " << source_ports1
            << ", \"source_component\": \"" << source_component(y) << "\""
            << ", \"bell_before\": " << bell0 << ", \"bell_after\": " << bell1
            << ", \"bell_from_joint_before\": " << bell_from_joint(joint0)
            << ", \"bell_from_joint_after\": " << bell_from_joint(joint1)
            << ", \"g16_before\": " << a16_0 << ", \"g16_after\": " << a16_1 << "},\n"
            << "  \"off_port\": {\"local_black_arm_contact_mask\": " << local_contact
            << ", \"source_touches_global_black_landing_component\": "
            << (global_contact ? "true" : "false")
            << ", \"black_landing_component_has_degree3_vertex_outside_z_ports\": "
            << (degree_branch ? "true" : "false") << "},\n"
            << "  \"nodes\": [\n"
            << "    {\"state\": 0, \"q\": " << state0.q << ", \"rank_index_q_plus_1\": " << state0.q+1
            << ", \"E\": " << state0.q*state0.q << ", \"bell\": " << bell0
            << ", \"joint_C\": " << joint0 << ", \"g16\": " << a16_0 << "},\n"
            << "    {\"state\": 1, \"q\": " << state1.q << ", \"rank_index_q_plus_1\": " << state1.q+1
            << ", \"E\": " << state1.q*state1.q << ", \"bell\": " << bell1
            << ", \"joint_C\": " << joint1 << ", \"g16\": " << a16_1 << "}\n"
            << "  ],\n"
            << "  \"edge\": {\"physical_move\": \"occupy z with every other site fixed\", "
               "\"changes_landing_rank\": true, \"changes_source_Bell\": true, "
               "\"changes_source_value_g16\": true, \"changes_joint_component_map\": true},\n"
            << "  \"sufficient_statistics\": {\"count\": 1, "
            << "\"sum_q0\": " << state0.q << ", \"sum_E0\": " << state0.q*state0.q
            << ", \"sum_a16_0\": " << a16_0 << ", \"sum_q0_a16_0\": " << state0.q*a16_0
            << ", \"sum_E0_a16_0\": " << state0.q*state0.q*a16_0
            << ", \"sum_q1\": " << state1.q << ", \"sum_E1\": " << state1.q*state1.q
            << ", \"sum_a16_1\": " << a16_1 << ", \"sum_q1_a16_1\": " << state1.q*a16_1
            << ", \"sum_E1_a16_1\": " << state1.q*state1.q*a16_1 << "}\n"
            << "}\n";
        if (!out) throw std::runtime_error("witness output write failed");
    }

    bool consume(std::uint64_t mask, const std::string& output) {
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
            const bool observable_diagonal_edge =
                !absent && state0.q != state1.q && bell0 != bell1 && a16_0 != a16_1;
            if (observable_diagonal_edge) {
                write_witness(output,mask,y,k_minus,state0,state1,arm_mask,corner_mask,
                              outer_occupied_join,outer_vacant_join,bell0,bell1,a16_0,a16_1);
                return true;
            }
        }
        return false;
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
    }

    std::uint64_t run(const std::string& output) {
        if (std::ifstream(output).good()) throw std::runtime_error("output already exists");
        for (std::uint64_t mask = 0; mask < BACKGROUNDS; ++mask)
            if (consume(mask,output)) return mask+1;
        throw std::runtime_error("complete frozen scan found no diagonal edge");
    }
};
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5)
            throw std::invalid_argument("usage: one_defect_witness_exact 5 0 kernel.tsv witness.json");
        const int a = int(integer(argv[1])), b = int(integer(argv[2]));
        if (a != 5 || b != 0) throw std::invalid_argument("the frozen first-witness scan is axis N25 (5,0)");
        const auto start = std::chrono::steady_clock::now();
        const Kernel kernel = read_kernel(argv[3]);
        const std::uint64_t scanned = Enumerator(a,b,kernel).run(argv[4]);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"first_observable_diagonal_edge_stop\",\"schema\":\"matching-one/p537-one-defect-witness-run/v1\","
                  << "\"N\":25,\"a\":" << a << ",\"b\":" << b
                  << ",\"fixed_x\":0,\"fixed_z_direction\":\"E\",\"backgrounds_scanned\":" << scanned
                  << ",\"kernel_rows\":" << kernel.rows << ",\"collar\":\"B_inf(z,1)\","
                  << "\"a_scale_denominator\":16,\"fixed_x_pair_scale_denominator\":25,"
                  << "\"c4_fixed_z_orbit_multiplier\":4,\"source_absent_sentinel\":" << KEY_SPACE
                  << ",\"elapsed_seconds\":" << elapsed << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
