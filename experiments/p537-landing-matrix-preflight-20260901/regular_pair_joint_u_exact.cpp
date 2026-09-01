// Exact fixed-origin two-insertion moments for canonical Kreg=K2bar+K0bar.
// N=25 only; origin is vacant, all other occupations are enumerated once.
// Kernel ports are (xN,xE,xS,xW,yN,yE,yS,yW). The original occupation
// determines q/E; the virtual colour closure never changes its rank.
// No root, U functional, conditional normalization, or Monte Carlo is used.
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
#include <utility>
#include <vector>

namespace {
constexpr int n = 25;
constexpr int physical_edges = 2*n;
constexpr std::uint64_t configurations = 1ULL << (n-1);
constexpr std::uint32_t key_space = 1U << 24;

// A 64 MiB direct-address table removes hashing from the inner enumeration.
// All stored values and accumulated moments remain exact integers.
struct Kernel {
    std::vector<std::int32_t> g16 = std::vector<std::int32_t>(key_space, 0);
    std::size_t rows = 0;
    std::int64_t max_abs_g16 = 0;
};

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
                if (fields[i] == "key" || fields[i] == "packed_key") key_column = static_cast<int>(i);
                if (fields[i] == "g16") g_column = static_cast<int>(i);
            }
            if (key_column < 0 || g_column < 0)
                throw std::runtime_error("kernel needs key (or packed_key) and g16 TSV columns");
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max(key_column, g_column)))
            throw std::runtime_error("incomplete kernel row");
        const auto wide_key = integer(fields[key_column]);
        if (wide_key < 0 || wide_key >= key_space)
            throw std::runtime_error("kernel key exceeds eight 3-bit labels");
        const auto key = static_cast<std::uint32_t>(wide_key);
        int max_label = -1;
        for (int i = 0; i < 8; ++i) {
            const int label = (key >> (3*i)) & 7;
            if (label > max_label + 1)
                throw std::runtime_error("kernel key is not canonical restricted-growth order");
            max_label = std::max(max_label, label);
        }
        const auto value = integer(fields[g_column]);
        if (value < std::numeric_limits<std::int32_t>::min() ||
            value > std::numeric_limits<std::int32_t>::max())
            throw std::runtime_error("kernel value does not fit exact int32 storage");
        kernel.g16[key] = static_cast<std::int32_t>(value);
        kernel.max_abs_g16 = std::max(kernel.max_abs_g16, value < 0 ? -value : value);
        seen.push_back(key);
    }
    if (!in.eof()) throw std::runtime_error("kernel read failed");
    if (seen.empty()) throw std::runtime_error("empty kernel");
    std::sort(seen.begin(), seen.end());
    if (std::adjacent_find(seen.begin(), seen.end()) != seen.end())
        throw std::runtime_error("duplicate kernel key");
    kernel.rows = seen.size();
    // Even every pair term across every configuration at int32's maximal
    // absolute magnitude is below INT64_MAX. q,E are respectively +/-1,0
    // and 0,1, so the same bound applies to both weighted sums.
    static_assert((n-1)*configurations*(1ULL << 31) <
                  static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()),
                  "the declared N=25 integer moment bound must fit int64");
    return kernel;
}

class RollbackComponents {
    std::array<int, n> parent{}, size{};
    std::vector<std::pair<int, int>> changed;
public:
    int components = 0;
    void activate(int v) { parent[v] = v; size[v] = 1; ++components; }
    int root(int v) const { while (parent[v] != v) v = parent[v]; return v; }
    std::size_t mark() const { return changed.size(); }
    void join(int a, int b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a, b);
        changed.emplace_back(b, size[a]);
        parent[b] = a; size[a] += size[b]; --components;
    }
    // Undo all joins since mark, then the single vertex activated at mark.
    void undo(std::size_t mark) {
        while (changed.size() > mark) {
            const auto [b, old_size] = changed.back(); changed.pop_back();
            size[parent[b]] = old_size; parent[b] = b; ++components;
        }
        --components;
    }
};

struct Totals {
    std::uint64_t count = 0, nonzero_G16_count = 0;
    std::int64_t sum_q = 0, sum_e = 0;
    std::int64_t sum_G16 = 0, sum_G16_q = 0, sum_G16_E = 0;
    std::int64_t sum_G16_contact = 0, sum_G16_q_contact = 0, sum_G16_E_contact = 0;
};

struct JointSum {
    std::int64_t total = 0, contact = 0;
};

enum class SourceMode { all, clean_same, clean_reversed, clean_total };

SourceMode source_mode(const std::string& text) {
    if (text == "all") return SourceMode::all;
    if (text == "clean_same") return SourceMode::clean_same;
    if (text == "clean_reversed") return SourceMode::clean_reversed;
    if (text == "clean_total") return SourceMode::clean_total;
    throw std::invalid_argument("mode must be all, clean_same, clean_reversed, or clean_total");
}

// Return 1 or -1 for the two nonzero branch-free C4 orbits, and zero for
// every other Bell-8 landing.  A branch-free two-bridge partition has six
// exterior blocks: two blocks occur once at each marked site and all four
// remaining port blocks are singletons.  The two shared ports must be
// adjacent at each mark; the sign records whether their cyclic orders agree.
int clean_orbit(std::uint32_t key) {
    std::array<int, 8> label{};
    std::array<std::array<int, 2>, 8> count{};
    int maximum = -1;
    for (int i = 0; i < 8; ++i) {
        label[i] = (key >> (3*i)) & 7;
        maximum = std::max(maximum, label[i]);
        ++count[label[i]][i >= 4];
    }
    if (maximum != 5) return 0;
    std::array<int, 2> shared{};
    int shared_count = 0;
    for (int value = 0; value <= maximum; ++value) {
        if (count[value][0] && count[value][1]) {
            if (count[value][0] != 1 || count[value][1] != 1 || shared_count == 2)
                return 0;
            shared[shared_count++] = value;
        }
    }
    if (shared_count != 2) return 0;
    std::array<std::array<int, 2>, 2> position{};
    for (int side = 0; side < 2; ++side) {
        for (int which = 0; which < 2; ++which) {
            position[side][which] = -1;
            for (int i = 0; i < 4; ++i)
                if (label[4*side+i] == shared[which]) position[side][which] = i;
        }
    }
    const int first = (position[0][1]-position[0][0]+4)%4;
    const int second = (position[1][1]-position[1][0]+4)%4;
    if (!((first == 1 || first == 3) && (second == 1 || second == 3))) return 0;
    return first == second ? 1 : -1;
}

class Enumerator {
    const Kernel& kernel;
    SourceMode mode;
    int a, b;
    std::array<Totals, n+1> histogram{};
    std::array<int, n> colour{}, occupied_root{};
    std::array<bool, n> origin_contact{};
    std::array<std::array<int, 4>, n> neighbors{}, edge_id{}; // N,E,S,W
    std::array<std::vector<int>, n> black_previous, white_previous;
    std::array<std::vector<std::array<int, 4>>, n> closing_faces;
    // IDs 0..24 denote occupied roots; 25..74 denote physical empty edges.
    // Only at most eight of these IDs are labelled for a particular pair.
    std::array<int, n+physical_edges> label_of{};
    RollbackComponents black, white;

    static int mod(int x) { x %= n; return x < 0 ? x+n : x; }
    int quotient_key(int x, int y) const { return n*mod(a*x+b*y) + mod(-b*x+a*y); }

    int outside_id(int center, int direction) const {
        const int neighbor = neighbors[center][direction];
        // Never ask for the black root of an inactive/vacant vertex.
        // Two vacant centers sharing an edge see exactly the same edge ID.
        return colour[neighbor] == 1 ? occupied_root[neighbor] : n+edge_id[center][direction];
    }

    JointSum joint_sum() {
        for (int v = 1; v < n; ++v)
            if (colour[v] == 1) occupied_root[v] = black.root(v);

        std::array<int, 4> origin_ids{};
        std::uint32_t origin_key = 0;
        int origin_labels = 0;
        for (int direction = 0; direction < 4; ++direction) {
            const int id = outside_id(0, direction);
            origin_ids[direction] = id;
            int& label = label_of[id];
            if (label < 0) label = origin_labels++;
            origin_key |= static_cast<std::uint32_t>(label) << (3*direction);
        }

        JointSum sum;
        for (int y = 1; y < n; ++y) if (colour[y] == 0) {
            // Keep the fixed origin's four canonical labels; only the new
            // labels introduced by y need to be cleared before the next y.
            std::array<int, 4> introduced{};
            int introduced_count = 0, next_label = origin_labels;
            std::uint32_t key = origin_key;
            for (int direction = 0; direction < 4; ++direction) {
                const int id = outside_id(y, direction);
                int& label = label_of[id];
                if (label < 0) {
                    label = next_label++;
                    introduced[introduced_count++] = id;
                }
                key |= static_cast<std::uint32_t>(label) << (3*(4+direction));
            }
            // Absent canonical keys in the delivered sparse kernel are zero.
            std::int64_t value = kernel.g16[key];
            const int orbit = clean_orbit(key);
            if (orbit && value != 4)
                throw std::logic_error("branch-free nonzero C4 orbit must have g16=4");
            if (mode == SourceMode::clean_same && orbit != 1) value = 0;
            if (mode == SourceMode::clean_reversed && orbit != -1) value = 0;
            if (mode == SourceMode::clean_total && orbit == 0) value = 0;
            sum.total += value;
            if (origin_contact[y]) sum.contact += value;
            for (int i = 0; i < introduced_count; ++i) label_of[introduced[i]] = -1;
        }
        for (int id : origin_ids) label_of[id] = -1;
        return sum;
    }

    void visit(int v, int k, int edges, int faces) {
        if (v == n) {
            // Digital Alexander: original NN black, matching white, and
            // occupied square-cell Euler characteristic. No virtual joins.
            const int q = black.components-white.components-(k-edges+faces);
            if (q < -1 || q > 1) throw std::logic_error("invalid digital rank value");
            const int e = q*q;
            const JointSum g16 = joint_sum();
            auto& row = histogram[k];
            ++row.count;
            row.sum_q += q; row.sum_e += e;
            row.sum_G16 += g16.total;
            row.sum_G16_q += q*g16.total;
            row.sum_G16_E += e*g16.total;
            row.sum_G16_contact += g16.contact;
            row.sum_G16_q_contact += q*g16.contact;
            row.sum_G16_E_contact += e*g16.contact;
            row.nonzero_G16_count += g16.total != 0;
            return;
        }

        colour[v] = 0;
        auto mark = white.mark(); white.activate(v);
        for (int u : white_previous[v]) if (colour[u] == 0) white.join(v, u);
        visit(v+1, k, edges, faces);
        white.undo(mark);

        colour[v] = 1;
        mark = black.mark(); black.activate(v);
        int extra_edges = 0, extra_faces = 0;
        for (int u : black_previous[v]) if (colour[u] == 1) {
            black.join(v, u); ++extra_edges;
        }
        for (const auto& face : closing_faces[v]) {
            bool full = true;
            for (int u : face) if (colour[u] != 1) full = false;
            extra_faces += full;
        }
        visit(v+1, k+1, edges+extra_edges, faces+extra_faces);
        black.undo(mark);
        colour[v] = -1;
    }

public:
    Enumerator(int aa, int bb, const Kernel& table, SourceMode selected)
        : kernel(table), mode(selected), a(aa), b(bb) {
        if (!((a == 5 && b == 0) || (a == 4 && b == 3)))
            throw std::invalid_argument("fixed N=25 pair is (5,0),(4,3)");
        colour.fill(-1);
        label_of.fill(-1);
        for (auto& edges : edge_id) edges.fill(-1);

        std::array<int, n*n> index; index.fill(-1);
        std::vector<std::pair<int, int>> representatives;
        representatives.emplace_back(0, 0); index[quotient_key(0, 0)] = 0;
        for (std::size_t i = 0; i < representatives.size(); ++i) {
            const auto [x, y] = representatives[i];
            for (const auto& step : std::array<std::pair<int, int>, 2>{{{1, 0}, {0, 1}}}) {
                const int xx = x+step.first, yy = y+step.second;
                const int h = quotient_key(xx, yy);
                if (index[h] < 0) {
                    index[h] = static_cast<int>(representatives.size());
                    representatives.emplace_back(xx, yy);
                }
            }
        }
        if (representatives.size() != n) throw std::logic_error("wrong quotient area");

        for (int v = 0; v < n; ++v) {
            const auto [x, y] = representatives[v];
            neighbors[v] = {index[quotient_key(x, y+1)], index[quotient_key(x+1, y)],
                            index[quotient_key(x, y-1)], index[quotient_key(x-1, y)]};
            for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) {
                if (dx == 0 && dy == 0) continue;
                const int u = index[quotient_key(x+dx, y+dy)];
                if (u < v) {
                    white_previous[v].push_back(u);
                    if (dx == 0 || dy == 0) black_previous[v].push_back(u);
                }
            }
            const std::array<int, 4> face{v, index[quotient_key(x+1, y)],
                                         index[quotient_key(x, y+1)], index[quotient_key(x+1, y+1)]};
            closing_faces[*std::max_element(face.begin(), face.end())].push_back(face);
        }
        for (int u : neighbors[0]) origin_contact[u] = true;

        int next_edge = 0;
        for (int v = 0; v < n; ++v) for (int direction : {0, 1}) {
            const int u = neighbors[v][direction], reverse = direction+2;
            if (neighbors[u][reverse] != v || edge_id[v][direction] >= 0 || edge_id[u][reverse] >= 0)
                throw std::logic_error("invalid reciprocal physical edge");
            edge_id[v][direction] = edge_id[u][reverse] = next_edge++;
        }
        if (next_edge != physical_edges) throw std::logic_error("wrong physical edge count");
        for (const auto& edges : edge_id) for (int id : edges)
            if (id < 0) throw std::logic_error("unassigned physical edge");
    }

    void run(const std::string& output) {
        if (std::ifstream(output).good()) throw std::runtime_error("output already exists");
        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot create output: " + output);
        out << "k,count,sum_q,sum_e,sum_G16,sum_G16_q,sum_G16_E,nonzero_G16_count,"
               "sum_G16_contact,sum_G16_q_contact,sum_G16_E_contact\n";

        // The origin is permanently vacant. Its matching component must be
        // present throughout the recursion over the remaining 24 vertices.
        colour[0] = 0;
        white.activate(0);
        visit(1, 0, 0, 0);

        std::uint64_t observed = 0;
        for (const auto& row : histogram) observed += row.count;
        if (observed != configurations) throw std::logic_error("wrong exact configuration count");
        for (int k = 0; k <= n; ++k) {
            const auto& row = histogram[k];
            out << k << ',' << row.count << ',' << row.sum_q << ',' << row.sum_e << ','
                << row.sum_G16 << ',' << row.sum_G16_q << ','
                << row.sum_G16_E << ',' << row.nonzero_G16_count << ',' << row.sum_G16_contact << ','
                << row.sum_G16_q_contact << ',' << row.sum_G16_E_contact << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("output write failed");
    }
};
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6)
            throw std::invalid_argument("usage: regular_pair_joint_u_exact a b kernel.tsv output.csv mode");
        const auto aa = integer(argv[1]), bb = integer(argv[2]);
        if (!((aa == 5 && bb == 0) || (aa == 4 && bb == 3)))
            throw std::invalid_argument("fixed N=25 pair is (5,0),(4,3)");
        const auto start = std::chrono::steady_clock::now();
        const Kernel kernel = read_kernel(argv[3]);
        const auto mode = source_mode(argv[5]);
        Enumerator(static_cast<int>(aa), static_cast<int>(bb), kernel, mode).run(argv[4]);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"completed\",\"a\":" << aa << ",\"b\":" << bb
                  << ",\"N\":25,\"fixed_origin\":0,\"origin_occupied\":false,\"configurations\":"
                  << configurations << ",\"kernel_rows\":" << kernel.rows
                  << ",\"max_abs_g16\":" << kernel.max_abs_g16
                  << ",\"elapsed_seconds\":" << seconds << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
