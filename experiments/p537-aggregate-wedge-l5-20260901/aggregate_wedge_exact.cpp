// Exact sharded finite-population coefficients for the P537 landing wedge.
//
// Fix the thermal site z=0, enumerate the other L^2-1 occupations, and retain
// the same radius-one axis-minus-diagonal landing character used by PR #545.
// For every background we compute the complete canonical pair source in the
// z=0 and z=1 states.  The output is a compact polynomial table; no root or
// floating-point decision is taken in this producer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::uint32_t KEY_SPACE = 1U << 24;

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) out.push_back(field);
    return out;
}

struct Kernel {
    std::vector<std::int32_t> g16 = std::vector<std::int32_t>(KEY_SPACE, 0);
    std::size_t rows = 0;
};

Kernel read_kernel(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read kernel: " + path);
    Kernel kernel;
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
            if (key_column < 0 || g_column < 0) throw std::runtime_error("kernel columns missing");
            continue;
        }
        const auto key = std::stoll(fields.at(key_column));
        const auto value = std::stoll(fields.at(g_column));
        if (key < 0 || key >= KEY_SPACE) throw std::runtime_error("kernel key out of range");
        if (value < std::numeric_limits<std::int32_t>::min() ||
            value > std::numeric_limits<std::int32_t>::max())
            throw std::runtime_error("kernel value out of range");
        kernel.g16[std::size_t(key)] = std::int32_t(value);
        ++kernel.rows;
    }
    if (!kernel.rows) throw std::runtime_error("empty kernel");
    return kernel;
}

struct DSU {
    std::vector<int> parent, size;
    explicit DSU(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int root(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void join(int a, int b) {
        a = root(a); b = root(b);
        if (a == b) return;
        if (size[a] < size[b]) std::swap(a, b);
        parent[b] = a; size[a] += size[b];
    }
};

struct State {
    int q = 0;
    std::vector<int> occupied_root;
};

struct GlobalRow {
    std::uint64_t count = 0;
    std::int64_t sum_q = 0;
    std::int64_t sum_source16 = 0; // ordered-pair g16 sum
};

struct LandingRow {
    std::int64_t signed_count = 0;
    std::int64_t signed_source_mid16 = 0; // h*(source16(z=0)+source16(z=1))
    std::uint64_t unsigned_count = 0;
};

class Torus {
    int L_ = 0, N_ = 0;
    const Kernel& kernel_;
    std::vector<std::array<int, 4>> neighbor_; // N,E,S,W
    std::vector<std::array<int, 4>> edge_id_;

    int vertex(int x, int y) const {
        x %= L_; y %= L_;
        if (x < 0) x += L_;
        if (y < 0) y += L_;
        return x + L_ * y;
    }

    static int sector(int dx, int dy) {
        if (dx == 1 && dy == 0) return 0;
        if (dx == 1 && dy == 1) return 1;
        if (dx == 0 && dy == 1) return 2;
        if (dx == -1 && dy == 1) return 3;
        if (dx == -1 && dy == 0) return 4;
        if (dx == -1 && dy == -1) return 5;
        if (dx == 0 && dy == -1) return 6;
        if (dx == 1 && dy == -1) return 7;
        throw std::logic_error("invalid collar offset");
    }

    static bool pair_distinct(const std::vector<int>& masks, int a, int b) {
        for (std::size_t i = 0; i < masks.size(); ++i)
            for (std::size_t j = 0; j < masks.size(); ++j)
                if (i != j && (masks[i] & (1 << a)) && (masks[j] & (1 << b))) return true;
        return false;
    }

public:
    Torus(int L, const Kernel& kernel) : L_(L), N_(L * L), kernel_(kernel),
        neighbor_(N_), edge_id_(N_) {
        if (L < 4 || N_ > 128) throw std::invalid_argument("require L>=4 and L^2<=128");
        const std::array<std::pair<int,int>,4> step{{{0,1},{1,0},{0,-1},{-1,0}}};
        for (int y = 0; y < L_; ++y) for (int x = 0; x < L_; ++x) {
            const int v = vertex(x,y);
            for (int d = 0; d < 4; ++d) neighbor_[v][d] = vertex(x+step[d].first,y+step[d].second);
        }
        std::unordered_map<std::uint64_t,int> ids;
        int next = 0;
        for (int v = 0; v < N_; ++v) for (int d = 0; d < 4; ++d) {
            const int u = neighbor_[v][d];
            const int lo = std::min(u,v), hi = std::max(u,v);
            const std::uint64_t key = std::uint64_t(lo) * std::uint64_t(N_) + std::uint64_t(hi);
            auto [it, inserted] = ids.emplace(key, next);
            if (inserted) ++next;
            edge_id_[v][d] = it->second;
        }
        if (next != 2 * N_) throw std::logic_error("physical edge count mismatch");
    }

    int size() const { return N_; }

    State evaluate(const std::vector<unsigned char>& occupied) const {
        DSU black(N_), white(N_);
        int k = 0, edges = 0, faces = 0;
        for (int v = 0; v < N_; ++v) k += occupied[v];
        for (int v = 0; v < N_; ++v) {
            if (occupied[v]) {
                for (int d : {0,1}) if (occupied[neighbor_[v][d]]) {
                    ++edges; black.join(v, neighbor_[v][d]);
                }
                const int north = neighbor_[v][0], east = neighbor_[v][1];
                const int northeast = neighbor_[north][1];
                faces += occupied[north] && occupied[east] && occupied[northeast];
            } else {
                // Four forward representatives of the matching-lattice edges.
                for (int d : {0,1}) if (!occupied[neighbor_[v][d]]) white.join(v, neighbor_[v][d]);
                const int north = neighbor_[v][0];
                const int south = neighbor_[v][2];
                const int northeast = neighbor_[north][1];
                const int southeast = neighbor_[south][1];
                if (!occupied[northeast]) white.join(v, northeast);
                if (!occupied[southeast]) white.join(v, southeast);
            }
        }
        std::vector<unsigned char> black_seen(N_,0), white_seen(N_,0);
        int cb = 0, cw = 0;
        State state;
        state.occupied_root.assign(N_,-1);
        for (int v = 0; v < N_; ++v) {
            if (occupied[v]) {
                const int r = black.root(v);
                state.occupied_root[v] = r;
                if (!black_seen[r]) { black_seen[r] = 1; ++cb; }
            } else {
                const int r = white.root(v);
                if (!white_seen[r]) { white_seen[r] = 1; ++cw; }
            }
        }
        state.q = cb - cw - (k - edges + faces);
        if (state.q < -1 || state.q > 1) throw std::logic_error("digital rank outside {-1,0,1}");
        return state;
    }

    int landing_h4(const std::vector<unsigned char>& occupied) const {
        struct Node { int dx, dy, v; };
        std::vector<Node> nodes;
        for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx)
            if (dx || dy) nodes.push_back({dx,dy,vertex(dx,dy)});
        auto component_masks = [&](bool on, bool matching) {
            std::vector<unsigned char> seen(nodes.size(),0);
            std::vector<int> masks;
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                if (seen[i] || bool(occupied[nodes[i].v]) != on) continue;
                std::vector<std::size_t> stack{i}; seen[i] = 1;
                int mask = 0; bool touches_port = false;
                while (!stack.empty()) {
                    const auto j = stack.back(); stack.pop_back();
                    const auto& a = nodes[j];
                    mask |= 1 << sector(a.dx,a.dy);
                    touches_port |= matching || (std::abs(a.dx)+std::abs(a.dy)==1);
                    for (std::size_t k = 0; k < nodes.size(); ++k) {
                        if (seen[k] || bool(occupied[nodes[k].v]) != on) continue;
                        const int ax = std::abs(a.dx-nodes[k].dx), ay = std::abs(a.dy-nodes[k].dy);
                        const bool adjacent = matching ? (std::max(ax,ay)==1) : (ax+ay==1);
                        if (adjacent) { seen[k] = 1; stack.push_back(k); }
                    }
                }
                if (touches_port) masks.push_back(mask);
            }
            return masks;
        };
        const auto open = component_masks(true,false);
        const auto closed = component_masks(false,true);
        const bool axis = (pair_distinct(open,0,4) && pair_distinct(closed,2,6)) ||
                          (pair_distinct(open,2,6) && pair_distinct(closed,0,4));
        const bool diagonal = (pair_distinct(open,1,5) && pair_distinct(closed,3,7)) ||
                              (pair_distinct(open,3,7) && pair_distinct(closed,1,5));
        return int(axis) - int(diagonal);
    }

    std::int64_t source16(const std::vector<unsigned char>& occupied, const State& state) const {
        std::int64_t unordered = 0;
        std::array<int,8> ids{}, labels{};
        for (int x = 0; x < N_; ++x) if (!occupied[x]) {
            for (int y = x + 1; y < N_; ++y) if (!occupied[y]) {
                for (int d = 0; d < 4; ++d) {
                    const int ux = neighbor_[x][d], uy = neighbor_[y][d];
                    ids[d] = occupied[ux] ? state.occupied_root[ux] : N_ + edge_id_[x][d];
                    ids[4+d] = occupied[uy] ? state.occupied_root[uy] : N_ + edge_id_[y][d];
                }
                int next = 0;
                std::uint32_t key = 0;
                for (int i = 0; i < 8; ++i) {
                    int label = -1;
                    for (int j = 0; j < i; ++j) if (ids[j] == ids[i]) { label = labels[j]; break; }
                    if (label < 0) label = next++;
                    labels[i] = label;
                    key |= std::uint32_t(label) << (3*i);
                }
                unordered += kernel_.g16[key];
            }
        }
        return 2 * unordered;
    }

    // Ordered contribution from one uniformly sampled first endpoint.  N times
    // this quantity is an unbiased estimator of source16() for a fixed state.
    std::int64_t source16_origin(const std::vector<unsigned char>& occupied,
                                 const State& state, int x) const {
        if (x < 0 || x >= N_) throw std::out_of_range("source origin");
        if (occupied[x]) return 0;
        std::int64_t total = 0;
        std::array<int,8> ids{}, labels{};
        for (int y = 0; y < N_; ++y) if (y != x && !occupied[y]) {
            for (int d = 0; d < 4; ++d) {
                const int ux = neighbor_[x][d], uy = neighbor_[y][d];
                ids[d] = occupied[ux] ? state.occupied_root[ux] : N_ + edge_id_[x][d];
                ids[4+d] = occupied[uy] ? state.occupied_root[uy] : N_ + edge_id_[y][d];
            }
            int next = 0;
            std::uint32_t key = 0;
            for (int i = 0; i < 8; ++i) {
                int label = -1;
                for (int j = 0; j < i; ++j) if (ids[j] == ids[i]) { label = labels[j]; break; }
                if (label < 0) label = next++;
                labels[i] = label;
                key |= std::uint32_t(label) << (3*i);
            }
            total += kernel_.g16[key];
        }
        return total;
    }
};

void write_output(const std::string& path, int L, int shard_index, int shard_count,
                  std::uint64_t begin, std::uint64_t end,
                  const std::vector<GlobalRow>& global,
                  const std::array<std::vector<LandingRow>,2>& landing) {
    std::ifstream probe(path);
    if (probe.good()) throw std::runtime_error("refusing to overwrite: " + path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open output: " + path);
    out << "# schema=matching-one/p537-aggregate-wedge-coefficients/v1\n";
    out << "# L=" << L << "\n# shard_index=" << shard_index << "\n# shard_count=" << shard_count
        << "\n# begin=" << begin << "\n# end=" << end << "\n";
    out << "kind\ttransition\tk\tcount\tsum_q\tsum_source16\tsigned_count\tsigned_source_mid16\tunsigned_count\n";
    for (std::size_t k = 0; k < global.size(); ++k)
        out << "global\t-\t" << k << '\t' << global[k].count << '\t' << global[k].sum_q
            << '\t' << global[k].sum_source16 << "\t0\t0\t0\n";
    for (int tr = 0; tr < 2; ++tr) for (std::size_t k = 0; k < landing[tr].size(); ++k)
        out << "landing\t" << (tr ? "12" : "01") << '\t' << k << "\t0\t0\t0\t"
            << landing[tr][k].signed_count << '\t' << landing[tr][k].signed_source_mid16
            << '\t' << landing[tr][k].unsigned_count << '\n';
}

} // namespace

#ifndef P537_LIBRARY_ONLY
int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: aggregate_wedge_exact KERNEL OUTPUT L SHARD_INDEX SHARD_COUNT\n";
        return 2;
    }
    try {
        const auto kernel = read_kernel(argv[1]);
        const std::string output = argv[2];
        const int L = std::stoi(argv[3]);
        const int shard_index = std::stoi(argv[4]);
        const int shard_count = std::stoi(argv[5]);
        Torus torus(L,kernel);
        const int N = torus.size();
        if (L > 5) throw std::invalid_argument("exact enumeration is bounded to L=4 or L=5");
        if (shard_count <= 0 || shard_index < 0 || shard_index >= shard_count)
            throw std::invalid_argument("invalid shard");
        const std::uint64_t total = 1ULL << (N-1);
        const std::uint64_t begin = total * std::uint64_t(shard_index) / std::uint64_t(shard_count);
        const std::uint64_t end = total * std::uint64_t(shard_index+1) / std::uint64_t(shard_count);
        std::vector<GlobalRow> global(N+1);
        std::array<std::vector<LandingRow>,2> landing{{std::vector<LandingRow>(N),std::vector<LandingRow>(N)}};
        std::vector<unsigned char> occupied(N,0);
        for (std::uint64_t mask = begin; mask < end; ++mask) {
            int k = 0;
            for (int v = 1; v < N; ++v) { occupied[v] = (mask >> (v-1)) & 1ULL; k += occupied[v]; }
            occupied[0] = 0;
            const State state0 = torus.evaluate(occupied);
            const int h4 = torus.landing_h4(occupied);
            const auto source0 = torus.source16(occupied,state0);
            occupied[0] = 1;
            const State state1 = torus.evaluate(occupied);
            const auto source1 = torus.source16(occupied,state1);
            occupied[0] = 0;

            ++global[k].count; global[k].sum_q += state0.q; global[k].sum_source16 += source0;
            ++global[k+1].count; global[k+1].sum_q += state1.q; global[k+1].sum_source16 += source1;

            const int r0 = state0.q + 1, r1 = state1.q + 1;
            int tr = -1;
            if (r0 == 0 && r1 == 1) tr = 0;
            if (r0 == 1 && r1 == 2) tr = 1;
            if (tr >= 0 && h4) {
                landing[tr][k].signed_count += h4;
                landing[tr][k].signed_source_mid16 += std::int64_t(h4) * (source0 + source1);
                ++landing[tr][k].unsigned_count;
            }
        }
        write_output(output,L,shard_index,shard_count,begin,end,global,landing);
        std::cerr << "completed masks=" << (end-begin) << " kernel_rows=" << kernel.rows << '\n';
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
#endif
