// Two prescribed cross-moments for the fixed N50 one-hole insertion.
// Adapted from execution bc17b81d:scripts/p337_endpoint_defect_exact.cpp.
// Only alternating four-neighbor configurations can change ambient rank.
// Track intact vacant components on that 1/8 subset, then restore the occupied
// A star algebraically. Subtract the paired S_minus*O_minus sums to obtain the
// missing correction to the existing full-population cross-moments.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

class RollbackComponents {
    std::array<int, 50> parent{}, size{};
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
        changed.emplace_back(b, size[a]); parent[b] = a; size[a] += size[b];
        --components;
    }
    void undo(std::size_t mark) {
        while (changed.size() > mark) {
            const auto [b, old_size] = changed.back(); changed.pop_back();
            size[parent[b]] = old_size; parent[b] = b; ++components;
        }
        --components;
    }
};

class Enumerator {
    static constexpr int n = 50, free_n = 25, fixed_occupied = 24;
    using Row = std::array<std::int64_t, 11>;
    std::array<Row, free_n+1> sums{};
    std::array<int, n> colour{}, free_position{};
    std::array<std::vector<int>, n> black_neighbors, white_neighbors;
    std::array<std::vector<std::array<int, 4>>, free_n> closing_faces;
    std::vector<int> free_vertices;
    RollbackComponents black, white, intact_white;
    std::vector<std::array<int, 4>> origin_faces;
    int a, b;
    static int mod(int x) { x %= n; return x < 0 ? x+n : x; }
    int key(int x, int y) const { return n*mod(a*x+b*y)+mod(-b*x+a*y); }

    void visit(int position, int k, int edges, int faces) {
        if (position == 4) {
            const int n0 = colour[free_vertices[0]], n1 = colour[free_vertices[1]];
            if (n0 == n1 || n0 != colour[free_vertices[2]] || n1 != colour[free_vertices[3]])
                return;
        }
        if (position == free_n) {
            const int total_k = fixed_occupied+k;
            const int q = black.components-white.components-(total_k-edges+faces);
            if (q < -1 || q > 1) throw std::logic_error("invalid digital parent rank value");
            const int e = q*q, bv = 2*n-4*total_k+edges;
            const int sstar = black.components+white.components+faces+bv;
            if (bv < 0) throw std::logic_error("negative vacant-edge count");
            int d = 0, c = 0, full_faces = 0;
            std::array<int, 4> contacts{};
            for (int u : black_neighbors[0]) if (colour[u] == 1) {
                ++d;
                const int component = black.root(u);
                bool seen = false;
                for (int j=0; j<c; ++j) seen |= contacts[j] == component;
                if (!seen) contacts[c++] = component;
            }
            for (const auto& face : origin_faces) {
                bool full = true;
                for (int u : face) if (u != 0 && colour[u] != 1) full = false;
                full_faces += full;
            }
            const int cb_plus = black.components + 1 - c;
            const int faces_plus = faces + full_faces;
            const int qp = cb_plus - intact_white.components
                         - (total_k + 1 - (edges+d) + faces_plus);
            const int ep = qp*qp;
            const int sp = cb_plus + intact_white.components + faces_plus;
            if (qp < -1 || qp > 1 || qp < q || qp-q > 1)
                throw std::logic_error("invalid paired rank restoration");
            const Row row{1, q, e, sstar, q*sstar, e*sstar,
                          qp, ep, sp, sstar*qp, sstar*ep};
            for (int j=0; j<11; ++j) sums[k][j] += row[j];
            return;
        }
        const int v = free_vertices[position];
        colour[v] = 0;
        auto mark = white.mark(); white.activate(v);
        auto intact_mark = intact_white.mark(); intact_white.activate(v);
        // Every fixed A neighbor is already active, regardless of its index.
        // Unvisited B neighbors have colour=-1; visited ones alone can join.
        for (int u : white_neighbors[v]) if (colour[u] == 0) {
            white.join(v, u);
            if (u != 0) intact_white.join(v, u);
        }
        visit(position+1, k, edges, faces);
        white.undo(mark); intact_white.undo(intact_mark);

        colour[v] = 1;
        mark = black.mark(); black.activate(v);
        int extra_edges = 0, extra_faces = 0;
        for (int u : black_neighbors[v]) if (colour[u] == 1) {
            black.join(v, u); ++extra_edges;
        }
        // A face is counted at its last free B activation, not max vertex ID.
        for (const auto& face : closing_faces[position]) {
            bool full = true;
            for (int u : face) if (colour[u] != 1) full = false;
            extra_faces += full;
        }
        visit(position+1, k+1, edges+extra_edges, faces+extra_faces);
        black.undo(mark); colour[v] = -1;
    }
public:
    Enumerator(int aa, int bb) : a(aa), b(bb) {
        if (!((a==5 && b==5) || (a==1 && b==7)))
            throw std::invalid_argument("fixed parent pair is (5,5),(1,7)");
        colour.fill(-1); free_position.fill(-1);
        std::array<int, n*n> index; index.fill(-1);
        std::vector<std::pair<int,int>> reps;
        reps.emplace_back(0,0); index[key(0,0)] = 0;
        for (std::size_t i=0; i<reps.size(); ++i) {
            auto [x,y] = reps[i];
            for (const auto& step : std::array<std::pair<int,int>,2>{{{1,0},{0,1}}}) {
                int xx=x+step.first, yy=y+step.second, h=key(xx,yy);
                if (index[h]<0) {
                    index[h]=static_cast<int>(reps.size()); reps.emplace_back(xx,yy);
                }
            }
        }
        if (reps.size()!=n) throw std::logic_error("wrong parent quotient area");
        int fixed_count = 0;
        for (int v=0; v<n; ++v) {
            const auto [x,y] = reps[v];
            if ((x+y)%2 == 0) {
                colour[v] = (v==0 ? 0 : 1);
                if (colour[v]) { black.activate(v); ++fixed_count; }
                else white.activate(v);
            } else {
                free_position[v] = static_cast<int>(free_vertices.size());
                free_vertices.push_back(v);
            }
        }
        if (fixed_count!=fixed_occupied || free_vertices.size()!=free_n || colour[0]!=0)
            throw std::logic_error("incorrect fixed-A/free-B partition");
        const auto all_free = free_vertices;
        free_vertices = {index[key(0,1)], index[key(1,0)], index[key(0,-1)], index[key(-1,0)]};
        for (int v : all_free)
            if (std::find(free_vertices.begin(), free_vertices.end(), v) == free_vertices.end())
                free_vertices.push_back(v);
        for (int j=0; j<free_n; ++j) free_position[free_vertices[j]] = j;
        // Fixed occupied A sites have no mutual NN edges; the sole fixed white
        // site is origin. Thus fixed components need no joins, edges or faces.
        for (int v=0; v<n; ++v) {
            auto [x,y] = reps[v];
            for (int dx=-1; dx<=1; ++dx) for (int dy=-1; dy<=1; ++dy) {
                if (dx==0 && dy==0) continue;
                int u=index[key(x+dx,y+dy)];
                white_neighbors[v].push_back(u);
                if (dx==0 || dy==0) black_neighbors[v].push_back(u);
            }
            std::array<int,4> face{v,index[key(x+1,y)],index[key(x,y+1)],index[key(x+1,y+1)]};
            int last_free=-1, free_corners=0;
            for (int u : face) if (free_position[u]>=0) {
                last_free=std::max(last_free,free_position[u]); ++free_corners;
            }
            if (free_corners!=2) throw std::logic_error("unit face does not have two free B corners");
            closing_faces[last_free].push_back(face);
            if (std::find(face.begin(), face.end(), 0) != face.end())
                origin_faces.push_back(face);
        }
    }

    void run(const char* path) {
        std::ifstream exists(path);
        if (exists.good()) throw std::runtime_error("output already exists");
        visit(0,0,0,0);
        std::int64_t count=0;
        for (const auto& row : sums) count += row[0];
        if (count!=(std::int64_t(1)<<(free_n-3))) throw std::logic_error("incomplete alternating-face population");
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot create output");
        out << "k,count,sum_q,sum_e,sum_sstar,sum_qsstar,sum_esstar,"
               "sum_qplus,sum_eplus,sum_splus,sum_sminus_qplus,sum_sminus_eplus\n";
        for (int k=0; k<=free_n; ++k) {
            out << k; for (auto value:sums[k]) out << ',' << value; out << '\n';
        }
        out.close(); if (!out) throw std::runtime_error("output write failed");
    }
};

int main(int argc, char** argv) {
    try {
        if (argc!=4) throw std::invalid_argument("usage: defect-cross-moments a b output.csv");
        const auto start=std::chrono::steady_clock::now();
        Enumerator(std::stoi(argv[1]),std::stoi(argv[2])).run(argv[3]);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "{\"status\":\"completed\",\"parent_N\":50,\"free_B_degree\":25,"
                     "\"fixed_A_occupied\":24,\"fixed_A_vacancy\":\"origin\","
                     "\"configurations\":4194304,\"elapsed_seconds\":" << seconds << "}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
