// Paired intact/one-hole N50 observer on exactly the same 25 free B sites.
// Quotient and rollback backend retained from 359bde9b endpoint-defect producer.
// Save only the two missing S_defect * O_intact cross products, not old scores.
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
    using Row = std::array<std::int64_t, 3>;
    std::array<Row, free_n+1> sums{};
    std::array<int, n> colour{}, free_position{};
    std::array<std::vector<int>, n> black_neighbors, white_neighbors;
    std::array<std::vector<std::array<int, 4>>, free_n> closing_faces;
    std::vector<int> free_vertices;
    RollbackComponents black, white, black_intact, white_intact;
    int a, b;
    static int mod(int x) { x %= n; return x < 0 ? x+n : x; }
    int key(int x, int y) const { return n*mod(a*x+b*y)+mod(-b*x+a*y); }

    void visit(int position, int k, int edges, int faces, int edges_intact, int faces_intact) {
        if (position == free_n) {
            const int total_k = fixed_occupied+k;
            const int q_defect = black.components-white.components-(total_k-edges+faces);
            const int q_intact = black_intact.components-white_intact.components
                                 -(total_k+1-edges_intact+faces_intact);
            if (q_defect < -1 || q_defect > 1 || q_intact < -1 || q_intact > 1)
                throw std::logic_error("invalid paired digital parent rank value");
            const int bv = 2*n-4*total_k+edges;
            const int s_defect = black.components+white.components+faces+bv;
            if (bv < 0) throw std::logic_error("negative defect vacant-edge count");
            const Row row{1, s_defect*q_intact, s_defect*q_intact*q_intact};
            for (int j=0; j<3; ++j) sums[k][j] += row[j];
            return;
        }
        const int v = free_vertices[position];
        colour[v] = 0;
        auto mark = white.mark(), mark_intact = white_intact.mark();
        white.activate(v); white_intact.activate(v);
        // All fixed A neighbors participate independent of their vertex IDs.
        // The shared colour array uses defect colours: origin alone differs.
        for (int u : white_neighbors[v]) if (colour[u] == 0) {
            white.join(v, u);
            if (u != 0) white_intact.join(v, u);
        }
        visit(position+1, k, edges, faces, edges_intact, faces_intact);
        white.undo(mark); white_intact.undo(mark_intact);

        colour[v] = 1;
        mark = black.mark(); mark_intact = black_intact.mark();
        black.activate(v); black_intact.activate(v);
        int extra_edges = 0, extra_faces = 0;
        int extra_edges_intact = 0, extra_faces_intact = 0;
        for (int u : black_neighbors[v]) {
            if (colour[u] == 1) { black.join(v, u); ++extra_edges; }
            if (u == 0 || colour[u] == 1) {
                black_intact.join(v, u); ++extra_edges_intact;
            }
        }
        for (const auto& face : closing_faces[position]) {
            bool full = true, full_intact = true;
            for (int u : face) {
                if (colour[u] != 1) full = false;
                if (u != 0 && colour[u] != 1) full_intact = false;
            }
            extra_faces += full; extra_faces_intact += full_intact;
        }
        visit(position+1, k+1, edges+extra_edges, faces+extra_faces,
              edges_intact+extra_edges_intact, faces_intact+extra_faces_intact);
        black.undo(mark); black_intact.undo(mark_intact); colour[v] = -1;
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
                black_intact.activate(v);  // every intact A site, including origin
            } else {
                free_position[v] = static_cast<int>(free_vertices.size());
                free_vertices.push_back(v);
            }
        }
        if (fixed_count!=fixed_occupied || free_vertices.size()!=free_n || colour[0]!=0)
            throw std::logic_error("incorrect fixed-A/free-B partition");
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
        }
    }

    void run(const char* path) {
        std::ifstream exists(path);
        if (exists.good()) throw std::runtime_error("output already exists");
        visit(0,0,0,0,0,0);
        std::int64_t count=0;
        for (const auto& row : sums) count += row[0];
        if (count!=(std::int64_t(1)<<free_n)) throw std::logic_error("incomplete configuration population");
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot create output");
        out << "k,count,sum_Sdef_qintact,sum_Sdef_Eintact\n";
        for (int k=0; k<=free_n; ++k) {
            out << k; for (auto value:sums[k]) out << ',' << value; out << '\n';
        }
        out.close(); if (!out) throw std::runtime_error("output write failed");
    }
};

int main(int argc, char** argv) {
    try {
        if (argc!=4) throw std::invalid_argument("usage: defect-reweighting-cross a b output.csv");
        const auto start=std::chrono::steady_clock::now();
        Enumerator(std::stoi(argv[1]),std::stoi(argv[2])).run(argv[3]);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "{\"status\":\"completed\",\"parent_N\":50,\"free_B_degree\":25,"
                     "\"fixed_A_occupied\":24,\"fixed_A_vacancy\":\"origin\","
                     "\"configurations\":33554432,\"elapsed_seconds\":" << seconds << "}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
