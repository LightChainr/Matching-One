#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct DSU {
    std::vector<int> parent;
    std::vector<unsigned char> rank;

    explicit DSU(int n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    }
};

struct Geometry {
    const char* id;
    int lambda_num;
    int lambda_den;
    int base_span;
    int x2_num;
    int x2_den;
};

constexpr std::array<Geometry, 4> kGeometries{{
    {"lambda_1_4", 1, 4, 15, 2, 5},
    {"lambda_1_3", 1, 3, 14, 1, 2},
    {"lambda_2_3", 2, 3, 15, 4, 5},
    {"lambda_3_4", 3, 4, 14, 6, 7},
}};

struct Edge {
    int first;
    int second;
};

struct Graph {
    int nx;
    int ny;
    int vertices;
    std::vector<Edge> edges;
    std::vector<std::vector<std::pair<int, int>>> adjacency;
    std::array<int, 4> terminals;
};

uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

bool open_bond(uint64_t seed, uint64_t sample, uint64_t edge) {
    return (splitmix64(seed ^ splitmix64(sample) ^ splitmix64(edge)) >> 63) != 0;
}

Graph build_graph(const Geometry& geometry, int level) {
    const int span = geometry.base_span * level;
    if ((span * geometry.x2_num) % geometry.x2_den != 0) {
        throw std::runtime_error("marked point is not integral");
    }
    Graph graph;
    graph.nx = 6 * span + 1;
    graph.ny = 4 * span + 1;
    graph.vertices = graph.nx * graph.ny;
    graph.adjacency.resize(graph.vertices);
    auto add_edge = [&](int first, int second) {
        const int index = static_cast<int>(graph.edges.size());
        graph.edges.push_back({first, second});
        graph.adjacency[first].push_back({index, second});
        graph.adjacency[second].push_back({index, first});
    };
    for (int y = 0; y < graph.ny; ++y) {
        for (int x = 0; x + 1 < graph.nx; ++x) {
            add_edge(y * graph.nx + x, y * graph.nx + x + 1);
        }
    }
    for (int y = 0; y + 1 < graph.ny; ++y) {
        for (int x = 0; x < graph.nx; ++x) {
            add_edge(y * graph.nx + x, (y + 1) * graph.nx + x);
        }
    }
    const int padding = 2 * span;
    graph.terminals = {{
        padding,
        padding + span * geometry.x2_num / geometry.x2_den,
        padding + span,
        padding + 2 * span,
    }};
    return graph;
}

bool forbidden_high_merge(DSU& dsu, const std::array<int, 4>& terminals) {
    const std::array<int, 2> first_group{{terminals[0], terminals[3]}};
    const std::array<int, 2> second_group{{terminals[1], terminals[2]}};
    for (int first : first_group) {
        for (int second : second_group) {
            if (dsu.find(first) == dsu.find(second)) return true;
        }
    }
    return false;
}

bool is_high_partition(DSU& dsu, const std::array<int, 4>& terminals) {
    const int first = dsu.find(terminals[0]);
    const int second = dsu.find(terminals[1]);
    const int third = dsu.find(terminals[2]);
    const int fourth = dsu.find(terminals[3]);
    return first == fourth && second == third && first != second;
}

struct Transcript {
    std::vector<int> edges;
    std::vector<unsigned char> values;
    bool high = false;
};

Transcript explore_high(const Graph& graph, uint64_t seed, uint64_t sample) {
    Transcript transcript;
    std::vector<unsigned char> revealed(graph.edges.size(), 0);
    std::vector<unsigned char> discovered(graph.vertices, 0);
    std::deque<int> frontier;
    for (int terminal : graph.terminals) {
        if (!discovered[terminal]) {
            discovered[terminal] = 1;
            frontier.push_back(terminal);
        }
    }
    DSU dsu(graph.vertices);
    while (!frontier.empty()) {
        const int vertex = frontier.front();
        frontier.pop_front();
        for (const auto& [edge_index, neighbor] : graph.adjacency[vertex]) {
            if (revealed[edge_index]) continue;
            revealed[edge_index] = 1;
            const bool value = open_bond(seed, sample, edge_index);
            transcript.edges.push_back(edge_index);
            transcript.values.push_back(static_cast<unsigned char>(value));
            if (!value) continue;
            dsu.unite(vertex, neighbor);
            if (!discovered[neighbor]) {
                discovered[neighbor] = 1;
                frontier.push_back(neighbor);
            }
            if (forbidden_high_merge(dsu, graph.terminals)) {
                transcript.high = false;
                return transcript;
            }
        }
    }
    transcript.high = is_high_partition(dsu, graph.terminals);
    return transcript;
}

int cluster_score_j(const Graph& graph, const std::vector<unsigned char>& bits) {
    if (bits.size() != graph.edges.size()) {
        throw std::runtime_error("configuration size mismatch");
    }
    DSU dsu(graph.vertices);
    int bonds = 0;
    for (size_t edge_index = 0; edge_index < graph.edges.size(); ++edge_index) {
        if (!bits[edge_index]) continue;
        ++bonds;
        dsu.unite(graph.edges[edge_index].first, graph.edges[edge_index].second);
    }
    int clusters = 0;
    for (int vertex = 0; vertex < graph.vertices; ++vertex) {
        if (dsu.find(vertex) == vertex) ++clusters;
    }
    return 2 * clusters + bonds;
}

struct Accumulator {
    uint64_t samples = 0;
    uint64_t count_high = 0;
    int64_t sum_delta_j_high = 0;
    uint64_t sum_delta_j_inner_square_high = 0;
    uint64_t sum_delta_j2_individual_high = 0;
    uint64_t sum_revealed = 0;
    uint64_t sum_revealed2 = 0;
    uint64_t sum_revealed_high = 0;
    uint64_t max_revealed = 0;
};

Accumulator run_batch(const Graph& graph, uint64_t begin, uint64_t count,
                      uint64_t outer_seed, uint64_t completion_seed,
                      int inner_replicates) {
    Accumulator result;
    for (uint64_t offset = 0; offset < count; ++offset) {
        const uint64_t sample = begin + offset;
        const Transcript transcript = explore_high(graph, outer_seed, sample);
        const uint64_t revealed = transcript.edges.size();
        ++result.samples;
        result.sum_revealed += revealed;
        result.sum_revealed2 += revealed * revealed;
        result.max_revealed = std::max(result.max_revealed, revealed);
        if (!transcript.high) continue;
        ++result.count_high;
        result.sum_revealed_high += revealed;
        int64_t inner_delta_sum = 0;
        for (int inner = 0; inner < inner_replicates; ++inner) {
            const uint64_t completion = sample * static_cast<uint64_t>(inner_replicates)
                                      + static_cast<uint64_t>(inner);
            std::vector<unsigned char> base(graph.edges.size(), 0);
            for (size_t edge = 0; edge < graph.edges.size(); ++edge) {
                base[edge] = static_cast<unsigned char>(
                    open_bond(completion_seed, completion, edge)
                );
            }
            std::vector<unsigned char> overwritten = base;
            for (size_t index = 0; index < transcript.edges.size(); ++index) {
                overwritten[transcript.edges[index]] = transcript.values[index];
            }
            const int delta_j = cluster_score_j(graph, overwritten)
                              - cluster_score_j(graph, base);
            if (std::abs(delta_j) > static_cast<int>(revealed)) {
                throw std::runtime_error("local score bound failed");
            }
            inner_delta_sum += delta_j;
            result.sum_delta_j2_individual_high += static_cast<uint64_t>(
                static_cast<int64_t>(delta_j) * delta_j
            );
        }
        result.sum_delta_j_high += inner_delta_sum;
        result.sum_delta_j_inner_square_high += static_cast<uint64_t>(
            inner_delta_sum * inner_delta_sum
        );
    }
    return result;
}

void write_header(std::ofstream& output) {
    output << "geometry_id,lambda_num,lambda_den,level,span_L,nx,ny,vertices,edges,"
              "batch,outer_seed,completion_seed,inner_replicates,sample_begin,samples,"
              "count_14_23,sum_delta_J_14_23,sum_delta_J_inner_square_14_23,"
              "sum_delta_J2_individual_14_23,sum_revealed,sum_revealed2,"
              "sum_revealed_14_23,max_revealed\n";
}

void write_row(std::ofstream& output, const Geometry& geometry, const Graph& graph,
               int level, int batch, uint64_t outer_seed, uint64_t completion_seed,
               int inner_replicates, uint64_t sample_begin,
               const Accumulator& row) {
    const int span = geometry.base_span * level;
    output << geometry.id << ',' << geometry.lambda_num << ',' << geometry.lambda_den
           << ',' << level << ',' << span << ',' << graph.nx << ',' << graph.ny
           << ',' << graph.vertices << ',' << graph.edges.size() << ',' << batch
           << ',' << outer_seed << ',' << completion_seed << ',' << inner_replicates
           << ',' << sample_begin << ',' << row.samples << ',' << row.count_high
           << ',' << row.sum_delta_j_high << ','
           << row.sum_delta_j_inner_square_high << ','
           << row.sum_delta_j2_individual_high << ',' << row.sum_revealed << ','
           << row.sum_revealed2 << ','
           << row.sum_revealed_high << ','
           << row.max_revealed << '\n';
}

int value_after(int argc, char** argv, const std::string& option, int fallback) {
    for (int index = 1; index + 1 < argc; ++index) {
        if (argv[index] == option) return std::stoi(argv[index + 1]);
    }
    return fallback;
}

uint64_t u64_after(int argc, char** argv, const std::string& option,
                   uint64_t fallback) {
    for (int index = 1; index + 1 < argc; ++index) {
        if (argv[index] == option) return std::stoull(argv[index + 1]);
    }
    return fallback;
}

std::string string_after(int argc, char** argv, const std::string& option) {
    for (int index = 1; index + 1 < argc; ++index) {
        if (argv[index] == option) return argv[index + 1];
    }
    throw std::runtime_error("missing required option " + option);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int level = value_after(argc, argv, "--level", 1);
        const int samples = value_after(argc, argv, "--samples", 0);
        const int batches = value_after(argc, argv, "--batches", 0);
        const int inner = value_after(argc, argv, "--inner", 0);
        const uint64_t outer_seed = u64_after(argc, argv, "--outer-seed", 0);
        const uint64_t completion_seed = u64_after(argc, argv, "--completion-seed", 0);
        const std::string output_path = string_after(argc, argv, "--output");
        if (level < 1 || samples < 1 || batches < 2 || inner < 1
            || samples % batches != 0 || outer_seed == completion_seed) {
            throw std::runtime_error(
                "require level>=1, samples>=1, batches>=2, inner>=1, "
                "samples%batches=0, and distinct seeds"
            );
        }
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("cannot open output");
        write_header(output);
        const uint64_t per_batch = static_cast<uint64_t>(samples / batches);
        for (const Geometry& geometry : kGeometries) {
            const Graph graph = build_graph(geometry, level);
            for (int batch = 0; batch < batches; ++batch) {
                const uint64_t begin = static_cast<uint64_t>(batch) * per_batch;
                const Accumulator row = run_batch(
                    graph, begin, per_batch, outer_seed, completion_seed, inner
                );
                write_row(
                    output, geometry, graph, level, batch, outer_seed,
                    completion_seed, inner, begin, row
                );
            }
        }
        std::cerr << "wrote " << output_path << " level=" << level
                  << " samples_per_geometry=" << samples
                  << " inner_replicates=" << inner << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
