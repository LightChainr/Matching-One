#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
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
    int k_num;
    int k_den;
};

constexpr std::array<Geometry, 4> kGeometries{{
    {"lambda_1_4", 1, 4, 15, 2, 5, 10, 3},
    {"lambda_1_3", 1, 3, 14, 1, 2, 3, 1},
    {"lambda_2_3", 2, 3, 15, 4, 5, 15, 4},
    {"lambda_3_4", 3, 4, 14, 6, 7, 14, 3},
}};

uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

bool open_bond(uint64_t seed, uint64_t sample, uint64_t edge) {
    // Geometry is intentionally absent: equal sample/edge counters are the
    // declared common-random-number coupling across the four rectangles.
    return (splitmix64(seed ^ splitmix64(sample) ^ splitmix64(edge)) >> 63) != 0;
}

struct Accumulator {
    uint64_t samples = 0;
    int64_t sum_j = 0;
    int64_t sum_j2 = 0;
    std::array<int64_t, 3> count{{0, 0, 0}};
    std::array<int64_t, 3> sum_j_channel{{0, 0, 0}};
    std::array<int64_t, 3> sum_j2_channel{{0, 0, 0}};
};

int channel(DSU& dsu, const std::array<int, 4>& terminals) {
    std::array<int, 4> root{};
    for (int i = 0; i < 4; ++i) root[i] = dsu.find(terminals[i]);
    if (root[0] == root[1] && root[0] == root[2] && root[0] == root[3]) return 0;
    if (root[0] == root[1] && root[2] == root[3] && root[0] != root[2]) return 1;
    if (root[0] == root[3] && root[1] == root[2] && root[0] != root[1]) return 2;
    return -1;
}

Accumulator run_batch(const Geometry& geometry, int level, uint64_t begin,
                      uint64_t count, uint64_t seed) {
    const int span = geometry.base_span * level;
    if ((span * geometry.x2_num) % geometry.x2_den != 0) {
        throw std::runtime_error("marked point is not integral");
    }
    const int nx = 6 * span + 1;
    const int ny = 4 * span + 1;
    const int vertices = nx * ny;
    const int horizontal_edges = (nx - 1) * ny;
    const int total_edges = horizontal_edges + nx * (ny - 1);
    const int padding = 2 * span;
    const std::array<int, 4> terminals{{
        padding,
        padding + span * geometry.x2_num / geometry.x2_den,
        padding + span,
        padding + 2 * span,
    }};

    Accumulator result;
    for (uint64_t offset = 0; offset < count; ++offset) {
        const uint64_t sample = begin + offset;
        DSU dsu(vertices);
        int bonds = 0;
        uint64_t edge_index = 0;
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x + 1 < nx; ++x, ++edge_index) {
                if (open_bond(seed, sample, edge_index)) {
                    ++bonds;
                    dsu.unite(y * nx + x, y * nx + x + 1);
                }
            }
        }
        for (int y = 0; y + 1 < ny; ++y) {
            for (int x = 0; x < nx; ++x, ++edge_index) {
                if (open_bond(seed, sample, edge_index)) {
                    ++bonds;
                    dsu.unite(y * nx + x, (y + 1) * nx + x);
                }
            }
        }
        if (static_cast<int>(edge_index) != total_edges) {
            throw std::runtime_error("edge enumeration mismatch");
        }
        int clusters = 0;
        for (int vertex = 0; vertex < vertices; ++vertex) {
            if (dsu.find(vertex) == vertex) ++clusters;
        }
        const int j = 2 * clusters + bonds;
        const int selected = channel(dsu, terminals);
        ++result.samples;
        result.sum_j += j;
        result.sum_j2 += static_cast<int64_t>(j) * j;
        if (selected >= 0) {
            ++result.count[selected];
            result.sum_j_channel[selected] += j;
            result.sum_j2_channel[selected] += static_cast<int64_t>(j) * j;
        }
    }
    return result;
}

void write_header(std::ofstream& output) {
    output << "geometry_id,lambda_num,lambda_den,level,span_L,nx,ny,vertices,edges,"
              "batch,seed,sample_begin,samples,sum_J,sum_J2,"
              "count_1234,sum_J_1234,sum_J2_1234,"
              "count_12_34,sum_J_12_34,sum_J2_12_34,"
              "count_14_23,sum_J_14_23,sum_J2_14_23\n";
}

void write_row(std::ofstream& output, const Geometry& geometry, int level,
               int batch, uint64_t seed, uint64_t sample_begin,
               const Accumulator& row) {
    const int span = geometry.base_span * level;
    const int nx = 6 * span + 1;
    const int ny = 4 * span + 1;
    const int vertices = nx * ny;
    const int edges = (nx - 1) * ny + nx * (ny - 1);
    output << geometry.id << ',' << geometry.lambda_num << ',' << geometry.lambda_den
           << ',' << level << ',' << span << ',' << nx << ',' << ny << ','
           << vertices << ',' << edges << ',' << batch << ',' << seed << ','
           << sample_begin << ',' << row.samples << ','
           << row.sum_j << ',' << row.sum_j2;
    for (int index = 0; index < 3; ++index) {
        output << ',' << row.count[index] << ',' << row.sum_j_channel[index] << ','
               << row.sum_j2_channel[index];
    }
    output << '\n';
}

int value_after(int argc, char** argv, const std::string& option, int fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == option) return std::stoi(argv[i + 1]);
    }
    return fallback;
}

uint64_t u64_after(int argc, char** argv, const std::string& option, uint64_t fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == option) return std::stoull(argv[i + 1]);
    }
    return fallback;
}

std::string string_after(int argc, char** argv, const std::string& option) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (argv[i] == option) return argv[i + 1];
    }
    throw std::runtime_error("missing required option " + option);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int level = value_after(argc, argv, "--level", 1);
        const int samples = value_after(argc, argv, "--samples", 0);
        const int batches = value_after(argc, argv, "--batches", 0);
        const uint64_t seed = u64_after(argc, argv, "--seed", 202608290263ULL);
        const std::string output_path = string_after(argc, argv, "--output");
        if (level < 1 || samples < 1 || batches < 1 || samples % batches != 0) {
            throw std::runtime_error("require level>=1, samples>=1, batches>=1 and samples%batches=0");
        }
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("cannot open output");
        write_header(output);
        const uint64_t per_batch = static_cast<uint64_t>(samples / batches);
        for (const Geometry& geometry : kGeometries) {
            for (int batch = 0; batch < batches; ++batch) {
                const uint64_t begin = static_cast<uint64_t>(batch) * per_batch;
                const Accumulator row = run_batch(
                    geometry, level, begin, per_batch, seed
                );
                write_row(output, geometry, level, batch, seed, begin, row);
            }
        }
        std::cerr << "wrote " << output_path << " for level=" << level
                  << " samples_per_geometry=" << samples << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
