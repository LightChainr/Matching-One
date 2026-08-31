// Fresh iid spatial readout of the prescribed regular-pair Q-jet lookup.
// Preparation only: execution requires the root's frozen contract and GO.
// No homology/rank, q/E, single-source response, or old archive is evaluated.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Options {
    int L = 0, batches = 0, samples_per_batch = 0;
    double p = -1;
    std::uint64_t seed = 0;
    bool seed_given = false;
    std::string lookup, output;
};

Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 == argc) throw std::invalid_argument("every option needs a value");
        const std::string key = argv[i], value = argv[i+1];
        if (key == "--L") o.L = std::stoi(value);
        else if (key == "--p") o.p = std::stod(value);
        else if (key == "--seed") { o.seed = std::stoull(value); o.seed_given = true; }
        else if (key == "--batches") o.batches = std::stoi(value);
        else if (key == "--samples-per-batch") o.samples_per_batch = std::stoi(value);
        else if (key == "--lookup") o.lookup = value;
        else if (key == "--output") o.output = value;
        else throw std::invalid_argument("unknown option: " + key);
    }
    if ((o.L != 16 && o.L != 32) || !(o.p > 0 && o.p < 1) || !o.seed_given ||
        o.batches <= 0 || o.samples_per_batch <= 0 || o.lookup.empty() || o.output.empty())
        throw std::invalid_argument("required: --L {16|32} --p P --seed U64 --batches B "
                                    "--samples-per-batch M --lookup kernel.tsv --output batches.csv");
    return o;
}

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

std::unordered_map<std::uint32_t, std::int64_t> read_lookup(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read lookup: " + path);
    std::unordered_map<std::uint32_t, std::int64_t> table;
    std::string line;
    int key_column = -1, g_column = -1;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tsv(line);
        if (key_column < 0) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "key") key_column = static_cast<int>(i);
                if (fields[i] == "g16") g_column = static_cast<int>(i);
            }
            if (key_column < 0 || g_column < 0)
                throw std::runtime_error("lookup TSV needs named key and g16 columns");
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max(key_column, g_column)))
            throw std::runtime_error("incomplete lookup row");
        const auto wide_key = std::stoull(fields[key_column]);
        if (wide_key >= (1ULL << 24)) throw std::runtime_error("lookup key exceeds eight 3-bit labels");
        const auto key = static_cast<std::uint32_t>(wide_key);
        int max_label = -1;
        for (int i = 0; i < 8; ++i) {
            const int label = (key >> (3*i)) & 7;
            if (label > max_label + 1) throw std::runtime_error("lookup key is not canonical restricted-growth order");
            max_label = std::max(max_label, label);
        }
        if (!table.emplace(key, std::stoll(fields[g_column])).second)
            throw std::runtime_error("duplicate lookup key");
    }
    if (table.empty()) throw std::runtime_error("empty lookup");
    return table;
}

class Components {
    std::vector<int> parent, size;
public:
    explicit Components(int n) : parent(n), size(n) {}
    void reset() {
        std::iota(parent.begin(), parent.end(), 0);
        std::fill(size.begin(), size.end(), 1);
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

struct Batch {
    std::int64_t sum_g16 = 0;
    std::array<std::int64_t,5> by_shared{};
    std::array<std::uint64_t,5> pairs_by_shared{};
    std::uint64_t eligible_pairs = 0, nonzero_pairs = 0;
};

class Sampler {
    const Options& o;
    const std::unordered_map<std::uint32_t,std::int64_t>& lookup;
    const int n;
    std::vector<std::array<int,4>> neighbors; // N,E,S,W
    std::array<std::pair<int,int>,32> pairs{};
    std::vector<unsigned char> occupied;
    Components components;
    std::mt19937_64 rng;
    const std::uint64_t threshold;
    int site(int x, int y) const { return ((y+o.L)%o.L)*o.L + (x+o.L)%o.L; }

    void configuration(Batch& batch) {
        // One standardized RNG word per site, row-major x-fast order.
        for (int v = 0; v < n; ++v) occupied[v] = (rng() >> 11) < threshold;
        components.reset();
        for (int v = 0; v < n; ++v) if (occupied[v]) {
            for (int direction : {0,1}) {
                const int u = neighbors[v][direction];
                if (occupied[u]) components.join(v,u);
            }
        }
        for (const auto& pair : pairs) {
            if (occupied[pair.first] || occupied[pair.second]) continue;
            std::array<int,8> outside{}, labels{};
            std::uint32_t key = 0;
            int next_label = 0;
            for (int i = 0; i < 8; ++i) {
                const int center = i < 4 ? pair.first : pair.second;
                const int neighbor = neighbors[center][i%4];
                // All eight incident physical edges are distinct (r>=4).
                // A vacant neighbor leaves this edge-node its own singleton.
                outside[i] = occupied[neighbor] ? components.root(neighbor) : n+i;
                int previous = 0;
                while (previous < i && outside[previous] != outside[i]) ++previous;
                labels[i] = previous < i ? labels[previous] : next_label++;
                key |= static_cast<std::uint32_t>(labels[i]) << (3*i);
            }
            int shared = 0;
            for (int label = 0; label < next_label; ++label) {
                const bool left = std::find(labels.begin(),labels.begin()+4,label) != labels.begin()+4;
                const bool right = std::find(labels.begin()+4,labels.end(),label) != labels.end();
                shared += left && right;
            }
            const auto found = lookup.find(key);
            if (found == lookup.end())
                throw std::runtime_error("missing exact g16 lookup key: " + std::to_string(key));
            const auto g16 = found->second;
            batch.sum_g16 += g16;
            batch.by_shared[shared] += g16;
            ++batch.pairs_by_shared[shared];
            ++batch.eligible_pairs;
            batch.nonzero_pairs += g16 != 0;
        }
    }
public:
    Sampler(const Options& options, const std::unordered_map<std::uint32_t,std::int64_t>& table)
        : o(options), lookup(table), n(o.L*o.L), neighbors(n), occupied(n),
          components(n), rng(o.seed),
          threshold(static_cast<std::uint64_t>(o.p*9007199254740992.0)) {
        for (int y = 0; y < o.L; ++y) for (int x = 0; x < o.L; ++x)
            neighbors[site(x,y)] = {site(x,y+1),site(x+1,y),site(x,y-1),site(x-1,y)};
        const int r = o.L/4;
        int jpair = 0;
        for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
            const int x = i*r, y = j*r;
            pairs[jpair++] = {site(x,y),site(x+r,y)};
            pairs[jpair++] = {site(x,y),site(x,y+r)};
        }
    }
    void run() {
        const auto start = std::chrono::steady_clock::now();
        const std::string metadata = o.output + ".metadata.json";
        if (std::ifstream(o.output).good() || std::ifstream(metadata).good())
            throw std::runtime_error("output/metadata exists; refusing overwrite");
        std::ofstream out(o.output);
        if (!out) throw std::runtime_error("cannot create output");
        out << "L,batch,samples,pairs,eligible_pairs,nonzero_pairs,sum_g16";
        for (int c = 0; c <= 4; ++c) out << ",sum_g16_shared" << c;
        for (int c = 0; c <= 4; ++c) out << ",pairs_shared" << c;
        out << '\n';
        for (int b = 0; b < o.batches; ++b) {
            Batch batch;
            for (int s = 0; s < o.samples_per_batch; ++s) configuration(batch);
            out << o.L << ',' << b << ',' << o.samples_per_batch << ','
                << 32ULL*o.samples_per_batch << ',' << batch.eligible_pairs << ','
                << batch.nonzero_pairs << ',' << batch.sum_g16;
            for (auto value : batch.by_shared) out << ',' << value;
            for (auto value : batch.pairs_by_shared) out << ',' << value;
            out << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("output write failed");
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream meta(metadata);
        if (!meta) throw std::runtime_error("cannot create metadata");
        meta << std::setprecision(17)
             << "{\n  \"status\": \"completed\",\n  \"schema\": \"p337.regular-pair-spatial-batches.v1\",\n"
             << "  \"L\": " << o.L << ",\n  \"N\": " << n << ",\n  \"r\": " << o.L/4
             << ",\n  \"p_requested\": " << o.p << ",\n  \"bernoulli_threshold_2pow53\": " << threshold
             << ",\n  \"p_implemented\": " << static_cast<double>(threshold)/9007199254740992.0
             << ",\n  \"seed\": " << o.seed << ",\n  \"rng\": \"std::mt19937_64; (word>>11)<threshold; one word/site\",\n"
             << "  \"site_order\": \"x-fast row-major; N=(0,+1), E=(+1,0), S=(0,-1), W=(-1,0)\",\n"
             << "  \"batches\": " << o.batches << ",\n  \"samples_per_batch\": " << o.samples_per_batch
             << ",\n  \"samples\": " << static_cast<std::uint64_t>(o.batches)*o.samples_per_batch
             << ",\n  \"pairs_per_configuration\": 32,\n  \"lookup_rows\": " << lookup.size()
             << ",\n  \"elapsed_seconds\": " << seconds
             << ",\n  \"inference_unit\": \"iid configuration; all 32 pairs are correlated readouts\",\n"
             << "  \"mean_g_denominator\": \"16*32*samples\",\n"
             << "  \"shared_counts_condition\": \"both endpoints vacant; occupied-endpoint pairs contribute zero\",\n"
             << "  \"rank_or_q_or_E_evaluated\": false\n}\n";
        meta.close();
        if (!meta) throw std::runtime_error("metadata write failed");
        std::cout << "{\"status\":\"completed\",\"samples\":"
                  << static_cast<std::uint64_t>(o.batches)*o.samples_per_batch
                  << ",\"elapsed_seconds\":" << seconds << "}\n";
    }
};

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc,argv);
        const auto lookup = read_lookup(options.lookup);
        Sampler(options,lookup).run();
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
