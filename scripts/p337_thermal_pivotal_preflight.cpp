// Deterministic thermal/pivotal preflight for the already-produced P337
// L32/r8 and L64/r16 occupation streams.
//
// This program deliberately replays only counters 0..31 from each frozen
// seed.  For every replayed configuration and every site z it evaluates the
// two forced states z=0 and z=1, then records pair/site midpoint primitives
// for the same 32 anchor/direction pairs used by the original producer.
// Pair/site absolute values are taken before aggregation.  No RNG seed or
// scientific sample count is accepted on the command line.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sys/resource.h>

namespace {

using Kernel = std::vector<std::int16_t>;

constexpr std::uint64_t bernoulli_threshold = 10934234699625173385ULL;
constexpr std::uint64_t seed_L32 = 2026083123593201ULL;
constexpr std::uint64_t seed_L64 = 2026083123596401ULL;
constexpr int frozen_configurations = 32;
constexpr int frozen_pairs = 32;
constexpr std::size_t kernel_size = 1U << 24;

enum CarrierBit : int {
    two_bridge_persistent_bit = 1,
    shared_transition_or_merger_bit = 2,
    kernel_preserving_topological_bit = 4,
    kernel_changed_bit = 8
};

struct Options {
    std::string kernel_path;
    std::string output_prefix;
    int threads = 0;
};

Options parse_options(int argc, char** argv) {
    Options result;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc) throw std::invalid_argument("every option needs a value");
        const std::string key = argv[i], value = argv[i + 1];
        if (key == "--kernel") result.kernel_path = value;
        else if (key == "--output-prefix") result.output_prefix = value;
        else if (key == "--threads") result.threads = std::stoi(value);
        else throw std::invalid_argument("unknown option: " + key);
    }
    if (result.kernel_path.empty() || result.output_prefix.empty() ||
        result.threads < 1 || result.threads > 64) {
        throw std::invalid_argument(
            "usage: p337_thermal_pivotal_preflight --kernel kernel.tsv "
            "--output-prefix PATH --threads 1..64");
    }
    return result;
}

std::uint64_t unsigned_integer(const std::string& text) {
    if (text.empty() || text.front() == '-')
        throw std::invalid_argument("invalid unsigned integer: " + text);
    std::size_t end = 0;
    const auto value = std::stoull(text, &end);
    if (end != text.size()) throw std::invalid_argument("invalid unsigned integer: " + text);
    return value;
}

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

Kernel read_kernel(const std::string& path, std::size_t& rows) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot read kernel: " + path);
    Kernel kernel(kernel_size, 0); // Sparse omissions are exact zeros.
    std::unordered_set<std::uint32_t> seen;
    std::string line;
    int key_column = -1, value_column = -1;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tsv(line);
        if (key_column < 0) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "key" || fields[i] == "packed_key") key_column = static_cast<int>(i);
                if (fields[i] == "g16") value_column = static_cast<int>(i);
            }
            if (key_column < 0 || value_column < 0)
                throw std::runtime_error("kernel TSV needs key (or packed_key) and g16 columns");
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max(key_column, value_column)))
            throw std::runtime_error("incomplete kernel row");
        const auto wide_key = unsigned_integer(fields[key_column]);
        std::size_t end = 0;
        const auto wide_value = std::stoll(fields[value_column], &end);
        if (end != fields[value_column].size() || wide_key >= kernel.size() ||
            wide_value < std::numeric_limits<std::int16_t>::min() ||
            wide_value > std::numeric_limits<std::int16_t>::max())
            throw std::runtime_error("invalid kernel entry");
        const auto key = static_cast<std::uint32_t>(wide_key);
        int maximum_label = -1;
        for (int i = 0; i < 8; ++i) {
            const int label = (key >> (3 * i)) & 7;
            if (label > maximum_label + 1)
                throw std::runtime_error("kernel key is not in restricted-growth order");
            maximum_label = std::max(maximum_label, label);
        }
        if (!seen.insert(key).second) throw std::runtime_error("duplicate kernel key");
        kernel[key] = static_cast<std::int16_t>(wide_value);
    }
    if (!input.eof() || seen.empty()) throw std::runtime_error("empty or failed kernel read");
    rows = seen.size();
    if (rows != 1874) throw std::runtime_error("frozen sparse kernel must contain 1874 rows");
    return kernel;
}

class Components {
    std::vector<int> parent_, size_;
public:
    explicit Components(int n) : parent_(n), size_(n) {}
    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
    }
    int root(int value) {
        while (parent_[value] != value) {
            parent_[value] = parent_[parent_[value]];
            value = parent_[value];
        }
        return value;
    }
    void join(int left, int right) {
        left = root(left); right = root(right);
        if (left == right) return;
        if (size_[left] < size_[right]) std::swap(left, right);
        parent_[right] = left;
        size_[left] += size_[right];
    }
};

struct PairDefinition {
    int x = -1, y = -1;
    int anchor_x = -1, anchor_y = -1;
    char orientation = '?';
};

struct PairEvaluation {
    std::int16_t g16 = 0;
    std::int8_t shared = -1; // -1 means an occupied marked endpoint.
    std::uint32_t key = 0;
};

struct StateSnapshot {
    int occupied_count = 0;
    int q = 0, e = 0;
    std::array<PairEvaluation, frozen_pairs> pair{};
};

class SquareGeometry {
    int L_, n_;
    const Kernel& kernel_;
    std::vector<std::array<int, 4>> neighbors_; // N,E,S,W
    std::vector<std::array<int, 4>> port_edges_;
public:
    SquareGeometry(int L, const Kernel& kernel)
        : L_(L), n_(L * L), kernel_(kernel), neighbors_(n_), port_edges_(n_) {
        for (int y = 0; y < L_; ++y) for (int x = 0; x < L_; ++x) {
            const int value = site(x, y);
            neighbors_[value] = {site(x, y + 1), site(x + 1, y),
                                 site(x, y - 1), site(x - 1, y)};
        }
        for (int value = 0; value < n_; ++value) {
            port_edges_[value] = {2 * value, 2 * value + 1,
                                  2 * neighbors_[value][2],
                                  2 * neighbors_[value][3] + 1};
        }
    }
    int L() const { return L_; }
    int n() const { return n_; }
    int site(int x, int y) const {
        x %= L_; y %= L_;
        if (x < 0) x += L_;
        if (y < 0) y += L_;
        return y * L_ + x;
    }
    std::pair<int, int> coordinates(int value) const { return {value % L_, value / L_}; }
    const std::array<int, 4>& neighbors(int value) const { return neighbors_[value]; }
    const std::array<int, 4>& port_edges(int value) const { return port_edges_[value]; }
    const Kernel& kernel() const { return kernel_; }

    std::array<PairDefinition, frozen_pairs> prescribed_pairs() const {
        if (L_ != 32 && L_ != 64)
            throw std::logic_error("prescribed pairs are defined only for L32/L64");
        std::array<PairDefinition, frozen_pairs> result{};
        const int r = L_ / 4;
        int ordinal = 0;
        for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
            const int x = i * r, y = j * r, anchor = site(x, y);
            result[ordinal++] = {anchor, site(x + r, y), x, y, 'H'};
            result[ordinal++] = {anchor, site(x, y + r), x, y, 'V'};
        }
        return result;
    }

    int torus_linf(int left, int right) const {
        const auto [lx, ly] = coordinates(left);
        const auto [rx, ry] = coordinates(right);
        const int dx0 = std::abs(lx - rx), dy0 = std::abs(ly - ry);
        return std::max(std::min(dx0, L_ - dx0), std::min(dy0, L_ - dy0));
    }
};

class StateEvaluator {
    const SquareGeometry& geometry_;
    Components black_, white_;
    std::vector<int> occupied_root_;

    bool occupied(const std::vector<unsigned char>& bits, int value,
                  int forced_site, bool forced_value) const {
        return value == forced_site ? forced_value : bits[value] != 0;
    }

    int component_count(Components& components, const std::vector<unsigned char>& bits,
                        bool black_sites, int forced_site, bool forced_value) {
        int result = 0;
        for (int value = 0; value < geometry_.n(); ++value) {
            const bool is_black = occupied(bits, value, forced_site, forced_value);
            if (is_black == black_sites && components.root(value) == value) ++result;
        }
        return result;
    }

    PairEvaluation evaluate_pair(const std::vector<unsigned char>& bits,
                                 int forced_site, bool forced_value,
                                 const PairDefinition& pair) {
        if (occupied(bits, pair.x, forced_site, forced_value) ||
            occupied(bits, pair.y, forced_site, forced_value)) return {};
        std::array<int, 8> outside{}, labels{};
        std::uint32_t key = 0;
        int next_label = 0;
        for (int i = 0; i < 8; ++i) {
            const int center = i < 4 ? pair.x : pair.y;
            const int direction = i % 4;
            const int neighbor = geometry_.neighbors(center)[direction];
            outside[i] = occupied(bits, neighbor, forced_site, forced_value)
                ? occupied_root_[neighbor]
                : geometry_.n() + geometry_.port_edges(center)[direction];
            int previous = 0;
            while (previous < i && outside[previous] != outside[i]) ++previous;
            labels[i] = previous < i ? labels[previous] : next_label++;
            key |= static_cast<std::uint32_t>(labels[i]) << (3 * i);
        }
        std::array<bool, 8> left{}, right{};
        for (int i = 0; i < 4; ++i) left[labels[i]] = true;
        for (int i = 4; i < 8; ++i) right[labels[i]] = true;
        int shared = 0;
        for (int label = 0; label < next_label; ++label) shared += left[label] && right[label];
        const auto g16 = geometry_.kernel()[key];
        if (shared <= 1 && g16 != 0)
            throw std::logic_error("s<=1 nonzero Bell8 kernel control failed");
        PairEvaluation result;
        result.g16 = g16;
        result.shared = static_cast<std::int8_t>(shared);
        result.key = key;
        return result;
    }

public:
    explicit StateEvaluator(const SquareGeometry& geometry)
        : geometry_(geometry), black_(geometry.n()), white_(geometry.n()),
          occupied_root_(geometry.n(), -1) {}

    StateSnapshot evaluate(const std::vector<unsigned char>& bits,
                           int forced_site, bool forced_value,
                           const std::array<PairDefinition, frozen_pairs>& pairs,
                           int pair_count = frozen_pairs) {
        black_.reset(); white_.reset();
        int occupied_count = 0, edges = 0, faces = 0;
        for (int value = 0; value < geometry_.n(); ++value) {
            const bool is_occupied = occupied(bits, value, forced_site, forced_value);
            occupied_count += is_occupied;
            if (is_occupied) {
                for (int direction : {0, 1}) {
                    const int other = geometry_.neighbors(value)[direction];
                    if (occupied(bits, other, forced_site, forced_value)) {
                        black_.join(value, other);
                        ++edges;
                    }
                }
                const int east = geometry_.neighbors(value)[1];
                const int north = geometry_.neighbors(value)[0];
                const int northeast = geometry_.neighbors(east)[0];
                faces += occupied(bits, east, forced_site, forced_value) &&
                         occupied(bits, north, forced_site, forced_value) &&
                         occupied(bits, northeast, forced_site, forced_value);
            } else {
                for (int direction : {0, 1}) {
                    const int other = geometry_.neighbors(value)[direction];
                    if (!occupied(bits, other, forced_site, forced_value)) white_.join(value, other);
                }
                const int northeast = geometry_.neighbors(geometry_.neighbors(value)[1])[0];
                const int northwest = geometry_.neighbors(geometry_.neighbors(value)[3])[0];
                if (!occupied(bits, northeast, forced_site, forced_value)) white_.join(value, northeast);
                if (!occupied(bits, northwest, forced_site, forced_value)) white_.join(value, northwest);
            }
        }
        const int black_components = component_count(
            black_, bits, true, forced_site, forced_value);
        const int white_components = component_count(
            white_, bits, false, forced_site, forced_value);
        const int q = black_components - white_components - (occupied_count - edges + faces);
        if (q < -1 || q > 1) throw std::logic_error("digital-Alexander q outside {-1,0,1}");
        std::fill(occupied_root_.begin(), occupied_root_.end(), -1);
        for (int value = 0; value < geometry_.n(); ++value)
            if (occupied(bits, value, forced_site, forced_value))
                occupied_root_[value] = black_.root(value);

        StateSnapshot result;
        result.occupied_count = occupied_count;
        result.q = q;
        result.e = q * q;
        for (int ordinal = 0; ordinal < pair_count; ++ordinal)
            result.pair[ordinal] = evaluate_pair(bits, forced_site, forced_value, pairs[ordinal]);
        return result;
    }
};

struct SiteFlip {
    int q0 = 0, q1 = 0, e0 = 0, e1 = 0;
    std::array<PairEvaluation, frozen_pairs> off{}, on{};
};

struct AggregateKey {
    int pair = 0, shell = 0, relation_mask = 0, carrier_mask = 0;
    bool operator<(const AggregateKey& other) const {
        return std::tie(pair, shell, relation_mask, carrier_mask) <
               std::tie(other.pair, other.shell, other.relation_mask, other.carrier_mask);
    }
};

std::int64_t absolute(std::int64_t value) { return value < 0 ? -value : value; }

struct Aggregate {
    std::uint64_t sites = 0, eligible_both = 0;
    std::uint64_t persistent_s2 = 0, shared_transition = 0;
    std::uint64_t topological = 0, kernel_only = 0, joint = 0;
    std::int64_t sum_g16_0 = 0, sum_g16_1 = 0;
    std::int64_t sum_abs_g16_0 = 0, sum_abs_g16_1 = 0;
    std::int64_t sum_delta_g16 = 0, sum_abs_delta_g16 = 0;
    std::int64_t sum_q0 = 0, sum_q1 = 0, sum_delta_q = 0, sum_abs_delta_q = 0;
    std::int64_t sum_e0 = 0, sum_e1 = 0, sum_delta_e = 0, sum_abs_delta_e = 0;
    std::int64_t sum_q_observable_num32 = 0, sum_abs_q_observable_num32 = 0;
    std::int64_t sum_q_kernel_num32 = 0, sum_abs_q_kernel_num32 = 0;
    std::int64_t sum_q_product_delta_num16 = 0, sum_abs_q_product_delta_num16 = 0;
    std::int64_t sum_e_observable_num32 = 0, sum_abs_e_observable_num32 = 0;
    std::int64_t sum_e_kernel_num32 = 0, sum_abs_e_kernel_num32 = 0;
    std::int64_t sum_e_product_delta_num16 = 0, sum_abs_e_product_delta_num16 = 0;
    std::int64_t q_identity_residual_num32 = 0, e_identity_residual_num32 = 0;
    std::int64_t q_identity_max_abs = 0, e_identity_max_abs = 0;
    std::int64_t sum_shared0 = 0, sum_shared1 = 0;

    void add(const SiteFlip& flip, const PairEvaluation& g0, const PairEvaluation& g1,
             int carrier_mask) {
        ++sites;
        eligible_both += g0.shared >= 0 && g1.shared >= 0;
        persistent_s2 += (carrier_mask & two_bridge_persistent_bit) != 0;
        shared_transition += (carrier_mask & shared_transition_or_merger_bit) != 0;
        topological += (carrier_mask & kernel_preserving_topological_bit) != 0;
        const std::int64_t delta_g16 = static_cast<std::int64_t>(g1.g16) - g0.g16;
        const int delta_q = flip.q1 - flip.q0, delta_e = flip.e1 - flip.e0;
        kernel_only += delta_g16 != 0 && delta_q == 0 && delta_e == 0;
        joint += delta_g16 != 0 && (delta_q != 0 || delta_e != 0);

        sum_g16_0 += g0.g16; sum_g16_1 += g1.g16;
        sum_abs_g16_0 += absolute(g0.g16); sum_abs_g16_1 += absolute(g1.g16);
        sum_delta_g16 += delta_g16; sum_abs_delta_g16 += absolute(delta_g16);
        sum_q0 += flip.q0; sum_q1 += flip.q1;
        sum_delta_q += delta_q; sum_abs_delta_q += absolute(delta_q);
        sum_e0 += flip.e0; sum_e1 += flip.e1;
        sum_delta_e += delta_e; sum_abs_delta_e += absolute(delta_e);
        sum_shared0 += g0.shared; sum_shared1 += g1.shared;

        const std::int64_t g_sum = static_cast<std::int64_t>(g0.g16) + g1.g16;
        const std::int64_t g_delta = delta_g16;
        const std::int64_t q_sum = flip.q0 + flip.q1;
        const std::int64_t e_sum = flip.e0 + flip.e1;
        const std::int64_t q_observable = g_sum * delta_q;
        const std::int64_t q_kernel = q_sum * g_delta;
        const std::int64_t q_product = static_cast<std::int64_t>(flip.q1) * g1.g16 -
                                       static_cast<std::int64_t>(flip.q0) * g0.g16;
        const std::int64_t e_observable = g_sum * delta_e;
        const std::int64_t e_kernel = e_sum * g_delta;
        const std::int64_t e_product = static_cast<std::int64_t>(flip.e1) * g1.g16 -
                                       static_cast<std::int64_t>(flip.e0) * g0.g16;
        const std::int64_t q_residual = q_observable + q_kernel - 2 * q_product;
        const std::int64_t e_residual = e_observable + e_kernel - 2 * e_product;
        if (q_residual != 0 || e_residual != 0)
            throw std::logic_error("pair/site midpoint product identity failed");

        sum_q_observable_num32 += q_observable;
        sum_abs_q_observable_num32 += absolute(q_observable);
        sum_q_kernel_num32 += q_kernel;
        sum_abs_q_kernel_num32 += absolute(q_kernel);
        sum_q_product_delta_num16 += q_product;
        sum_abs_q_product_delta_num16 += absolute(q_product);
        sum_e_observable_num32 += e_observable;
        sum_abs_e_observable_num32 += absolute(e_observable);
        sum_e_kernel_num32 += e_kernel;
        sum_abs_e_kernel_num32 += absolute(e_kernel);
        sum_e_product_delta_num16 += e_product;
        sum_abs_e_product_delta_num16 += absolute(e_product);
        q_identity_residual_num32 += q_residual;
        e_identity_residual_num32 += e_residual;
        q_identity_max_abs = std::max(q_identity_max_abs, absolute(q_residual));
        e_identity_max_abs = std::max(e_identity_max_abs, absolute(e_residual));
    }
};

struct SizeStats {
    int L = 0;
    std::uint64_t seed = 0;
    double wall_seconds = 0, cpu_seconds = 0;
    std::uint64_t configuration_rows = 0, shell_rows = 0;
};

int shell_index(int distance) {
    if (distance == 0) return 0;
    int result = 1;
    while ((1 << result) <= distance) ++result;
    return result;
}

std::pair<int, int> shell_bounds(int shell, int L) {
    if (shell == 0) return {0, 0};
    const int lower = 1 << (shell - 1);
    const int upper = std::min(L / 2, (1 << shell) - 1);
    return {lower, upper};
}

int relation_mask(const SquareGeometry& geometry, int z, const PairDefinition& pair) {
    int result = 0;
    if (z == pair.x) result |= 1;
    if (z == pair.y) result |= 2;
    if (std::find(geometry.neighbors(pair.x).begin(), geometry.neighbors(pair.x).end(), z) !=
        geometry.neighbors(pair.x).end()) result |= 4;
    if (std::find(geometry.neighbors(pair.y).begin(), geometry.neighbors(pair.y).end(), z) !=
        geometry.neighbors(pair.y).end()) result |= 8;
    return result;
}

std::uint64_t occupation_fingerprint(const std::vector<unsigned char>& bits) {
    std::uint64_t value = 1469598103934665603ULL;
    for (const auto bit : bits) {
        value ^= bit;
        value *= 1099511628211ULL;
    }
    return value;
}

std::uint32_t pack_labels(const std::array<int, 8>& labels) {
    std::uint32_t key = 0;
    for (int i = 0; i < 8; ++i) key |= static_cast<std::uint32_t>(labels[i]) << (3 * i);
    return key;
}

void run_tiny_controls(const Kernel& kernel) {
    SquareGeometry geometry(4, kernel);
    std::array<PairDefinition, frozen_pairs> pairs{};
    pairs[0] = {0, 2, 0, 0, 'H'};
    StateEvaluator evaluator(geometry);
    const int z = 5;

    std::vector<unsigned char> first(16, 0);
    for (int value : {3, 4, 6, 8, 10, 11}) first[value] = 1;
    const auto first0 = evaluator.evaluate(first, z, false, pairs, 1);
    const auto first1 = evaluator.evaluate(first, z, true, pairs, 1);
    const auto expected_key = pack_labels({0, 1, 2, 3, 0, 3, 4, 5});
    if (first0.q != -1 || first1.q != 0 || first0.pair[0].g16 != 4 ||
        first1.pair[0].g16 != 4 || first0.pair[0].key != expected_key ||
        first1.pair[0].key != expected_key) {
        std::ostringstream message;
        message << "first L4 kernel-preserving topological control failed: q="
                << first0.q << "->" << first1.q << ",g16="
                << first0.pair[0].g16 << "->" << first1.pair[0].g16
                << ",shared=" << static_cast<int>(first0.pair[0].shared) << "->"
                << static_cast<int>(first1.pair[0].shared) << ",key="
                << first0.pair[0].key << "->" << first1.pair[0].key
                << ",expected_key=" << expected_key;
        throw std::logic_error(message.str());
    }

    std::vector<unsigned char> second(16, 0);
    for (int value : {3, 4, 6}) second[value] = 1;
    const auto second0 = evaluator.evaluate(second, z, false, pairs, 1);
    const auto second1 = evaluator.evaluate(second, z, true, pairs, 1);
    if (second0.q != -1 || second1.q != -1 || second0.pair[0].g16 != 0 ||
        second1.pair[0].g16 != 4) {
        std::ostringstream message;
        message << "second L4 kernel-only control failed: q="
                << second0.q << "->" << second1.q << ",g16="
                << second0.pair[0].g16 << "->" << second1.pair[0].g16
                << ",shared=" << static_cast<int>(second0.pair[0].shared) << "->"
                << static_cast<int>(second1.pair[0].shared) << ",key="
                << second0.pair[0].key << "->" << second1.pair[0].key;
        throw std::logic_error(message.str());
    }
}

std::uint64_t peak_rss_bytes() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

void write_aggregate(std::ostream& output, int L, std::uint64_t seed, int counter,
                     const PairDefinition& pair, const AggregateKey& key,
                     const Aggregate& value) {
    const auto [lower, upper] = shell_bounds(key.shell, L);
    output << L << ',' << seed << ',' << counter << ',' << key.pair << ','
           << pair.x << ',' << pair.y << ',' << pair.anchor_x << ',' << pair.anchor_y << ','
           << pair.orientation << ',' << key.shell << ',' << lower << ',' << upper << ','
           << key.relation_mask << ',' << key.carrier_mask << ',' << value.sites << ','
           << value.eligible_both << ',' << value.persistent_s2 << ','
           << value.shared_transition << ',' << value.topological << ','
           << value.kernel_only << ',' << value.joint << ','
           << value.sum_g16_0 << ',' << value.sum_g16_1 << ','
           << value.sum_abs_g16_0 << ',' << value.sum_abs_g16_1 << ','
           << value.sum_delta_g16 << ',' << value.sum_abs_delta_g16 << ','
           << value.sum_q0 << ',' << value.sum_q1 << ',' << value.sum_delta_q << ','
           << value.sum_abs_delta_q << ',' << value.sum_e0 << ',' << value.sum_e1 << ','
           << value.sum_delta_e << ',' << value.sum_abs_delta_e << ','
           << value.sum_q_observable_num32 << ',' << value.sum_abs_q_observable_num32 << ','
           << value.sum_q_kernel_num32 << ',' << value.sum_abs_q_kernel_num32 << ','
           << value.sum_q_product_delta_num16 << ',' << value.sum_abs_q_product_delta_num16 << ','
           << value.sum_e_observable_num32 << ',' << value.sum_abs_e_observable_num32 << ','
           << value.sum_e_kernel_num32 << ',' << value.sum_abs_e_kernel_num32 << ','
           << value.sum_e_product_delta_num16 << ',' << value.sum_abs_e_product_delta_num16 << ','
           << value.q_identity_residual_num32 << ',' << value.e_identity_residual_num32 << ','
           << value.q_identity_max_abs << ',' << value.e_identity_max_abs << ','
           << value.sum_shared0 << ',' << value.sum_shared1 << '\n';
}

SizeStats run_size(int L, std::uint64_t seed, int threads, const Kernel& kernel,
                   std::ostream& configuration_output, std::ostream& shell_output) {
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    SquareGeometry geometry(L, kernel);
    const auto pairs = geometry.prescribed_pairs();
    std::mt19937_64 rng(seed);
    std::vector<unsigned char> bits(geometry.n());
    SizeStats stats; stats.L = L; stats.seed = seed;

    for (int counter = 0; counter < frozen_configurations; ++counter) {
        int occupied_count = 0;
        for (int value = 0; value < geometry.n(); ++value) {
            bits[value] = rng() < bernoulli_threshold;
            occupied_count += bits[value] != 0;
        }
        StateEvaluator original_evaluator(geometry);
        const auto original = original_evaluator.evaluate(bits, -1, false, pairs);
        if (original.occupied_count != occupied_count)
            throw std::logic_error("original occupied count mismatch");

        std::vector<SiteFlip> flips(geometry.n());
        std::atomic<int> next_site{0};
        std::atomic<bool> worker_failed{false};
        std::exception_ptr worker_error;
        std::mutex worker_error_mutex;
        const int worker_count = std::min(threads, geometry.n());
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (int worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&]() {
                try {
                    StateEvaluator evaluator(geometry);
                    while (!worker_failed.load(std::memory_order_relaxed)) {
                        const int z = next_site.fetch_add(1, std::memory_order_relaxed);
                        if (z >= geometry.n()) break;
                        const auto off = evaluator.evaluate(bits, z, false, pairs);
                        const auto on = evaluator.evaluate(bits, z, true, pairs);
                        auto& target = flips[z];
                        target.q0 = off.q; target.q1 = on.q;
                        target.e0 = off.e; target.e1 = on.e;
                        target.off = off.pair; target.on = on.pair;
                    }
                } catch (...) {
                    worker_failed.store(true, std::memory_order_relaxed);
                    std::lock_guard<std::mutex> lock(worker_error_mutex);
                    if (!worker_error) worker_error = std::current_exception();
                }
            });
        }
        for (auto& worker : workers) worker.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::int64_t original_sum_g16 = 0;
        std::array<std::int64_t, 5> original_sum_by_shared{};
        std::array<std::uint64_t, 5> original_pairs_by_shared{};
        std::uint64_t original_eligible = 0, original_nonzero = 0;
        for (const auto& pair : original.pair) {
            original_sum_g16 += pair.g16;
            if (pair.shared >= 0) {
                ++original_eligible;
                ++original_pairs_by_shared[pair.shared];
                original_sum_by_shared[pair.shared] += pair.g16;
                original_nonzero += pair.g16 != 0;
            }
        }
        configuration_output << L << ',' << seed << ',' << counter << ','
            << occupation_fingerprint(bits) << ',' << occupied_count << ',' << original.q << ','
            << original.e << ',' << original_sum_g16 << ',' << original_eligible << ','
            << original_nonzero;
        for (const auto value : original_sum_by_shared) configuration_output << ',' << value;
        for (const auto value : original_pairs_by_shared) configuration_output << ',' << value;
        configuration_output << '\n';
        ++stats.configuration_rows;

        std::map<AggregateKey, Aggregate> aggregates;
        for (int z = 0; z < geometry.n(); ++z) {
            const auto& flip = flips[z];
            const bool selected_on = bits[z] != 0;
            if ((selected_on ? flip.q1 : flip.q0) != original.q ||
                (selected_on ? flip.e1 : flip.e0) != original.e)
                throw std::logic_error("forced state does not reproduce original q/E");
            for (int ordinal = 0; ordinal < frozen_pairs; ++ordinal) {
                const auto& g0 = flip.off[ordinal];
                const auto& g1 = flip.on[ordinal];
                const auto& selected = selected_on ? g1 : g0;
                if (selected.g16 != original.pair[ordinal].g16 ||
                    selected.shared != original.pair[ordinal].shared ||
                    selected.key != original.pair[ordinal].key)
                    throw std::logic_error("forced state does not reproduce original pair kernel");
                const auto& pair = pairs[ordinal];
                const int distance = std::min(
                    geometry.torus_linf(z, pair.x), geometry.torus_linf(z, pair.y));
                const int shell = shell_index(distance);
                const int relation = relation_mask(geometry, z, pair);
                const int delta_g16 = static_cast<int>(g1.g16) - g0.g16;
                const int delta_q = flip.q1 - flip.q0, delta_e = flip.e1 - flip.e0;
                const bool endpoints_valid = g0.shared >= 0 && g1.shared >= 0;
                const bool active = delta_g16 != 0 || delta_q != 0 || delta_e != 0;
                int carrier = 0;
                if (active && endpoints_valid && g0.shared == 2 && g1.shared == 2 &&
                    (g0.g16 != 0 || g1.g16 != 0)) carrier |= two_bridge_persistent_bit;
                if (active && endpoints_valid && g0.shared != g1.shared &&
                    (g0.g16 != 0 || g1.g16 != 0)) carrier |= shared_transition_or_merger_bit;
                if (delta_g16 == 0 && g0.g16 != 0 && (delta_q != 0 || delta_e != 0))
                    carrier |= kernel_preserving_topological_bit;
                if (delta_g16 != 0) carrier |= kernel_changed_bit;
                AggregateKey key{ordinal, shell, relation, carrier};
                aggregates[key].add(flip, g0, g1, carrier);
            }
        }
        for (const auto& item : aggregates) {
            write_aggregate(shell_output, L, seed, counter, pairs[item.first.pair],
                            item.first, item.second);
            ++stats.shell_rows;
        }
    }
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    stats.cpu_seconds = static_cast<double>(std::clock() - cpu_start) / CLOCKS_PER_SEC;
    return stats;
}

std::string json_escape(const std::string& value) {
    std::string result;
    for (const char character : value) {
        if (character == '\\' || character == '"') result.push_back('\\');
        result.push_back(character);
    }
    return result;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto total_start = std::chrono::steady_clock::now();
        const auto total_cpu_start = std::clock();
        const auto options = parse_options(argc, argv);
        const std::string config_path = options.output_prefix + ".config.csv";
        const std::string shell_path = options.output_prefix + ".shell.csv";
        const std::string metadata_path = options.output_prefix + ".metadata.json";
        if (std::ifstream(config_path).good() || std::ifstream(shell_path).good() ||
            std::ifstream(metadata_path).good())
            throw std::runtime_error("an output sidecar exists; refusing overwrite");

        std::size_t kernel_rows = 0;
        const auto kernel = read_kernel(options.kernel_path, kernel_rows);
        run_tiny_controls(kernel);

        std::ofstream configuration_output(config_path), shell_output(shell_path);
        if (!configuration_output || !shell_output)
            throw std::runtime_error("cannot create output CSVs");
        configuration_output << "L,seed,configuration_counter,occupation_fnv1a64,K,q,E,"
            "sum_g16,eligible_pairs,nonzero_pairs";
        for (int shared = 0; shared <= 4; ++shared)
            configuration_output << ",sum_g16_shared" << shared;
        for (int shared = 0; shared <= 4; ++shared)
            configuration_output << ",pairs_shared" << shared;
        configuration_output << '\n';

        shell_output << "L,seed,configuration_counter,pair_ordinal,x,y,anchor_x,anchor_y,"
            "orientation,shell,shell_lower,shell_upper,relation_mask,carrier_mask,sites,"
            "eligible_both,persistent_s2_count,shared_transition_count,"
            "kernel_preserving_topological_count,kernel_only_count,joint_count,"
            "sum_g16_0,sum_g16_1,sum_abs_g16_0,sum_abs_g16_1,"
            "sum_delta_g16,sum_abs_delta_g16,sum_q0,sum_q1,sum_delta_q,sum_abs_delta_q,"
            "sum_E0,sum_E1,sum_delta_E,sum_abs_delta_E,"
            "sum_q_observable_num32,sum_abs_q_observable_num32,"
            "sum_q_kernel_num32,sum_abs_q_kernel_num32,"
            "sum_q_product_delta_num16,sum_abs_q_product_delta_num16,"
            "sum_E_observable_num32,sum_abs_E_observable_num32,"
            "sum_E_kernel_num32,sum_abs_E_kernel_num32,"
            "sum_E_product_delta_num16,sum_abs_E_product_delta_num16,"
            "q_identity_residual_num32,E_identity_residual_num32,"
            "q_identity_max_abs,E_identity_max_abs,sum_shared0,sum_shared1\n";

        const auto stats32 = run_size(32, seed_L32, options.threads, kernel,
                                      configuration_output, shell_output);
        const auto stats64 = run_size(64, seed_L64, options.threads, kernel,
                                      configuration_output, shell_output);
        configuration_output.close(); shell_output.close();
        if (!configuration_output || !shell_output)
            throw std::runtime_error("output CSV write failed");

        const double total_wall = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - total_start).count();
        const double total_cpu = static_cast<double>(std::clock() - total_cpu_start) / CLOCKS_PER_SEC;
        const auto rss = peak_rss_bytes();
        std::ofstream metadata(metadata_path);
        if (!metadata) throw std::runtime_error("cannot create metadata JSON");
        metadata << std::setprecision(17)
            << "{\n  \"schema\": \"matching-one.p337-thermal-pivotal-preflight.run.v1\",\n"
            << "  \"status\": \"completed_frozen_64_configuration_replay\",\n"
            << "  \"kernel_path\": \"" << json_escape(options.kernel_path) << "\",\n"
            << "  \"kernel_rows\": " << kernel_rows << ",\n"
            << "  \"p_requested_decimal\": \"0.592746050790\",\n"
            << "  \"bernoulli_threshold_2pow64\": " << bernoulli_threshold << ",\n"
            << "  \"rng\": \"std::mt19937_64; one word per site; x-fast row-major\",\n"
            << "  \"configurations_per_size\": " << frozen_configurations << ",\n"
            << "  \"pairs_per_configuration\": " << frozen_pairs << ",\n"
            << "  \"threads\": " << options.threads << ",\n"
            << "  \"tiny_controls\": {\"kernel_preserving_topological\": true, "
               "\"kernel_only\": true},\n"
            << "  \"carrier_bits\": {\"1\": \"two_bridge_persistent\", "
               "\"2\": \"shared_transition_or_merger\", "
               "\"4\": \"kernel_preserving_topological\", "
               "\"8\": \"kernel_changed_kernel_only_or_joint\"},\n"
            << "  \"absolute_value_rule\": \"take abs for each pair/site primitive before shell aggregation\",\n"
            << "  \"midpoint_units\": \"observable and kernel numerators divide by 32; product-delta numerators divide by 16\",\n"
            << "  \"L32\": {\"seed\": " << stats32.seed << ", \"configuration_rows\": "
            << stats32.configuration_rows << ", \"shell_rows\": " << stats32.shell_rows
            << ", \"wall_seconds\": " << stats32.wall_seconds
            << ", \"cpu_seconds\": " << stats32.cpu_seconds << "},\n"
            << "  \"L64\": {\"seed\": " << stats64.seed << ", \"configuration_rows\": "
            << stats64.configuration_rows << ", \"shell_rows\": " << stats64.shell_rows
            << ", \"wall_seconds\": " << stats64.wall_seconds
            << ", \"cpu_seconds\": " << stats64.cpu_seconds << "},\n"
            << "  \"total_wall_seconds\": " << total_wall << ",\n"
            << "  \"total_cpu_seconds\": " << total_cpu << ",\n"
            << "  \"peak_RSS_bytes\": " << rss << ",\n"
            << "  \"configuration_output\": \"" << json_escape(config_path) << "\",\n"
            << "  \"shell_output\": \"" << json_escape(shell_path) << "\",\n"
            << "  \"sha256_boundary\": \"execution wrapper records binary, kernel and output SHA256; this standard-library producer does not shell out\",\n"
            << "  \"evidence_boundary\": \"deterministic replay and cost/semantic preflight; not a significance estimate or independent evidence\"\n}\n";
        metadata.close();
        if (!metadata) throw std::runtime_error("metadata write failed");

        std::cout << std::setprecision(17)
                  << "{\"status\":\"completed_frozen_64_configuration_replay\","
                  << "\"wall_seconds\":" << total_wall << ",\"cpu_seconds\":" << total_cpu
                  << ",\"peak_RSS_bytes\":" << rss << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
