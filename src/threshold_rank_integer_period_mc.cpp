// Production threshold-rank Newman--Ziff engine for arbitrary 2x2 integer
// period matrices.
//
// The finite graph is Z^2 / P Z^2, with the columns of P interpreted as the
// two lifted period vectors.  Vertices are enumerated once through a column
// Hermite normal form
//
//     H = P V = ((h11,h12),(0,h22)),  det(V)=+/-1,
//
// using label = rx + h11*ry for 0<=rx<h11, 0<=ry<h22.  Closed lifted cycle
// displacements are converted to period-basis windings exactly using
// adj(P)*displacement/det(P).  No floating-point geometry arithmetic is used.
//
// The first production designs are the nonprimitive norm-4 Gaussian children
// N=260 and N=340, whose Smith invariants are (2,130) and (2,170).

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Int = std::int64_t;

struct Vector {
    Int x = 0;
    Int y = 0;
};

struct Matrix {
    Int a = 0;
    Int b = 0;
    Int c = 0;
    Int d = 0;
};

struct Edge {
    int i;
    int j;
    int dx;
    int dy;
};

Int checked_int(__int128 value, const char* context) {
    if (value < std::numeric_limits<Int>::min() ||
        value > std::numeric_limits<Int>::max()) {
        throw std::overflow_error(std::string(context) + " overflows int64");
    }
    return static_cast<Int>(value);
}

Int safe_abs(Int value, const char* context) {
    if (value == std::numeric_limits<Int>::min()) {
        throw std::overflow_error(std::string(context) + " absolute value overflows int64");
    }
    return std::llabs(value);
}

Int determinant(const Matrix& matrix) {
    return checked_int(static_cast<__int128>(matrix.a) * matrix.d -
                       static_cast<__int128>(matrix.b) * matrix.c,
                       "period determinant");
}

Int positive_mod(Int value, Int modulus) {
    value %= modulus;
    return value < 0 ? value + modulus : value;
}

Int positive_mod128(__int128 value, Int modulus) {
    value %= modulus;
    if (value < 0) value += modulus;
    return static_cast<Int>(value);
}

struct Bezout {
    Int gcd;
    Int x;
    Int y;
};

Bezout extended_gcd_nonnegative(Int a, Int b) {
    if (a < 0 || b < 0 || (a == 0 && b == 0)) {
        throw std::invalid_argument("extended gcd requires nonnegative nonzero input");
    }
    Int old_r = a, r = b;
    Int old_s = 1, s = 0;
    Int old_t = 0, t = 1;
    while (r != 0) {
        const Int quotient = old_r / r;
        const Int next_r = old_r - quotient * r;
        const Int next_s = old_s - quotient * s;
        const Int next_t = old_t - quotient * t;
        old_r = r; r = next_r;
        old_s = s; s = next_s;
        old_t = t; t = next_t;
    }
    return {old_r, old_s, old_t};
}

Bezout extended_gcd(Int a, Int b) {
    const Bezout positive = extended_gcd_nonnegative(
        safe_abs(a, "matrix entry"), safe_abs(b, "matrix entry"));
    return {positive.gcd, a < 0 ? -positive.x : positive.x,
            b < 0 ? -positive.y : positive.y};
}

struct QuotientCoordinates {
    Matrix periods;
    Int det;
    Int order;
    Int h11;
    Int h12;
    Int h22;
    Int smith1;
    Int smith2;

    explicit QuotientCoordinates(Matrix input) : periods(input), det(determinant(input)) {
        for (const Int entry : std::array<Int, 4>{{
                 periods.a, periods.b, periods.c, periods.d}}) {
            if (entry < std::numeric_limits<int>::min() ||
                entry > std::numeric_limits<int>::max()) {
                throw std::invalid_argument(
                    "period matrix entries must fit signed 32-bit production range");
            }
        }
        if (det == 0) throw std::invalid_argument("period matrix must be nonsingular");
        if (det == std::numeric_limits<Int>::min()) {
            throw std::overflow_error("absolute determinant overflows int64");
        }
        order = std::llabs(det);
        if (order > std::numeric_limits<int>::max()) {
            throw std::invalid_argument("quotient order exceeds production int indexing");
        }

        // If u*c+v*d=g, the unimodular column transform with columns
        // (d/g,-c/g) and (u,v) makes the lower row (0,g).  A sign flip and
        // one column reduction give the unique upper column-HNF convention.
        const Bezout lower = extended_gcd(periods.c, periods.d);
        h22 = lower.gcd;
        h11 = order / h22;
        const __int128 raw_h12 = static_cast<__int128>(periods.a) * lower.x +
                                 static_cast<__int128>(periods.b) * lower.y;
        h12 = positive_mod128(raw_h12, h11);

        smith1 = std::gcd(
            std::gcd(safe_abs(periods.a, "matrix entry"),
                     safe_abs(periods.b, "matrix entry")),
            std::gcd(safe_abs(periods.c, "matrix entry"),
                     safe_abs(periods.d, "matrix entry")));
        smith2 = order / smith1;

        if (h11 <= 0 || h22 <= 0 || h11 * h22 != order ||
            h12 < 0 || h12 >= h11 || smith1 <= 0 || smith1 * smith2 != order) {
            throw std::logic_error("invalid HNF/Smith quotient construction");
        }
    }

    int label(Vector point) const {
        const Int quotient_y = (point.y - positive_mod(point.y, h22)) / h22;
        const Int ry = point.y - quotient_y * h22;
        const Int rx = positive_mod(point.x - quotient_y * h12, h11);
        const Int value = rx + h11 * ry;
        if (value < 0 || value >= order) throw std::logic_error("invalid quotient label");
        return static_cast<int>(value);
    }

    Vector representative(int label_value) const {
        if (label_value < 0 || label_value >= order) {
            throw std::out_of_range("quotient label outside range");
        }
        return {label_value % h11, label_value / h11};
    }

    Vector winding(Int dx, Int dy) const {
        const __int128 numerator0 = static_cast<__int128>(periods.d) * dx -
                                    static_cast<__int128>(periods.b) * dy;
        const __int128 numerator1 = -static_cast<__int128>(periods.c) * dx +
                                    static_cast<__int128>(periods.a) * dy;
        if (numerator0 % det != 0 || numerator1 % det != 0) {
            throw std::logic_error("cycle displacement is outside the period lattice");
        }
        return {checked_int(numerator0 / det, "winding coordinate"),
                checked_int(numerator1 / det, "winding coordinate")};
    }

    Vector period_vector(Vector winding_value) const {
        return {
            checked_int(static_cast<__int128>(periods.a) * winding_value.x +
                        static_cast<__int128>(periods.b) * winding_value.y,
                        "period vector"),
            checked_int(static_cast<__int128>(periods.c) * winding_value.x +
                        static_cast<__int128>(periods.d) * winding_value.y,
                        "period vector")
        };
    }
};

struct Geometry {
    QuotientCoordinates quotient;
    int n;
    std::vector<Edge> primal_edges;
    std::vector<Edge> matching_edges;
    std::vector<std::vector<int>> primal_incident;
    std::vector<std::vector<int>> matching_incident;

    explicit Geometry(Matrix periods)
        : quotient(periods), n(static_cast<int>(quotient.order)) {}
};

struct PairDesign {
    int n;
    int a1;
    int b1;
    Matrix first;
    int a2;
    int b2;
    Matrix second;
    std::string id;
};

const std::vector<PairDesign> kDesigns = {
    {260, 16, 2, {16, -2, 2, 16}, 14, 8, {14, -8, 8, 14}, "N65_to_N260_q4"},
    {340, 18, 4, {18, -4, 4, 18}, 14, 12, {14, -12, 12, 14}, "N85_to_N340_q4"},
};

Vector primitive(Vector value) {
    const Int divisor = std::gcd(safe_abs(value.x, "winding"),
                                 safe_abs(value.y, "winding"));
    if (divisor == 0) return value;
    value.x /= divisor;
    value.y /= divisor;
    if (value.x < 0 || (value.x == 0 && value.y < 0)) {
        value.x = -value.x;
        value.y = -value.y;
    }
    return value;
}

class HomologyUnionFind {
  public:
    explicit HomologyUnionFind(const QuotientCoordinates& quotient)
        : quotient_(quotient), parent_(quotient.order), size_(quotient.order),
          delta_x_(quotient.order), delta_y_(quotient.order), rank_(quotient.order),
          basis_(quotient.order) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
    }

    struct FindResult {
        int root;
        Int dx;
        Int dy;
    };

    FindResult find(int vertex) {
        if (parent_[vertex] == vertex) return {vertex, 0, 0};
        const int old_parent = parent_[vertex];
        const FindResult above = find(old_parent);
        delta_x_[vertex] += above.dx;
        delta_y_[vertex] += above.dy;
        parent_[vertex] = above.root;
        return {above.root, delta_x_[vertex], delta_y_[vertex]};
    }

    void extend(int root, Vector value) {
        if ((value.x == 0 && value.y == 0) || rank_[root] == 2) return;
        value = primitive(value);
        if (rank_[root] == 0) {
            basis_[root][0] = value;
            rank_[root] = 1;
            return;
        }
        const Vector first = basis_[root][0];
        if (static_cast<__int128>(first.x) * value.y !=
            static_cast<__int128>(first.y) * value.x) {
            basis_[root][1] = value;
            rank_[root] = 2;
        }
    }

    void add_edge(const Edge& edge) {
        FindResult first = find(edge.i);
        FindResult second = find(edge.j);
        Int root_dx = first.dx + edge.dx - second.dx;
        Int root_dy = first.dy + edge.dy - second.dy;
        if (first.root == second.root) {
            extend(first.root, quotient_.winding(root_dx, root_dy));
            return;
        }
        if (size_[first.root] < size_[second.root]) {
            std::swap(first, second);
            root_dx = -root_dx;
            root_dy = -root_dy;
        }
        parent_[second.root] = first.root;
        delta_x_[second.root] = root_dx;
        delta_y_[second.root] = root_dy;
        size_[first.root] += size_[second.root];
        for (std::uint8_t index = 0; index < rank_[second.root]; ++index) {
            extend(first.root, basis_[second.root][index]);
        }
        rank_[second.root] = 0;
    }

    bool component_crosses(int vertex) { return rank_[find(vertex).root] == 2; }

    int component_rank(int vertex) { return rank_[find(vertex).root]; }

    Vector component_line(int vertex) {
        const int root = find(vertex).root;
        if (rank_[root] != 1) {
            throw std::logic_error("projective line requested outside rank one");
        }
        return primitive(basis_[root][0]);
    }

  private:
    const QuotientCoordinates& quotient_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<Int> delta_x_;
    std::vector<Int> delta_y_;
    std::vector<std::uint8_t> rank_;
    std::vector<std::array<Vector, 2>> basis_;
};

std::vector<std::vector<int>> make_incident(int n, const std::vector<Edge>& edges) {
    std::vector<std::vector<int>> incident(n);
    for (int index = 0; index < static_cast<int>(edges.size()); ++index) {
        const Edge& edge = edges[index];
        incident[edge.i].push_back(index);
        if (edge.j != edge.i) incident[edge.j].push_back(index);
    }
    return incident;
}

Geometry make_geometry(Matrix periods) {
    Geometry geometry(periods);
    geometry.primal_edges.reserve(2 * geometry.n);
    geometry.matching_edges.reserve(4 * geometry.n);
    const std::array<Vector, 4> steps = {{{1, 0}, {0, 1}, {1, 1}, {1, -1}}};
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const Vector source = geometry.quotient.representative(vertex);
        for (std::size_t index = 0; index < steps.size(); ++index) {
            const Vector step = steps[index];
            const Edge edge{vertex,
                            geometry.quotient.label({source.x + step.x,
                                                     source.y + step.y}),
                            static_cast<int>(step.x), static_cast<int>(step.y)};
            if (index < 2) geometry.primal_edges.push_back(edge);
            geometry.matching_edges.push_back(edge);
        }
    }
    geometry.primal_incident = make_incident(geometry.n, geometry.primal_edges);
    geometry.matching_incident = make_incident(geometry.n, geometry.matching_edges);
    return geometry;
}

class ThresholdEngine {
  public:
    explicit ThresholdEngine(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n), union_find_(geometry.quotient) {}

    int first_cross(const std::vector<int>& permutation, bool matching, bool reverse) {
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        const std::vector<Edge>& edges = matching ? geometry_.matching_edges
                                                  : geometry_.primal_edges;
        const std::vector<std::vector<int>>& incident = matching
            ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) union_find_.add_edge(edge);
            }
            if (union_find_.component_crosses(vertex)) return offset + 1;
        }
        throw std::logic_error("fully occupied graph did not cross wrap");
    }

    std::pair<int, int> ranks(const std::vector<int>& permutation) {
        const int k_plus = first_cross(permutation, false, false);
        const int reverse_white = first_cross(permutation, true, true);
        const int k_minus = geometry_.n - reverse_white + 1;
        if (k_minus > k_plus) throw std::logic_error("K_minus exceeds K_plus");
        return {k_minus, k_plus};
    }

    struct ProjectiveBirth {
        int tau1 = 0;
        int tau2 = 0;
        Vector ell;
        bool direct_rank2 = false;
    };

    ProjectiveBirth projective_birth(const std::vector<int>& permutation) {
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        ProjectiveBirth result;
        int ambient_rank = 0;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[offset];
            active_[vertex] = 1;
            for (const int edge_index : geometry_.primal_incident[vertex]) {
                const Edge& edge = geometry_.primal_edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) union_find_.add_edge(edge);
            }
            const int rank = union_find_.component_rank(vertex);
            if (ambient_rank == 0 && rank == 1) {
                result.tau1 = offset + 1;
                result.ell = union_find_.component_line(vertex);
                ambient_rank = 1;
            } else if (ambient_rank == 0 && rank == 2) {
                result.tau1 = result.tau2 = offset + 1;
                result.direct_rank2 = true;
                return result;
            } else if (ambient_rank == 1 && rank == 1) {
                const Vector current = union_find_.component_line(vertex);
                if (current.x != result.ell.x || current.y != result.ell.y) {
                    throw std::logic_error("rank-one projective line changed along filtration");
                }
            } else if (ambient_rank == 1 && rank == 2) {
                result.tau2 = offset + 1;
                return result;
            }
        }
        throw std::logic_error("fully occupied graph did not reach ambient rank two");
    }

    int reverse_kminus(const std::vector<int>& permutation) {
        const int reverse_white = first_cross(permutation, true, true);
        return geometry_.n - reverse_white + 1;
    }

  private:
    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    HomologyUnionFind union_find_;
};

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

class SplitMixStream {
  public:
    explicit SplitMixStream(std::uint64_t state) : state_(state) {}
    std::uint64_t next() {
        state_ += 0x9e3779b97f4a7c15ULL;
        std::uint64_t value = state_;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }
    std::uint64_t below(std::uint64_t bound) {
        const std::uint64_t remainder =
            (std::numeric_limits<std::uint64_t>::max() % bound + 1) % bound;
        const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max() - remainder;
        while (true) {
            const std::uint64_t value = next();
            if (remainder == 0 || value <= maximum) return value % bound;
        }
    }
  private:
    std::uint64_t state_;
};

void counter_permutation(int n, std::uint64_t seed, std::uint64_t replica,
                         std::vector<int>& permutation) {
    permutation.resize(n);
    std::iota(permutation.begin(), permutation.end(), 0);
    const std::uint64_t stream_key = splitmix64(
        seed ^ splitmix64(replica + 0xd1b54a32d192ed03ULL));
    SplitMixStream generator(stream_key);
    for (int stop = n - 1; stop > 0; --stop) {
        const int other = static_cast<int>(generator.below(static_cast<std::uint64_t>(stop + 1)));
        std::swap(permutation[stop], permutation[other]);
    }
}

struct RankCounts {
    std::vector<std::uint64_t> minus;
    std::vector<std::uint64_t> plus;
    std::uint64_t samples = 0;
    std::uint64_t sum_minus = 0;
    std::uint64_t sum_plus = 0;
    std::uint64_t sum_minus2 = 0;
    std::uint64_t sum_plus2 = 0;
    std::uint64_t sum_product = 0;
    std::uint64_t sum_gap = 0;
    std::uint64_t sum_gap2 = 0;

    explicit RankCounts(int n = 0) : minus(n + 1), plus(n + 1) {}

    void add(int k_minus, int k_plus) {
        if (!(1 <= k_minus && k_minus <= k_plus &&
              k_plus < static_cast<int>(plus.size()))) {
            throw std::logic_error("invalid threshold rank pair");
        }
        ++minus[k_minus];
        ++plus[k_plus];
        ++samples;
        sum_minus += k_minus;
        sum_plus += k_plus;
        sum_minus2 += static_cast<std::uint64_t>(k_minus) * k_minus;
        sum_plus2 += static_cast<std::uint64_t>(k_plus) * k_plus;
        sum_product += static_cast<std::uint64_t>(k_minus) * k_plus;
        const std::uint64_t gap = static_cast<std::uint64_t>(k_plus - k_minus);
        sum_gap += gap;
        sum_gap2 += gap * gap;
    }
};

struct PairBatch {
    RankCounts first;
    RankCounts second;
    explicit PairBatch(int n = 0) : first(n), second(n) {}
};

using BirthKey = std::tuple<int, int, Int, Int, bool>;

struct ProjectiveCounts {
    RankCounts ranks;
    std::map<BirthKey, std::uint64_t> births;

    explicit ProjectiveCounts(int n = 0) : ranks(n) {}

    void add(const ThresholdEngine::ProjectiveBirth& birth, int reverse_kminus) {
        if (birth.tau1 != reverse_kminus) {
            throw std::logic_error("projective tau1 disagrees with Alexander K_minus");
        }
        if (birth.direct_rank2) {
            if (birth.tau1 != birth.tau2 || birth.ell.x != 0 || birth.ell.y != 0) {
                throw std::logic_error("invalid DIRECT_RANK2 record");
            }
        } else {
            if (!(birth.tau1 < birth.tau2) ||
                std::gcd(safe_abs(birth.ell.x, "projective line"),
                         safe_abs(birth.ell.y, "projective line")) != 1 ||
                birth.ell.x < 0 || (birth.ell.x == 0 && birth.ell.y < 0)) {
                throw std::logic_error("invalid primitive projective line record");
            }
        }
        ranks.add(birth.tau1, birth.tau2);
        ++births[{birth.tau1, birth.tau2, birth.ell.x, birth.ell.y,
                  birth.direct_rank2}];
    }
};

struct ProjectivePairBatch {
    ProjectiveCounts first;
    ProjectiveCounts second;
    explicit ProjectivePairBatch(int n = 0) : first(n), second(n) {}
};

std::vector<int> remap_permutation(const Geometry& source, const Geometry& target,
                                   const std::vector<int>& permutation) {
    std::vector<int> mapped;
    mapped.reserve(permutation.size());
    for (const int label : permutation) {
        mapped.push_back(target.quotient.label(source.quotient.representative(label)));
    }
    return mapped;
}

void self_test() {
    const Geometry gaussian = make_geometry({2, -1, 1, 2});
    ThresholdEngine gaussian_engine(gaussian);
    RankCounts counts(gaussian.n);
    ProjectiveCounts marked_counts(gaussian.n);
    std::vector<int> permutation(gaussian.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto ranks = gaussian_engine.ranks(permutation);
        counts.add(ranks.first, ranks.second);
        const auto birth = gaussian_engine.projective_birth(permutation);
        marked_counts.add(birth, gaussian_engine.reverse_kminus(permutation));
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (counts.samples != 120 || counts.minus[3] != 120 || counts.plus[4] != 120) {
        throw std::runtime_error("N=5 all-permutation rank histogram regression failed");
    }
    if (marked_counts.ranks.minus != counts.minus || marked_counts.ranks.plus != counts.plus ||
        std::any_of(marked_counts.births.begin(), marked_counts.births.end(),
                    [](const auto& item) { return std::get<4>(item.first); })) {
        throw std::runtime_error("N=5 projective birth reconstruction failed");
    }

    const Geometry axis = make_geometry({2, 0, 0, 2});
    ThresholdEngine axis_engine(axis);
    int direct_paths = 0;
    permutation.resize(axis.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto birth = axis_engine.projective_birth(permutation);
        if (birth.direct_rank2) ++direct_paths;
        if (birth.tau1 != axis_engine.reverse_kminus(permutation)) {
            throw std::runtime_error("axis projective birth disagrees with K_minus");
        }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (direct_paths != 8) {
        throw std::runtime_error("axis L=2 DIRECT_RANK2 regression failed");
    }
    permutation.resize(gaussian.n);

    const QuotientCoordinates arbitrary({3, 1, 1, 2});
    if (arbitrary.order != 5 || arbitrary.smith1 != 1 || arbitrary.smith2 != 5) {
        throw std::runtime_error("arbitrary-matrix Smith regression failed");
    }
    for (const Vector winding : std::array<Vector, 4>{{{1, 0}, {0, 1}, {2, -3}, {-4, 7}}}) {
        const Vector displacement = arbitrary.period_vector(winding);
        const Vector recovered = arbitrary.winding(displacement.x, displacement.y);
        if (recovered.x != winding.x || recovered.y != winding.y ||
            arbitrary.label(displacement) != arbitrary.label({0, 0})) {
            throw std::runtime_error("exact winding/quotient regression failed");
        }
    }
    for (const Matrix matrix : std::array<Matrix, 6>{{
             {3, 1, 1, 2}, {2, -1, 1, 3}, {-2, 1, 1, 2},
             {4, 0, 0, 3}, {1, 4, 2, -1}, {0, 3, -2, 1}}}) {
        const QuotientCoordinates quotient(matrix);
        const Vector period0{matrix.a, matrix.c};
        const Vector period1{matrix.b, matrix.d};
        for (int label = 0; label < quotient.order; ++label) {
            const Vector point = quotient.representative(label);
            if (quotient.label(point) != label ||
                quotient.label({point.x + period0.x, point.y + period0.y}) != label ||
                quotient.label({point.x + period1.x, point.y + period1.y}) != label) {
                throw std::runtime_error("general HNF quotient-label regression failed");
            }
        }
    }

    // Same lattice after the unimodular basis change P -> P*((1,1),(0,1)).
    const Geometry first = make_geometry({3, 4, 1, 3});
    const Geometry second = make_geometry({3, 1, 1, 2});
    ThresholdEngine first_engine(first);
    ThresholdEngine second_engine(second);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto first_ranks = first_engine.ranks(permutation);
        const auto mapped = remap_permutation(first, second, permutation);
        if (first_ranks != second_engine.ranks(mapped)) {
            throw std::runtime_error("unimodular period-basis rank regression failed");
        }
    } while (std::next_permutation(permutation.begin(), permutation.end()));

    const QuotientCoordinates q260({16, -2, 2, 16});
    const QuotientCoordinates q340({18, -4, 4, 18});
    if (q260.order != 260 || q260.smith1 != 2 || q260.smith2 != 130 ||
        q340.order != 340 || q340.smith1 != 2 || q340.smith2 != 170) {
        throw std::runtime_error("norm-4 Smith regression failed");
    }

    counter_permutation(5, 17, 0, permutation);
    const std::vector<int> expected = {4, 3, 1, 0, 2};
    if (permutation != expected || gaussian_engine.ranks(permutation) != std::make_pair(3, 4)) {
        throw std::runtime_error("Python-compatible counter/permutation regression failed");
    }
    std::cout << "self-test passed: arbitrary integer periods, exact HNF quotient/winding, "
                 "basis invariance, projective birth reconstruction, DIRECT_RANK2, "
                 "Smith(2,130)/(2,170), N=5 all permutations\n";
}

struct Options {
    std::uint64_t samples = 1000000;
    int batches = 100;
    std::uint64_t seed = 20260828;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    int only_n = 0;
    std::string git_commit = "unknown";
    std::filesystem::path output_prefix;
    bool self_test = false;
    bool custom = false;
    Matrix first_matrix;
    Matrix second_matrix;
    int first_a = 0;
    int first_b = 0;
    int second_a = 0;
    int second_b = 0;
    bool first_rep_set = false;
    bool second_rep_set = false;
    bool projective_births = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N          replicas per period pair (default 1000000)\n"
        << "  --batches B          equal batches (default 100)\n"
        << "  --seed S             unsigned 64-bit seed (default 20260828)\n"
        << "  --replica-offset K   first sample counter (default 0)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --n N                predefined N=260 or N=340 norm-4 pair\n"
        << "  --first-matrix A B C D   custom first row-major period matrix\n"
        << "  --second-matrix A B C D  custom second row-major period matrix\n"
        << "  --first-rep A B      optional Gaussian lineage label in CSV\n"
        << "  --second-rep A B     optional Gaussian lineage label in CSV\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix PATH writes .hist.csv, .moments.csv, .metadata.json\n"
        << "  --projective-births  also writes sparse .births.csv with tau1,ell1,tau2\n"
        << "  --self-test          exact tiny and Smith regressions, then exit\n";
    std::exit(status);
}

template <typename T>
T parse_number(const std::string& text, const std::string& option) {
    if constexpr (std::is_unsigned_v<T>) {
        if (!text.empty() && text.front() == '-') {
            throw std::invalid_argument("negative value for " + option);
        }
    }
    std::istringstream input(text);
    T value{};
    input >> value;
    if (!input || !input.eof()) throw std::invalid_argument("invalid value for " + option);
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    bool first_matrix_set = false;
    bool second_matrix_set = false;
    auto need = [&](int& index, const std::string& option) -> std::string {
        if (++index >= argc) throw std::invalid_argument(option + " needs more values");
        return argv[index];
    };
    auto parse_matrix = [&](int& index, const std::string& option) -> Matrix {
        Matrix value;
        value.a = parse_number<Int>(need(index, option), option);
        value.b = parse_number<Int>(need(index, option), option);
        value.c = parse_number<Int>(need(index, option), option);
        value.d = parse_number<Int>(need(index, option), option);
        return value;
    };
    auto parse_rep = [&](int& index, const std::string& option) -> std::pair<int, int> {
        const int a = parse_number<int>(need(index, option), option);
        const int b = parse_number<int>(need(index, option), option);
        return {a, b};
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--samples") options.samples = parse_number<std::uint64_t>(need(i, arg), arg);
        else if (arg == "--batches") options.batches = parse_number<int>(need(i, arg), arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(need(i, arg), arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(need(i, arg), arg);
        else if (arg == "--threads") options.threads = parse_number<int>(need(i, arg), arg);
        else if (arg == "--n") options.only_n = parse_number<int>(need(i, arg), arg);
        else if (arg == "--git-commit") options.git_commit = need(i, arg);
        else if (arg == "--output-prefix") options.output_prefix = need(i, arg);
        else if (arg == "--first-matrix") { options.first_matrix = parse_matrix(i, arg); first_matrix_set = true; }
        else if (arg == "--second-matrix") { options.second_matrix = parse_matrix(i, arg); second_matrix_set = true; }
        else if (arg == "--first-rep") {
            const auto value = parse_rep(i, arg); options.first_a = value.first;
            options.first_b = value.second; options.first_rep_set = true;
        }
        else if (arg == "--second-rep") {
            const auto value = parse_rep(i, arg); options.second_a = value.first;
            options.second_b = value.second; options.second_rep_set = true;
        }
        else if (arg == "--projective-births") options.projective_births = true;
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") usage(argv[0], 0);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (options.self_test) return options;
    if (options.output_prefix.empty()) throw std::invalid_argument("--output-prefix required");
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be positive and divisible by batches>=2");
    }
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (first_matrix_set != second_matrix_set) {
        throw std::invalid_argument("custom runs require both period matrices");
    }
    options.custom = first_matrix_set;
    if (options.custom && options.only_n != 0) {
        throw std::invalid_argument("--n cannot be combined with custom matrices");
    }
    if (!options.custom && options.only_n != 260 && options.only_n != 340) {
        throw std::invalid_argument("choose predefined --n 260 or --n 340, or custom matrices");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("replica counter range overflows uint64");
    }
    return options;
}

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (const char ch : value) {
        if (ch == '\\' || ch == '"') output << '\\' << ch;
        else if (ch == '\n') output << "\\n";
        else output << ch;
    }
    return output.str();
}

std::string utc_now() {
    const std::time_t value = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
    std::tm time{};
#ifdef _WIN32
    gmtime_s(&time, &value);
#else
    gmtime_r(&value, &time);
#endif
    std::ostringstream output;
    output << std::put_time(&time, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

std::string matrix_json(const Matrix& value) {
    std::ostringstream output;
    output << "[[" << value.a << ',' << value.b << "],[" << value.c << ',' << value.d << "]]";
    return output.str();
}

PairDesign custom_design(const Options& options) {
    const QuotientCoordinates first(options.first_matrix);
    const QuotientCoordinates second(options.second_matrix);
    if (first.order != second.order) {
        throw std::invalid_argument("paired matrices must have equal absolute determinant");
    }
    auto default_rep = [](Int value) -> int {
        if (value < std::numeric_limits<int>::min() ||
            value > std::numeric_limits<int>::max()) {
            throw std::invalid_argument(
                "matrix column is outside CSV lineage-label range; pass --first-rep/--second-rep");
        }
        return static_cast<int>(value);
    };
    const int a1 = options.first_rep_set ? options.first_a : default_rep(options.first_matrix.a);
    const int b1 = options.first_rep_set ? options.first_b : default_rep(options.first_matrix.c);
    const int a2 = options.second_rep_set ? options.second_a : default_rep(options.second_matrix.a);
    const int b2 = options.second_rep_set ? options.second_b : default_rep(options.second_matrix.c);
    return {static_cast<int>(first.order), a1, b1, options.first_matrix,
            a2, b2, options.second_matrix, "custom"};
}

void run_design(const PairDesign& design, const Options& options,
                std::ofstream& histogram, std::ofstream& moments,
                std::ofstream* births) {
    const Geometry first_geometry = make_geometry(design.first);
    const Geometry second_geometry = make_geometry(design.second);
    if (first_geometry.n != design.n || second_geometry.n != design.n) {
        throw std::logic_error("design N does not match period determinant");
    }
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<ProjectivePairBatch> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) output.emplace_back(design.n);

#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        ProjectivePairBatch local(design.n);
        ThresholdEngine first_engine(first_geometry);
        ThresholdEngine second_engine(second_geometry);
        std::vector<int> permutation;
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(design.n, options.seed, replica, permutation);
            if (options.projective_births) {
                const auto first = first_engine.projective_birth(permutation);
                const auto second = second_engine.projective_birth(permutation);
                local.first.add(first, first_engine.reverse_kminus(permutation));
                local.second.add(second, second_engine.reverse_kminus(permutation));
            } else {
                const auto first = first_engine.ranks(permutation);
                const auto second = second_engine.ranks(permutation);
                local.first.ranks.add(first.first, first.second);
                local.second.ranks.add(second.first, second.second);
            }
        }
        output[batch] = std::move(local);
    }

    auto write_orientation = [&](int batch, const char* orientation, int a, int b,
                                 const ProjectiveCounts& projective) {
        const RankCounts& counts = projective.ranks;
        for (int rank = 1; rank <= design.n; ++rank) {
            if (counts.minus[rank]) {
                histogram << design.n << ',' << a << ',' << b << ',' << orientation << ','
                          << batch << ',' << counts.samples << ",minus," << rank << ','
                          << counts.minus[rank] << '\n';
            }
            if (counts.plus[rank]) {
                histogram << design.n << ',' << a << ',' << b << ',' << orientation << ','
                          << batch << ',' << counts.samples << ",plus," << rank << ','
                          << counts.plus[rank] << '\n';
            }
        }
        moments << design.n << ',' << a << ',' << b << ',' << orientation << ',' << batch
                << ',' << counts.samples << ',' << counts.sum_minus << ',' << counts.sum_plus
                << ',' << counts.sum_minus2 << ',' << counts.sum_plus2 << ','
                << counts.sum_product << ',' << counts.sum_gap << ',' << counts.sum_gap2 << '\n';
        if (births != nullptr) {
            for (const auto& item : projective.births) {
                const auto& [key, count] = item;
                const auto& [tau1, tau2, ell_x, ell_y, direct] = key;
                *births << design.n << ',' << a << ',' << b << ',' << orientation << ','
                        << batch << ',' << counts.samples << ',' << tau1 << ',' << tau2 << ','
                        << (direct ? "DIRECT_RANK2" : "LINE") << ',' << ell_x << ','
                        << ell_y << ',' << count << '\n';
            }
        }
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write_orientation(batch, "first", design.a1, design.b1, output[batch].first);
        write_orientation(batch, "second", design.a2, design.b2, output[batch].second);
    }
    std::cout << "completed " << design.id << " N=" << design.n << " samples="
              << options.samples << '\n';
}

int run(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    if (options.self_test) { self_test(); return 0; }
    PairDesign design = options.custom ? custom_design(options) : kDesigns.front();
    if (!options.custom) {
        design = *std::find_if(kDesigns.begin(), kDesigns.end(),
            [&](const PairDesign& candidate) { return candidate.n == options.only_n; });
    }

    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path histogram_path = options.output_prefix.string() + ".hist.csv";
    const std::filesystem::path moments_path = options.output_prefix.string() + ".moments.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string() + ".metadata.json";
    const std::filesystem::path births_path = options.output_prefix.string() + ".births.csv";
    std::ofstream histogram(histogram_path), moments(moments_path);
    if (!histogram || !moments) throw std::runtime_error("cannot open output files");
    histogram << "n,a,b,orientation,batch,samples,kind,k,count\n";
    moments << "n,a,b,orientation,batch,samples,sum_kminus,sum_kplus,sum_kminus2,"
               "sum_kplus2,sum_product,sum_gap,sum_gap2\n";
    std::ofstream births;
    if (options.projective_births) {
        births.open(births_path);
        if (!births) throw std::runtime_error("cannot open projective birth output");
        births << "n,a,b,orientation,batch,samples,tau1,tau2,kind,ell_x,ell_y,count\n";
    }
    const auto started = std::chrono::steady_clock::now();
    run_design(design, options, histogram, moments,
               options.projective_births ? &births : nullptr);
    histogram.close();
    moments.close();
    if (births) births.close();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    const QuotientCoordinates first(design.first), second(design.second);
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    std::ofstream metadata(metadata_path);
    if (!metadata) throw std::runtime_error("cannot open metadata output");
    metadata << "{\n"
             << "  \"engine\": \"general integer-period threshold-rank Newman-Ziff\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"command\": \"" << json_escape(command.str()) << "\",\n"
             << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
             << "  \"openmp\": "
#ifdef _OPENMP
             << "true,\n"
#else
             << "false,\n"
#endif
             << "  \"threads_requested\": " << options.threads << ",\n"
             << "  \"samples_per_pair\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": "
             << options.replica_offset + options.samples << ",\n"
             << "  \"rng\": \"counter-derived SplitMix64 stream plus unbiased Fisher-Yates\",\n"
             << "  \"coupling\": \"same HNF mixed-radix label permutation shared by the pair\",\n"
             << "  \"quotient_coordinates\": \"column-HNF representatives (rx,ry), label=rx+h11*ry\",\n"
             << "  \"period_matrix_convention\": \"row-major matrix; columns are lifted period vectors\",\n"
             << "  \"channel\": \"rank-2 cross wrapping\",\n"
             << "  \"K_plus\": \"first black primal cross rank, 1-based\",\n"
             << "  \"K_minus\": \"first black rank after white matching cross is lost; N-r+1\",\n"
             << "  \"sparse_joint_histogram\": false,\n"
             << "  \"per_batch_joint_moments\": true,\n"
             << "  \"projective_births\": " << (options.projective_births ? "true" : "false") << ",\n"
             << "  \"projective_line\": \"primitive period-basis vector canonical up to sign\",\n"
             << "  \"DIRECT_RANK2\": \"typed 0-to-2 birth with no projective line\",\n"
             << "  \"integral_saturation\": \"iota=1 by c1a72e5; no varying channel recorded\",\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"designs\": [\n"
             << "    {\"id\": \"" << json_escape(design.id) << "\", \"N\": " << design.n
             << ", \"first\": [" << design.a1 << ',' << design.b1 << "]"
             << ", \"second\": [" << design.a2 << ',' << design.b2 << "]"
             << ", \"first_period_matrix\": " << matrix_json(design.first)
             << ", \"second_period_matrix\": " << matrix_json(design.second)
             << ", \"first_HNF\": [[" << first.h11 << ',' << first.h12 << "],[0," << first.h22 << "]]"
             << ", \"second_HNF\": [[" << second.h11 << ',' << second.h12 << "],[0," << second.h22 << "]]"
             << ", \"first_smith_invariants\": [" << first.smith1 << ',' << first.smith2 << "]"
             << ", \"second_smith_invariants\": [" << second.smith1 << ',' << second.smith2 << "]}\n"
             << "  ],\n"
             << "  \"histogram_csv\": \"" << json_escape(histogram_path.string()) << "\",\n"
             << "  \"moments_csv\": \"" << json_escape(moments_path.string()) << "\",\n"
             << "  \"births_csv\": "
             << (options.projective_births ? "\"" + json_escape(births_path.string()) + "\"" : "null")
             << "\n"
             << "}\n";
    std::cout << "wrote " << histogram_path << "\nwrote " << moments_path
              << "\nwrote " << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
