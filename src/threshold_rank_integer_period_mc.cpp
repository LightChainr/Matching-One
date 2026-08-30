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
#include <cmath>
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

bool same_vector(const Vector& first, const Vector& second) {
    return first.x == second.x && first.y == second.y;
}

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

Vector scaled(Vector value, Int factor) {
    return {checked_int(static_cast<__int128>(value.x) * factor, "scaled winding"),
            checked_int(static_cast<__int128>(value.y) * factor, "scaled winding")};
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

struct LocalMark {
    bool valid = false;
    int axis = 0;
    int diagonal = 0;
    int landed = 0;
    int h4 = 0;
};

const std::array<Vector, 8> kLocalPoints = {{{1, 0}, {1, 1}, {0, 1}, {-1, 1},
                                             {-1, 0}, {-1, -1}, {0, -1}, {1, -1}}};

bool local_adjacent(Vector first, Vector second, bool matching) {
    const Int dx = safe_abs(first.x - second.x, "local dx");
    const Int dy = safe_abs(first.y - second.y, "local dy");
    return matching ? (std::max(dx, dy) == 1) : (dx + dy == 1);
}

std::vector<int> local_component_masks(const Geometry& geometry,
                                       const std::vector<std::uint8_t>& active,
                                       int root, bool matching, bool enabled_value,
                                       bool& injective) {
    const Vector origin = geometry.quotient.representative(root);
    std::array<int, 8> vertices{};
    for (int index = 0; index < 8; ++index) {
        vertices[index] = geometry.quotient.label(
            {origin.x + kLocalPoints[index].x, origin.y + kLocalPoints[index].y});
        if (vertices[index] == root) injective = false;
        for (int prior = 0; prior < index; ++prior) {
            if (vertices[index] == vertices[prior]) injective = false;
        }
    }
    if (!injective) return {};
    std::array<bool, 8> unseen{};
    for (int index = 0; index < 8; ++index) {
        unseen[index] = static_cast<bool>(active[vertices[index]]) == enabled_value;
    }
    std::vector<int> masks;
    for (int start = 0; start < 8; ++start) {
        if (!unseen[start]) continue;
        unseen[start] = false;
        std::vector<int> stack{start};
        int mask = 0;
        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            mask |= 1 << current;
            for (int other = 0; other < 8; ++other) {
                if (unseen[other] &&
                    local_adjacent(kLocalPoints[current], kLocalPoints[other], matching)) {
                    unseen[other] = false;
                    stack.push_back(other);
                }
            }
        }
        masks.push_back(mask);
    }
    return masks;
}

bool distinct_pair(const std::vector<int>& masks, int first, int second) {
    for (int i = 0; i < static_cast<int>(masks.size()); ++i) {
        for (int j = 0; j < static_cast<int>(masks.size()); ++j) {
            if (i != j && (masks[i] & (1 << first)) && (masks[j] & (1 << second))) {
                return true;
            }
        }
    }
    return false;
}

LocalMark local_landing_mark(const Geometry& geometry,
                             const std::vector<std::uint8_t>& active,
                             int root, bool open_matching) {
    bool injective = true;
    const auto opened = local_component_masks(
        geometry, active, root, open_matching, true, injective);
    const auto closed = local_component_masks(
        geometry, active, root, !open_matching, false, injective);
    if (!injective) return {};
    const bool axis =
        (distinct_pair(opened, 0, 4) && distinct_pair(closed, 2, 6)) ||
        (distinct_pair(opened, 2, 6) && distinct_pair(closed, 0, 4));
    const bool diagonal =
        (distinct_pair(opened, 1, 5) && distinct_pair(closed, 3, 7)) ||
        (distinct_pair(opened, 3, 7) && distinct_pair(closed, 1, 5));
    return {true, static_cast<int>(axis), static_cast<int>(diagonal),
            static_cast<int>(axis || diagonal), static_cast<int>(axis) - static_cast<int>(diagonal)};
}

struct PathInsertion {
    int k_before = 0;
    int q_before = -1;
    int site = -1;
    int gate01 = 0;
    int gate12 = 0;
    int components_before = 0;
    int components_after = 0;
    int euler_near = 0;
    Vector line{0, 0};
    Int index = 0;
    LocalMark local;
    LocalMark far_axis;
    LocalMark far_diagonal;
    int far_rotation = -1;
};

struct BirthTrace {
    int k1 = 0;
    int k2 = 0;
    int site1 = -1;
    int site2 = -1;
    Vector line{0, 0};
    Int index1 = 0;
    Int index2 = 0;
    LocalMark mark1;
    LocalMark mark2;
    std::vector<PathInsertion> path;
};

struct GeometryPilotRecord {
    bool at_risk = false;
    int k1 = 0;
    int k2 = 0;
    Vector line{0, 0};
    int essential_size = 0;
    int essential_carriers = 0;
    int occupied_frontier = 0;
    int vacant_frontier = 0;
    int boundary_cut_edges = 0;
    int boundary_multicontact_sites = 0;
    int boundary_contact_pairs = 0;
    int core_vertices = 0;
    int core_edges = 0;
    int articulation_vertices = 0;
    int bridges = 0;
    int boundary_axis_imbalance = 0;
    int boundary_corner_balance = 0;
    int frontier_components = 0;
    int largest_frontier_component = 0;
    int frontier_component_sumsq = 0;
    int h2 = 0;
    int h2_theta = 0;
    int h2_figure8 = 0;
    int h2_separate = 0;
    int h2_direction_positive = 0;
    int h2_direction_negative = 0;
    int h2_direction_mixed = 0;
    int next_site = -1;
    int next_exit = 0;
    int checkpoint_b1_safe_count = -1;
    int branch_common_safe = 0;
    int branch_suffix_site1 = -1;
    int branch_suffix_site2 = -1;
    int branch_clone1_survives = 0;
    int branch_clone2_survives = 0;
    int branch_both_survive = 0;
};

struct CarrierShapeMetrics {
    int boundary_cut_edges = 0;
    int boundary_multicontact_sites = 0;
    int boundary_contact_pairs = 0;
    int core_vertices = 0;
    int core_edges = 0;
    int articulation_vertices = 0;
    int bridges = 0;
    int boundary_axis_imbalance = 0;
    int boundary_corner_balance = 0;
    int frontier_components = 0;
    int largest_frontier_component = 0;
    int frontier_component_sumsq = 0;
};

CarrierShapeMetrics carrier_shape_metrics(
    int n, const std::vector<std::uint8_t>& active,
    const std::vector<std::uint8_t>& essential, const std::vector<Edge>& edges) {
    if (static_cast<int>(active.size()) != n ||
        static_cast<int>(essential.size()) != n) {
        throw std::invalid_argument("carrier shape mask length differs from N");
    }
    CarrierShapeMetrics output;
    std::vector<int> vacant_contacts(n, 0);
    std::vector<std::uint8_t> vacant_contact_directions(n, 0);
    std::vector<std::vector<std::pair<int, int>>> adjacency(n);
    std::vector<int> degree(n, 0);
    for (int edge_index = 0; edge_index < static_cast<int>(edges.size()); ++edge_index) {
        const Edge& edge = edges[edge_index];
        if (essential[edge.i] && !active[edge.j]) {
            ++output.boundary_cut_edges;
            ++vacant_contacts[edge.j];
            vacant_contact_directions[edge.j] |= edge.dx != 0 ? (1U << 2) : (1U << 3);
            output.boundary_axis_imbalance += edge.dx != 0 ? 1 : -1;
        }
        if (essential[edge.j] && !active[edge.i]) {
            ++output.boundary_cut_edges;
            ++vacant_contacts[edge.i];
            vacant_contact_directions[edge.i] |= edge.dx != 0 ? (1U << 0) : (1U << 1);
            output.boundary_axis_imbalance += edge.dx != 0 ? 1 : -1;
        }
        if (!(essential[edge.i] && essential[edge.j])) continue;
        adjacency[edge.i].push_back({edge.j, edge_index});
        adjacency[edge.j].push_back({edge.i, edge_index});
        ++degree[edge.i];
        ++degree[edge.j];
    }
    for (int vertex = 0; vertex < n; ++vertex) {
        if (!active[vertex] && vacant_contacts[vertex] >= 2) {
            ++output.boundary_multicontact_sites;
            output.boundary_contact_pairs +=
                vacant_contacts[vertex] * (vacant_contacts[vertex] - 1) / 2;
        }
        if (!active[vertex] && vacant_contacts[vertex] > 0) {
            const std::uint8_t mask = vacant_contact_directions[vertex];
            int adjacent = 0;
            for (int direction = 0; direction < 4; ++direction) {
                adjacent += static_cast<int>((mask & (1U << direction)) &&
                    (mask & (1U << ((direction + 1) % 4))));
            }
            const int opposite = static_cast<int>((mask & 1U) && (mask & 4U)) +
                                 static_cast<int>((mask & 2U) && (mask & 8U));
            output.boundary_corner_balance += adjacent - opposite;
        }
    }

    // Organization of the active boundary: connected vacant-frontier arcs in
    // the same nearest-neighbour quotient graph.  Only aggregate sizes leave
    // this traversal; no configuration or component labels are persisted.
    std::vector<std::uint8_t> frontier(n, 0), frontier_seen(n, 0);
    std::vector<std::vector<int>> frontier_adjacency(n);
    for (int vertex = 0; vertex < n; ++vertex) {
        frontier[vertex] = static_cast<std::uint8_t>(
            !active[vertex] && vacant_contacts[vertex] > 0);
    }
    for (const Edge& edge : edges) {
        if (frontier[edge.i] && frontier[edge.j] && edge.i != edge.j) {
            frontier_adjacency[edge.i].push_back(edge.j);
            frontier_adjacency[edge.j].push_back(edge.i);
        }
    }
    for (int start = 0; start < n; ++start) {
        if (!frontier[start] || frontier_seen[start]) continue;
        ++output.frontier_components;
        frontier_seen[start] = 1;
        std::vector<int> stack{start};
        int component_size = 0;
        while (!stack.empty()) {
            const int vertex = stack.back();
            stack.pop_back();
            ++component_size;
            for (const int other : frontier_adjacency[vertex]) {
                if (!frontier_seen[other]) {
                    frontier_seen[other] = 1;
                    stack.push_back(other);
                }
            }
        }
        output.largest_frontier_component =
            std::max(output.largest_frontier_component, component_size);
        output.frontier_component_sumsq += component_size * component_size;
    }

    std::vector<std::uint8_t> core = essential;
    std::vector<int> queue;
    queue.reserve(n);
    for (int vertex = 0; vertex < n; ++vertex) {
        if (core[vertex] && degree[vertex] < 2) queue.push_back(vertex);
    }
    for (std::size_t offset = 0; offset < queue.size(); ++offset) {
        const int vertex = queue[offset];
        if (!core[vertex]) continue;
        core[vertex] = 0;
        for (const auto& item : adjacency[vertex]) {
            const int other = item.first;
            if (other != vertex && core[other] && --degree[other] < 2) {
                queue.push_back(other);
            }
        }
    }
    for (const std::uint8_t value : core) output.core_vertices += value;
    for (const Edge& edge : edges) {
        if (core[edge.i] && core[edge.j]) ++output.core_edges;
    }

    std::vector<int> discovery(n, -1), low(n, -1);
    std::vector<std::uint8_t> articulation(n, 0);
    int clock = 0;
    auto visit = [&](auto&& self, int vertex, int parent_edge) -> void {
        discovery[vertex] = low[vertex] = clock++;
        int children = 0;
        for (const auto& item : adjacency[vertex]) {
            const int other = item.first;
            const int edge_index = item.second;
            if (discovery[other] < 0) {
                ++children;
                self(self, other, edge_index);
                low[vertex] = std::min(low[vertex], low[other]);
                if (parent_edge >= 0 && low[other] >= discovery[vertex]) {
                    articulation[vertex] = 1;
                }
                if (low[other] > discovery[vertex]) ++output.bridges;
            } else if (edge_index != parent_edge) {
                low[vertex] = std::min(low[vertex], discovery[other]);
            }
        }
        if (parent_edge < 0 && children > 1) articulation[vertex] = 1;
    };
    for (int vertex = 0; vertex < n; ++vertex) {
        if (essential[vertex] && discovery[vertex] < 0) visit(visit, vertex, -1);
    }
    for (const std::uint8_t value : articulation) output.articulation_vertices += value;
    return output;
}

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
    {65, 8, 1, {8, -1, 1, 8}, 7, 4, {7, -4, 4, 7}, "N65_q2_parent"},
    {130, 9, 7, {9, -7, 7, 9}, 11, 3, {11, -3, 3, 11}, "N130_q2_child"},
    {145, 12, 1, {12, -1, 1, 12}, 9, 8, {9, -8, 8, 9}, "N145_max_leverage"},
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
    struct ComponentMark {
        int rank = 0;
        Vector line{0, 0};
        Int index = 0;
    };

    explicit HomologyUnionFind(const QuotientCoordinates& quotient)
        : quotient_(quotient), parent_(quotient.order), size_(quotient.order),
          delta_x_(quotient.order), delta_y_(quotient.order), rank_(quotient.order),
          basis_(quotient.order), index_(quotient.order) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
        std::fill(index_.begin(), index_.end(), 0);
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
        if (rank_[root] == 0) {
            const Int divisor = std::gcd(safe_abs(value.x, "winding"),
                                         safe_abs(value.y, "winding"));
            basis_[root][0] = primitive(value);
            index_[root] = divisor;
            rank_[root] = 1;
            return;
        }
        const Vector first = basis_[root][0];
        if (static_cast<__int128>(first.x) * value.y ==
            static_cast<__int128>(first.y) * value.x) {
            Int coefficient = 0;
            if (first.x != 0) {
                if (value.x % first.x != 0 ||
                    value.y != (value.x / first.x) * first.y) {
                    throw std::logic_error("collinear winding is not integral on primitive line");
                }
                coefficient = value.x / first.x;
            } else {
                if (first.y == 0 || value.y % first.y != 0 || value.x != 0) {
                    throw std::logic_error("vertical winding is not integral on primitive line");
                }
                coefficient = value.y / first.y;
            }
            index_[root] = std::gcd(index_[root], safe_abs(coefficient, "line coefficient"));
        } else {
            basis_[root][1] = primitive(value);
            rank_[root] = 2;
            index_[root] = 0;
        }
    }

    bool add_edge(const Edge& edge) {
        FindResult first = find(edge.i);
        FindResult second = find(edge.j);
        Int root_dx = first.dx + edge.dx - second.dx;
        Int root_dy = first.dy + edge.dy - second.dy;
        if (first.root == second.root) {
            extend(first.root, quotient_.winding(root_dx, root_dy));
            return false;
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
            const Vector generator = index == 0 && rank_[second.root] == 1
                ? scaled(basis_[second.root][0], index_[second.root])
                : basis_[second.root][index];
            extend(first.root, generator);
        }
        rank_[second.root] = 0;
        index_[second.root] = 0;
        return true;
    }

    bool component_crosses(int vertex) { return rank_[find(vertex).root] == 2; }

    int component_size(int vertex) { return size_[find(vertex).root]; }

    ComponentMark component_mark(int vertex) {
        const int root = find(vertex).root;
        if (rank_[root] == 1) return {1, basis_[root][0], index_[root]};
        return {static_cast<int>(rank_[root]), {0, 0}, 0};
    }

  private:
    const QuotientCoordinates& quotient_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<Int> delta_x_;
    std::vector<Int> delta_y_;
    std::vector<std::uint8_t> rank_;
    std::vector<std::array<Vector, 2>> basis_;
    std::vector<Int> index_;
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

int euler_cell_residue(const Geometry& geometry,
                       const std::vector<std::uint8_t>& active) {
    if (static_cast<int>(active.size()) != geometry.n) {
        throw std::invalid_argument("Euler mask length differs from N");
    }
    int vertices = 0;
    int edges = 0;
    int faces = 0;
    for (const std::uint8_t value : active) vertices += value;
    for (const Edge& edge : geometry.primal_edges) {
        if (active[edge.i] && active[edge.j]) ++edges;
    }
    for (int root = 0; root < geometry.n; ++root) {
        const Vector point = geometry.quotient.representative(root);
        const int east = geometry.quotient.label({point.x + 1, point.y});
        const int north = geometry.quotient.label({point.x, point.y + 1});
        const int northeast = geometry.quotient.label({point.x + 1, point.y + 1});
        if (active[root] && active[east] && active[north] && active[northeast]) {
            ++faces;
        }
    }
    return vertices - edges + faces;
}

int euler_local_residue_r2(const Geometry& geometry,
                           const std::vector<std::uint8_t>& active, int root) {
    if (static_cast<int>(active.size()) != geometry.n) {
        throw std::invalid_argument("local Euler mask length differs from N");
    }
    const Vector origin = geometry.quotient.representative(root);
    int residue = 0;
    // Attribute V-E+F0 to the southwest cell anchor.  The Chebyshev-R2
    // window is translation invariant and D4 symmetric; repeated anchors on
    // tiny quotients deliberately retain multiplicity.
    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            const Vector point{origin.x + dx, origin.y + dy};
            const int anchor = geometry.quotient.label(point);
            const int east = geometry.quotient.label({point.x + 1, point.y});
            const int north = geometry.quotient.label({point.x, point.y + 1});
            const int northeast = geometry.quotient.label({point.x + 1, point.y + 1});
            const int occupied = active[anchor];
            const int east_edge = occupied && active[east];
            const int north_edge = occupied && active[north];
            const int empty_face = !occupied && !active[east] &&
                                   !active[north] && !active[northeast];
            residue += occupied - east_edge - north_edge + empty_face;
        }
    }
    return residue;
}

class ThresholdEngine {
  public:
    explicit ThresholdEngine(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n), union_find_(geometry.quotient) {}

    int first_cross(const std::vector<int>& permutation, bool matching, bool reverse) {
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        graph_components_ = 0;
        const std::vector<Edge>& edges = matching ? geometry_.matching_edges
                                                  : geometry_.primal_edges;
        const std::vector<std::vector<int>>& incident = matching
            ? geometry_.matching_incident : geometry_.primal_incident;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            active_[vertex] = 1;
            ++graph_components_;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j] && union_find_.add_edge(edge)) {
                    --graph_components_;
                }
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

    BirthTrace trace(const std::vector<int>& permutation, bool matching, bool reverse,
                     int far_radius = 0,
                     const std::vector<int>* far_rotations = nullptr) {
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        graph_components_ = 0;
        const std::vector<Edge>& edges = matching ? geometry_.matching_edges
                                                  : geometry_.primal_edges;
        const std::vector<std::vector<int>>& incident = matching
            ? geometry_.matching_incident : geometry_.primal_incident;
        BirthTrace trace;
        trace.path.reserve(geometry_.n);
        int global_rank = 0;
        Vector plateau_line{0, 0};
        Int plateau_index = 0;
        for (int offset = 0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[reverse ? geometry_.n - 1 - offset : offset];
            const LocalMark local = local_landing_mark(
                geometry_, active_, vertex, matching);
            LocalMark far_axis;
            LocalMark far_diagonal;
            int far_rotation = -1;
            if (far_radius > 0) {
                if (!far_rotations || static_cast<int>(far_rotations->size()) != geometry_.n) {
                    throw std::logic_error("separated observer lacks its C4 rotation stream");
                }
                const int black_k = reverse ? geometry_.n - 1 - offset : offset;
                far_rotation = (*far_rotations)[black_k];
                Vector axis{far_radius, 0};
                Vector diagonal{far_radius, far_radius};
                for (int turn = 0; turn < far_rotation; ++turn) {
                    axis = {-axis.y, axis.x};
                    diagonal = {-diagonal.y, diagonal.x};
                }
                const Vector root = geometry_.quotient.representative(vertex);
                const int axis_root = geometry_.quotient.label(
                    {root.x + axis.x, root.y + axis.y});
                const int diagonal_root = geometry_.quotient.label(
                    {root.x + diagonal.x, root.y + diagonal.y});
                far_axis = local_landing_mark(
                    geometry_, active_, axis_root, matching);
                far_diagonal = local_landing_mark(
                    geometry_, active_, diagonal_root, matching);
            }
            const int before_rank = global_rank;
            const int components_before = graph_components_;
            const int euler_near = euler_local_residue_r2(geometry_, active_, vertex);
            active_[vertex] = 1;
            ++graph_components_;
            for (const int edge_index : incident[vertex]) {
                const Edge& edge = edges[edge_index];
                if (active_[edge.i] && active_[edge.j] && union_find_.add_edge(edge)) {
                    --graph_components_;
                }
            }
            const HomologyUnionFind::ComponentMark component =
                union_find_.component_mark(vertex);
            if (component.rank == 2) global_rank = 2;
            else if (component.rank == 1 && global_rank == 0) global_rank = 1;

            const int gate01 = before_rank == 0 && global_rank >= 1;
            const int gate12 = before_rank <= 1 && global_rank == 2;
            if (global_rank - before_rank != gate01 + gate12) {
                throw std::logic_error("rank increment did not split into birth gates");
            }
            PathInsertion insertion;
            insertion.k_before = offset;
            insertion.q_before = before_rank - 1;
            insertion.site = vertex;
            insertion.gate01 = gate01;
            insertion.gate12 = gate12;
            insertion.components_before = components_before;
            insertion.components_after = graph_components_;
            insertion.euler_near = euler_near;
            insertion.local = local;
            insertion.far_axis = far_axis;
            insertion.far_diagonal = far_diagonal;
            insertion.far_rotation = far_rotation;

            if (gate01 && gate12) {
                trace.k1 = trace.k2 = offset + 1;
                trace.site1 = trace.site2 = vertex;
                trace.mark1 = trace.mark2 = local;
            } else if (gate01) {
                if (component.rank != 1 || component.index < 1) {
                    throw std::logic_error("strict first birth lacks a rank-one mark");
                }
                plateau_line = component.line;
                plateau_index = component.index;
                insertion.line = plateau_line;
                insertion.index = plateau_index;
                trace.k1 = offset + 1;
                trace.site1 = vertex;
                trace.line = plateau_line;
                trace.index1 = plateau_index;
                trace.mark1 = local;
            } else if (gate12) {
                if (plateau_index < 1) {
                    throw std::logic_error("strict second birth lacks the plateau mark");
                }
                insertion.line = plateau_line;
                insertion.index = plateau_index;
                trace.k2 = offset + 1;
                trace.site2 = vertex;
                trace.index2 = plateau_index;
                trace.mark2 = local;
            }
            trace.path.push_back(insertion);

            if (global_rank == 1 && component.rank == 1) {
                if (plateau_index == 0) {
                    plateau_line = component.line;
                    plateau_index = component.index;
                } else {
                    if (!same_vector(plateau_line, component.line)) {
                        throw std::logic_error(
                            "disconnected rank-one components have different ambient lines");
                    }
                    plateau_index = std::gcd(plateau_index, component.index);
                }
            }
        }
        if (!(1 <= trace.k1 && trace.k1 <= trace.k2 && trace.k2 <= geometry_.n)) {
            throw std::logic_error("marked trace did not contain both essential births");
        }
        return trace;
    }

    GeometryPilotRecord geometry_pilot(const std::vector<int>& permutation, int k0,
                                        int branch_suffix_offset1 = -1,
                                        int branch_suffix_offset2 = -1) {
        if (static_cast<int>(permutation.size()) != geometry_.n ||
            k0 <= 0 || k0 >= geometry_.n) {
            throw std::invalid_argument("invalid current-geometry pilot layer");
        }
        std::fill(active_.begin(), active_.end(), 0);
        union_find_.reset();
        graph_components_ = 0;
        int global_rank = 0;
        int k1 = 0;
        Vector plateau_line{0, 0};
        for (int offset = 0; offset < k0; ++offset) {
            const int vertex = permutation[offset];
            active_[vertex] = 1;
            ++graph_components_;
            for (const int edge_index : geometry_.primal_incident[vertex]) {
                const Edge& edge = geometry_.primal_edges[edge_index];
                if (active_[edge.i] && active_[edge.j] && union_find_.add_edge(edge)) {
                    --graph_components_;
                }
            }
            const HomologyUnionFind::ComponentMark component =
                union_find_.component_mark(vertex);
            if (component.rank == 2) {
                global_rank = 2;
            } else if (component.rank == 1) {
                if (global_rank == 0) {
                    global_rank = 1;
                    k1 = offset + 1;
                    plateau_line = component.line;
                } else if (global_rank == 1 &&
                           !same_vector(plateau_line, component.line)) {
                    // This branch should be forbidden by disjoint-carrier
                    // geometry on T2, but the ambient image definition is
                    // explicit here rather than relying on that theorem.
                    global_rank = 2;
                }
            }
        }

        GeometryPilotRecord record;
        if (global_rank != 1) return record;
        record.at_risk = true;
        record.k1 = k1;
        record.line = plateau_line;
        record.next_site = permutation[k0];

        std::map<int, HomologyUnionFind::ComponentMark> roots;
        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (!active_[vertex]) continue;
            const int root = union_find_.find(vertex).root;
            if (!roots.count(root)) roots[root] = union_find_.component_mark(vertex);
        }
        for (const auto& item : roots) {
            const HomologyUnionFind::ComponentMark& mark = item.second;
            if (mark.rank == 1) {
                if (!same_vector(mark.line, plateau_line)) {
                    throw std::logic_error("rank-one current carriers do not share ell");
                }
                ++record.essential_carriers;
                record.essential_size += union_find_.component_size(item.first);
            }
        }

        std::vector<std::uint8_t> essential(geometry_.n, 0);
        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (active_[vertex] && union_find_.component_mark(vertex).rank == 1) {
                essential[vertex] = 1;
            }
        }

        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (active_[vertex]) {
                if (!essential[vertex]) continue;
                bool frontier = false;
                for (const int edge_index : geometry_.primal_incident[vertex]) {
                    const Edge& edge = geometry_.primal_edges[edge_index];
                    const int other = edge.i == vertex ? edge.j : edge.i;
                    if (!active_[other]) frontier = true;
                }
                record.occupied_frontier += static_cast<int>(frontier);
            } else {
                bool frontier = false;
                for (const int edge_index : geometry_.primal_incident[vertex]) {
                    const Edge& edge = geometry_.primal_edges[edge_index];
                    const int other = edge.i == vertex ? edge.j : edge.i;
                    if (essential[other]) frontier = true;
                }
                record.vacant_frontier += static_cast<int>(frontier);
            }
        }

        const CarrierShapeMetrics shape = carrier_shape_metrics(
            geometry_.n, active_, essential, geometry_.primal_edges);
        record.boundary_cut_edges = shape.boundary_cut_edges;
        record.boundary_multicontact_sites = shape.boundary_multicontact_sites;
        record.boundary_contact_pairs = shape.boundary_contact_pairs;
        record.core_vertices = shape.core_vertices;
        record.core_edges = shape.core_edges;
        record.articulation_vertices = shape.articulation_vertices;
        record.bridges = shape.bridges;
        record.boundary_axis_imbalance = shape.boundary_axis_imbalance;
        record.boundary_corner_balance = shape.boundary_corner_balance;
        record.frontier_components = shape.frontier_components;
        record.largest_frontier_component = shape.largest_frontier_component;
        record.frontier_component_sumsq = shape.frontier_component_sumsq;

        for (int vertex = 0; vertex < geometry_.n; ++vertex) {
            if (active_[vertex]) continue;
            const Completion completion = one_step_completion(vertex, plateau_line);
            if (!completion.crosses) continue;
            ++record.h2;
            if (completion.type == 0) ++record.h2_theta;
            else if (completion.type == 1) ++record.h2_figure8;
            else ++record.h2_separate;
            if (completion.direction == 1) ++record.h2_direction_positive;
            else if (completion.direction == -1) ++record.h2_direction_negative;
            else ++record.h2_direction_mixed;
            if (vertex == record.next_site) record.next_exit = 1;
        }
        if (record.h2 != record.h2_theta + record.h2_figure8 + record.h2_separate ||
            record.h2 != record.h2_direction_positive +
                         record.h2_direction_negative + record.h2_direction_mixed) {
            throw std::logic_error("H2 trigger decompositions do not sum to H2");
        }
        record.checkpoint_b1_safe_count = geometry_.n - k0 - record.h2;

        const bool branching = branch_suffix_offset1 >= k0 + 1 &&
                               branch_suffix_offset1 < geometry_.n &&
                               branch_suffix_offset2 >= k0 + 1 &&
                               branch_suffix_offset2 < geometry_.n;
        if (branching) {
            record.branch_suffix_site1 = permutation[branch_suffix_offset1];
            record.branch_suffix_site2 = permutation[branch_suffix_offset2];
            record.branch_common_safe = 1 - record.next_exit;
            if (record.branch_common_safe) {
                // The common update is applied once, then this exact successor
                // state is cloned.  The two suffix sites come from independent
                // counter streams; equality across clones is allowed, as it is
                // under two independent draws from the same remaining pool.
                ThresholdEngine successor = *this;
                const HomologyUnionFind::ComponentMark after_common =
                    successor.insert_occupied(record.next_site);
                if (after_common.rank == 2 ||
                    (after_common.rank == 1 &&
                     !same_vector(after_common.line, plateau_line))) {
                    throw std::logic_error(
                        "safe common update did not preserve the plateau line");
                }
                record.branch_clone1_survives = static_cast<int>(
                    !successor.one_step_completion(
                        record.branch_suffix_site1, plateau_line).crosses);
                record.branch_clone2_survives = static_cast<int>(
                    !successor.one_step_completion(
                        record.branch_suffix_site2, plateau_line).crosses);
                record.branch_both_survive =
                    record.branch_clone1_survives * record.branch_clone2_survives;
            }
        }
        for (int offset = k0; offset < geometry_.n; ++offset) {
            const int vertex = permutation[offset];
            active_[vertex] = 1;
            ++graph_components_;
            for (const int edge_index : geometry_.primal_incident[vertex]) {
                const Edge& edge = geometry_.primal_edges[edge_index];
                if (active_[edge.i] && active_[edge.j] && union_find_.add_edge(edge)) {
                    --graph_components_;
                }
            }
            const HomologyUnionFind::ComponentMark component =
                union_find_.component_mark(vertex);
            const bool crossed = component.rank == 2 ||
                (component.rank == 1 && !same_vector(component.line, plateau_line));
            if (crossed) {
                record.k2 = offset + 1;
                break;
            }
        }
        if (record.k2 == 0 || record.next_exit != static_cast<int>(record.k2 == k0 + 1)) {
            throw std::logic_error("H2 next-site membership disagrees with realized K2");
        }
        return record;
    }

  private:
    struct Completion {
        bool crosses = false;
        int type = 2;       // 0 theta, 1 joined figure-eight, 2 separate carrier
        int direction = 0;  // +1/-1, or 0 for mixed completing signs
    };

    HomologyUnionFind::ComponentMark insert_occupied(int vertex) {
        if (active_[vertex]) throw std::logic_error("attempted duplicate insertion");
        active_[vertex] = 1;
        ++graph_components_;
        for (const int edge_index : geometry_.primal_incident[vertex]) {
            const Edge& edge = geometry_.primal_edges[edge_index];
            if (active_[edge.i] && active_[edge.j] && union_find_.add_edge(edge)) {
                --graph_components_;
            }
        }
        return union_find_.component_mark(vertex);
    }

    Completion one_step_completion(int vertex, Vector plateau_line) {
        struct RootContacts {
            HomologyUnionFind::ComponentMark mark;
            std::vector<Vector> root_positions;
        };
        std::map<int, RootContacts> contacts;
        bool touches_rank_one = false;
        for (const int edge_index : geometry_.primal_incident[vertex]) {
            const Edge& edge = geometry_.primal_edges[edge_index];
            const int other = edge.i == vertex ? edge.j : edge.i;
            if (!active_[other]) continue;
            const Vector step = edge.i == vertex
                ? Vector{edge.dx, edge.dy} : Vector{-edge.dx, -edge.dy};
            const HomologyUnionFind::FindResult found = union_find_.find(other);
            const Vector root_position{step.x - found.dx, step.y - found.dy};
            auto inserted = contacts.emplace(
                found.root,
                RootContacts{union_find_.component_mark(other), {}});
            inserted.first->second.root_positions.push_back(root_position);
            if (inserted.first->second.mark.rank == 1) touches_rank_one = true;
        }

        bool theta = false;
        bool figure8 = false;
        bool separate = false;
        bool positive = false;
        bool negative = false;
        for (const auto& item : contacts) {
            const RootContacts& root = item.second;
            if (root.root_positions.size() < 2) continue;
            const Vector anchor = root.root_positions.front();
            for (std::size_t index = 1; index < root.root_positions.size(); ++index) {
                const Vector displacement{
                    root.root_positions[index].x - anchor.x,
                    root.root_positions[index].y - anchor.y};
                const Vector winding = geometry_.quotient.winding(
                    displacement.x, displacement.y);
                const __int128 cross = static_cast<__int128>(plateau_line.x) * winding.y -
                                       static_cast<__int128>(plateau_line.y) * winding.x;
                if (cross == 0) continue;
                if (root.mark.rank == 1) theta = true;
                else if (touches_rank_one) figure8 = true;
                else separate = true;
                const Vector canonical = primitive(winding);
                const __int128 signed_cross =
                    static_cast<__int128>(plateau_line.x) * canonical.y -
                    static_cast<__int128>(plateau_line.y) * canonical.x;
                positive = positive || signed_cross > 0;
                negative = negative || signed_cross < 0;
            }
        }
        Completion result;
        result.crosses = theta || figure8 || separate;
        result.type = theta ? 0 : (figure8 ? 1 : 2);
        result.direction = positive == negative ? 0 : (positive ? 1 : -1);
        return result;
    }

    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    HomologyUnionFind union_find_;
    int graph_components_ = 0;
};

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

void separated_rotations(int n, std::uint64_t seed, std::uint64_t replica,
                         std::vector<int>& rotations) {
    rotations.resize(n);
    const std::uint64_t stream = splitmix64(
        seed ^ splitmix64(replica + 0x4f1bbcdc6765d7f5ULL));
    for (int k = 0; k < n; ++k) {
        rotations[k] = static_cast<int>(
            splitmix64(stream ^ splitmix64(static_cast<std::uint64_t>(k))) & 3ULL);
    }
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

struct Spin4 {
    long double real = 0;
    long double imaginary = 0;
    Vector physical{0, 0};
};

Spin4 spin4_mark(const QuotientCoordinates& quotient, Vector line) {
    if (line.x == 0 && line.y == 0) return {};
    const Vector physical = primitive(quotient.period_vector(line));
    const long double x = static_cast<long double>(physical.x);
    const long double y = static_cast<long double>(physical.y);
    const long double radius2 = x * x + y * y;
    const long double denominator = radius2 * radius2;
    return {(x * x * x * x - 6 * x * x * y * y + y * y * y * y) / denominator,
            (4 * x * y * (x * x - y * y)) / denominator,
            physical};
}

bool same_local(const LocalMark& first, const LocalMark& second) {
    return first.valid == second.valid && first.axis == second.axis &&
           first.diagonal == second.diagonal && first.landed == second.landed &&
           first.h4 == second.h4;
}

struct MarkedKey {
    int k1 = 0;
    int k2 = 0;
    int site1 = -1;
    int site2 = -1;
    Vector line{0, 0};
    Int index1 = 0;
    Int index2 = 0;
    Vector physical{0, 0};
    LocalMark mark1;
    LocalMark mark2;

    bool operator<(const MarkedKey& other) const {
        return std::tie(k1, k2, site1, site2, line.x, line.y, index1, index2,
                        physical.x, physical.y, mark1.valid, mark1.axis,
                        mark1.diagonal, mark1.landed, mark1.h4, mark2.valid,
                        mark2.axis, mark2.diagonal, mark2.landed, mark2.h4) <
               std::tie(other.k1, other.k2, other.site1, other.site2,
                        other.line.x, other.line.y, other.index1, other.index2,
                        other.physical.x, other.physical.y, other.mark1.valid,
                        other.mark1.axis, other.mark1.diagonal, other.mark1.landed,
                        other.mark1.h4, other.mark2.valid, other.mark2.axis,
                        other.mark2.diagonal, other.mark2.landed, other.mark2.h4);
    }
};

struct PathMoments {
    std::uint64_t samples = 0;
    std::int64_t sum_q = 0;
    std::uint64_t sum_q2 = 0;
    std::uint64_t sum_gate01 = 0;
    std::uint64_t sum_gate12 = 0;
    std::uint64_t sum_inactive_gate01 = 0;
    std::uint64_t sum_inactive_gate12 = 0;
    std::uint64_t sum_active_S = 0;
    std::int64_t sum_active_D = 0;
    std::uint64_t sum_inactive_S = 0;
    std::int64_t sum_inactive_D = 0;
    std::uint64_t sum_site_S = 0;
    std::int64_t sum_site_D = 0;
    long double sum_J_S_real = 0;
    long double sum_J_S_imaginary = 0;
    long double sum_J_D_real = 0;
    long double sum_J_D_imaginary = 0;
    long double sum_q_J_D_real = 0;
    long double sum_q_J_D_imaginary = 0;
    std::int64_t sum_O_ext = 0;
    std::uint64_t sum_O_ext2 = 0;
    std::int64_t sum_O_near = 0;
    std::uint64_t sum_O_near2 = 0;
    std::int64_t sum_O_ext_O_near = 0;
    long double sum_O_ext_J_S_real = 0;
    long double sum_O_ext_J_S_imaginary = 0;
    long double sum_O_ext_J_D_real = 0;
    long double sum_O_ext_J_D_imaginary = 0;
    long double sum_O_near_J_S_real = 0;
    long double sum_O_near_J_S_imaginary = 0;
    long double sum_O_near_J_D_real = 0;
    long double sum_O_near_J_D_imaginary = 0;
    long double sum_J_D_conj_J_S_real = 0;
    long double sum_J_D_conj_J_S_imaginary = 0;
    long double sum_abs_J_S2 = 0;
    std::int64_t sum_local_S = 0;
    std::int64_t sum_local_D = 0;
    long double sum_O_sep_axis = 0;
    long double sum_O_sep_diagonal = 0;
    long double sum_O_sep4 = 0;
    long double sum_O_sep4_2 = 0;
    long double sum_O_sep_axis_internal_h4 = 0;
    long double sum_O_sep_diagonal_internal_h4 = 0;
    long double sum_O_sep4_J_S_real = 0;
    long double sum_O_sep4_J_S_imaginary = 0;
    long double sum_O_sep4_J_D_real = 0;
    long double sum_O_sep4_J_D_imaginary = 0;

    void add(const PathInsertion& active, const PathInsertion& inactive,
             const QuotientCoordinates& quotient, int n) {
        if (active.site != inactive.site ||
            active.gate01 != inactive.gate12 ||
            active.gate12 != inactive.gate01 ||
            !same_vector(active.line, inactive.line)) {
            throw std::logic_error("active/inactive insertion source did not complement-pair");
        }
        if (active.k_before + inactive.k_before != n - 1) {
            throw std::logic_error("active/inactive insertion ranks did not reverse-pair");
        }
        const std::uint64_t absent = static_cast<std::uint64_t>(n - active.k_before);
        const int active_s = active.gate01 + active.gate12;
        const int active_d = active.gate12 - active.gate01;
        const int inactive_s = inactive.gate01 + inactive.gate12;
        const int inactive_d = inactive.gate12 - inactive.gate01;
        // q=(r_primal-r_matching)/2.  Under black insertion the matching
        // insertion is traversed backwards, hence the minus sign in full D.
        const int full_s = (active_s + inactive_s) / 2;
        const int full_d = (active_d - inactive_d) / 2;
        if (2 * full_s != active_s + inactive_s ||
            2 * full_d != active_d - inactive_d) {
            throw std::logic_error("full active/inactive source was not integral");
        }
        const Spin4 active_chi = spin4_mark(quotient, active.line);
        const Spin4 inactive_chi = spin4_mark(quotient, inactive.line);
        const long double site_s = static_cast<long double>(absent) * full_s;
        const long double site_d = static_cast<long double>(absent) * full_d;
        const long double chi_real = (active_chi.real + inactive_chi.real) / 2;
        const long double chi_imaginary = (active_chi.imaginary + inactive_chi.imaginary) / 2;
        // This is a configuration observer, not a function of ambient rank:
        //
        //   O_ext = C_black^NN - C_white^matching - q.
        //
        // The black component count is read before inserting the next site.
        // The reverse matching trace is read after inserting that same site,
        // which is exactly the white complement of the black configuration.
        // On square cellulations O_ext=V-E+F0 configuration by configuration.
        const std::int64_t o_ext =
            static_cast<std::int64_t>(active.components_before) -
            static_cast<std::int64_t>(inactive.components_after) - active.q_before;
        const std::int64_t o_near = active.euler_near;
        const long double j_s_real = site_s * chi_real;
        const long double j_s_imaginary = site_s * chi_imaginary;
        const long double j_d_real = site_d * chi_real;
        const long double j_d_imaginary = site_d * chi_imaginary;
        ++samples;
        sum_q += active.q_before;
        sum_q2 += static_cast<std::uint64_t>(active.q_before * active.q_before);
        sum_gate01 += active.gate01;
        sum_gate12 += active.gate12;
        sum_inactive_gate01 += inactive.gate01;
        sum_inactive_gate12 += inactive.gate12;
        sum_active_S += absent * static_cast<std::uint64_t>(active_s);
        sum_active_D += static_cast<std::int64_t>(absent) * active_d;
        sum_inactive_S += absent * static_cast<std::uint64_t>(inactive_s);
        sum_inactive_D += static_cast<std::int64_t>(absent) * inactive_d;
        sum_site_S += absent * static_cast<std::uint64_t>(full_s);
        sum_site_D += static_cast<std::int64_t>(absent) * full_d;
        sum_J_S_real += j_s_real;
        sum_J_S_imaginary += j_s_imaginary;
        sum_J_D_real += j_d_real;
        sum_J_D_imaginary += j_d_imaginary;
        sum_q_J_D_real += active.q_before * j_d_real;
        sum_q_J_D_imaginary += active.q_before * j_d_imaginary;
        sum_O_ext += o_ext;
        sum_O_ext2 += static_cast<std::uint64_t>(o_ext * o_ext);
        sum_O_near += o_near;
        sum_O_near2 += static_cast<std::uint64_t>(o_near * o_near);
        sum_O_ext_O_near += o_ext * o_near;
        sum_O_ext_J_S_real += o_ext * j_s_real;
        sum_O_ext_J_S_imaginary += o_ext * j_s_imaginary;
        sum_O_ext_J_D_real += o_ext * j_d_real;
        sum_O_ext_J_D_imaginary += o_ext * j_d_imaginary;
        sum_O_near_J_S_real += o_near * j_s_real;
        sum_O_near_J_S_imaginary += o_near * j_s_imaginary;
        sum_O_near_J_D_real += o_near * j_d_real;
        sum_O_near_J_D_imaginary += o_near * j_d_imaginary;
        sum_J_D_conj_J_S_real +=
            j_d_real * j_s_real + j_d_imaginary * j_s_imaginary;
        sum_J_D_conj_J_S_imaginary +=
            j_d_imaginary * j_s_real - j_d_real * j_s_imaginary;
        sum_abs_J_S2 += j_s_real * j_s_real + j_s_imaginary * j_s_imaginary;
        if (active.local.valid && inactive.local.valid) {
            const int local_s = (active_s * active.local.h4 +
                                 inactive_s * inactive.local.h4) / 2;
            const int local_d = (active_d * active.local.h4 -
                                 inactive_d * inactive.local.h4) / 2;
            sum_local_S += static_cast<std::int64_t>(absent) * local_s;
            sum_local_D += static_cast<std::int64_t>(absent) * local_d;
        }
        if (active.far_axis.valid && active.far_diagonal.valid &&
            inactive.far_axis.valid && inactive.far_diagonal.valid) {
            // The two direction orbits are deliberately retained rather than
            // collapsed. Packed as axis+i*diagonal they form the exact rank-2
            // scalar/spin-4 response basis. O_sep4 is twice the normalized H4
            // projection, so it stays integral on every complement-paired path.
            const long double axis =
                (active.far_axis.landed + inactive.far_axis.landed) / 2.0L;
            const long double diagonal =
                (active.far_diagonal.landed + inactive.far_diagonal.landed) / 2.0L;
            const long double o_sep4 = axis - diagonal;
            sum_O_sep_axis += axis;
            sum_O_sep_diagonal += diagonal;
            sum_O_sep4 += o_sep4;
            sum_O_sep4_2 += o_sep4 * o_sep4;
            sum_O_sep_axis_internal_h4 +=
                (active.far_axis.h4 + inactive.far_axis.h4) / 2.0L;
            sum_O_sep_diagonal_internal_h4 +=
                (active.far_diagonal.h4 + inactive.far_diagonal.h4) / 2.0L;
            sum_O_sep4_J_S_real += o_sep4 * j_s_real;
            sum_O_sep4_J_S_imaginary += o_sep4 * j_s_imaginary;
            sum_O_sep4_J_D_real += o_sep4 * j_d_real;
            sum_O_sep4_J_D_imaginary += o_sep4 * j_d_imaginary;
        }
    }
};

struct ComplementAudit {
    std::uint64_t endpoint_failures = 0;
    std::uint64_t site_failures = 0;
    std::uint64_t line_failures = 0;
    std::uint64_t local_mark_failures = 0;
    std::uint64_t index_mismatches = 0;
    std::uint64_t separated_mark_failures = 0;

    void add(const BirthTrace& primal, const BirthTrace& matching, int n) {
        if (primal.k1 + matching.k2 != n + 1 ||
            primal.k2 + matching.k1 != n + 1) ++endpoint_failures;
        if (primal.site1 != matching.site2 || primal.site2 != matching.site1) {
            ++site_failures;
        }
        if (!same_vector(primal.line, matching.line)) ++line_failures;
        if (!same_local(primal.mark1, matching.mark2) ||
            !same_local(primal.mark2, matching.mark1)) ++local_mark_failures;
        if (primal.index1 != matching.index2 || primal.index2 != matching.index1) {
            ++index_mismatches;
        }
        if (primal.path.size() == matching.path.size()) {
            for (const PathInsertion& active : primal.path) {
                const PathInsertion& inactive = matching.path[n - 1 - active.k_before];
                if (active.far_rotation != inactive.far_rotation ||
                    !same_local(active.far_axis, inactive.far_axis) ||
                    !same_local(active.far_diagonal, inactive.far_diagonal)) {
                    ++separated_mark_failures;
                }
            }
        }
    }
};

struct MarkedOrientation {
    std::vector<PathMoments> path;
    std::map<MarkedKey, std::uint64_t> sparse;
    ComplementAudit complement;

    explicit MarkedOrientation(int n = 0) : path(n) {}

    void add(const BirthTrace& primal, const BirthTrace& matching,
             const Geometry& geometry) {
        if (static_cast<int>(primal.path.size()) != geometry.n) {
            throw std::logic_error("marked path length differs from N");
        }
        const Spin4 chi = spin4_mark(geometry.quotient, primal.line);
        MarkedKey key;
        key.k1 = primal.k1;
        key.k2 = primal.k2;
        key.site1 = primal.site1;
        key.site2 = primal.site2;
        key.line = primal.line;
        key.index1 = primal.index1;
        key.index2 = primal.index2;
        key.physical = chi.physical;
        key.mark1 = primal.mark1;
        key.mark2 = primal.mark2;
        ++sparse[key];
        for (const PathInsertion& insertion : primal.path) {
            const PathInsertion& inactive = matching.path[geometry.n - 1 - insertion.k_before];
            path[insertion.k_before].add(
                insertion, inactive, geometry.quotient, geometry.n);
        }
        complement.add(primal, matching, geometry.n);
    }
};

struct MarkedPairBatch {
    MarkedOrientation first;
    MarkedOrientation second;
    explicit MarkedPairBatch(int n = 0) : first(n), second(n) {}
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
    std::vector<int> permutation(gaussian.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const auto ranks = gaussian_engine.ranks(permutation);
        counts.add(ranks.first, ranks.second);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (counts.samples != 120 || counts.minus[3] != 120 || counts.plus[4] != 120) {
        throw std::runtime_error("N=5 all-permutation rank histogram regression failed");
    }

    // The virtual H2 oracle must equal explicit one-step replay.  N=5 has a
    // rank-one plateau after three sites for every permutation, so all 120
    // orders exercise the current-state path without a sampling caveat.
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        const GeometryPilotRecord pilot = gaussian_engine.geometry_pilot(permutation, 3);
        if (!pilot.at_risk || pilot.k1 != 3) {
            throw std::runtime_error("N=5 geometry-pilot risk-set regression failed");
        }
        int brute_h2 = 0;
        for (int candidate = 3; candidate < gaussian.n; ++candidate) {
            std::vector<int> replay = permutation;
            std::swap(replay[3], replay[candidate]);
            const BirthTrace trace = gaussian_engine.trace(replay, false, false);
            brute_h2 += static_cast<int>(trace.k2 == 4);
        }
        const BirthTrace actual = gaussian_engine.trace(permutation, false, false);
        if (pilot.h2 != brute_h2 || pilot.k2 != actual.k2 ||
            pilot.next_exit != static_cast<int>(actual.k2 == 4)) {
            throw std::runtime_error("virtual H2 differs from one-step replay");
        }
    } while (std::next_permutation(permutation.begin(), permutation.end()));

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

    // Saturation index must use raw winding coefficients before primitive reduction.
    const QuotientCoordinates index_quotient({4, 0, 0, 4});
    HomologyUnionFind index_union(index_quotient);
    index_union.add_edge({0, 0, 8, 0});  // winding (2,0)
    auto index_mark = index_union.component_mark(0);
    if (index_mark.rank != 1 || !same_vector(index_mark.line, {1, 0}) ||
        index_mark.index != 2) {
        throw std::runtime_error("raw-winding saturation-index regression failed");
    }
    index_union.add_edge({0, 0, 16, 0});  // winding (4,0), gcd remains two
    index_mark = index_union.component_mark(0);
    if (index_mark.index != 2) {
        throw std::runtime_error("same-line saturation gcd regression failed");
    }
    index_union.add_edge({0, 0, 4, 0});  // winding (1,0), saturation becomes one
    index_mark = index_union.component_mark(0);
    if (index_mark.index != 1) {
        throw std::runtime_error("saturation-to-primitive regression failed");
    }

    // chi4 is evaluated in the lifted Euclidean P*ell frame, not period coordinates.
    const Spin4 gaussian_chi = spin4_mark(gaussian.quotient, {1, 0});
    if (!same_vector(gaussian_chi.physical, {2, 1}) ||
        std::fabs(gaussian_chi.real + static_cast<long double>(7) / 25) > 1e-18L ||
        std::fabs(gaussian_chi.imaginary - static_cast<long double>(24) / 25) > 1e-18L) {
        throw std::runtime_error("lifted-Euclidean chi4 frame regression failed");
    }

    // The degenerate N=4 control contains direct 0->2 births.  They have two
    // gates, no canonical line, D=0, and reverse-complement endpoint/site exchange.
    const Geometry axis_l2 = make_geometry({2, 0, 0, 2});
    ThresholdEngine axis_l2_engine(axis_l2);
    std::vector<int> axis_permutation(axis_l2.n);
    std::iota(axis_permutation.begin(), axis_permutation.end(), 0);
    std::uint64_t direct_births = 0;
    do {
        const BirthTrace primal = axis_l2_engine.trace(axis_permutation, false, false);
        const BirthTrace matching = axis_l2_engine.trace(axis_permutation, true, true);
        ComplementAudit audit;
        audit.add(primal, matching, axis_l2.n);
        if (audit.endpoint_failures || audit.site_failures || audit.line_failures) {
            throw std::runtime_error("direct-birth complement regression failed");
        }
        if (primal.k1 == primal.k2) {
            ++direct_births;
            const PathInsertion& insertion = primal.path[primal.k1 - 1];
            if (insertion.gate01 != 1 || insertion.gate12 != 1 ||
                insertion.line.x != 0 || insertion.line.y != 0 ||
                insertion.index != 0) {
                throw std::runtime_error("direct 0->2 marked schema regression failed");
            }
        }
        std::vector<std::uint8_t> active(axis_l2.n, 0);
        for (int k = 0; k < axis_l2.n; ++k) {
            PathMoments one;
            one.add(primal.path[k], matching.path[axis_l2.n - 1 - k],
                    axis_l2.quotient, axis_l2.n);
            if (one.sum_O_ext != euler_cell_residue(axis_l2, active) ||
                one.sum_O_near != euler_local_residue_r2(
                    axis_l2, active, axis_permutation[k]) ||
                std::fabs(one.sum_J_D_conj_J_S_imaginary) > 1e-18L) {
                throw std::runtime_error("external/local Euler/Gram path regression failed");
            }
            active[axis_permutation[k]] = 1;
        }
    } while (std::next_permutation(axis_permutation.begin(), axis_permutation.end()));
    if (direct_births == 0) {
        throw std::runtime_error("N=4 direct-birth control found no simultaneous births");
    }

    // One-pass bottleneck summaries are checked against a direct tiny graph:
    // a four-cycle with one leaf and one vacant site touching two cycle sites.
    const std::vector<Edge> shape_edges = {
        {0, 1, 0, 0}, {1, 2, 0, 0}, {2, 3, 0, 0}, {3, 0, 0, 0},
        {3, 4, 0, 0}, {0, 5, 1, 0}, {1, 5, 0, 1},
        {2, 6, 1, 0}, {5, 6, 0, 1}};
    const std::vector<std::uint8_t> shape_active = {1, 1, 1, 1, 1, 0, 0};
    const std::vector<std::uint8_t> shape_essential = {1, 1, 1, 1, 1, 0, 0};
    const CarrierShapeMetrics shape = carrier_shape_metrics(
        7, shape_active, shape_essential, shape_edges);
    auto components_after = [&](int removed_vertex, int removed_edge) {
        std::vector<std::uint8_t> seen(7, 0);
        int components = 0;
        for (int start = 0; start < 7; ++start) {
            if (!shape_essential[start] || start == removed_vertex || seen[start]) continue;
            ++components;
            seen[start] = 1;
            std::vector<int> stack{start};
            while (!stack.empty()) {
                const int vertex = stack.back();
                stack.pop_back();
                for (int index = 0; index < static_cast<int>(shape_edges.size()); ++index) {
                    if (index == removed_edge) continue;
                    const Edge& edge = shape_edges[index];
                    if (!(shape_essential[edge.i] && shape_essential[edge.j])) continue;
                    int other = -1;
                    if (edge.i == vertex) other = edge.j;
                    else if (edge.j == vertex) other = edge.i;
                    if (other >= 0 && other != removed_vertex && !seen[other]) {
                        seen[other] = 1;
                        stack.push_back(other);
                    }
                }
            }
        }
        return components;
    };
    const int base_components = components_after(-1, -1);
    int brute_articulations = 0;
    for (int vertex = 0; vertex < 5; ++vertex) {
        brute_articulations += components_after(vertex, -1) > base_components;
    }
    int brute_bridges = 0;
    for (int edge = 0; edge < 5; ++edge) {
        brute_bridges += components_after(-1, edge) > base_components;
    }
    if (shape.boundary_cut_edges != 3 || shape.boundary_multicontact_sites != 1 ||
        shape.boundary_contact_pairs != 1 || shape.core_vertices != 4 ||
        shape.core_edges != 4 || shape.articulation_vertices != brute_articulations ||
        shape.bridges != brute_bridges || brute_articulations != 1 || brute_bridges != 1 ||
        shape.boundary_axis_imbalance != 1 || shape.boundary_corner_balance != 1 ||
        shape.frontier_components != 1 || shape.largest_frontier_component != 2 ||
        shape.frontier_component_sumsq != 4) {
        throw std::runtime_error("tiny carrier bottleneck/2-core oracle failed");
    }
    std::cout << "self-test passed: arbitrary integer periods, exact HNF quotient/winding, "
                 "basis invariance, saturation gcd, lifted chi4, direct births, "
                 "Euler external observer/Gram, exact virtual H2, carrier bottleneck/core, "
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
    bool marked_births = false;
    int far_radius = 0;
    int geometry_pilot_k0 = 0;
    bool branching_clones = false;
    bool custom = false;
    Matrix first_matrix;
    Matrix second_matrix;
    int first_a = 0;
    int first_b = 0;
    int second_a = 0;
    int second_b = 0;
    bool first_rep_set = false;
    bool second_rep_set = false;
};

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " [options]\n"
        << "  --samples N          replicas per period pair (default 1000000)\n"
        << "  --batches B          equal batches (default 100)\n"
        << "  --seed S             unsigned 64-bit seed (default 20260828)\n"
        << "  --replica-offset K   first sample counter (default 0)\n"
        << "  --threads T          OpenMP threads; 0 uses runtime default\n"
        << "  --n N                predefined N=65,130,145 marked or N=260,340 norm-4 pair\n"
        << "  --first-matrix A B C D   custom first row-major period matrix\n"
        << "  --second-matrix A B C D  custom second row-major period matrix\n"
        << "  --first-rep A B      optional Gaussian lineage label in CSV\n"
        << "  --second-rep A B     optional Gaussian lineage label in CSV\n"
        << "  --git-commit SHA     provenance string\n"
        << "  --output-prefix PATH writes .hist.csv, .moments.csv, .metadata.json\n"
        << "  --marked-births      also writes sparse birth and microcanonical path streams\n"
        << "  --far-radius R       add paired axis/diagonal local-arm observer at distance R\n"
        << "  --geometry-pilot-k0 K  write current rank-one geometry rows after K sites\n"
        << "  --branching-clones   after the common next update, draw two independent suffix sites\n"
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
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--marked-births") options.marked_births = true;
        else if (arg == "--far-radius") options.far_radius = parse_number<int>(need(i, arg), arg);
        else if (arg == "--geometry-pilot-k0") {
            options.geometry_pilot_k0 = parse_number<int>(need(i, arg), arg);
        }
        else if (arg == "--branching-clones") options.branching_clones = true;
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
    if (options.far_radius < 0 || (options.far_radius > 0 && !options.marked_births)) {
        throw std::invalid_argument("far-radius must be nonnegative and requires --marked-births");
    }
    if (options.geometry_pilot_k0 < 0 ||
        (options.geometry_pilot_k0 > 0 && options.marked_births)) {
        throw std::invalid_argument(
            "geometry-pilot-k0 must be nonnegative and is exclusive with --marked-births");
    }
    if (options.branching_clones && options.geometry_pilot_k0 <= 0) {
        throw std::invalid_argument("branching-clones requires --geometry-pilot-k0");
    }
    if (first_matrix_set != second_matrix_set) {
        throw std::invalid_argument("custom runs require both period matrices");
    }
    options.custom = first_matrix_set;
    if (options.custom && options.only_n != 0) {
        throw std::invalid_argument("--n cannot be combined with custom matrices");
    }
    if (!options.custom && options.only_n != 65 && options.only_n != 130 &&
        options.only_n != 145 && options.only_n != 260 && options.only_n != 340) {
        throw std::invalid_argument(
            "choose predefined --n 65, 130, 145, 260 or 340, or custom matrices");
    }
    if (options.geometry_pilot_k0 >=
        static_cast<int>(options.custom
            ? QuotientCoordinates(options.first_matrix).order : options.only_n)) {
        throw std::invalid_argument("geometry-pilot-k0 must be strictly below N");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max() - options.samples) {
        throw std::invalid_argument("replica counter range overflows uint64");
    }
    return options;
}

struct PilotRow {
    std::uint64_t replica = 0;
    GeometryPilotRecord record;
};

void run_geometry_pilot_design(const PairDesign& design, const Options& options,
                               std::ofstream& output) {
    const Geometry first_geometry = make_geometry(design.first);
    const Geometry second_geometry = make_geometry(design.second);
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<std::vector<PilotRow>> first_rows(options.batches);
    std::vector<std::vector<PilotRow>> second_rows(options.batches);
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        ThresholdEngine first_engine(first_geometry);
        ThresholdEngine second_engine(second_geometry);
        std::vector<int> permutation;
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch) * per_batch;
        first_rows[batch].reserve(per_batch);
        second_rows[batch].reserve(per_batch);
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(design.n, options.seed, replica, permutation);
            int suffix_offset1 = -1;
            int suffix_offset2 = -1;
            if (options.branching_clones) {
                const std::uint64_t remaining = static_cast<std::uint64_t>(
                    design.n - options.geometry_pilot_k0 - 1);
                if (remaining == 0) {
                    throw std::logic_error("branching pilot has no suffix site");
                }
                SplitMixStream suffix1(splitmix64(
                    options.seed ^ splitmix64(replica ^ 0x429b2c101ULL)));
                SplitMixStream suffix2(splitmix64(
                    options.seed ^ splitmix64(replica ^ 0x429b2c202ULL)));
                suffix_offset1 = options.geometry_pilot_k0 + 1 +
                    static_cast<int>(suffix1.below(remaining));
                suffix_offset2 = options.geometry_pilot_k0 + 1 +
                    static_cast<int>(suffix2.below(remaining));
            }
            const GeometryPilotRecord first =
                first_engine.geometry_pilot(permutation, options.geometry_pilot_k0,
                                            suffix_offset1, suffix_offset2);
            const GeometryPilotRecord second =
                second_engine.geometry_pilot(permutation, options.geometry_pilot_k0,
                                             suffix_offset1, suffix_offset2);
            if (first.at_risk) first_rows[batch].push_back({replica, first});
            if (second.at_risk) second_rows[batch].push_back({replica, second});
        }
    }

    auto write_rows = [&](int batch, const char* orientation, int a, int b,
                          const std::vector<PilotRow>& rows) {
        for (const PilotRow& row : rows) {
            const GeometryPilotRecord& value = row.record;
            output << design.n << ',' << a << ',' << b << ',' << orientation << ','
                   << batch << ',' << row.replica << ',' << options.geometry_pilot_k0 << ','
                   << value.k1 << ',' << value.k2 << ','
                   << options.geometry_pilot_k0 - value.k1 << ','
                   << value.line.x << ',' << value.line.y << ','
                   << value.essential_size << ',' << value.essential_carriers << ','
                   << value.occupied_frontier << ',' << value.vacant_frontier << ','
                   << value.boundary_cut_edges << ','
                   << value.boundary_multicontact_sites << ','
                   << value.boundary_contact_pairs << ',' << value.core_vertices << ','
                   << value.core_edges << ',' << value.articulation_vertices << ','
                   << value.bridges << ',' << value.boundary_axis_imbalance << ','
                   << value.boundary_corner_balance << ',' << value.frontier_components << ','
                   << value.largest_frontier_component << ','
                   << value.frontier_component_sumsq << ','
                   << value.h2 << ',' << value.h2_theta << ',' << value.h2_figure8 << ','
                   << value.h2_separate << ',' << value.h2_direction_positive << ','
                   << value.h2_direction_negative << ',' << value.h2_direction_mixed << ','
                   << value.next_site << ',' << value.next_exit << ','
                   << value.checkpoint_b1_safe_count << ','
                   << value.branch_common_safe << ',' << value.branch_suffix_site1 << ','
                   << value.branch_suffix_site2 << ','
                   << value.branch_clone1_survives << ','
                   << value.branch_clone2_survives << ','
                   << value.branch_both_survive << '\n';
        }
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write_rows(batch, "first", design.a1, design.b1, first_rows[batch]);
        write_rows(batch, "second", design.a2, design.b2, second_rows[batch]);
    }
    std::cout << "completed current-k0 geometry pilot " << design.id
              << " N=" << design.n << " k0=" << options.geometry_pilot_k0
              << " samples=" << options.samples << '\n';
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
                std::ofstream* marked_births, std::ofstream* path_moments,
                std::ofstream* complement_audit) {
    const Geometry first_geometry = make_geometry(design.first);
    const Geometry second_geometry = make_geometry(design.second);
    if (first_geometry.n != design.n || second_geometry.n != design.n) {
        throw std::logic_error("design N does not match period determinant");
    }
    const std::uint64_t per_batch = options.samples / options.batches;
    std::vector<PairBatch> output;
    output.reserve(options.batches);
    for (int batch = 0; batch < options.batches; ++batch) output.emplace_back(design.n);
    std::vector<MarkedPairBatch> marked_output;
    if (options.marked_births) {
        marked_output.reserve(options.batches);
        for (int batch = 0; batch < options.batches; ++batch) {
            marked_output.emplace_back(design.n);
        }
    }

#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        PairBatch local(design.n);
        ThresholdEngine first_engine(first_geometry);
        ThresholdEngine second_engine(second_geometry);
        std::vector<int> permutation;
        std::vector<int> far_rotations;
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch) * per_batch;
        for (std::uint64_t replica = begin; replica < begin + per_batch; ++replica) {
            counter_permutation(design.n, options.seed, replica, permutation);
            if (options.marked_births) {
                if (options.far_radius > 0) {
                    separated_rotations(design.n, options.seed, replica, far_rotations);
                }
                const auto* rotations = options.far_radius > 0 ? &far_rotations : nullptr;
                const BirthTrace first_primal = first_engine.trace(
                    permutation, false, false, options.far_radius, rotations);
                const BirthTrace first_matching = first_engine.trace(
                    permutation, true, true, options.far_radius, rotations);
                const BirthTrace second_primal = second_engine.trace(
                    permutation, false, false, options.far_radius, rotations);
                const BirthTrace second_matching = second_engine.trace(
                    permutation, true, true, options.far_radius, rotations);
                if (first_primal.k1 != design.n - first_matching.k2 + 1 ||
                    first_primal.k2 != design.n - first_matching.k1 + 1 ||
                    second_primal.k1 != design.n - second_matching.k2 + 1 ||
                    second_primal.k2 != design.n - second_matching.k1 + 1) {
                    throw std::logic_error("direct and reverse-complement birth ranks differ");
                }
                local.first.add(first_primal.k1, first_primal.k2);
                local.second.add(second_primal.k1, second_primal.k2);
                marked_output[batch].first.add(
                    first_primal, first_matching, first_geometry);
                marked_output[batch].second.add(
                    second_primal, second_matching, second_geometry);
            } else {
                const auto first = first_engine.ranks(permutation);
                const auto second = second_engine.ranks(permutation);
                local.first.add(first.first, first.second);
                local.second.add(second.first, second.second);
            }
        }
        output[batch] = std::move(local);
    }

    auto write_orientation = [&](int batch, const char* orientation, int a, int b,
                                 const RankCounts& counts) {
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
    };
    for (int batch = 0; batch < options.batches; ++batch) {
        write_orientation(batch, "first", design.a1, design.b1, output[batch].first);
        write_orientation(batch, "second", design.a2, design.b2, output[batch].second);
    }
    if (options.marked_births) {
        if (!marked_births || !path_moments || !complement_audit) {
            throw std::logic_error("marked output streams were not opened");
        }
        auto write_marked_orientation = [&](int batch, const char* orientation,
                                            int a, int b, const Geometry& geometry,
                                            const MarkedOrientation& marked) {
            for (const auto& item : marked.sparse) {
                const MarkedKey& key = item.first;
                const Spin4 chi = spin4_mark(geometry.quotient, key.line);
                *marked_births << design.n << ',' << a << ',' << b << ',' << orientation
                    << ',' << batch << ',' << per_batch << ',' << key.k1 << ',' << key.k2
                    << ',' << static_cast<int>(key.k1 == key.k2) << ',' << key.site1
                    << ',' << key.site2 << ','
                    << static_cast<int>(key.line.x == 0 && key.line.y == 0) << ','
                    << key.line.x << ',' << key.line.y << ',' << key.index1 << ','
                    << key.index2 << ',' << key.physical.x << ',' << key.physical.y << ','
                    << std::setprecision(21) << chi.real << ',' << chi.imaginary << ','
                    << static_cast<int>(key.mark1.valid) << ',' << key.mark1.axis << ','
                    << key.mark1.diagonal << ',' << key.mark1.landed << ',' << key.mark1.h4
                    << ',' << static_cast<int>(key.mark2.valid) << ',' << key.mark2.axis
                    << ',' << key.mark2.diagonal << ',' << key.mark2.landed << ','
                    << key.mark2.h4 << ',' << item.second << '\n';
            }
            for (int k = 0; k < design.n; ++k) {
                const PathMoments& row = marked.path[k];
                *path_moments << design.n << ',' << a << ',' << b << ',' << orientation
                    << ',' << batch << ',' << row.samples << ',' << k << ',' << row.sum_q
                    << ',' << row.sum_q2 << ',' << row.sum_gate01 << ',' << row.sum_gate12
                    << ',' << row.sum_inactive_gate01 << ',' << row.sum_inactive_gate12
                    << ',' << row.sum_active_S << ',' << row.sum_active_D
                    << ',' << row.sum_inactive_S << ',' << row.sum_inactive_D
                    << ',' << row.sum_site_S << ',' << row.sum_site_D << ','
                    << std::setprecision(21) << row.sum_J_S_real << ','
                    << row.sum_J_S_imaginary << ',' << row.sum_J_D_real << ','
                    << row.sum_J_D_imaginary << ',' << row.sum_q_J_D_real << ','
                    << row.sum_q_J_D_imaginary << ',' << row.sum_O_ext << ','
                    << row.sum_O_ext2 << ',' << row.sum_O_near << ','
                    << row.sum_O_near2 << ',' << row.sum_O_ext_O_near << ','
                    << row.sum_O_ext_J_S_real << ','
                    << row.sum_O_ext_J_S_imaginary << ',' << row.sum_O_ext_J_D_real << ','
                    << row.sum_O_ext_J_D_imaginary << ',' << row.sum_O_near_J_S_real << ','
                    << row.sum_O_near_J_S_imaginary << ',' << row.sum_O_near_J_D_real << ','
                    << row.sum_O_near_J_D_imaginary << ','
                    << row.sum_J_D_conj_J_S_real << ','
                    << row.sum_J_D_conj_J_S_imaginary << ',' << row.sum_abs_J_S2 << ','
                    << row.sum_local_S << ',' << row.sum_local_D << ','
                    << row.sum_O_sep_axis << ',' << row.sum_O_sep_diagonal << ','
                    << row.sum_O_sep4 << ',' << row.sum_O_sep4_2 << ','
                    << row.sum_O_sep_axis_internal_h4 << ','
                    << row.sum_O_sep_diagonal_internal_h4 << ','
                    << row.sum_O_sep4_J_S_real << ','
                    << row.sum_O_sep4_J_S_imaginary << ','
                    << row.sum_O_sep4_J_D_real << ','
                    << row.sum_O_sep4_J_D_imaginary << '\n';
            }
            const ComplementAudit& audit = marked.complement;
            *complement_audit << design.n << ',' << a << ',' << b << ',' << orientation
                << ',' << batch << ',' << per_batch << ',' << audit.endpoint_failures << ','
                << audit.site_failures << ',' << audit.line_failures << ','
                << audit.local_mark_failures << ',' << audit.index_mismatches << ','
                << audit.separated_mark_failures << '\n';
        };
        for (int batch = 0; batch < options.batches; ++batch) {
            write_marked_orientation(batch, "first", design.a1, design.b1,
                                     first_geometry, marked_output[batch].first);
            write_marked_orientation(batch, "second", design.a2, design.b2,
                                     second_geometry, marked_output[batch].second);
        }
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
    const std::filesystem::path marked_path = options.output_prefix.string() + ".marked_births.csv";
    const std::filesystem::path path_moments_path = options.output_prefix.string() + ".path.csv";
    const std::filesystem::path audit_path = options.output_prefix.string() + ".complement_audit.csv";
    const std::filesystem::path geometry_pilot_path =
        options.output_prefix.string() + ".geometry_pilot.csv";
    std::ofstream histogram(histogram_path), moments(moments_path);
    if (!histogram || !moments) throw std::runtime_error("cannot open output files");
    std::ofstream marked_births, path_moments, complement_audit;
    if (options.marked_births) {
        marked_births.open(marked_path);
        path_moments.open(path_moments_path);
        complement_audit.open(audit_path);
        if (!marked_births || !path_moments || !complement_audit) {
            throw std::runtime_error("cannot open marked-birth output files");
        }
    }
    std::ofstream geometry_pilot;
    if (options.geometry_pilot_k0 > 0) {
        geometry_pilot.open(geometry_pilot_path);
        if (!geometry_pilot) throw std::runtime_error("cannot open geometry-pilot output");
        geometry_pilot
            << "n,a,b,orientation,batch,replica,k0,k1,k2,age_steps,ell_u,ell_v,"
               "essential_size,essential_carriers,occupied_frontier,vacant_frontier,"
               "boundary_cut_edges,boundary_multicontact_sites,boundary_contact_pairs,"
               "core_vertices,core_edges,articulation_vertices,bridges,"
               "boundary_axis_imbalance,boundary_corner_balance,frontier_components,"
               "largest_frontier_component,frontier_component_sumsq,"
               "H2,H2_theta,H2_figure8,H2_separate,H2_direction_positive,"
               "H2_direction_negative,H2_direction_mixed,next_site,next_exit,"
               "checkpoint_b1_safe_count,branch_common_safe,branch_suffix_site1,"
               "branch_suffix_site2,branch_clone1_survives,"
               "branch_clone2_survives,branch_both_survive\n";
    }
    histogram << "n,a,b,orientation,batch,samples,kind,k,count\n";
    moments << "n,a,b,orientation,batch,samples,sum_kminus,sum_kplus,sum_kminus2,"
               "sum_kplus2,sum_product,sum_gap,sum_gap2\n";
    if (options.marked_births) {
        marked_births
            << "n,a,b,orientation,batch,samples,k1,k2,direct_0_to_2,site01,site12,"
               "line_null,ell_u,ell_v,iota01,iota12,physical_x,physical_y,chi4_re,"
               "chi4_im,mark01_valid,mark01_axis,mark01_diagonal,mark01_landed,"
               "mark01_h4,mark12_valid,mark12_axis,mark12_diagonal,mark12_landed,"
               "mark12_h4,count\n";
        path_moments
            << "n,a,b,orientation,batch,samples,k,sum_q,sum_q2,sum_gate01,sum_gate12,"
               "sum_inactive_gate01,sum_inactive_gate12,sum_active_S,sum_active_D,"
               "sum_inactive_S,sum_inactive_D,"
               "sum_site_S,sum_site_D,sum_J_S_re,sum_J_S_im,sum_J_D_re,sum_J_D_im,"
               "sum_q_J_D_re,sum_q_J_D_im,sum_O_ext,sum_O_ext2,"
               "sum_O_near,sum_O_near2,sum_O_ext_O_near,"
               "sum_O_ext_J_S_re,sum_O_ext_J_S_im,sum_O_ext_J_D_re,"
               "sum_O_ext_J_D_im,sum_O_near_J_S_re,sum_O_near_J_S_im,"
               "sum_O_near_J_D_re,sum_O_near_J_D_im,"
               "sum_J_D_conj_J_S_re,sum_J_D_conj_J_S_im,"
               "sum_abs_J_S2,sum_local_S,sum_local_D,"
               "sum_O_sep_axis,sum_O_sep_diagonal,sum_O_sep4,sum_O_sep4_2,"
               "sum_O_sep_axis_internal_h4,sum_O_sep_diagonal_internal_h4,"
               "sum_O_sep4_J_S_re,sum_O_sep4_J_S_im,"
               "sum_O_sep4_J_D_re,sum_O_sep4_J_D_im\n";
        complement_audit
            << "n,a,b,orientation,batch,samples,endpoint_failures,site_failures,"
               "line_failures,local_mark_failures,index_mismatches,"
               "separated_mark_failures\n";
    }
    const auto started = std::chrono::steady_clock::now();
    if (options.geometry_pilot_k0 > 0) {
        run_geometry_pilot_design(design, options, geometry_pilot);
    } else {
        run_design(design, options, histogram, moments,
                   options.marked_births ? &marked_births : nullptr,
                   options.marked_births ? &path_moments : nullptr,
                   options.marked_births ? &complement_audit : nullptr);
    }
    histogram.close();
    moments.close();
    if (options.marked_births) {
        marked_births.close();
        path_moments.close();
        complement_audit.close();
    }
    if (geometry_pilot) geometry_pilot.close();
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
             << "  \"marked_birth_schema\": " << (options.marked_births ? "true" : "false") << ",\n"
             << "  \"geometry_pilot_k0\": " << options.geometry_pilot_k0 << ",\n"
             << "  \"geometry_pilot_semantics\": \"rank-one state after k0 occupied sites; H2 is the exact number of vacant sites whose one-step insertion gives ambient rank two\",\n"
             << "  \"branching_clones\": " << (options.branching_clones ? "true" : "false") << ",\n"
             << "  \"branching_semantics\": \"per rank-one checkpoint: permutation[k0] is one common uniform update; two independently tagged counter streams draw one uniform remaining site each from the identical successor clone; common absorption scores both clone survivals zero\",\n"
             << "  \"branching_suffix_stream_tags\": [\"0x429b2c101\", \"0x429b2c202\"],\n"
             << "  \"marked_birth_semantics\": \"strict 0->1 uses post-line; strict 1->2 uses pre-line; direct 0->2 has null line, D=0, S=2\",\n"
             << "  \"line_coordinates\": \"ell_u,ell_v are primitive period-basis winding coordinates\",\n"
             << "  \"chi4_frame\": \"physical lifted Euclidean direction primitive(P*ell), never raw period coordinates\",\n"
             << "  \"saturation_index\": \"gcd of raw winding coefficients on the primitive rational line before primitive reduction\",\n"
             << "  \"path_horvitz\": \"at pre-insertion k multiply the next-site gate by N-k; canonical Russo scorer later multiplies by N/(N-k) under Bin(N-1,k)\",\n"
             << "  \"full_source\": \"sum_site_S=(active_S+inactive_S)/2 and sum_site_D=(active_D-inactive_D)/2 on the paired primal/matching reverse insertion; raw sides are retained\",\n"
             << "  \"external_observer\": \"O_ext=C_black_NN-C_white_matching-q=V-E+F0; evaluated on the pre-insertion configuration and outside the q-only algebra\",\n"
             << "  \"external_contact_split\": \"O_near is the Chebyshev-radius-2 D4-symmetric sum of southwest-anchored local V-E+F0 densities around the next insertion site; O_far=O_ext-O_near\",\n"
             << "  \"external_products\": \"O_ext/O_near first and second moments, O_ext*O_near, both times J_S4/J_D4, plus J_D4*conj(J_S4) and |J_S4|^2 same-path Gram rows; q*J_D4 retained only as contact control\",\n"
             << "  \"separated_observer_radius\": " << options.far_radius << ",\n"
             << "  \"separated_observer\": \"at every pre-insertion root a counter-random common C4 rotation samples one axis anchor R*(1,0) and one diagonal anchor R*(1,1); their local arm landing values are retained as the typed complex pair axis+i*diagonal and O_sep4=axis-diagonal is the twice-normalized spatial H4 projection\",\n"
             << "  \"separated_products\": \"O_sep axis/diagonal means, O_sep4 first/second moments, internal local-H4 type controls, and O_sep4 times J_S4/J_D4 in the same path batch\",\n"
             << "  \"sparse_joint_histogram\": " << (options.marked_births ? "true" : "false") << ",\n"
             << "  \"per_batch_joint_moments\": true,\n"
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
             << "  \"marked_births_csv\": "
             << (options.marked_births ? "\"" + json_escape(marked_path.string()) + "\"" : "null") << ",\n"
             << "  \"path_csv\": "
             << (options.marked_births ? "\"" + json_escape(path_moments_path.string()) + "\"" : "null") << ",\n"
             << "  \"complement_audit_csv\": "
             << (options.marked_births ? "\"" + json_escape(audit_path.string()) + "\"" : "null") << ",\n"
             << "  \"geometry_pilot_csv\": "
             << (options.geometry_pilot_k0 > 0
                    ? "\"" + json_escape(geometry_pilot_path.string()) + "\"" : "null") << "\n"
             << "}\n";
    std::cout << "wrote " << histogram_path << "\nwrote " << moments_path
              << "\nwrote " << metadata_path << '\n';
    if (options.marked_births) {
        std::cout << "wrote " << marked_path << "\nwrote " << path_moments_path
                  << "\nwrote " << audit_path << '\n';
    }
    if (options.geometry_pilot_k0 > 0) {
        std::cout << "wrote " << geometry_pilot_path << '\n';
    }
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
