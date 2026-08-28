#pragma once

// C00 general 2x2 integer-period homology union-find, ported from
// scripts/torus_homology.py.  Closed cover displacements convert to generator
// windings by the exact identity w = adj(P) d / det(P) with integer arithmetic
// only.  This is the single DSU used by pell_matching_mc.cpp and
// gaussian_orientation_mc.cpp; do not add a second topology implementation.

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace matching {

struct PeriodMatrix {
    // Row-major [[a00, a01], [a10, a11]] with columns equal to the two
    // declared generators, matching scripts/torus_homology.py.
    std::int64_t a00 = 1;
    std::int64_t a01 = 0;
    std::int64_t a10 = 0;
    std::int64_t a11 = 1;

    static PeriodMatrix diagonal(std::int64_t px, std::int64_t py) {
        if (px <= 0 || py <= 0) {
            throw std::invalid_argument("diagonal periods must be positive");
        }
        return PeriodMatrix{px, 0, 0, py};
    }

    static PeriodMatrix gaussian(std::int64_t a, std::int64_t b) {
        if (a == 0 && b == 0) {
            throw std::invalid_argument("Gaussian generators must be nonzero");
        }
        return PeriodMatrix{a, -b, b, a};
    }

    std::int64_t det() const { return a00 * a11 - a01 * a10; }

    PeriodMatrix adjugate() const {
        return PeriodMatrix{a11, -a01, -a10, a00};
    }
};

struct Winding {
    std::int64_t x = 0;
    std::int64_t y = 0;
};

inline Winding apply_matrix(const PeriodMatrix& m, Winding v) {
    return Winding{m.a00 * v.x + m.a01 * v.y, m.a10 * v.x + m.a11 * v.y};
}

inline Winding primitive_winding(Winding v) {
    const std::int64_t ax = v.x < 0 ? -v.x : v.x;
    const std::int64_t ay = v.y < 0 ? -v.y : v.y;
    std::int64_t g = std::gcd(ax, ay);
    if (g == 0) {
        throw std::invalid_argument("the zero vector has no winding direction");
    }
    v.x /= g;
    v.y /= g;
    if (v.x < 0 || (v.x == 0 && v.y < 0)) {
        v.x = -v.x;
        v.y = -v.y;
    }
    return v;
}

struct WrappingChannels {
    int max_rank = 0;
    bool direction_0 = false;
    bool direction_1 = false;
    bool either = false;
    bool both = false;
    bool cross = false;

    int flag(int channel) const {
        switch (channel) {
            case 0: return static_cast<int>(cross);
            case 1: return static_cast<int>(both);
            case 2: return static_cast<int>(either);
            case 3: return static_cast<int>(direction_0);
            case 4: return static_cast<int>(direction_1);
            default: throw std::invalid_argument("channel index must be 0..4");
        }
    }
};

class HomologyUnionFind {
  public:
    HomologyUnionFind(int n, PeriodMatrix period)
        : n_(n),
          period_(period),
          det_(period.det()),
          adj_(period.adjugate()),
          parent_(static_cast<std::size_t>(n)),
          size_(static_cast<std::size_t>(n)),
          delta_x_(static_cast<std::size_t>(n)),
          delta_y_(static_cast<std::size_t>(n)),
          rank_(static_cast<std::size_t>(n)),
          basis0_(static_cast<std::size_t>(n)),
          basis1_(static_cast<std::size_t>(n)) {
        if (n < 0) throw std::invalid_argument("n must be nonnegative");
        if (det_ == 0) throw std::invalid_argument("period matrix must have nonzero determinant");
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(rank_.begin(), rank_.end(), 0);
        max_rank_ = 0;
        any_d0_ = false;
        any_d1_ = false;
        any_cross_ = false;
    }

    struct FindResult {
        int root;
        std::int64_t dx;
        std::int64_t dy;
    };

    FindResult find(int x) {
        int node = x;
        std::int64_t total_x = 0;
        std::int64_t total_y = 0;
        while (parent_[static_cast<std::size_t>(node)] != node) {
            total_x += delta_x_[static_cast<std::size_t>(node)];
            total_y += delta_y_[static_cast<std::size_t>(node)];
            node = parent_[static_cast<std::size_t>(node)];
        }
        const int root = node;

        node = x;
        std::int64_t remaining_x = total_x;
        std::int64_t remaining_y = total_y;
        while (parent_[static_cast<std::size_t>(node)] != node) {
            const int next = parent_[static_cast<std::size_t>(node)];
            const std::int64_t step_x = delta_x_[static_cast<std::size_t>(node)];
            const std::int64_t step_y = delta_y_[static_cast<std::size_t>(node)];
            parent_[static_cast<std::size_t>(node)] = root;
            delta_x_[static_cast<std::size_t>(node)] = remaining_x;
            delta_y_[static_cast<std::size_t>(node)] = remaining_y;
            remaining_x -= step_x;
            remaining_y -= step_y;
            node = next;
        }
        return {root, total_x, total_y};
    }

    Winding winding(std::int64_t dx, std::int64_t dy) const {
        const Winding numerator = apply_matrix(adj_, Winding{dx, dy});
        if (numerator.x % det_ != 0 || numerator.y % det_ != 0) {
            throw std::runtime_error(
                "closed-cycle displacement is not in the quotient period lattice");
        }
        return Winding{numerator.x / det_, numerator.y / det_};
    }

    void add_edge(int i, int j, std::int64_t edge_dx, std::int64_t edge_dy) {
        const FindResult fi = find(i);
        const FindResult fj = find(j);
        std::int64_t root_dx = fi.dx + edge_dx - fj.dx;
        std::int64_t root_dy = fi.dy + edge_dy - fj.dy;

        if (fi.root == fj.root) {
            extend_basis(fi.root, winding(root_dx, root_dy));
            return;
        }

        int root_i = fi.root;
        int root_j = fj.root;
        if (size_[static_cast<std::size_t>(root_i)] < size_[static_cast<std::size_t>(root_j)]) {
            std::swap(root_i, root_j);
            root_dx = -root_dx;
            root_dy = -root_dy;
        }

        parent_[static_cast<std::size_t>(root_j)] = root_i;
        delta_x_[static_cast<std::size_t>(root_j)] = root_dx;
        delta_y_[static_cast<std::size_t>(root_j)] = root_dy;
        size_[static_cast<std::size_t>(root_i)] += size_[static_cast<std::size_t>(root_j)];
        const int rank_j = rank_[static_cast<std::size_t>(root_j)];
        if (rank_j >= 1) {
            extend_basis(root_i, basis0_[static_cast<std::size_t>(root_j)]);
        }
        if (rank_j >= 2) {
            extend_basis(root_i, basis1_[static_cast<std::size_t>(root_j)]);
        }
        rank_[static_cast<std::size_t>(root_j)] = 0;
    }

    WrappingChannels channels() const {
        WrappingChannels out;
        out.max_rank = max_rank_;
        out.direction_0 = any_d0_;
        out.direction_1 = any_d1_;
        out.either = max_rank_ > 0;
        out.both = any_d0_ && any_d1_;
        out.cross = any_cross_;
        return out;
    }

    bool either() const { return max_rank_ > 0; }

    int n() const { return n_; }
    const PeriodMatrix& period_matrix() const { return period_; }

  private:
    void note_vector(Winding v) {
        if (v.x != 0) any_d0_ = true;
        if (v.y != 0) any_d1_ = true;
    }

    void extend_basis(int root, Winding vector) {
        if (vector.x == 0 && vector.y == 0) return;
        auto& rank = rank_[static_cast<std::size_t>(root)];
        if (rank == 2) return;
        vector = primitive_winding(vector);
        if (rank == 0) {
            basis0_[static_cast<std::size_t>(root)] = vector;
            rank = 1;
            note_vector(vector);
        } else {
            const Winding b0 = basis0_[static_cast<std::size_t>(root)];
            if (b0.x * vector.y == b0.y * vector.x) return;
            basis1_[static_cast<std::size_t>(root)] = vector;
            rank = 2;
            note_vector(vector);
            any_cross_ = true;
        }
        if (rank > max_rank_) max_rank_ = rank;
    }

    int n_;
    PeriodMatrix period_;
    std::int64_t det_;
    PeriodMatrix adj_;
    std::vector<int> parent_;
    std::vector<int> size_;
    std::vector<std::int64_t> delta_x_;
    std::vector<std::int64_t> delta_y_;
    std::vector<int> rank_;
    std::vector<Winding> basis0_;
    std::vector<Winding> basis1_;
    int max_rank_ = 0;
    bool any_d0_ = false;
    bool any_d1_ = false;
    bool any_cross_ = false;
};

inline const char* channel_name(int channel) {
    switch (channel) {
        case 0: return "cross";
        case 1: return "both";
        case 2: return "either";
        case 3: return "direction_0";
        case 4: return "direction_1";
        default: return "unknown";
    }
}

}  // namespace matching
