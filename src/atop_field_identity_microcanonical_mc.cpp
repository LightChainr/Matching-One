// Microcanonical A_top--rank-birth H4 field-identity stream.
//
// For every counter-keyed site permutation, B_k is the first k sites.  The
// outgoing insertion is uniform among the N-k inactive roots and the incoming
// insertion is uniform among the k active roots.  Therefore
//
//   J(B_k) = (N-k) j(next | B_k) + k j(last | B_k\last)
//
// is an unbiased estimator of the full site sum used by 078bd61.  Saving this
// at every k permits an exact binomial transform at a finite matching root
// chosen inside every delete-one replicate; no Bernoulli p is selected after
// seeing the field-identity response.

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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Int = std::int64_t;

struct Vector { Int x = 0; Int y = 0; };
struct Matrix { Int a = 0; Int b = 0; Int c = 0; Int d = 0; };
struct Edge { int i; int j; int dx; int dy; };

Int checked_int(__int128 value, const char* context) {
    if (value < std::numeric_limits<Int>::min() ||
        value > std::numeric_limits<Int>::max()) {
        throw std::overflow_error(std::string(context) + " overflows int64");
    }
    return static_cast<Int>(value);
}

Int safe_abs(Int value) {
    if (value == std::numeric_limits<Int>::min()) {
        throw std::overflow_error("absolute value overflows int64");
    }
    return std::llabs(value);
}

Int determinant(const Matrix& value) {
    return checked_int(static_cast<__int128>(value.a) * value.d -
                       static_cast<__int128>(value.b) * value.c,
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

struct Bezout { Int gcd; Int x; Int y; };

Bezout extended_gcd_nonnegative(Int a, Int b) {
    if (a < 0 || b < 0 || (a == 0 && b == 0)) {
        throw std::invalid_argument("invalid extended-gcd input");
    }
    Int old_r = a, r = b, old_s = 1, s = 0, old_t = 0, t = 1;
    while (r != 0) {
        const Int q = old_r / r;
        const Int nr = old_r - q * r, ns = old_s - q * s, nt = old_t - q * t;
        old_r = r; r = nr; old_s = s; s = ns; old_t = t; t = nt;
    }
    return {old_r, old_s, old_t};
}

Bezout extended_gcd(Int a, Int b) {
    const Bezout p = extended_gcd_nonnegative(safe_abs(a), safe_abs(b));
    return {p.gcd, a < 0 ? -p.x : p.x, b < 0 ? -p.y : p.y};
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

    explicit QuotientCoordinates(Matrix input)
        : periods(input), det(determinant(input)) {
        if (det == 0) throw std::invalid_argument("period matrix is singular");
        order = safe_abs(det);
        if (order > std::numeric_limits<int>::max()) {
            throw std::invalid_argument("quotient too large");
        }
        const Bezout lower = extended_gcd(periods.c, periods.d);
        h22 = lower.gcd;
        h11 = order / h22;
        h12 = positive_mod128(static_cast<__int128>(periods.a) * lower.x +
                              static_cast<__int128>(periods.b) * lower.y, h11);
        smith1 = std::gcd(std::gcd(safe_abs(periods.a), safe_abs(periods.b)),
                          std::gcd(safe_abs(periods.c), safe_abs(periods.d)));
        smith2 = order / smith1;
    }

    int label(Vector point) const {
        const Int qy = (point.y - positive_mod(point.y, h22)) / h22;
        const Int ry = point.y - qy * h22;
        const Int rx = positive_mod(point.x - qy * h12, h11);
        return static_cast<int>(rx + h11 * ry);
    }

    Vector representative(int label_value) const {
        return {label_value % h11, label_value / h11};
    }

    Vector winding(Int dx, Int dy) const {
        const __int128 n0 = static_cast<__int128>(periods.d) * dx -
                            static_cast<__int128>(periods.b) * dy;
        const __int128 n1 = -static_cast<__int128>(periods.c) * dx +
                            static_cast<__int128>(periods.a) * dy;
        if (n0 % det != 0 || n1 % det != 0) {
            throw std::logic_error("cycle is outside the period lattice");
        }
        return {checked_int(n0 / det, "winding"), checked_int(n1 / det, "winding")};
    }

    Vector period_vector(Vector winding_value) const {
        return {
            checked_int(static_cast<__int128>(periods.a) * winding_value.x +
                        static_cast<__int128>(periods.b) * winding_value.y,
                        "physical line"),
            checked_int(static_cast<__int128>(periods.c) * winding_value.x +
                        static_cast<__int128>(periods.d) * winding_value.y,
                        "physical line")
        };
    }
};

Vector primitive(Vector value) {
    const Int divisor = std::gcd(safe_abs(value.x), safe_abs(value.y));
    if (divisor == 0) return value;
    value.x /= divisor;
    value.y /= divisor;
    if (value.x < 0 || (value.x == 0 && value.y < 0)) {
        value.x = -value.x;
        value.y = -value.y;
    }
    return value;
}

struct Geometry {
    QuotientCoordinates quotient;
    int n;
    std::vector<Edge> edges;
    std::vector<std::vector<int>> incident;
    explicit Geometry(Matrix periods)
        : quotient(periods), n(static_cast<int>(quotient.order)) {}
};

Geometry make_geometry(Matrix periods) {
    Geometry geometry(periods);
    geometry.edges.reserve(2 * geometry.n);
    geometry.incident.resize(geometry.n);
    const std::array<Vector, 2> steps = {{{1, 0}, {0, 1}}};
    for (int vertex = 0; vertex < geometry.n; ++vertex) {
        const Vector source = geometry.quotient.representative(vertex);
        for (const Vector step : steps) {
            const Edge edge{vertex,
                            geometry.quotient.label({source.x + step.x, source.y + step.y}),
                            static_cast<int>(step.x), static_cast<int>(step.y)};
            const int index = static_cast<int>(geometry.edges.size());
            geometry.edges.push_back(edge);
            geometry.incident[edge.i].push_back(index);
            if (edge.j != edge.i) geometry.incident[edge.j].push_back(index);
        }
    }
    return geometry;
}

class AmbientHomologyUnionFind {
  public:
    explicit AmbientHomologyUnionFind(const QuotientCoordinates& quotient)
        : quotient_(quotient), parent_(quotient.order), size_(quotient.order),
          delta_x_(quotient.order), delta_y_(quotient.order), component_rank_(quotient.order),
          component_basis_(quotient.order) {
        reset();
    }

    void reset() {
        std::iota(parent_.begin(), parent_.end(), 0);
        std::fill(size_.begin(), size_.end(), 1);
        std::fill(delta_x_.begin(), delta_x_.end(), 0);
        std::fill(delta_y_.begin(), delta_y_.end(), 0);
        std::fill(component_rank_.begin(), component_rank_.end(), 0);
        ambient_rank_ = 0;
        ambient_basis_ = {};
    }

    struct FindResult { int root; Int dx; Int dy; };

    FindResult find(int vertex) {
        if (parent_[vertex] == vertex) return {vertex, 0, 0};
        const int old_parent = parent_[vertex];
        const FindResult above = find(old_parent);
        delta_x_[vertex] += above.dx;
        delta_y_[vertex] += above.dy;
        parent_[vertex] = above.root;
        return {above.root, delta_x_[vertex], delta_y_[vertex]};
    }

    static void extend_basis(std::uint8_t& rank, std::array<Vector, 2>& basis,
                             Vector value) {
        if ((value.x == 0 && value.y == 0) || rank == 2) return;
        value = primitive(value);
        if (rank == 0) {
            basis[0] = value;
            rank = 1;
            return;
        }
        const Vector first = basis[0];
        if (static_cast<__int128>(first.x) * value.y !=
            static_cast<__int128>(first.y) * value.x) {
            basis[1] = value;
            rank = 2;
        }
    }

    void add_cycle(int root, Vector winding) {
        const std::uint8_t old_rank = component_rank_[root];
        extend_basis(component_rank_[root], component_basis_[root], winding);
        if (component_rank_[root] > old_rank) {
            extend_basis(ambient_rank_, ambient_basis_, winding);
        }
    }

    void add_edge(const Edge& edge) {
        FindResult first = find(edge.i), second = find(edge.j);
        Int dx = first.dx + edge.dx - second.dx;
        Int dy = first.dy + edge.dy - second.dy;
        if (first.root == second.root) {
            add_cycle(first.root, quotient_.winding(dx, dy));
            return;
        }
        if (size_[first.root] < size_[second.root]) {
            std::swap(first, second);
            dx = -dx; dy = -dy;
        }
        parent_[second.root] = first.root;
        delta_x_[second.root] = dx;
        delta_y_[second.root] = dy;
        size_[first.root] += size_[second.root];
        for (std::uint8_t i = 0; i < component_rank_[second.root]; ++i) {
            extend_basis(component_rank_[first.root], component_basis_[first.root],
                         component_basis_[second.root][i]);
        }
        component_rank_[second.root] = 0;
    }

    int ambient_rank() const { return ambient_rank_; }
    Vector ambient_line() const {
        if (ambient_rank_ != 1) throw std::logic_error("ambient line requested outside rank one");
        return ambient_basis_[0];
    }

  private:
    const QuotientCoordinates& quotient_;
    std::vector<int> parent_, size_;
    std::vector<Int> delta_x_, delta_y_;
    std::vector<std::uint8_t> component_rank_;
    std::vector<std::array<Vector, 2>> component_basis_;
    std::uint8_t ambient_rank_ = 0;
    std::array<Vector, 2> ambient_basis_{};
};

struct Character { long double re = 0; long double im = 0; };

Character spin4(const QuotientCoordinates& quotient, Vector line) {
    const Vector physical = quotient.period_vector(line);
    const long double x = physical.x, y = physical.y;
    const long double r2 = x * x + y * y;
    const long double denominator = r2 * r2;
    return {(x*x*x*x - 6*x*x*y*y + y*y*y*y) / denominator,
            (4*x*y*(x*x-y*y)) / denominator};
}

struct Mark {
    int i01 = 0;
    int i12 = 0;
    int i02 = 0;
    Character s4;
    Character d4;
};

Mark insertion_mark(const QuotientCoordinates& quotient, int before_rank,
                    Vector before_line, int after_rank, Vector after_line) {
    if (!(0 <= before_rank && before_rank <= after_rank && after_rank <= 2)) {
        throw std::logic_error("ambient rank insertion is not monotone");
    }
    Mark mark;
    mark.i01 = before_rank == 0 && after_rank >= 1;
    mark.i12 = before_rank <= 1 && after_rank == 2;
    mark.i02 = before_rank == 0 && after_rank == 2;
    if (after_rank - before_rank != mark.i01 + mark.i12) {
        throw std::logic_error("rank jump did not split into I01+I12");
    }
    if (!mark.i02 && (mark.i01 || mark.i12)) {
        const Vector line = (before_rank == 0) ? after_line : before_line;
        const Character chi = spin4(quotient, line);
        const int s = mark.i01 + mark.i12;
        const int d = mark.i12 - mark.i01;
        mark.s4 = {s * chi.re, s * chi.im};
        mark.d4 = {d * chi.re, d * chi.im};
    }
    return mark;
}

struct LevelStats {
    long double samples = 0;
    long double sum_q = 0, sum_q2 = 0;
    long double sum_i01 = 0, sum_i12 = 0, sum_i02 = 0;
    long double sum_re_js4 = 0, sum_im_js4 = 0;
    long double sum_re_jd4 = 0, sum_im_jd4 = 0;
    long double sum_q_re_js4 = 0, sum_q_im_js4 = 0;
    long double sum_q_re_jd4 = 0, sum_q_im_jd4 = 0;
    long double sum_birth_mass = 0;

    void add(int q, int k, int n, const Mark* incoming, const Mark* outgoing) {
        const long double wi = static_cast<long double>(k);
        const long double wo = static_cast<long double>(n-k);
        auto scalar = [&](auto member) {
            return (incoming ? wi * (incoming->*member) : 0.0L) +
                   (outgoing ? wo * (outgoing->*member) : 0.0L);
        };
        auto character = [&](bool d, bool imag) {
            long double value = 0;
            if (incoming) {
                const Character& c = d ? incoming->d4 : incoming->s4;
                value += wi * (imag ? c.im : c.re);
            }
            if (outgoing) {
                const Character& c = d ? outgoing->d4 : outgoing->s4;
                value += wo * (imag ? c.im : c.re);
            }
            return value;
        };
        const long double i01 = scalar(&Mark::i01);
        const long double i12 = scalar(&Mark::i12);
        const long double i02 = scalar(&Mark::i02);
        const long double rs = character(false, false), is = character(false, true);
        const long double rd = character(true, false), id = character(true, true);
        ++samples;
        sum_q += q; sum_q2 += q*q;
        sum_i01 += i01; sum_i12 += i12; sum_i02 += i02;
        sum_re_js4 += rs; sum_im_js4 += is;
        sum_re_jd4 += rd; sum_im_jd4 += id;
        sum_q_re_js4 += q*rs; sum_q_im_js4 += q*is;
        sum_q_re_jd4 += q*rd; sum_q_im_jd4 += q*id;
        sum_birth_mass += i01+i12;
    }
};

class SweepEngine {
  public:
    explicit SweepEngine(const Geometry& geometry)
        : geometry_(geometry), active_(geometry.n), uf_(geometry.quotient) {}

    void accumulate(const std::vector<int>& permutation, std::vector<LevelStats>& levels) {
        std::fill(active_.begin(), active_.end(), 0);
        uf_.reset();
        Mark incoming;
        bool has_incoming = false;
        for (int k = 0; k < geometry_.n; ++k) {
            const int before_rank = uf_.ambient_rank();
            const Vector before_line = before_rank == 1 ? uf_.ambient_line() : Vector{};
            const int vertex = permutation[k];
            active_[vertex] = 1;
            for (const int edge_index : geometry_.incident[vertex]) {
                const Edge& edge = geometry_.edges[edge_index];
                if (active_[edge.i] && active_[edge.j]) uf_.add_edge(edge);
            }
            const int after_rank = uf_.ambient_rank();
            const Vector after_line = after_rank == 1 ? uf_.ambient_line() : Vector{};
            const Mark outgoing = insertion_mark(
                geometry_.quotient, before_rank, before_line, after_rank, after_line);
            levels[k].add(before_rank-1, k, geometry_.n,
                          has_incoming ? &incoming : nullptr, &outgoing);
            incoming = outgoing;
            has_incoming = true;
        }
        levels[geometry_.n].add(uf_.ambient_rank()-1, geometry_.n, geometry_.n,
                                &incoming, nullptr);
        if (uf_.ambient_rank() != 2) {
            throw std::logic_error("fully occupied graph has ambient rank below two");
        }
    }

  private:
    const Geometry& geometry_;
    std::vector<std::uint8_t> active_;
    AmbientHomologyUnionFind uf_;
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
    const std::uint64_t key = splitmix64(
        seed ^ splitmix64(replica + 0xd1b54a32d192ed03ULL));
    SplitMixStream generator(key);
    for (int stop = n-1; stop > 0; --stop) {
        const int other = static_cast<int>(generator.below(stop+1));
        std::swap(permutation[stop], permutation[other]);
    }
}

std::uint64_t permutation_digest(std::uint64_t replica,
                                 const std::vector<int>& permutation) {
    std::uint64_t value = splitmix64(replica ^ 0xa0761d6478bd642fULL);
    for (const int site : permutation) value = splitmix64(value ^ static_cast<std::uint64_t>(site));
    return value;
}

struct Options {
    std::uint64_t samples = 1000000;
    int batches = 100;
    std::uint64_t seed = 0;
    std::uint64_t replica_offset = 0;
    int threads = 0;
    Matrix matrix;
    bool matrix_set = false;
    std::string modulus;
    Int z_a = 0, z_b = 0;
    std::string git_commit = "unknown";
    std::string binary_sha256 = "unknown";
    std::string phase = "P275_PHASE1_MATCHING_ROOT";
    std::filesystem::path output_prefix;
    bool self_test = false;
};

template <typename T>
T parse_number(const std::string& text, const std::string& option) {
    if constexpr (std::is_unsigned_v<T>) {
        if (!text.empty() && text.front() == '-') throw std::invalid_argument("negative " + option);
    }
    std::istringstream input(text); T value{}; input >> value;
    if (!input || !input.eof()) throw std::invalid_argument("invalid " + option);
    return value;
}

[[noreturn]] void usage(const char* program, int status) {
    std::ostream& out = status == 0 ? std::cout : std::cerr;
    out << "Usage: " << program << " --matrix A B C D --modulus ID --z A B "
        << "--samples N --batches B --seed S --replica-offset K "
        << "--git-commit SHA --binary-sha256 SHA --output-prefix PATH [--threads T]\n"
        << "       " << program << " --self-test\n";
    std::exit(status);
}

Options parse_options(int argc, char** argv) {
    Options options;
    auto need = [&](int& i, const std::string& option) -> std::string {
        if (++i >= argc) throw std::invalid_argument(option + " needs a value");
        return argv[i];
    };
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--samples") options.samples = parse_number<std::uint64_t>(need(i,arg),arg);
        else if (arg == "--batches") options.batches = parse_number<int>(need(i,arg),arg);
        else if (arg == "--seed") options.seed = parse_number<std::uint64_t>(need(i,arg),arg);
        else if (arg == "--replica-offset") options.replica_offset = parse_number<std::uint64_t>(need(i,arg),arg);
        else if (arg == "--threads") options.threads = parse_number<int>(need(i,arg),arg);
        else if (arg == "--modulus") options.modulus = need(i,arg);
        else if (arg == "--git-commit") options.git_commit = need(i,arg);
        else if (arg == "--binary-sha256") options.binary_sha256 = need(i,arg);
        else if (arg == "--phase") options.phase = need(i,arg);
        else if (arg == "--output-prefix") options.output_prefix = need(i,arg);
        else if (arg == "--matrix") {
            options.matrix.a = parse_number<Int>(need(i,arg),arg);
            options.matrix.b = parse_number<Int>(need(i,arg),arg);
            options.matrix.c = parse_number<Int>(need(i,arg),arg);
            options.matrix.d = parse_number<Int>(need(i,arg),arg);
            options.matrix_set = true;
        } else if (arg == "--z") {
            options.z_a = parse_number<Int>(need(i,arg),arg);
            options.z_b = parse_number<Int>(need(i,arg),arg);
        } else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") usage(argv[0],0);
        else throw std::invalid_argument("unknown option " + arg);
    }
    if (options.self_test) return options;
    if (!options.matrix_set || options.modulus.empty() || options.output_prefix.empty() ||
        options.git_commit == "unknown" || options.binary_sha256 == "unknown") {
        throw std::invalid_argument("matrix/modulus/provenance/output are required");
    }
    if (options.samples == 0 || options.batches < 2 ||
        options.samples % static_cast<std::uint64_t>(options.batches) != 0) {
        throw std::invalid_argument("samples must be divisible by batches>=2");
    }
    if (options.replica_offset > std::numeric_limits<std::uint64_t>::max()-options.samples) {
        throw std::invalid_argument("counter range overflow");
    }
    const QuotientCoordinates quotient(options.matrix);
    if (quotient.smith1 != 1 || quotient.smith2 != quotient.order) {
        throw std::invalid_argument("P275 Phase 1 requires Smith invariants (1,N)");
    }
    return options;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        if (ch == '\\' || ch == '"') out << '\\' << ch;
        else if (ch == '\n') out << "\\n";
        else out << ch;
    }
    return out.str();
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
    std::ostringstream out; out << std::put_time(&time, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

void self_test() {
    // This exact all-permutation transform reproduces the committed 078bd61
    // Gaussian-(2,1) p=1/2 connected coupling, not merely endpoint ranks.
    const Geometry geometry = make_geometry({2,-1,1,2});
    std::vector<std::vector<LevelStats>> paths;
    std::vector<int> permutation(geometry.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    do {
        std::vector<LevelStats> levels(geometry.n+1);
        SweepEngine engine(geometry);
        engine.accumulate(permutation, levels);
        paths.push_back(std::move(levels));
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    if (paths.size() != 120) throw std::runtime_error("N=5 permutation count failed");
    long double mean_q = 0, mean_j_re = 0, mean_j_im = 0;
    long double mean_qj_re = 0, mean_qj_im = 0, mean_birth = 0;
    for (int k = 0; k <= geometry.n; ++k) {
        long double q = 0, birth = 0, d_re = 0, d_im = 0;
        long double qd_re = 0, qd_im = 0;
        for (const auto& path : paths) {
            q += path[k].sum_q;
            birth += path[k].sum_birth_mass;
            d_re += path[k].sum_re_jd4;
            d_im += path[k].sum_im_jd4;
            qd_re += path[k].sum_q_re_jd4;
            qd_im += path[k].sum_q_im_jd4;
        }
        q /= paths.size(); birth /= paths.size(); d_re /= paths.size();
        d_im /= paths.size(); qd_re /= paths.size(); qd_im /= paths.size();
        if ((k == 0 && q != -1) || (k == geometry.n && q != 1) ||
            !std::isfinite(static_cast<double>(birth+d_re+d_im))) {
            throw std::runtime_error("tiny microcanonical moment regression failed");
        }
        // Bin(5,1/2) weights.
        const std::array<int,6> choose{{1,5,10,10,5,1}};
        const long double weight = static_cast<long double>(choose[k])/32.0L;
        mean_q += weight*q;
        mean_j_re += weight*d_re; mean_j_im += weight*d_im;
        mean_qj_re += weight*qd_re; mean_qj_im += weight*qd_im;
        mean_birth += weight*birth;
    }
    const long double cov_re = mean_qj_re-mean_q*mean_j_re;
    const long double cov_im = mean_qj_im-mean_q*mean_j_im;
    if (std::fabs(cov_re + 49.0L/128.0L) > 1e-12L ||
        std::fabs(cov_im - 21.0L/16.0L) > 1e-12L ||
        std::fabs(mean_birth - 25.0L/8.0L) > 1e-12L) {
        std::ostringstream error;
        error << std::setprecision(20)
              << "078bd61 exact connected-coupling regression failed: cov=("
              << cov_re << ',' << cov_im << ") B=" << mean_birth
              << " means q/J=(" << mean_q << ';' << mean_j_re << ',' << mean_j_im << ')';
        throw std::runtime_error(error.str());
    }
    // The same period-basis line must acquire a different physical H4 phase
    // at the three moduli; using winding coordinates directly would fail.
    const Character ci = spin4(QuotientCoordinates({7,-1,1,7}), {1,0});
    const Character c2 = spin4(QuotientCoordinates({4,-6,3,8}), {1,0});
    const Character c25 = spin4(QuotientCoordinates({4,-5,2,10}), {1,0});
    const auto distance = [](Character x, Character y) {
        return std::fabs(x.re-y.re)+std::fabs(x.im-y.im);
    };
    if (distance(ci,c2) < 0.1L || distance(ci,c25) < 0.1L ||
        distance(c2,c25) < 0.1L) {
        throw std::runtime_error("modulus-sensitive physical-line H4 regression failed");
    }
    std::vector<int> expected{0,1,2,3,4};
    counter_permutation(5, 17, 0, permutation);
    if (permutation == expected || permutation_digest(0, permutation) == 0) {
        throw std::runtime_error("counter field/digest regression failed");
    }
    std::cout << "self-test passed: N=5 exact 078bd61 covariance, "
                 "incoming+outgoing root estimator, modulus-sensitive physical H4\n";
}

int run(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    if (options.self_test) { self_test(); return 0; }
    const Geometry geometry = make_geometry(options.matrix);
    const std::uint64_t per_batch = options.samples/options.batches;
    std::vector<std::vector<LevelStats>> output(
        options.batches, std::vector<LevelStats>(geometry.n+1));
    std::vector<std::uint64_t> digests(options.batches, 0);
    const auto started = std::chrono::steady_clock::now();
#ifdef _OPENMP
    if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
#pragma omp parallel for schedule(static)
    for (int batch = 0; batch < options.batches; ++batch) {
        SweepEngine engine(geometry);
        std::vector<int> permutation;
        std::uint64_t digest = 0;
        const std::uint64_t begin = options.replica_offset +
                                    static_cast<std::uint64_t>(batch)*per_batch;
        for (std::uint64_t replica = begin; replica < begin+per_batch; ++replica) {
            counter_permutation(geometry.n, options.seed, replica, permutation);
            engine.accumulate(permutation, output[batch]);
            digest ^= permutation_digest(replica, permutation);
        }
        digests[batch] = digest;
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();

    const auto parent = options.output_prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path csv_path = options.output_prefix.string()+".batches.csv";
    const std::filesystem::path metadata_path = options.output_prefix.string()+".metadata.json";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open batch CSV");
    csv << "n,modulus,batch,k,counter_first,counter_last_exclusive,samples,"
           "sum_q,sum_q2,sum_I01,sum_I12,sum_I02,sum_Re_J_S4,sum_Im_J_S4,"
           "sum_Re_J_D4,sum_Im_J_D4,sum_q_Re_J_S4,sum_q_Im_J_S4,"
           "sum_q_Re_J_D4,sum_q_Im_J_D4,sum_birth_mass,priority_field_digest\n";
    csv << std::setprecision(21);
    for (int batch = 0; batch < options.batches; ++batch) {
        const std::uint64_t first = options.replica_offset +
            static_cast<std::uint64_t>(batch)*per_batch;
        for (int k = 0; k <= geometry.n; ++k) {
            const LevelStats& s = output[batch][k];
            csv << geometry.n << ',' << options.modulus << ',' << batch << ',' << k << ','
                << first << ',' << first+per_batch << ',' << s.samples << ','
                << s.sum_q << ',' << s.sum_q2 << ',' << s.sum_i01 << ',' << s.sum_i12
                << ',' << s.sum_i02 << ',' << s.sum_re_js4 << ',' << s.sum_im_js4
                << ',' << s.sum_re_jd4 << ',' << s.sum_im_jd4 << ','
                << s.sum_q_re_js4 << ',' << s.sum_q_im_js4 << ','
                << s.sum_q_re_jd4 << ',' << s.sum_q_im_jd4 << ','
                << s.sum_birth_mass << ',' << std::hex << digests[batch] << std::dec << '\n';
        }
    }
    csv.close();

    const QuotientCoordinates quotient(options.matrix);
    std::ofstream metadata(metadata_path);
    if (!metadata) throw std::runtime_error("cannot open metadata");
    metadata << "{\n"
             << "  \"schema\": \"matching-one/p275-atop-field-identity-microcanonical/v1\",\n"
             << "  \"phase\": \"" << json_escape(options.phase) << "\",\n"
             << "  \"generated_utc\": \"" << utc_now() << "\",\n"
             << "  \"git_commit\": \"" << json_escape(options.git_commit) << "\",\n"
             << "  \"binary_sha256\": \"" << json_escape(options.binary_sha256) << "\",\n"
             << "  \"compiler\": \"" << json_escape(__VERSION__) << "\",\n"
             << "  \"openmp\": "
#ifdef _OPENMP
             << "true,\n"
#else
             << "false,\n"
#endif
             << "  \"threads_requested\": " << options.threads << ",\n"
             << "  \"N\": " << geometry.n << ",\n"
             << "  \"modulus\": \"" << json_escape(options.modulus) << "\",\n"
             << "  \"z\": [" << options.z_a << ',' << options.z_b << "],\n"
             << "  \"period_matrix\": [[" << options.matrix.a << ',' << options.matrix.b
             << "],[" << options.matrix.c << ',' << options.matrix.d << "]],\n"
             << "  \"smith_invariants\": [" << quotient.smith1 << ',' << quotient.smith2 << "],\n"
             << "  \"samples\": " << options.samples << ",\n"
             << "  \"batches\": " << options.batches << ",\n"
             << "  \"seed\": " << options.seed << ",\n"
             << "  \"replica_counter_first\": " << options.replica_offset << ",\n"
             << "  \"replica_counter_last_exclusive\": " << options.replica_offset+options.samples << ",\n"
             << "  \"rng\": \"counter-derived SplitMix64 plus unbiased Fisher-Yates\",\n"
             << "  \"field_contract\": \"same N/seed/counter gives byte-identical label priority field\",\n"
             << "  \"root_estimator\": \"k*last_active_mark+(N-k)*next_inactive_mark\",\n"
             << "  \"evaluation_p\": \"solve M_N(p)=E_p[A_top]=0 inside every delete-one replicate\",\n"
             << "  \"A_top\": \"ambient_H1_image_rank_minus_one\",\n"
             << "  \"J_D4_sign\": \"(I12-I01)*(x+iy)^4/abs(x+iy)^4\",\n"
             << "  \"birth_mass\": \"I01+I12 including direct 0-to-2 twice\",\n"
             << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
             << "  \"batch_csv\": \"" << json_escape(csv_path.string()) << "\"\n"
             << "}\n";
    std::cout << "completed N=" << geometry.n << " modulus=" << options.modulus
              << " samples=" << options.samples << " elapsed=" << elapsed << "s\n"
              << "wrote " << csv_path << "\nwrote " << metadata_path << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try { return run(argc, argv); }
    catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
