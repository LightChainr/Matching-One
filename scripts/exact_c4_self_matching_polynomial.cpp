// Exact brute-force polynomial for primitive odd C4 self-matching quotients.
//
// The quotient (a,b) has N=a^2+b^2 cyclic labels j=a*x+b*y mod N.
// The checkerboard triangulation contains square NN edges and, from even j,
// the two diagonals joining even-parity corners.  A displacement-aware DSU
// retains up to two independent winding generators for each component and
// accumulates all five repository wrapping channels by occupation count.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int kMaxSites = 26;
constexpr int kChannels = 5;
constexpr std::array<const char *, kChannels> kChannelNames{
    "cross", "both", "either", "direction_0", "direction_1"};

struct Edge {
  std::uint8_t i;
  std::uint8_t j;
  std::int8_t dx;
  std::int8_t dy;
};

struct Classification {
  std::array<bool, kChannels> value{};
};

struct HomologyDsu {
  std::array<std::uint8_t, kMaxSites> parent{};
  std::array<std::uint8_t, kMaxSites> size{};
  std::array<std::int16_t, kMaxSites> delta_x{};
  std::array<std::int16_t, kMaxSites> delta_y{};
  std::array<std::uint8_t, kMaxSites> rank{};
  std::array<std::int16_t, kMaxSites> g0_x{};
  std::array<std::int16_t, kMaxSites> g0_y{};
  std::array<std::int16_t, kMaxSites> g1_x{};
  std::array<std::int16_t, kMaxSites> g1_y{};
  int a = 0;
  int b = 0;
  int n = 0;

  void initialize(std::uint32_t mask, int input_a, int input_b, int input_n) {
    a = input_a;
    b = input_b;
    n = input_n;
    for (int i = 0; i < n; ++i) {
      if ((mask >> i) & 1U) {
        parent[i] = static_cast<std::uint8_t>(i);
        size[i] = 1;
        delta_x[i] = 0;
        delta_y[i] = 0;
        rank[i] = 0;
        g0_x[i] = g0_y[i] = g1_x[i] = g1_y[i] = 0;
      }
    }
  }

  int find(int x, int &dx, int &dy) {
    if (parent[x] == x) {
      dx = 0;
      dy = 0;
      return x;
    }
    const int old_parent = parent[x];
    int px = 0;
    int py = 0;
    const int root = find(old_parent, px, py);
    dx = static_cast<int>(delta_x[x]) + px;
    dy = static_cast<int>(delta_y[x]) + py;
    parent[x] = static_cast<std::uint8_t>(root);
    delta_x[x] = static_cast<std::int16_t>(dx);
    delta_y[x] = static_cast<std::int16_t>(dy);
    return root;
  }

  void add_generator(int root, int wx, int wy) {
    if ((wx == 0 && wy == 0) || rank[root] == 2) return;
    if (rank[root] == 0) {
      g0_x[root] = static_cast<std::int16_t>(wx);
      g0_y[root] = static_cast<std::int16_t>(wy);
      rank[root] = 1;
      return;
    }
    const int determinant = static_cast<int>(g0_x[root]) * wy -
                            static_cast<int>(g0_y[root]) * wx;
    if (determinant != 0) {
      g1_x[root] = static_cast<std::int16_t>(wx);
      g1_y[root] = static_cast<std::int16_t>(wy);
      rank[root] = 2;
    }
  }

  void add_edge(const Edge &edge) {
    int ix = 0, iy = 0, jx = 0, jy = 0;
    int ri = find(edge.i, ix, iy);
    int rj = find(edge.j, jx, jy);
    int root_dx = ix + edge.dx - jx;
    int root_dy = iy + edge.dy - jy;
    if (ri == rj) {
      // P=[[a,-b],[b,a]], so P^-1 d = adj(P)d/N.
      const int numerator_x = a * root_dx + b * root_dy;
      const int numerator_y = -b * root_dx + a * root_dy;
      if (numerator_x % n != 0 || numerator_y % n != 0) {
        throw std::runtime_error("cycle displacement is not in period lattice");
      }
      add_generator(ri, numerator_x / n, numerator_y / n);
      return;
    }
    if (size[ri] < size[rj]) {
      std::swap(ri, rj);
      root_dx = -root_dx;
      root_dy = -root_dy;
    }
    parent[rj] = static_cast<std::uint8_t>(ri);
    delta_x[rj] = static_cast<std::int16_t>(root_dx);
    delta_y[rj] = static_cast<std::int16_t>(root_dy);
    size[ri] = static_cast<std::uint8_t>(size[ri] + size[rj]);
    if (rank[rj] >= 1) add_generator(ri, g0_x[rj], g0_y[rj]);
    if (rank[rj] >= 2) add_generator(ri, g1_x[rj], g1_y[rj]);
  }
};

std::vector<Edge> make_edges(int a, int b) {
  const int n = a * a + b * b;
  std::vector<Edge> edges;
  edges.reserve(3 * n);
  const auto target = [n](int value) {
    value %= n;
    return value < 0 ? value + n : value;
  };
  for (int j = 0; j < n; ++j) {
    edges.push_back(Edge{static_cast<std::uint8_t>(j),
                         static_cast<std::uint8_t>(target(j + a)), 1, 0});
    edges.push_back(Edge{static_cast<std::uint8_t>(j),
                         static_cast<std::uint8_t>(target(j + b)), 0, 1});
    if ((j & 1) == 0) {
      edges.push_back(Edge{static_cast<std::uint8_t>(j),
                           static_cast<std::uint8_t>(target(j + a + b)), 1, 1});
      edges.push_back(Edge{static_cast<std::uint8_t>(j),
                           static_cast<std::uint8_t>(target(j + a - b)), 1, -1});
    }
  }
  if (static_cast<int>(edges.size()) != 3 * n) {
    throw std::runtime_error("checkerboard triangulation must have 3N edges");
  }
  return edges;
}

Classification classify(std::uint32_t mask, int a, int b,
                        const std::vector<Edge> &edges, HomologyDsu &dsu) {
  const int n = a * a + b * b;
  dsu.initialize(mask, a, b, n);
  for (const Edge &edge : edges) {
    if (((mask >> edge.i) & 1U) && ((mask >> edge.j) & 1U)) {
      dsu.add_edge(edge);
    }
  }

  bool direction_0 = false;
  bool direction_1 = false;
  bool cross = false;
  for (int vertex = 0; vertex < n; ++vertex) {
    if (!((mask >> vertex) & 1U) || dsu.parent[vertex] != vertex) continue;
    if (dsu.rank[vertex] >= 1) {
      direction_0 = direction_0 || dsu.g0_x[vertex] != 0;
      direction_1 = direction_1 || dsu.g0_y[vertex] != 0;
    }
    if (dsu.rank[vertex] >= 2) {
      direction_0 = direction_0 || dsu.g1_x[vertex] != 0;
      direction_1 = direction_1 || dsu.g1_y[vertex] != 0;
      cross = true;
    }
  }
  Classification result;
  result.value = {cross, direction_0 && direction_1,
                  direction_0 || direction_1, direction_0, direction_1};
  return result;
}

std::string integer_string(__int128 value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value > 0) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

std::uint64_t binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t value = 1;
  for (int i = 1; i <= k; ++i) value = value * (n - k + i) / i;
  return value;
}

template <class T>
std::string json_array(const std::vector<T> &values) {
  std::ostringstream out;
  out << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ',';
    out << values[i];
  }
  out << ']';
  return out.str();
}

std::string json_array_i128(const std::vector<__int128> &values) {
  std::ostringstream out;
  out << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ',';
    out << integer_string(values[i]);
  }
  out << ']';
  return out.str();
}

}  // namespace

int main(int argc, char **argv) {
  int a = 0;
  int b = 0;
  unsigned threads = std::max(1U, std::thread::hardware_concurrency());
  std::string json_path;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--a" && i + 1 < argc) a = std::stoi(argv[++i]);
    else if (arg == "--b" && i + 1 < argc) b = std::stoi(argv[++i]);
    else if (arg == "--threads" && i + 1 < argc) threads = std::stoul(argv[++i]);
    else if (arg == "--json" && i + 1 < argc) json_path = argv[++i];
    else if (arg == "--help") {
      std::cout << "usage: exact_c4_self_matching_polynomial --a 5 --b 1 "
                   "[--threads N] [--json PATH]\n";
      return 0;
    } else {
      throw std::invalid_argument("unknown or incomplete argument: " + arg);
    }
  }
  if (a <= 0 || b <= 0 || (a & 1) == 0 || (b & 1) == 0 ||
      std::gcd(a, b) != 1) {
    throw std::invalid_argument("require positive, odd, primitive --a/--b");
  }
  const int n = a * a + b * b;
  if (n > kMaxSites) throw std::invalid_argument("N must not exceed 26");
  if (threads < 1) throw std::invalid_argument("--threads must be positive");

  const std::uint64_t total = std::uint64_t{1} << n;
  threads = std::min<std::uint64_t>(threads, total);
  const auto edges = make_edges(a, b);
  using Counts = std::array<std::vector<std::uint64_t>, kChannels>;
  std::vector<Counts> counts_by_thread(threads);
  for (auto &counts : counts_by_thread) {
    for (auto &channel : counts) channel.assign(n + 1, 0);
  }

  const auto started = std::chrono::steady_clock::now();
  std::vector<std::thread> workers;
  for (unsigned tid = 0; tid < threads; ++tid) {
    workers.emplace_back([&, tid]() {
      HomologyDsu dsu;
      const std::uint64_t begin = total * tid / threads;
      const std::uint64_t end = total * (tid + 1) / threads;
      for (std::uint64_t raw_mask = begin; raw_mask < end; ++raw_mask) {
        const auto mask = static_cast<std::uint32_t>(raw_mask);
        const int k = __builtin_popcount(mask);
        const auto result = classify(mask, a, b, edges, dsu);
        for (int channel = 0; channel < kChannels; ++channel) {
          if (result.value[channel]) ++counts_by_thread[tid][channel][k];
        }
      }
    });
  }
  for (auto &worker : workers) worker.join();
  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();

  Counts wrapping;
  for (auto &channel : wrapping) channel.assign(n + 1, 0);
  for (unsigned tid = 0; tid < threads; ++tid) {
    for (int channel = 0; channel < kChannels; ++channel) {
      for (int k = 0; k <= n; ++k) {
        wrapping[channel][k] += counts_by_thread[tid][channel][k];
      }
    }
  }

  std::array<std::vector<__int128>, kChannels> matching;
  std::array<std::vector<__int128>, kChannels> power;
  for (int channel = 0; channel < kChannels; ++channel) {
    matching[channel].assign(n + 1, 0);
    power[channel].assign(n + 1, 0);
    for (int k = 0; k <= n; ++k) {
      matching[channel][k] =
          static_cast<__int128>(wrapping[channel][k]) -
          static_cast<__int128>(wrapping[channel][n - k]);
    }
    for (int k = 0; k <= n; ++k) {
      for (int degree = k; degree <= n; ++degree) {
        const __int128 term = matching[channel][k] *
                              binomial(n - k, degree - k);
        power[channel][degree] += ((degree - k) & 1) ? -term : term;
      }
    }
    while (power[channel].size() > 1 && power[channel].back() == 0) {
      power[channel].pop_back();
    }
  }

  std::ostringstream payload;
  payload.precision(17);
  payload << "{\n"
          << "  \"schema\": \"matching-one/c4-self-matching-exact-cpp/v1\",\n"
          << "  \"geometry\": {\"a\": " << a << ", \"b\": " << b
          << ", \"N\": " << n << ", \"edges\": " << edges.size()
          << ", \"period_matrix\": [[" << a << "," << -b << "],[" << b
          << "," << a << "]]},\n"
          << "  \"configurations\": " << total << ",\n"
          << "  \"threads\": " << threads << ",\n"
          << "  \"elapsed_seconds\": " << elapsed << ",\n"
          << "  \"channels\": {\n";
  for (int channel = 0; channel < kChannels; ++channel) {
    payload << "    \"" << kChannelNames[channel] << "\": {\n"
            << "      \"R_bernstein_integer_coefficients\": "
            << json_array(wrapping[channel]) << ",\n"
            << "      \"M_bernstein_integer_coefficients\": "
            << json_array_i128(matching[channel]) << ",\n"
            << "      \"M_power_coefficients_ascending\": "
            << json_array_i128(power[channel]) << "\n"
            << "    }" << (channel + 1 == kChannels ? "\n" : ",\n");
  }
  payload << "  }\n}\n";
  std::cout << payload.str();
  if (!json_path.empty()) {
    std::ofstream output(json_path);
    if (!output) throw std::runtime_error("cannot open output path: " + json_path);
    output << payload.str();
  }
  return 0;
}
