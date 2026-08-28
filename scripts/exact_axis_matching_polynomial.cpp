// Exact brute-force matching polynomial for an axis-aligned L x L square torus.
//
// This is deliberately a small frontier kernel, not a transfer-matrix framework.
// It evaluates all 2^(L^2) occupation masks, classifies wrapping with a
// displacement-aware union-find, and accumulates exact Bernstein coefficients.
// The implementation is portable C++17 and uses std::thread rather than OpenMP.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int kMaxSites = 25;

struct Edge {
  std::uint8_t i;
  std::uint8_t j;
  std::int8_t dx;
  std::int8_t dy;
};

struct WrapDsu {
  std::array<std::uint8_t, kMaxSites> parent{};
  std::array<std::uint8_t, kMaxSites> size{};
  std::array<std::int8_t, kMaxSites> delta_x{};
  std::array<std::int8_t, kMaxSites> delta_y{};
  std::array<bool, kMaxSites> wrap{};

  void initialize(std::uint32_t mask, int n) {
    for (int i = 0; i < n; ++i) {
      if ((mask >> i) & 1U) {
        parent[i] = static_cast<std::uint8_t>(i);
        size[i] = 1;
        delta_x[i] = 0;
        delta_y[i] = 0;
        wrap[i] = false;
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
    delta_x[x] = static_cast<std::int8_t>(dx);
    delta_y[x] = static_cast<std::int8_t>(dy);
    return root;
  }

  bool add_edge(const Edge &edge) {
    int ix = 0, iy = 0, jx = 0, jy = 0;
    int ri = find(edge.i, ix, iy);
    int rj = find(edge.j, jx, jy);
    const int root_dx = ix + edge.dx - jx;
    const int root_dy = iy + edge.dy - jy;
    if (ri == rj) {
      if (root_dx != 0 || root_dy != 0) wrap[ri] = true;
      return wrap[ri];
    }
    if (size[ri] >= size[rj]) {
      parent[rj] = static_cast<std::uint8_t>(ri);
      delta_x[rj] = static_cast<std::int8_t>(root_dx);
      delta_y[rj] = static_cast<std::int8_t>(root_dy);
      size[ri] = static_cast<std::uint8_t>(size[ri] + size[rj]);
      wrap[ri] = wrap[ri] || wrap[rj];
      return wrap[ri];
    }
    parent[ri] = static_cast<std::uint8_t>(rj);
    delta_x[ri] = static_cast<std::int8_t>(-root_dx);
    delta_y[ri] = static_cast<std::int8_t>(-root_dy);
    size[rj] = static_cast<std::uint8_t>(size[rj] + size[ri]);
    wrap[rj] = wrap[ri] || wrap[rj];
    return wrap[rj];
  }
};

std::vector<Edge> make_edges(int length, bool matching) {
  const std::array<std::array<int, 2>, 4> vectors{{
      {{1, 0}}, {{0, 1}}, {{1, 1}}, {{1, -1}},
  }};
  const int vector_count = matching ? 4 : 2;
  std::vector<Edge> edges;
  edges.reserve(length * length * vector_count);
  for (int y = 0; y < length; ++y) {
    for (int x = 0; x < length; ++x) {
      const int i = y * length + x;
      for (int v = 0; v < vector_count; ++v) {
        const int dx = vectors[v][0];
        const int dy = vectors[v][1];
        const int tx = (x + dx + length) % length;
        const int ty = (y + dy + length) % length;
        edges.push_back(Edge{static_cast<std::uint8_t>(i),
                             static_cast<std::uint8_t>(ty * length + tx),
                             static_cast<std::int8_t>(dx),
                             static_cast<std::int8_t>(dy)});
      }
    }
  }
  return edges;
}

bool wraps(std::uint32_t mask, int n, const std::vector<Edge> &edges,
           WrapDsu &dsu) {
  dsu.initialize(mask, n);
  for (const Edge &edge : edges) {
    if (((mask >> edge.i) & 1U) && ((mask >> edge.j) & 1U) &&
        dsu.add_edge(edge)) {
      return true;
    }
  }
  return false;
}

std::string integer_string(__int128 value) {
  if (value == 0) return "0";
  bool negative = value < 0;
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
  int length = 0;
  unsigned threads = std::max(1U, std::thread::hardware_concurrency());
  std::string json_path;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--L" && i + 1 < argc) length = std::stoi(argv[++i]);
    else if (arg == "--threads" && i + 1 < argc) threads = std::stoul(argv[++i]);
    else if (arg == "--json" && i + 1 < argc) json_path = argv[++i];
    else if (arg == "--help") {
      std::cout << "usage: exact_axis_matching_polynomial --L 5 [--threads N] [--json PATH]\n";
      return 0;
    } else {
      throw std::invalid_argument("unknown or incomplete argument: " + arg);
    }
  }
  if (length < 1 || length > 5) throw std::invalid_argument("--L must be in 1..5");
  if (threads < 1) throw std::invalid_argument("--threads must be positive");

  const int n = length * length;
  const std::uint64_t total = std::uint64_t{1} << n;
  threads = std::min<std::uint64_t>(threads, total);
  const auto primal_edges = make_edges(length, false);
  const auto matching_edges = make_edges(length, true);
  std::vector<std::vector<std::uint64_t>> primal_by_thread(
      threads, std::vector<std::uint64_t>(n + 1));
  std::vector<std::vector<std::uint64_t>> matching_by_thread(
      threads, std::vector<std::uint64_t>(n + 1));

  const auto started = std::chrono::steady_clock::now();
  std::vector<std::thread> workers;
  for (unsigned tid = 0; tid < threads; ++tid) {
    workers.emplace_back([&, tid]() {
      WrapDsu primal_dsu;
      WrapDsu matching_dsu;
      const std::uint64_t begin = total * tid / threads;
      const std::uint64_t end = total * (tid + 1) / threads;
      for (std::uint64_t raw_mask = begin; raw_mask < end; ++raw_mask) {
        const auto mask = static_cast<std::uint32_t>(raw_mask);
        const int k = __builtin_popcount(mask);
        if (wraps(mask, n, primal_edges, primal_dsu)) ++primal_by_thread[tid][k];
        if (wraps(mask, n, matching_edges, matching_dsu)) ++matching_by_thread[tid][k];
      }
    });
  }
  for (auto &worker : workers) worker.join();
  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();

  std::vector<std::uint64_t> primal(n + 1), matching(n + 1);
  for (unsigned tid = 0; tid < threads; ++tid) {
    for (int k = 0; k <= n; ++k) {
      primal[k] += primal_by_thread[tid][k];
      matching[k] += matching_by_thread[tid][k];
    }
  }
  std::vector<__int128> bernstein(n + 1);
  for (int k = 0; k <= n; ++k) {
    bernstein[k] = static_cast<__int128>(primal[k]) - matching[n - k];
  }
  std::vector<__int128> power(n + 1);
  for (int k = 0; k <= n; ++k) {
    if (bernstein[k] == 0) continue;
    for (int degree = k; degree <= n; ++degree) {
      const __int128 term = bernstein[k] * binomial(n - k, degree - k);
      power[degree] += ((degree - k) & 1) ? -term : term;
    }
  }
  while (power.size() > 1 && power.back() == 0) power.pop_back();

  std::ostringstream payload;
  payload.precision(17);
  payload << "{\n"
          << "  \"schema\": \"exact axis matching polynomial c++ v1\",\n"
          << "  \"geometry\": \"axis\",\n"
          << "  \"L\": " << length << ",\n"
          << "  \"N\": " << n << ",\n"
          << "  \"configurations\": " << total << ",\n"
          << "  \"threads\": " << threads << ",\n"
          << "  \"elapsed_seconds\": " << elapsed << ",\n"
          << "  \"primal_wrap_counts_by_occupancy\": " << json_array(primal) << ",\n"
          << "  \"matching_wrap_counts_by_occupancy\": " << json_array(matching) << ",\n"
          << "  \"bernstein_counts\": " << json_array_i128(bernstein) << ",\n"
          << "  \"power_coefficients_ascending\": " << json_array_i128(power) << "\n"
          << "}\n";
  std::cout << payload.str();
  if (!json_path.empty()) {
    std::ofstream output(json_path);
    if (!output) throw std::runtime_error("cannot open output path: " + json_path);
    output << payload.str();
  }
  return 0;
}
