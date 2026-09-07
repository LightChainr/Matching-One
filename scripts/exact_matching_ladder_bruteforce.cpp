// Exact brute-force Bernstein ladder for the finite Mertens-Ziff matching
// function on axis and diamond tori (issue #639).
//
// Reference definition (scripts/exact_matching_polynomial.py, geometry from
// scripts/matched_torus_reference.py):
//
//   D(C) = 1{black NN wraps} - 1{white NN+NNN wraps},  a_k = sum_{|C|=k} D(C)
//
// Two independent kernels:
//
//   K1 "rebuild": for every one of the 2^N site configurations, build a
//       displacement-aware union-find from scratch for the black NN graph and
//       for the white NN+NNN matching graph, and accumulate a_k. Direct port
//       of the Python reference (same directed-edge lists, same union-by-size
//       with path compression).
//
//   K2 "dfs": depth-first enumeration of the same 2^N configurations with an
//       incremental, fully journaled union-find (union by size, NO path
//       compression) maintained for BOTH colors simultaneously. Each
//       undirected edge is activated exactly once, when its higher-index
//       endpoint receives its color, and undone on backtracking. This is not
//       a transfer matrix and uses no symmetry assumptions: it visits the
//       identical 2^N configuration space.
//
// Both kernels must reproduce the committed L=2/3/4 integers bit-for-bit
// before any new rung is reported (GOVERNANCE.md section 2A).
//
// Build:  g++ -O2 -std=c++17 -pthread -o exact_ladder exact_matching_ladder_bruteforce.cpp
// Usage:  ./exact_ladder --geometry axis --L 5 --kernel k2 --threads 14
//         [--split B] [--json out.json] [--dump-edges]

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

struct Edge {
  int i;
  int j;
  int dx;
  int dy;
};

struct Geometry {
  std::string name;
  int L = 0;
  int n = 0;
  std::string physical_period;
  std::vector<Edge> primal;    // black NN, both directions, exactly as Python
  std::vector<Edge> matching;  // white NN+NNN, both directions, exactly as Python
};

int mod(int a, int p) {
  int r = a % p;
  return r < 0 ? r + p : r;
}

// Exact port of matched_torus_reference._make_edges: for every site (in
// enumeration order) and every vector, one directed edge i -> j with raw
// displacement (dx, dy), target taken modulo the period.
std::vector<Edge> make_edges(const std::vector<std::pair<int, int>>& coords,
                             const std::map<std::pair<int, int>, int>& ids,
                             int period,
                             const std::vector<std::pair<int, int>>& vectors) {
  std::vector<Edge> edges;
  for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
    const int x = coords[i].first;
    const int y = coords[i].second;
    for (const auto& v : vectors) {
      const int tx = mod(x + v.first, period);
      const int ty = mod(y + v.second, period);
      edges.push_back(Edge{i, ids.at({tx, ty}), v.first, v.second});
    }
  }
  return edges;
}

Geometry axis_geometry(int L) {
  Geometry g;
  g.name = "axis";
  g.L = L;
  g.physical_period = std::to_string(L);
  std::vector<std::pair<int, int>> coords;
  coords.reserve(L * L);
  for (int y = 0; y < L; ++y)
    for (int x = 0; x < L; ++x) coords.emplace_back(x, y);
  std::map<std::pair<int, int>, int> ids;
  for (int i = 0; i < static_cast<int>(coords.size()); ++i) ids[coords[i]] = i;
  g.n = static_cast<int>(coords.size());
  g.primal = make_edges(coords, ids, L, {{1, 0}, {0, 1}});
  g.matching = make_edges(coords, ids, L, {{1, 0}, {0, 1}, {1, 1}, {1, -1}});
  return g;
}

Geometry diamond_geometry(int L) {
  Geometry g;
  g.name = "diamond";
  g.L = L;
  g.physical_period = "sqrt(2)*" + std::to_string(L);
  const int period = 2 * L;
  std::vector<std::pair<int, int>> coords;
  for (int u = 0; u < period; ++u)
    for (int v = 0; v < period; ++v)
      if (mod(u - v, 2) == 0) coords.emplace_back(u, v);
  std::map<std::pair<int, int>, int> ids;
  for (int i = 0; i < static_cast<int>(coords.size()); ++i) ids[coords[i]] = i;
  g.n = static_cast<int>(coords.size());
  g.primal = make_edges(coords, ids, period, {{1, -1}, {1, 1}});
  g.matching = make_edges(coords, ids, period, {{1, -1}, {1, 1}, {2, 0}, {0, 2}});
  return g;
}

// ---------------------------------------------------------------------------
// K1: compressing union-find, exact port of matched_torus_reference.
// ---------------------------------------------------------------------------

class WrapDsu {
 public:
  explicit WrapDsu(int n) : n_(n) {
    parent_.resize(n);
    size_.assign(n, 1);
    dx_.assign(n, 0);
    dy_.assign(n, 0);
    wrap_.assign(n, 0);
  }

  void reset() {
    for (int i = 0; i < n_; ++i) {
      parent_[i] = i;
      size_[i] = 1;
      dx_[i] = 0;
      dy_[i] = 0;
      wrap_[i] = 0;
    }
  }

  // Returns root; dx/dy receive position(x) - position(root).  Mirrors the
  // Python reference exactly: after path compression the *updated* deltas
  // (dx + px, dy + py) are returned, not the pre-compression values.
  int find(int x, int& dx, int& dy) {
    if (parent_[x] == x) {
      dx = 0;
      dy = 0;
      return x;
    }
    const int p = parent_[x];
    int px = 0;
    int py = 0;
    const int root = find(p, px, py);
    const int old_dx = dx_[x];
    const int old_dy = dy_[x];
    parent_[x] = root;
    dx_[x] = old_dx + px;
    dy_[x] = old_dy + py;
    dx = dx_[x];
    dy = dy_[x];
    return root;
  }

  // Requires position(j) = position(i) + (edge_dx, edge_dy).
  void add_edge(int i, int j, int edge_dx, int edge_dy) {
    int ix = 0, iy = 0, jx = 0, jy = 0;
    const int ri = find(i, ix, iy);
    const int rj = find(j, jx, jy);
    const int root_dx = ix + edge_dx - jx;
    const int root_dy = iy + edge_dy - jy;
    if (ri == rj) {
      if (root_dx != 0 || root_dy != 0) wrap_[ri] = 1;
      return;
    }
    if (size_[ri] >= size_[rj]) {
      parent_[rj] = ri;
      dx_[rj] = root_dx;
      dy_[rj] = root_dy;
      size_[ri] += size_[rj];
      wrap_[ri] = wrap_[ri] | wrap_[rj];
    } else {
      parent_[ri] = rj;
      dx_[ri] = -root_dx;
      dy_[ri] = -root_dy;
      size_[rj] += size_[ri];
      wrap_[rj] = wrap_[ri] | wrap_[rj];
    }
  }

  // Port of cluster_stats' final pass: find() over all active vertices and
  // report whether any component carries non-zero displacement.
  bool any_wrap(const std::vector<char>& active) {
    for (int i = 0; i < n_; ++i) {
      if (!active[i]) continue;
      int dx = 0, dy = 0;
      const int root = find(i, dx, dy);
      if (wrap_[root]) return true;
    }
    return false;
  }

 private:
  int n_;
  std::vector<int> parent_, size_, dx_, dy_;
  std::vector<char> wrap_;
};

void kernel1_threaded(const Geometry& g, int threads,
                      std::vector<long long>& counts, long long& visited,
                      double& seconds) {
  const int n = g.n;
  counts.assign(n + 1, 0);
  visited = 0;
  const unsigned long long total = 1ULL << n;
  const unsigned long long chunk = (total + threads - 1) / threads;
  std::vector<std::vector<long long>> locals(threads);
  std::vector<long long> dones(threads, 0);
  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t]() {
      const unsigned long long begin = t * chunk;
      const unsigned long long end = std::min(total, begin + chunk);
      locals[t].assign(n + 1, 0);
      if (begin >= end) return;
      WrapDsu dsu(n);
      std::vector<char> black(n), white(n);
      std::vector<long long>& local = locals[t];
      for (unsigned long long mask = begin; mask < end; ++mask) {
        for (int i = 0; i < n; ++i) {
          const bool b = (mask >> i) & 1ULL;
          black[i] = b ? 1 : 0;
          white[i] = b ? 0 : 1;
        }
        const int k = __builtin_popcountll(mask);
        dsu.reset();
        for (const Edge& e : g.primal)
          if (black[e.i] && black[e.j]) dsu.add_edge(e.i, e.j, e.dx, e.dy);
        const int bw = dsu.any_wrap(black) ? 1 : 0;
        dsu.reset();
        for (const Edge& e : g.matching)
          if (white[e.i] && white[e.j]) dsu.add_edge(e.i, e.j, e.dx, e.dy);
        const int ww = dsu.any_wrap(white) ? 1 : 0;
        if (getenv("LADDER_TRACE") && static_cast<int>(mask) < 4096)
          std::printf("mask=%llu k=%d bw=%d ww=%d d=%d\n",
                      (unsigned long long)mask, k, bw, ww, bw - ww);
        local[k] += bw - ww;
      }
      dones[t] = static_cast<long long>(end - begin);
    });
  }
  for (auto& th : pool) th.join();
  auto t1 = std::chrono::steady_clock::now();
  seconds = std::chrono::duration<double>(t1 - t0).count();
  for (int t = 0; t < threads; ++t) {
    visited += dones[t];
    for (int k = 0; k <= n; ++k) counts[k] += locals[t][k];
  }
}

// ---------------------------------------------------------------------------
// K2: journaled union-find (union by size, no path compression) + DFS over
// all 2^N colorings. Both color DSUs live simultaneously; each undirected
// edge is added exactly once, when its higher-index endpoint receives its
// color, and undone on backtracking.
// ---------------------------------------------------------------------------

class JDsu {
 public:
  explicit JDsu(int n) : n_(n) {
    parent_.resize(n);
    size_.resize(n);
    dx_.resize(n);
    dy_.resize(n);
    wrap_.resize(n);
  }

  void reset() {
    for (int i = 0; i < n_; ++i) {
      parent_[i] = i;
      size_[i] = 1;
      dx_[i] = 0;
      dy_[i] = 0;
      wrap_[i] = 0;
    }
    journal_.clear();
    events_ = 0;
  }

  // Iterative find without path compression: sums displacements up to root.
  int find(int x, int& dx, int& dy) const {
    dx = 0;
    dy = 0;
    while (parent_[x] != x) {
      dx += dx_[x];
      dy += dy_[x];
      x = parent_[x];
    }
    return x;
  }

  void add_edge(int i, int j, int edge_dx, int edge_dy) {
    int ix = 0, iy = 0, jx = 0, jy = 0;
    const int ri = find(i, ix, iy);
    const int rj = find(j, jx, jy);
    const int root_dx = ix + edge_dx - jx;
    const int root_dy = iy + edge_dy - jy;
    if (ri == rj) {
      if (root_dx != 0 || root_dy != 0) {
        if (!wrap_[ri]) {
          journal_.push_back(Mut{0, ri, 0, 0, 0, 0});
          wrap_[ri] = 1;
          ++events_;
        }
      }
      return;
    }
    if (size_[ri] >= size_[rj]) {
      journal_.push_back(Mut{1, rj, ri, size_[ri], wrap_[ri], wrap_[rj]});
      parent_[rj] = ri;
      dx_[rj] = root_dx;
      dy_[rj] = root_dy;
      size_[ri] += size_[rj];
      wrap_[ri] = wrap_[ri] | wrap_[rj];
    } else {
      journal_.push_back(Mut{1, ri, rj, size_[rj], wrap_[rj], wrap_[ri]});
      parent_[ri] = rj;
      dx_[ri] = -root_dx;
      dy_[ri] = -root_dy;
      size_[rj] += size_[ri];
      wrap_[rj] = wrap_[ri] | wrap_[rj];
    }
  }

  void undo_to(size_t mark) {
    while (journal_.size() > mark) {
      const Mut& m = journal_.back();
      if (m.kind == 0) {
        wrap_[m.root] = 0;
        --events_;
      } else {
        parent_[m.root] = m.root;
        dx_[m.root] = 0;
        dy_[m.root] = 0;
        size_[m.newroot] = m.old_size_newroot;
        wrap_[m.newroot] = m.old_wrap_newroot;
        wrap_[m.root] = m.old_wrap_root;
      }
      journal_.pop_back();
    }
  }

  size_t mark() const { return journal_.size(); }
  bool any_wrap() const { return events_ > 0; }

 private:
  struct Mut {
    int kind;  // 0 = wrap-set on existing root; 1 = merge
    int root;  // kind 0: root whose wrap was set; kind 1: merged (child) root
    int newroot;  // kind 1: surviving root
    int old_size_newroot;
    uint8_t old_wrap_newroot;
    uint8_t old_wrap_root;
  };
  int n_;
  std::vector<int> parent_, size_, dx_, dy_;
  std::vector<uint8_t> wrap_;
  std::vector<Mut> journal_;
  long long events_ = 0;
};

// Undirected adjacency: for vertex v (the higher index of the pair), entries
// (u, dx, dy) with position(u) = position(v) + (dx, dy).
struct AdjEntry {
  int u;
  int dx;
  int dy;
};

std::vector<std::vector<AdjEntry>> build_adjacency(const std::vector<Edge>& edges,
                                                   int n) {
  std::vector<std::vector<AdjEntry>> adj(n);
  // NO deduplication by unordered pair: two directed edges between the same
  // site pair can carry different raw displacements (e.g. on an L=2 torus,
  // +1 and -1 are distinct raw displacements whose difference is the wrap).
  // The reference processes every directed edge, so we keep every one, in
  // the same relative order.  Each entry is stored at the higher-index
  // endpoint v, in the add_edge(v, u, d) convention pos(u) = pos(v) + d.
  for (const Edge& e : edges) {
    const int lo = std::min(e.i, e.j);
    const int hi = std::max(e.i, e.j);
    if (lo == hi) continue;
    if (e.i == hi) {
      // pos(j=lo) = pos(i=hi) + d  =>  entry at hi: (lo, d).
      adj[hi].push_back(AdjEntry{lo, e.dx, e.dy});
    } else {
      // e.j == hi: pos(lo) = pos(hi) + d  =>  entry at hi: (lo, -d).
      adj[hi].push_back(AdjEntry{lo, -e.dx, -e.dy});
    }
  }
  return adj;
}

struct K2Ctx {
  JDsu* black;
  JDsu* white;
  const std::vector<std::vector<AdjEntry>>* adj_black;
  const std::vector<std::vector<AdjEntry>>* adj_white;
  uint8_t* color;  // per-vertex: 1 = black, 0 = white, 2 = unassigned
  int n;
  std::vector<long long> counts;
};

void dfs(K2Ctx& ctx, int idx, int black_cnt) {
  if (idx == ctx.n) {
    if (const char* tr = getenv("LADDER_TRACE2")) {
      if (atoi(tr)) {
        unsigned long long m = 0;
        for (int i = 0; i < ctx.n; ++i)
          if (ctx.color[i] == 1) m |= (1ULL << i);
        std::printf("k2mask=%llu k=%d bw=%d ww=%d\n", m, black_cnt,
                    ctx.black->any_wrap() ? 1 : 0, ctx.white->any_wrap() ? 1 : 0);
      }
    }
    ctx.counts[black_cnt] +=
        (ctx.black->any_wrap() ? 1 : 0) - (ctx.white->any_wrap() ? 1 : 0);
    return;
  }
  // Branch white: vertex idx is white.
  {
    ctx.color[idx] = 0;
    const size_t m = ctx.white->mark();
    for (const AdjEntry& e : (*ctx.adj_white)[idx]) {
      if (e.u < idx && ctx.color[e.u] == 0) ctx.white->add_edge(idx, e.u, e.dx, e.dy);
    }
    dfs(ctx, idx + 1, black_cnt);
    ctx.white->undo_to(m);
  }
  // Branch black: vertex idx is black.
  {
    ctx.color[idx] = 1;
    const size_t m = ctx.black->mark();
    for (const AdjEntry& e : (*ctx.adj_black)[idx]) {
      if (e.u < idx && ctx.color[e.u] == 1) ctx.black->add_edge(idx, e.u, e.dx, e.dy);
    }
    dfs(ctx, idx + 1, black_cnt + 1);
    ctx.black->undo_to(m);
  }
}

void kernel2_threaded(const Geometry& g, int threads, int split,
                      std::vector<long long>& counts, long long& visited,
                      double& seconds) {
  const int n = g.n;
  counts.assign(n + 1, 0);
  visited = 0;
  auto adj_black = std::make_shared<std::vector<std::vector<AdjEntry>>>(
      build_adjacency(g.primal, n));
  auto adj_white = std::make_shared<std::vector<std::vector<AdjEntry>>>(
      build_adjacency(g.matching, n));
  if (split <= 0 || split >= n) split = std::min(n, 12);
  const int B = split;
  const unsigned long long tasks = 1ULL << B;

  std::vector<std::vector<long long>> locals(threads);
  std::vector<long long> leaves(threads, 0);
  std::atomic<unsigned long long> next_task{0};
  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t]() {
      JDsu black(n), white(n);
      std::vector<uint8_t> color(n, 0);
      K2Ctx ctx{&black, &white, &*adj_black, &*adj_white, color.data(), n, {}};
      ctx.counts.assign(n + 1, 0);
      for (;;) {
        const unsigned long long task = next_task.fetch_add(1);
        if (task >= tasks) break;
        black.reset();
        white.reset();
        std::fill(color.begin(), color.end(), 2);  // 2 = unassigned
        int black_cnt = 0;
        for (int v = 0; v < B; ++v) {
          const int c = static_cast<int>((task >> v) & 1ULL);
          color[v] = static_cast<uint8_t>(c);
          if (c) ++black_cnt;
          // Add each undirected edge exactly once, when its higher-index
          // endpoint is colored: adjacency is stored at the higher index.
          // The endpoint v itself must carry the matching color.
          if (c == 1) {
            for (const AdjEntry& e : (*adj_black)[v])
              if (color[e.u] == 1) black.add_edge(v, e.u, e.dx, e.dy);
          } else {
            for (const AdjEntry& e : (*adj_white)[v])
              if (color[e.u] == 0) white.add_edge(v, e.u, e.dx, e.dy);
          }
        }
        dfs(ctx, B, black_cnt);
        leaves[t] += (1LL << (n - B));
      }
      locals[t] = std::move(ctx.counts);
    });
  }
  for (auto& th : pool) th.join();
  auto t1 = std::chrono::steady_clock::now();
  seconds = std::chrono::duration<double>(t1 - t0).count();
  for (int t = 0; t < threads; ++t) {
    visited += leaves[t];
    for (int k = 0; k <= n; ++k) counts[k] += locals[t][k];
  }
}

void dump_edges(const Geometry& g) {
  auto dump = [](const char* label, const std::vector<Edge>& edges) {
    std::printf("%s\n", label);
    for (const Edge& e : edges)
      std::printf("  %d %d %d %d\n", e.i, e.j, e.dx, e.dy);
  };
  dump("primal", g.primal);
  dump("matching", g.matching);
}

void dump_adjacency(const Geometry& g) {
  const int n = g.n;
  auto adj_black = build_adjacency(g.primal, n);
  auto adj_white = build_adjacency(g.matching, n);
  auto dump = [](const char* label, const std::vector<std::vector<AdjEntry>>& adj) {
    std::printf("%s\n", label);
    for (int v = 0; v < static_cast<int>(adj.size()); ++v)
      for (const AdjEntry& e : adj[v])
        std::printf("  v=%d u=%d d=(%d,%d)\n", v, e.u, e.dx, e.dy);
  };
  dump("adj_black", adj_black);
  dump("adj_white", adj_white);
}

}  // namespace

int main(int argc, char** argv) {
  std::string geometry_name;
  int L = 0;
  std::string kernel = "k2";
  int threads = static_cast<int>(std::thread::hardware_concurrency());
  int split = 0;
  std::string json_path;
  bool dump = false;

  for (int i = 1; i < argc; ++i) {
    auto need = [&]() -> const char* { return argv[++i]; };
    if (!std::strcmp(argv[i], "--geometry")) geometry_name = need();
    else if (!std::strcmp(argv[i], "--L")) L = std::atoi(need());
    else if (!std::strcmp(argv[i], "--kernel")) kernel = need();
    else if (!std::strcmp(argv[i], "--threads")) threads = std::atoi(need());
    else if (!std::strcmp(argv[i], "--split")) split = std::atoi(need());
    else if (!std::strcmp(argv[i], "--json")) json_path = need();
    else if (!std::strcmp(argv[i], "--dump-edges")) dump = true;
    else {
      std::fprintf(stderr, "unknown arg %s\n", argv[i]);
      return 2;
    }
  }
  if (geometry_name.empty() || L <= 0) {
    std::fprintf(stderr,
                 "usage: exact_ladder --geometry axis|diamond --L N "
                 "[--kernel k1|k2] [--threads T] [--split B] [--json f] "
                 "[--dump-edges]\n");
    return 2;
  }
  if (threads <= 0) threads = 1;

  Geometry g = geometry_name == "axis" ? axis_geometry(L) : diamond_geometry(L);
  if (dump) {
    dump_edges(g);
    dump_adjacency(g);
    return 0;
  }

  std::vector<long long> counts;
  long long visited = 0;
  double seconds = 0.0;

  if (kernel == "k1") {
    if (g.n > 40) {
      std::fprintf(stderr, "K1 refuses N=%d > 40 (64-bit mask limit; the "
                           "Python reference refuses N>26)\n", g.n);
      return 3;
    }
    kernel1_threaded(g, threads, counts, visited, seconds);
  } else {
    kernel2_threaded(g, threads, split, counts, visited, seconds);
  }

  std::printf("geometry: %s\n", g.name.c_str());
  std::printf("L: %d\n", g.L);
  std::printf("N: %d\n", g.n);
  std::printf("physical_period: %s\n", g.physical_period.c_str());
  std::printf("kernel: %s\n", kernel.c_str());
  std::printf("threads: %d\n", threads);
  std::printf("bernstein_counts:");
  for (long long v : counts) std::printf(" %lld", v);
  std::printf("\n");
  std::printf("configs_visited: %lld (expect %lld)\n", visited, 1LL << g.n);
  std::printf("wall_seconds: %.3f\n", seconds);

  if (!json_path.empty()) {
    FILE* f = std::fopen(json_path.c_str(), "w");
    if (f) {
      std::fprintf(f, "{\n");
      std::fprintf(f, "  \"geometry\": \"%s\",\n", g.name.c_str());
      std::fprintf(f, "  \"L\": %d,\n", g.L);
      std::fprintf(f, "  \"N\": %d,\n", g.n);
      std::fprintf(f, "  \"physical_period\": \"%s\",\n",
                   g.physical_period.c_str());
      std::fprintf(f, "  \"kernel\": \"%s\",\n", kernel.c_str());
      std::fprintf(f, "  \"threads\": %d,\n", threads);
      std::fprintf(f, "  \"bernstein_counts\": [");
      for (size_t k = 0; k < counts.size(); ++k) {
        std::fprintf(f, "%lld%s", counts[k],
                     k + 1 == counts.size() ? "" : ", ");
      }
      std::fprintf(f, "],\n");
      std::fprintf(f, "  \"configs_visited\": %lld,\n", visited);
      std::fprintf(f, "  \"wall_seconds\": %.3f\n", seconds);
      std::fprintf(f, "}\n");
      std::fclose(f);
    }
  }
  return 0;
}
