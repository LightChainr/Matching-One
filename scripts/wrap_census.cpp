// Wrapping-type census with #646 winding labels (issues #640, #651).
//
// Base: scripts/exact_matching_ladder_bruteforce.cpp from PR #649 (branch
// exact-bernstein-ladder-639) — the two independent exact kernels K1/K2 that
// reproduce the committed axis L=2..5, diamond L=2..4 Bernstein integers.
// This file adds per-configuration wrap-type instrumentation ONLY; the DSU,
// edge lists, enumeration order and D(C) semantics are untouched.
//
// Wrap labels (#646 semantics, scripts/probe635_sector_decomposition.py):
//   wrap_homology() maintains per-cluster wrap_x / wrap_y flags from the raw
//   displacement deltas.  Label of a color class:
//     none       no cluster wraps in either direction
//     x / y      some cluster wraps exactly one axis (up to sign)
//     both-same  a single cluster wraps both axes
//     both-two   wraps both axes but only via two distinct clusters
//   The union path is identical to the base file's find/add_edge (same delta
//   convention), so the label computation shares the displacement bookkeeping
//   that D(C) uses.  For the merged-label path (K1 and K2 replay) we carry
//   wrap_x / wrap_y per root instead of the single wrap_ bit.
//
// Outputs per geometry+kernel:
//   - collapsed D per k: a_k = sum_{|C|=k} 1{black NN wraps} - 1{white NN+NNN
//     wraps}.  This is the tripwire: must equal the committed Bernstein
//     integers bit-for-bit.
//   - #646 joint label tables: per k, 5x5 counts [black label][white label].
//   - #640 coarse tables: per k, 4x4 counts [black hxv][white hxv] where each
//     side is 0=none,1=x,2=y,3=both (h/v/both/neither from the issue text).
//   - signed D attributed to (black label, white label) pairs: per k,
//     contribution sums.  Verdict data for #651: D mass in labels other than
//     (both-same, both-same) must be zero for A to continue.
//
// Build: g++ -O2 -std=c++17 -pthread -o wrap_census wrap_census.cpp
// Usage: ./wrap_census --geometry axis|diamond --L N [--kernel k1|k2|both]
//        [--threads T] [--split B] --out prefix.json
//        [--labels k]  print the 5x5 joint table for occupation k

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
  std::vector<Edge> primal;    // black NN, both directions
  std::vector<Edge> matching;  // white NN+NNN, both directions
};

int mod(int a, int p) {
  int r = a % p;
  return r < 0 ? r + p : r;
}

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
// Labels: 0=none 1=x 2=y 3=both-same 4=both-two
// ---------------------------------------------------------------------------
// Joint per-k counters: [black label][white label], plus signed-D tables.
struct Tables {
  int n = 0;
  // occurrence counts (regardless of D)
  std::vector<long long> joint;          // 5*5 per k
  std::vector<long long> coarse;         // 4*4 per k (#640: 0=none 1=x 2=y 3=both)
  // signed D attribution
  std::vector<long long> d_joint;        // +1 at (black label, white label) when db=1
  std::vector<long long> d_joint_neg;    // -1 when dw=1
  std::vector<long long> collapsed;      // per-k D
  long long visited = 0;

  void init(int nn) {
    n = nn;
    joint.assign((n + 1) * 25, 0);
    coarse.assign((n + 1) * 16, 0);
    d_joint.assign((n + 1) * 25, 0);
    d_joint_neg.assign((n + 1) * 25, 0);
    collapsed.assign(n + 1, 0);
    visited = 0;
  }

  void record(int k, int lb, int lw, int db, int dw) {
    joint[(size_t)k * 25 + lb * 5 + lw] += 1;
    const int cb = (lb == 4 ? 3 : lb);  // coarse: both-two folds into both
    const int cw = (lw == 4 ? 3 : lw);
    coarse[(size_t)k * 16 + cb * 4 + cw] += 1;
    if (db) d_joint[(size_t)k * 25 + lb * 5 + lw] += 1;
    if (dw) d_joint_neg[(size_t)k * 25 + lb * 5 + lw] -= 1;
    collapsed[k] += db - dw;
  }

  void merge(const Tables& o) {
    for (size_t i = 0; i < joint.size(); ++i) joint[i] += o.joint[i];
    for (size_t i = 0; i < coarse.size(); ++i) coarse[i] += o.coarse[i];
    for (size_t i = 0; i < d_joint.size(); ++i) d_joint[i] += o.d_joint[i];
    for (size_t i = 0; i < d_joint_neg.size(); ++i) d_joint_neg[i] += o.d_joint_neg[i];
    for (size_t i = 0; i < collapsed.size(); ++i) collapsed[i] += o.collapsed[i];
    visited += o.visited;
  }

  void add_collapsed_only(int k, int db, int dw) { collapsed[k] += db - dw; }
};

const char* LABEL_NAMES[5] = {"none", "x", "y", "both-same", "both-two"};

// ---------------------------------------------------------------------------
// K1: compressing union-find with x/y wrap flags (base WrapDsu + label split).
// The find/add_edge delta convention is identical to the base file.
// ---------------------------------------------------------------------------
class LabelDsu {
 public:
  explicit LabelDsu(int n)
      : n_(n), parent_(n), size_(n, 1), dx_(n, 0), dy_(n, 0), wx_(n, 0),
        wy_(n, 0) {}

  void reset() {
    for (int i = 0; i < n_; ++i) {
      parent_[i] = i;
      size_[i] = 1;
      dx_[i] = 0;
      dy_[i] = 0;
      wx_[i] = 0;
      wy_[i] = 0;
    }
  }

  int find(int x, int& dx, int& dy) {
    if (parent_[x] == x) {
      dx = 0;
      dy = 0;
      return x;
    }
    const int p = parent_[x];
    int px = 0, py = 0;
    const int root = find(p, px, py);
    const int odx = dx_[x], ody = dy_[x];
    parent_[x] = root;
    dx_[x] = odx + px;
    dy_[x] = ody + py;
    dx = dx_[x];
    dy = dy_[x];
    return root;
  }

  void add_edge(int i, int j, int edge_dx, int edge_dy) {
    int ix = 0, iy = 0, jx = 0, jy = 0;
    const int ri = find(i, ix, iy);
    const int rj = find(j, jx, jy);
    const int rdx = ix + edge_dx - jx;
    const int rdy = iy + edge_dy - jy;
    if (ri == rj) {
      if (rdx != 0) wx_[ri] = 1;
      if (rdy != 0) wy_[ri] = 1;
      return;
    }
    if (size_[ri] >= size_[rj]) {
      parent_[rj] = ri;
      dx_[rj] = rdx;
      dy_[rj] = rdy;
      size_[ri] += size_[rj];
      wx_[ri] |= wx_[rj];
      wy_[ri] |= wy_[rj];
    } else {
      parent_[ri] = rj;
      dx_[ri] = -rdx;
      dy_[ri] = -rdy;
      size_[rj] += size_[ri];
      wx_[rj] |= wx_[ri];
      wy_[rj] |= wy_[ri];
    }
  }

  // Label over all active vertices.  Returns (wraps_any, label).
  std::pair<bool, int> label(const std::vector<char>& active) {
    bool ax = false, ay = false, same = false;
    for (int i = 0; i < n_; ++i) {
      if (!active[i]) continue;
      int dx = 0, dy = 0;
      const int r = find(i, dx, dy);
      if (wx_[r]) ax = true;
      if (wy_[r]) ay = true;
    }
    if (ax && ay) {
      for (int i = 0; i < n_; ++i) {
        if (!active[i]) continue;
        int dx = 0, dy = 0;
        const int r = find(i, dx, dy);
        if (wx_[r] && wy_[r]) {
          same = true;
          break;
        }
      }
    }
    const bool any = ax || ay;
    int lab;
    if (!ax && !ay) lab = 0;
    else if (ax && !ay) lab = 1;
    else if (ay && !ax) lab = 2;
    else if (same) lab = 3;
    else lab = 4;
    return {any, lab};
  }

 private:
  int n_;
  std::vector<int> parent_, size_, dx_, dy_;
  std::vector<char> wx_, wy_;
};

void kernel1_threaded(const Geometry& g, int threads, Tables& out,
                      double& seconds) {
  const int n = g.n;
  out.init(n);
  const unsigned long long total = 1ULL << n;
  const unsigned long long chunk = (total + threads - 1) / threads;
  std::vector<Tables> locals(threads);
  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t]() {
      const unsigned long long begin = t * chunk;
      const unsigned long long end = std::min(total, begin + chunk);
      locals[t].init(n);
      if (begin >= end) return;
      LabelDsu dsu(n);
      std::vector<char> black(n), white(n);
      Tables& local = locals[t];
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
        const auto bw = dsu.label(black);
        dsu.reset();
        for (const Edge& e : g.matching)
          if (white[e.i] && white[e.j]) dsu.add_edge(e.i, e.j, e.dx, e.dy);
        const auto ww = dsu.label(white);
        local.record(k, bw.second, ww.second, bw.first ? 1 : 0, ww.first ? 1 : 0);
      }
      local.visited = end - begin;
    });
  }
  for (auto& th : pool) th.join();
  auto t1 = std::chrono::steady_clock::now();
  seconds = std::chrono::duration<double>(t1 - t0).count();
  for (int t = 0; t < threads; ++t) out.merge(locals[t]);
}

// ---------------------------------------------------------------------------
// K2: journaled union-find (no path compression) + DFS, x/y wrap flags.
// The journal/undo logic is identical to the base file; only the wrap flag is
// split into wx_/wy_.
// ---------------------------------------------------------------------------
class JLabelDsu {
 public:
  explicit JLabelDsu(int n)
      : n_(n), parent_(n), size_(n), dx_(n), dy_(n), wx_(n), wy_(n) {}

  void reset() {
    for (int i = 0; i < n_; ++i) {
      parent_[i] = i;
      size_[i] = 1;
      dx_[i] = 0;
      dy_[i] = 0;
      wx_[i] = 0;
      wy_[i] = 0;
    }
    journal_.clear();
  }

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
    const int rdx = ix + edge_dx - jx;
    const int rdy = iy + edge_dy - jy;
    if (ri == rj) {
      if (rdx != 0 && !wx_[ri]) {
        journal_.push_back(Mut{2, ri, 0, 0, 0, 0});
        wx_[ri] = 1;
      }
      if (rdy != 0 && !wy_[ri]) {
        journal_.push_back(Mut{3, ri, 0, 0, 0, 0});
        wy_[ri] = 1;
      }
      return;
    }
    if (size_[ri] >= size_[rj]) {
      journal_.push_back(Mut{1, rj, ri, size_[ri], wx_[ri], wy_[ri]});
      parent_[rj] = ri;
      dx_[rj] = rdx;
      dy_[rj] = rdy;
      size_[ri] += size_[rj];
      wx_[ri] |= wx_[rj];
      wy_[ri] |= wy_[rj];
    } else {
      journal_.push_back(Mut{1, ri, rj, size_[rj], wx_[rj], wy_[rj]});
      parent_[ri] = rj;
      dx_[ri] = -rdx;
      dy_[ri] = -rdy;
      size_[rj] += size_[ri];
      wx_[rj] |= wx_[ri];
      wy_[rj] |= wy_[ri];
    }
  }

  void undo_to(size_t mark) {
    while (journal_.size() > mark) {
      const Mut& m = journal_.back();
      if (m.kind == 1) {
        parent_[m.root] = m.root;
        dx_[m.root] = 0;
        dy_[m.root] = 0;
        size_[m.newroot] = m.old_size_newroot;
        wx_[m.newroot] = m.old_a;
        wy_[m.newroot] = m.old_b;
      } else if (m.kind == 2) {
        wx_[m.root] = 0;
      } else {
        wy_[m.root] = 0;
      }
      journal_.pop_back();
    }
  }

  size_t mark() const { return journal_.size(); }

  // Label over all vertices (active == color != unassigned handled by caller:
  // every vertex is colored at the leaf).
  std::pair<bool, int> label(int n, const uint8_t* color, uint8_t active_val) {
    bool ax = false, ay = false, same = false;
    for (int i = 0; i < n; ++i) {
      if (color[i] != active_val) continue;
      int dx = 0, dy = 0;
      const int r = find(i, dx, dy);
      if (wx_[r]) ax = true;
      if (wy_[r]) ay = true;
    }
    if (ax && ay) {
      for (int i = 0; i < n; ++i) {
        if (color[i] != active_val) continue;
        int dx = 0, dy = 0;
        const int r = find(i, dx, dy);
        if (wx_[r] && wy_[r]) {
          same = true;
          break;
        }
      }
    }
    const bool any = ax || ay;
    int lab;
    if (!ax && !ay) lab = 0;
    else if (ax && !ay) lab = 1;
    else if (ay && !ax) lab = 2;
    else if (same) lab = 3;
    else lab = 4;
    return {any, lab};
  }

 private:
  struct Mut {
    int kind;      // 1=merge, 2=wx set, 3=wy set
    int root;      // child root (merge) or root whose flag was set
    int newroot;   // merge: surviving root
    int old_size_newroot;
    char old_a;
    char old_b;
  };
  int n_;
  std::vector<int> parent_, size_, dx_, dy_;
  std::vector<char> wx_, wy_;
  std::vector<Mut> journal_;
};

struct AdjEntry {
  int u;
  int dx;
  int dy;
};

std::vector<std::vector<AdjEntry>> build_adjacency(const std::vector<Edge>& edges,
                                                   int n) {
  std::vector<std::vector<AdjEntry>> adj(n);
  // No dedup: identical to the base file (directed edges kept individually).
  for (const Edge& e : edges) {
    const int lo = std::min(e.i, e.j);
    const int hi = std::max(e.i, e.j);
    if (lo == hi) continue;
    if (e.i == hi) adj[hi].push_back(AdjEntry{lo, e.dx, e.dy});
    else adj[hi].push_back(AdjEntry{lo, -e.dx, -e.dy});
  }
  return adj;
}

struct K2Ctx {
  JLabelDsu* black;
  JLabelDsu* white;
  const std::vector<std::vector<AdjEntry>>* adj_black;
  const std::vector<std::vector<AdjEntry>>* adj_white;
  uint8_t* color;
  int n;
  Tables* local;
};

void dfs(K2Ctx& ctx, int idx, int black_cnt) {
  if (idx == ctx.n) {
    const int k = black_cnt;
    const auto bw = ctx.black->label(ctx.n, ctx.color, 1);
    const auto ww = ctx.white->label(ctx.n, ctx.color, 0);
    ctx.local->record(k, bw.second, ww.second, bw.first ? 1 : 0,
                      ww.first ? 1 : 0);
    return;
  }
  {
    ctx.color[idx] = 0;
    const size_t m = ctx.white->mark();
    for (const AdjEntry& e : (*ctx.adj_white)[idx])
      if (e.u < idx && ctx.color[e.u] == 0)
        ctx.white->add_edge(idx, e.u, e.dx, e.dy);
    dfs(ctx, idx + 1, black_cnt);
    ctx.white->undo_to(m);
  }
  {
    ctx.color[idx] = 1;
    const size_t m = ctx.black->mark();
    for (const AdjEntry& e : (*ctx.adj_black)[idx])
      if (e.u < idx && ctx.color[e.u] == 1)
        ctx.black->add_edge(idx, e.u, e.dx, e.dy);
    dfs(ctx, idx + 1, black_cnt + 1);
    ctx.black->undo_to(m);
  }
}

void kernel2_threaded(const Geometry& g, int threads, int split, Tables& out,
                      double& seconds) {
  const int n = g.n;
  out.init(n);
  auto adj_black = std::make_shared<std::vector<std::vector<AdjEntry>>>(
      build_adjacency(g.primal, n));
  auto adj_white = std::make_shared<std::vector<std::vector<AdjEntry>>>(
      build_adjacency(g.matching, n));
  if (split <= 0 || split >= n) split = std::min(n, 12);
  const int B = split;
  const unsigned long long tasks = 1ULL << B;

  std::vector<Tables> locals(threads);
  std::atomic<unsigned long long> next_task{0};
  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t]() {
      JLabelDsu black(n), white(n);
      std::vector<uint8_t> color(n, 0);
      locals[t].init(n);
      K2Ctx ctx{&black, &white, &*adj_black, &*adj_white, color.data(), n,
                &locals[t]};
      for (;;) {
        const unsigned long long task = next_task.fetch_add(1);
        if (task >= tasks) break;
        black.reset();
        white.reset();
        std::fill(color.begin(), color.end(), 2);
        int black_cnt = 0;
        for (int v = 0; v < B; ++v) {
          const int c = static_cast<int>((task >> v) & 1ULL);
          color[v] = static_cast<uint8_t>(c);
          if (c) ++black_cnt;
          if (c == 1) {
            for (const AdjEntry& e : (*adj_black)[v])
              if (color[e.u] == 1) black.add_edge(v, e.u, e.dx, e.dy);
          } else {
            for (const AdjEntry& e : (*adj_white)[v])
              if (color[e.u] == 0) white.add_edge(v, e.u, e.dx, e.dy);
          }
        }
        dfs(ctx, B, black_cnt);
      }
      locals[t].visited = 1LL << (n - B);  // per task; merged below as tasks*2^(n-B)
    });
  }
  for (auto& th : pool) th.join();
  auto t1 = std::chrono::steady_clock::now();
  seconds = std::chrono::duration<double>(t1 - t0).count();
  out.init(n);
  for (int t = 0; t < threads; ++t) out.merge(locals[t]);
  out.visited = 1LL << n;
}

void dump_tables_json(FILE* f, const char* key, const std::vector<long long>& v,
                      int n, int stride, const char** names, bool comma) {
  std::fprintf(f, "  \"%s\": [\n", key);
  for (int k = 0; k <= n; ++k) {
    std::fprintf(f, "   [");
    for (int a = 0; a < stride; ++a) {
      for (int b = 0; b < stride; ++b) {
        std::fprintf(f, "%lld%s", v[(size_t)k * stride * stride + a * stride + b],
                     (a + 1 == stride && b + 1 == stride) ? "" : ",");
      }
    }
    std::fprintf(f, "]%s\n", k == n ? "" : ",");
  }
  std::fprintf(f, "  ]%s\n", comma ? "," : "");
}

}  // namespace

int main(int argc, char** argv) {
  std::string geometry_name;
  int L = 0;
  std::string kernel = "both";
  int threads = static_cast<int>(std::thread::hardware_concurrency());
  int split = 0;
  std::string out_path;
  int labels_k = -1;

  for (int i = 1; i < argc; ++i) {
    auto need = [&]() -> const char* { return argv[++i]; };
    if (!std::strcmp(argv[i], "--geometry")) geometry_name = need();
    else if (!std::strcmp(argv[i], "--L")) L = std::atoi(need());
    else if (!std::strcmp(argv[i], "--kernel")) kernel = need();
    else if (!std::strcmp(argv[i], "--threads")) threads = std::atoi(need());
    else if (!std::strcmp(argv[i], "--split")) split = std::atoi(need());
    else if (!std::strcmp(argv[i], "--out")) out_path = need();
    else if (!std::strcmp(argv[i], "--labels-k")) labels_k = std::atoi(need());
    else {
      std::fprintf(stderr, "unknown arg %s\n", argv[i]);
      return 2;
    }
  }
  if (geometry_name.empty() || L <= 0 || out_path.empty()) {
    std::fprintf(stderr,
                 "usage: wrap_census --geometry axis|diamond --L N --out f "
                 "[--kernel k1|k2|both] [--threads T] [--split B] "
                 "[--labels-k K]\n");
    return 2;
  }
  if (threads <= 0) threads = 1;

  Geometry g = geometry_name == "axis" ? axis_geometry(L) : diamond_geometry(L);
  const int n = g.n;

  FILE* f = std::fopen(out_path.c_str(), "w");
  if (!f) {
    std::fprintf(stderr, "cannot open %s\n", out_path.c_str());
    return 3;
  }
  std::fprintf(f, "{\n");
  std::fprintf(f, " \"geometry\": \"%s\",\n", g.name.c_str());
  std::fprintf(f, " \"L\": %d,\n", g.L);
  std::fprintf(f, " \"N\": %d,\n", n);
  std::fprintf(f, " \"physical_period\": \"%s\",\n", g.physical_period.c_str());
  std::fprintf(f, " \"threads\": %d,\n", threads);
  std::fprintf(f, " \"configs_expect\": %lld,\n", 1LL << n);

  auto run_kernel = [&](const char* kn, bool comma_after_kernel) {
    Tables t;
    double seconds = 0.0;
    if (!std::strcmp(kn, "k1")) {
      if (n > 40) {
        std::fprintf(stderr, "K1 refuses N=%d > 40\n", n);
        std::exit(3);
      }
      kernel1_threaded(g, threads, t, seconds);
    } else {
      kernel2_threaded(g, threads, split, t, seconds);
    }
    std::fprintf(f, " \"%s\": {\n", kn);
    std::fprintf(f, "  \"configs_visited\": %lld,\n", t.visited);
    std::fprintf(f, "  \"wall_seconds\": %.3f,\n", seconds);
    std::fprintf(f, "  \"collapsed_D_per_k\": [");
    for (int k = 0; k <= n; ++k)
      std::fprintf(f, "%lld%s", t.collapsed[k], k == n ? "" : ",");
    std::fprintf(f, "],\n");
    dump_tables_json(f, "joint_label_counts_per_k", t.joint, n, 5, LABEL_NAMES, true);
    dump_tables_json(f, "coarse_4x4_counts_per_k", t.coarse, n, 4, nullptr, true);
    dump_tables_json(f, "D_positive_attribution_per_k", t.d_joint, n, 5, LABEL_NAMES, true);
    dump_tables_json(f, "D_negative_attribution_per_k", t.d_joint_neg, n, 5, LABEL_NAMES, false);
    std::fprintf(f, " }%s\n", comma_after_kernel ? "," : "");

    // human-readable 5x5 for one k
    if (labels_k >= 0 && labels_k <= n) {
      std::printf("joint labels at k=%d (%s):\n", labels_k, kn);
      for (int a = 0; a < 5; ++a) {
        std::printf("  %-10s", LABEL_NAMES[a]);
        for (int b = 0; b < 5; ++b)
          std::printf(" %12lld",
                      t.joint[(size_t)labels_k * 25 + a * 5 + b]);
        std::printf("\n");
      }
      std::printf("  columns: none/x/y/both-same/both-two\n");
    }
    std::fflush(f);
    return t;
  };

  Tables t1res, t2res;
  bool have1 = false, have2 = false;
  if (kernel == "k1" || kernel == "both") {
    t1res = run_kernel("k1", kernel == "both");
    have1 = true;
  }
  if (kernel == "k2" || kernel == "both") {
    t2res = run_kernel("k2", false);
    have2 = true;
  }
  if (have1 && have2) {
    bool agree = t1res.collapsed == t2res.collapsed;
    std::fprintf(f, " ,\"k1_k2_collapsed_agree\": %s\n", agree ? "true" : "false");
  }
  std::fprintf(f, "}\n");
  std::fclose(f);
  return 0;
}
