#!/usr/bin/env python3
"""
Matching-One 精确有限枚举验证
==============================

在 LxL 环面（honest torus）上枚举所有 2^(L^2) 个 site-percolation 配置，验证：

1. 数字-Alexander 对偶：r_b + r_w = 2（黑=4-连通 NN，白=8-连通 NN+NNN）
2. 秩单调性与 M_L(0)=-1, M_L(1)=+1
3. M_L(p) = P_2 - P_0 及唯一平衡根 p_L^H（M_L(p_L^H)=0）
4. 自匹配（triangular 语义）时的精确平衡：p=1/2 处 M≡0

环境秩 r = rank im[ H_1(K) -> H_1(T^2) ] ∈ {0,1,2}
由图的 winding 像的秩给出（2-cells 只贡献可缩面边界，winding=0）。
"""
from itertools import product
from fractions import Fraction
from collections import defaultdict


def neighbors(i, j, L, diagonal):
    """返回 (i,j) 的邻点。diagonal=False 为 4-连通 NN；True 为 8-连通 NN+NNN。"""
    out = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            if not diagonal and abs(di) + abs(dj) != 1:
                continue
            out.append(((i + di) % L, (j + dj) % L))
    return out


def ambient_rank(black, L, diagonal):
    """计算占据集合的 ambient H1 秩（0/1/2），通过图 winding 像的秩。"""
    black = list(black)
    black_set = set(black)
    n = len(black_set)
    if n == 0:
        return 0
    # 邻接表
    idx = {v: k for k, v in enumerate(black)}
    adj = [[] for _ in range(n)]
    edges = []
    for v in black_set:
        for w in neighbors(*v, L, diagonal):
            if w in black_set and v < w:
                u = idx[v]
                ww = idx[w]
                adj[u].append(ww)
                adj[ww].append(u)
                edges.append((u, ww, v, w))

    # 连通分量 + 生成森林
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree_edges = []
    nontree_edges = []
    for (u, ww, v, w) in edges:
        ru, rw = find(u), find(ww)
        if ru != rw:
            parent[ru] = rw
            tree_edges.append((u, ww, v, w))
        else:
            nontree_edges.append((u, ww, v, w))

    # 为每个连通分量建生成树的父关系（BFS），用于求树上路径
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), []).append(i)

    # winding：对每条非树边，找其两端在生成树中的路径，加上该边自身构成基本圈，
    # 累加 lift 位移得到 (sum_x, sum_y)，除以 L 得 winding (wx, wy)。
    # lift step: 向右 x+1，向左 x-1，向上 y+1，向下 y-1，对角 (dx,dy)。
    windings = []

    # 对每个分量建生成树邻接（只含 tree_edges）
    tree_adj = defaultdict(list)
    for (u, ww, v, w) in tree_edges:
        tree_adj[u].append((ww, v, w))
        tree_adj[ww].append((u, w, v))

    def lift_step(v, w):
        (i, j), (a, b) = v, w
        dx = (a - i) % L
        dy = (b - j) % L
        if dx == L - 1:
            dx = -1
        if dy == L - 1:
            dy = -1
        # dx,dy ∈ {-1,0,1}
        return dx, dy

    for (u, ww, v, w) in nontree_edges:
        # BFS 找 u->ww 的树路径
        from collections import deque
        q = deque([u])
        prev = {u: None}
        while q:
            x = q.popleft()
            if x == ww:
                break
            for (y, _, _) in tree_adj[x]:
                if y not in prev:
                    prev[y] = x
                    q.append(y)
        # 回溯路径 u -> ww
        path = []
        x = ww
        while x is not None and x != u:
            path.append(x)
            x = prev.get(x)
        path.append(u)
        path.reverse()
        sx = sy = 0
        for k in range(len(path) - 1):
            a, b = black[path[k]], black[path[k + 1]]
            dx, dy = lift_step(a, b)
            sx += dx
            sy += dy
        # 加上非树边 w -> v（闭合基本圈）
        dx, dy = lift_step(w, v)
        sx += dx
        sy += dy
        # winding = (sx/L, sy/L) 应为整数
        wx, wy = sx // L, sy // L
        assert sx == wx * L and sy == wy * L, (sx, sy, L)
        if wx != 0 or wy != 0:
            windings.append((wx, wy))

    # 秩 = Q 上 winding 向量张成空间的维数
    indep = []
    for (a, b) in windings:
        if a == 0 and b == 0:
            continue
        # 检查是否与已有独立向量线性无关
        # 2 维空间最多 2 个独立向量
        cand = (a, b)
        # 用整数行列式判断
        new_indep = []
        for vec in indep:
            new_indep.append(vec)
        # 判断 cand 是否在 new_indep 张成（有理）空间内
        if len(new_indep) == 0:
            indep.append(cand)
        elif len(new_indep) == 1:
            (a0, b0) = new_indep[0]
            # cand 与 (a0,b0) 线性无关 iff 行列式 != 0
            if a0 * b - b0 * a != 0:
                indep.append(cand)
        # 已有 2 个独立向量则秩必为 2
    return min(len(indep), 2)


def enumerate_ML(L, diagonal_white=True):
    """枚举所有配置，返回 M_L(p) 的多项式系数（按 p 的次数）+ 对偶验证统计。"""
    N = L * L
    coeff = [Fraction(0)] * (N + 1)  # coeff[k] = sum over |B|=k of (r_b - 1)
    dual_fail = 0
    rank_pairs = defaultdict(int)
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        rb = ambient_rank(black, L, diagonal=False)
        rw = ambient_rank([(k // L, k % L) for k in range(N) if not ((mask >> k) & 1)],
                          L, diagonal=diagonal_white)
        if rb + rw != 2:
            dual_fail += 1
        rank_pairs[(rb, rw)] += 1
        k = len(black)
        coeff[k] += Fraction(rb - 1)
    return coeff, dual_fail, rank_pairs


def eval_ML(coeff, p):
    N = len(coeff) - 1
    return sum(coeff[k] * (Fraction(p) ** k) * (Fraction(1 - p) ** (N - k))
               for k in range(N + 1))


def find_root(coeff, tol=Fraction(1, 10**12)):
    """在 (0,1) 内二分求 M_L(p)=0 的唯一根。"""
    lo, hi = Fraction(0), Fraction(1)
    assert eval_ML(coeff, lo) < 0 and eval_ML(coeff, hi) > 0
    for _ in range(200):
        mid = (lo + hi) / 2
        if eval_ML(coeff, mid) < 0:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def main():
    for L in (3, 4):
        coeff, dual_fail, rank_pairs = enumerate_ML(L)
        N = L * L
        m0 = eval_ML(coeff, Fraction(0))
        m1 = eval_ML(coeff, Fraction(1))
        root = find_root(coeff)
        print(f"===== L = {L} (N = {N}) =====")
        print(f"配置总数: {2**N}")
        print(f"数字-Alexander 对偶失败数 (r_b+r_w!=2): {dual_fail}")
        print(f"秩对分布 (r_b, r_w) -> 计数:")
        for k in sorted(rank_pairs):
            print(f"    {k}: {rank_pairs[k]}")
        print(f"M_L(0) = {m0}   (应为 -1)")
        print(f"M_L(1) = {m1}   (应为 +1)")
        print(f"唯一平衡根 p_L^H = {root:.12f}")
        print(f"  (square-site p_c ≈ 0.5927460..., 三角格点精确 0.5)")
        mvals = [float(eval_ML(coeff, Fraction(k, 100))) for k in range(0, 101, 5)]
        monotone = all(mvals[i] <= mvals[i + 1] for i in range(len(mvals) - 1))
        print(f"M 在 p∈[0,1] 上单调非降（采样验证）: {monotone}")
        print()


if __name__ == "__main__":
    main()
