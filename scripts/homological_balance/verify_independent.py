#!/usr/bin/env python3
"""
独立交叉验证 exact_torus_enum.py 的 ambient_rank 计算。

用「单 BFS 提升坐标」方法（与基本圈/生成树方法在实现上独立）：
  - 对每个连通分量做一次 BFS，记录每个顶点的提升坐标 p(u) ∈ Z^2；
  - 对每条边 (u,v)，winding = (p(u) + step(u->v) - p(v)) / L；
  - 所有非零 winding 的 Q-秩 = 该分量（及总和）的环境 H1 像秩。

同时验证若干手工可判定的配置。
"""
from collections import deque, defaultdict
from itertools import product


def lift_step(v, w, L):
    (i, j), (a, b) = v, w
    dx = (a - i) % L
    dy = (b - j) % L
    if dx == L - 1:
        dx = -1
    if dy == L - 1:
        dy = -1
    return dx, dy


def neighbors(i, j, L, diagonal):
    out = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            if not diagonal and abs(di) + abs(dj) != 1:
                continue
            out.append(((i + di) % L, (j + dj) % L))
    return out


def rank_z2(vecs):
    indep = []
    for (a, b) in vecs:
        if a == 0 and b == 0:
            continue
        # 判断 (a,b) 是否与 indep 有理线性无关
        free = True
        if len(indep) >= 2:
            break
        if len(indep) == 1:
            (a0, b0) = indep[0]
            if a0 * b - b0 * a == 0:
                free = False
        if free:
            indep.append((a, b))
            if len(indep) == 2:
                break
    return len(indep)


def ambient_rank_v2(black, L, diagonal):
    """独立实现：单 BFS 提升坐标。"""
    black = list(black)
    S = set(black)
    if not S:
        return 0
    idx = {v: i for i, v in enumerate(black)}
    n = len(black)
    seen = [False] * n
    total_windings = []
    for s in range(n):
        if seen[s]:
            continue
        # BFS 求分量内每个顶点的提升坐标 p
        p = {s: black[s]}  # vertex idx -> (X, Y)
        q = deque([s])
        seen[s] = True
        order = [s]
        while q:
            u = q.popleft()
            (i, j) = black[u]
            (X, Y) = p[u]
            for (a, b) in neighbors(i, j, L, diagonal):
                if (a, b) not in S:
                    continue
                v = idx[(a, b)]
                dx, dy = lift_step((i, j), (a, b), L)
                if v not in p:
                    p[v] = (X + dx, Y + dy)
                    seen[v] = True
                    q.append(v)
                    order.append(v)
        # 该分量的边 winding
        windings = []
        for u in order:
            (i, j) = black[u]
            (X, Y) = p[u]
            for (a, b) in neighbors(i, j, L, diagonal):
                if (a, b) not in S:
                    continue
                v = idx[(a, b)]
                if u > v:  # 每条无向边只算一次
                    continue
                dx, dy = lift_step((i, j), (a, b), L)
                wx_num = X + dx - p[v][0]
                wy_num = Y + dy - p[v][1]
                assert wx_num % L == 0 and wy_num % L == 0, (black[u], black[v], L)
                wx, wy = wx_num // L, wy_num // L
                if wx != 0 or wy != 0:
                    windings.append((wx, wy))
        r_comp = rank_z2(windings)
        # 各分量的 winding 子群之和的秩 <= 2
        total_windings.extend(windings)
    return rank_z2(total_windings)


def manual_checks():
    """手工可判定配置。"""
    from exact_torus_enum import ambient_rank as ar1, ambient_rank as ar_orig
    L = 3
    allpts = list(product(range(L), range(L)))
    cases = {
        "全黑(9点)": set(allpts),
        "空(0点)": set(),
        "单点(0,0)": {(0, 0)},
        "第0行3点(水平环)": {(0, 0), (0, 1), (0, 2)},
        "第0列3点(垂直环)": {(0, 0), (1, 0), (2, 0)},
        "对角(0,0),(1,1),(2,2)": {(0, 0), (1, 1), (2, 2)},
        "两条交叉环(0行+0列)": {(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)},
    }
    print("=== 手工配置验证 (L=3, 黑=4连通 / 白=8连通) ===")
    for name, B in cases.items():
        rb = ambient_rank_v2(B, L, False)
        W = set(allpts) - B
        rw = ambient_rank_v2(W, L, True)
        print(f"{name:28s}  r_b={rb}  r_w={rw}  r_b+r_w={rb+rw}  (应=2)")
    return


def cross_check(L):
    from exact_torus_enum import ambient_rank as ar1
    N = L * L
    mismatch = 0
    for mask in range(1 << N):
        B = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        r1 = ar1(B, L, False)
        r2 = ambient_rank_v2(B, L, False)
        if r1 != r2:
            mismatch += 1
            if mismatch <= 5:
                print(f"  MISMATCH mask={mask} B={B} r1={r1} r2={r2}")
    print(f"L={L}: 两套独立算法对 {2**N} 个配置的黑秩结果不一致数 = {mismatch}")
    return mismatch


if __name__ == "__main__":
    manual_checks()
    print()
    total_bad = 0
    for L in (3, 4):
        total_bad += cross_check(L)
    print()
    print("交叉验证总不一致数:", total_bad, "(0 = 两套算法完全一致)")
