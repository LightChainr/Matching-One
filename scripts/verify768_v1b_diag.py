# -*- coding: utf-8 -*-
"""
v768_v1b_diag.py --- 判定 e_n = dim(L_{-1} 在不可约商 level n 的像) 到底是多少

theory768 (g2_virasoro.py) 用 proj_quot 得到 e_4 = 2 (=d_3)；
我用「complement + mod rad_n」得到 e_4 = 1。两者若都对就矛盾，必须判定。

这里给 5 条互相独立的路线（n=3,4,5）：
  A. dim((L_{-1}P_{n-1} + rad_n)/rad_n)                       —— 我的路线
  B. dim((L_{-1}V_{n-1}  + rad_n)/rad_n)                      —— theory768 的路线(等价形式)
  C. rank(诱导映射) 用显式商坐标（rref 主元坐标做代表）
  D. d_n - (dim W - dim rad_n),  W={v∈V_n : L_1 v ∈ rad_{n-1}}  —— 用 L_1 的核（伴随）
  E. 显式找 0 != u ∈ V_{n-1}, u ∉ rad_{n-1}, 且 L_{-1}u ∈ rad_n，
     并给出证书：<L_{-1}u , x> = 0 对所有 x ∈ V_n 的基。
  F. 分解证书：dim(L_{-1}V_{n-1} ∩ rad_n) 显式算出来。

PY39COMPAT_MARKER=1
"""
import json
from fractions import Fraction

import v768_v1_verma as V

NMAX = 6
CV, HV = "0", "5/8"


def rref(rows, ncols):
    return V.frac_rref(rows, ncols)


def rank(rows, ncols):
    return V.frac_rank(rows, ncols)


def rowspace_basis(rows, ncols):
    if not rows:
        return []
    R, piv = rref(rows, ncols)
    return [R[i] for i in range(len(piv))]


def mat_of(mode, nfrom, nto, B):
    """L_mode : V_nfrom -> V_nto ，形状 (len(nto), len(nfrom))"""
    if nto < 0:
        return []
    M = [[Fraction(0)] * len(B[nfrom]) for _ in range(len(B[nto]))]
    for j, w in enumerate(B[nfrom]):
        d = V.mul_left(CV, HV, mode, {w: Fraction(1)})
        for kw, kv in d.items():
            if kw in B[nto]:
                M[B[nto].index(kw)][j] = kv
    return M


def complement_indices(rad_rows, n):
    """返回标准基下标，与 rad_rows 组成 V_n 的一组基（即商空间代表）"""
    span = [list(r) for r in rad_rows]
    cur = rank(span, n) if span else 0
    comp = []
    for i in range(n):
        e = [Fraction(0)] * n
        e[i] = Fraction(1)
        if rank(span + [e], n) > cur:
            span.append(e)
            comp.append(i)
            cur += 1
            if cur == n:
                break
    return comp


def main():
    out = {"marker": "PY39COMPAT_MARKER=1", "c": CV, "h": HV, "nmax": NMAX}
    B = {n: V.level_basis(n) for n in range(NMAX + 1)}
    G = {n: V.gram_matrix(CV, HV, n, B) for n in range(NMAX + 1)}
    d = {n: rank(G[n], len(B[n])) for n in range(NMAX + 1)}
    rad = {n: V.frac_nullspace(G[n]) for n in range(NMAX + 1)}
    out["d"] = d
    out["dim_rad"] = {n: len(rad[n]) for n in range(NMAX + 1)}
    res = {}
    for n in range(2, NMAX + 1):
        bn, bp = len(B[n]), len(B[n - 1])
        A = mat_of(-1, n - 1, n, B)              # L_{-1}: V_{n-1} -> V_n  (bn x bp)
        Rn = [list(v) for v in rad[n]]
        Rp = [list(v) for v in rad[n - 1]]
        cA = complement_indices([list(r) for r in Rp], len(B[n - 1]))
        # A.
        colsA = [[A[r][i] for r in range(bn)] for i in cA]
        eA = rank(colsA + Rn, bn) - len(Rn) if (colsA or Rn) else 0
        # B.
        colsAll = [[A[r][i] for r in range(bn)] for i in range(bp)]
        eB = rank(colsAll + Rn, bn) - len(Rn) if (colsAll or Rn) else 0
        # C. 诱导映射的显式商坐标
        cN = complement_indices([list(r) for r in Rn], bn)
        # 把 V_n 的向量投影到 non-pivot 坐标（与 complement 等价的投影）
        proj = {}
        for i in cN:
            e = [Fraction(0)] * bn
            e[i] = Fraction(1)
            # 沿 rad 方向消掉主元
            outv = e[:]
            # 用 rad 的 rref 消主元
            Rb = rowspace_basis([list(r) for r in Rn], bn)
            if Rb:
                Rr, piv = rref(Rb, bn)
                acc = [Fraction(0)] * bn
                for ii, pc in enumerate(piv):
                    acc = [a + outv[pc] * Rr[ii][k] for k, a in enumerate(acc)]
                outv = [a - b for a, b in zip(outv, acc)]
            proj[i] = outv
        colsC = []
        for i in cA:
            v = [A[r][i] for r in range(bn)]
            # 投影
            Rb = rowspace_basis([list(r) for r in Rn], bn)
            outv = v[:]
            if Rb:
                Rr, piv = rref(Rb, bn)
                acc = [Fraction(0)] * bn
                for ii, pc in enumerate(piv):
                    acc = [a + outv[pc] * Rr[ii][k] for k, a in enumerate(acc)]
                outv = [a - b for a, b in zip(outv, acc)]
            colsC.append(outv)
        eC = rank(colsC, bn) if colsC else 0
        # D. 用 L_1
        A1 = mat_of(1, n, n - 1, B)              # L_1: V_n -> V_{n-1}
        Rb = rowspace_basis([list(r) for r in Rp], len(B[n - 1]))
        if Rb:
            Rr, piv = rref(Rb, len(B[n - 1]))
            PA1 = [[Fraction(0)] * bn for _ in range(len(B[n - 1]))]
            for r in range(len(B[n - 1])):
                for j in range(bn):
                    PA1[r][j] = A1[r][j]
            for ii, pc in enumerate(piv):
                for r in range(len(B[n - 1])):
                    f = PA1[pc][0] if False else None
            # 显式：W = ker(P·A1)，P 是投影到 rad_{n-1} 的补
            PA = []
            for j in range(bn):
                col = [A1[r][j] for r in range(len(B[n - 1]))]
                acc = [Fraction(0)] * len(B[n - 1])
                for ii, pc in enumerate(piv):
                    acc = [a + col[pc] * Rr[ii][k] for k, a in enumerate(acc)]
                PA.append([a - b for a, b in zip(col, acc)])
            # PA[j] 是长度 (n-1) 的列；组成矩阵 (n-1) x n，求核维数
            Mt = [[PA[j][r] for j in range(bn)] for r in range(len(B[n - 1]))]
            dimW = bn - rank(Mt, bn)
        else:
            dimW = bn - rank(A1, bn)
        eD = d[n] - (dimW - len(Rn))
        # F. dim(L_{-1}V_{n-1} ∩ rad_n)
        # L_{-1}V_{n-1} ∩ rad_n 的维数 = dim L_{-1}V_{n-1} + dim rad_n - dim(L_{-1}V_{n-1}+rad_n)
        dimImg = rank(colsAll, bn) if colsAll else 0
        inter = len(Rn) + dimImg - (eB + len(Rn))
        res[n] = {"d_n": d[n], "eA": eA, "eB": eB, "eC": eC, "eD": eD,
                  "dLm1_prev": d[n - 1], "dimW": dimW,
                  "dimLm1V_prev": dimImg, "dim_inter_with_rad_n": inter}
        print("n=%d d=%d | eA=%d eB=%d eC=%d eD=%d | d_{n-1}=%d dim(L_-1 V_prev)=%d "
              "inter=%d dimW=%d" % (n, d[n], eA, eB, eC, eD, d[n - 1], dimImg, inter, dimW),
              flush=True)
    out["per_level"] = res

    # E. 显式核向量证书（n=4）
    n = 4
    bn, bp = len(B[n]), len(B[n - 1])
    A = mat_of(-1, n - 1, n, B)
    Rn = [list(v) for v in rad[n]]
    Rp = [list(v) for v in rad[n - 1]]
    # 求 u ∈ V_3 使 L_{-1}u ∈ rad_4 且 u ∉ rad_3
    # 即 A u ∈ span(Rn)。构造投影 P4 消掉 rad_4 方向后求核。
    Rb = rowspace_basis(Rn, bn)
    Rr, piv = rref(Rb, bn)
    PAu = []
    for j in range(bp):
        col = [A[r][j] for r in range(bn)]
        acc = [Fraction(0)] * bn
        for ii, pc in enumerate(piv):
            acc = [a + col[pc] * Rr[ii][k] for k, a in enumerate(acc)]
        PAu.append([a - b for a, b in zip(col, acc)])
    Mt = [[PAu[j][r] for j in range(bp)] for r in range(bn)]
    ns = V.frac_nullspace(Mt)
    cert = []
    for u in ns:
        # 检查 u ∉ rad_3
        inrad3 = rank([list(u)] + [list(r) for r in Rp], bp) == len(Rp)
        # 检查 L_{-1}u ⊥ V_4
        Au = [sum(A[r][j] * u[j] for j in range(bp)) for r in range(bn)]
        # L_{-1}u 的范数/正交性：<L_{-1}u, x> = u^T A^T G_n x
        orth = all(sum(Au[r] * G[n][r][c] * x[c] for r in range(bn) for c in range(bn)) == 0
                   for x in [[Fraction(1) if k == i else Fraction(0) for k in range(bn)]
                             for i in range(bn)])
        selfnorm = sum(Au[r] * G[n][r][c] * Au[c] for r in range(bn) for c in range(bn))
        cert.append({"u": [str(x) for x in u], "u_in_rad3": inrad3,
                     "Lm1u_orthogonal_to_V4": orth, "Lm1u_norm": str(selfnorm),
                     "Lm1u_is_zero": all(x == 0 for x in Au)})
    out["kernel_certificate_n4"] = cert
    print("kernel certificate n=4:", json.dumps(cert, ensure_ascii=False), flush=True)

    with open("v768_v1b_diag.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
