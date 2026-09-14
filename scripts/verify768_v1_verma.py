# -*- coding: utf-8 -*-
"""
v768_v1_verma.py  ---  verify768 / V1

独立重算 (不用 theory768 的脚本、不用它的 Gram 实现)：
  (A) Verma 等级维数 p(n)（自写 DP）
  (B) c=0,h=5/8 的 Gram 矩阵 G_n 的 rank / nullity
      -> 判 rank(G_n) == p(n) - p(n-2)（这次用正确公式！theory768 的 JSON 里
         rank_equals_p_shift2=false 是脚本 bug：它把 rank 拿去和 p(n-2) 比）
  (C) 各级奇异向量维数
  (D) chi = -3 L_{-2}|h> + 2 L_{-1}^2|h> 的范数（两条独立路径）
  (E) level 1..4 模去 L_{-1} 后的维数（不可约商 & Verma 两种口径）
  (F) 把 0,0,1,1 是否只是「rank=p(n)-p(n-2) + L_{-1} 单射」的算术重述单独判出来

约化机制（与 theory768 不同、且结构上必然终止）：
  把 L_m 逐个从右到左乘到一个「已正规序」的 PBW 词（负模升序）上，
  用递归 insert/push 实现：
     L_m L_a = L_a L_m + (m-a) L_{m+a}   (+ (C/12)(m^3-m) delta_{m+a,0})
  每次递归 word 长度严格减少 => 终止。全部分数精确。

PY39COMPAT_MARKER=1
"""
import json
import os
import sys
from fractions import Fraction

assert (255).bit_count() == 8, "py39 compat shim not active"
MARK = "PY39COMPAT_MARKER=1"
NMAX = int(os.environ.get("V768_NMAX", "7"))


def partitions_counts(nmax):
    p = [0] * (nmax + 1)
    p[0] = 1
    for k in range(1, nmax + 1):
        for n in range(k, nmax + 1):
            p[n] += p[n - k]
    return p


# ---------- 核心：单模乘到正规序词上 ----------
def reduce_one(cval, hval, m, word, memo=None):
    """word: 正规序 PBW 词（负模升序，即模长降序）；返回 dict word->coeff。
    L_m 乘在 word 左端。"""
    C = Fraction(cval)
    if memo is None:
        memo = {}
    key = (m, word)
    if key in memo:
        return memo[key]
    if not word:
        if m < 0:
            r = {(m,): Fraction(1)}
        elif m == 0:
            r = {(): Fraction(hval)}
        else:
            r = {}
        memo[key] = r
        return r
    a = word[0]
    rest = word[1:]
    if m < 0 and m <= a:
        r = {(m,) + word: Fraction(1)}
        memo[key] = r
        return r
    # t1 = L_a · (L_m · |rest>)
    t1 = {}
    for w, c in reduce_one(cval, hval, m, rest, memo).items():
        for w2, c2 in reduce_one(cval, hval, a, w, memo).items():
            t1[w2] = t1.get(w2, Fraction(0)) + c * c2
    t2 = {}
    c = m + a
    if c == 0:
        for w, cc in reduce_one(cval, hval, 0, rest, memo).items():
            t2[w] = t2.get(w, Fraction(0)) + (m - a) * cc
        cent = C / 12 * (m ** 3 - m)
        t2[rest] = t2.get(rest, Fraction(0)) + cent
    else:
        for w, cc in reduce_one(cval, hval, c, rest, memo).items():
            t2[w] = t2.get(w, Fraction(0)) + (m - a) * cc
    r = {}
    for w, v in t1.items():
        r[w] = r.get(w, Fraction(0)) + v
    for w, v in t2.items():
        r[w] = r.get(w, Fraction(0)) + v
    r = {w: v for w, v in r.items() if v != 0}
    memo[key] = r
    return r


def act(cval, hval, state, m):
    """L_m · state"""
    acc = {}
    for word, coeff in state.items():
        for w, v in reduce_one(cval, hval, m, word).items():
            acc[w] = acc.get(w, Fraction(0)) + coeff * v
    return {w: v for w, v in acc.items() if v != 0}


def reduce_seq(cval, hval, seq):
    st = {(): Fraction(1)}
    for m in reversed(tuple(seq)):
        st = act(cval, hval, st, m)
    return st


def mul_left(cval, hval, m, state):
    return act(cval, hval, state, m)


def level_basis(n):
    res = []

    def rec(rem, biggest, cur):
        if rem == 0:
            res.append(tuple(sorted(cur, reverse=True)))
            return
        for k in range(min(rem, biggest), 0, -1):
            rec(rem - k, k, cur + [k])

    rec(n, n, [])
    out = sorted(set(tuple(sorted((-x for x in w))) for w in res))
    return out


def gram_entry(cval, hval, wi, wj):
    """<h| wi^dagger wj |h>；wi,wj 为负模升序 PBW 词"""
    As = [-m for m in wi]
    seq = list(reversed(As)) + list(wj)
    d = reduce_seq(cval, hval, tuple(seq))
    return d.get((), Fraction(0))


def gram_matrix(cval, hval, n, B):
    idx = {w: i for i, w in enumerate(B[n])}
    G = [[Fraction(0)] * len(B[n]) for _ in range(len(B[n]))]
    for wi in B[n]:
        for wj in B[n]:
            G[idx[wi]][idx[wj]] = gram_entry(cval, hval, wi, wj)
    return G


# ---------- 精确有理数线性代数 ----------
def frac_rref(rows, ncols):
    M = [[Fraction(x) for x in r] for r in rows]
    m = len(M)
    piv_cols = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, m):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
    return M, piv_cols


def frac_rank(rows, ncols):
    if not rows:
        return 0
    if rows and len(rows[0]) == 0:
        return 0
    M, piv = frac_rref(rows, ncols)
    return len(piv)


def frac_nullspace(G):
    m = len(G)
    n = len(G[0]) if m else 0
    if m == 0:
        return []
    M, piv_cols = frac_rref(G, n)
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for f in free:
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for ri, pc in enumerate(piv_cols):
            v[pc] = -M[ri][f]
        basis.append(v)
    return basis


def run(cval, hval, label):
    res = {"label": label, "c": str(cval), "h": str(hval)}
    p = partitions_counts(NMAX)
    res["p"] = p
    B = {n: level_basis(n) for n in range(NMAX + 1)}
    res["verma_dims"] = {n: len(B[n]) for n in range(NMAX + 1)}
    assert res["verma_dims"] == {n: p[n] for n in range(NMAX + 1)}, "PBW 计数 != p(n)"
    G = {n: gram_matrix(cval, hval, n, B) for n in range(NMAX + 1)}
    rank = {n: frac_rank(G[n], len(B[n])) for n in range(NMAX + 1)}
    res["gram_rank"] = rank
    res["gram_nullity"] = {n: len(B[n]) - rank[n] for n in range(NMAX + 1)}
    ok = {}
    for n in range(NMAX + 1):
        tgt = p[n] - (p[n - 2] if n >= 2 else 0)
        ok[n] = (rank[n] == tgt)
    res["rank_eq_p_minus_p_shift2"] = ok
    res["rank_eq_p_minus_p_shift2_all"] = all(ok.values())
    res["theory768_style_flag_p_shift2_only"] = all(
        rank[n] == (p[n - 2] if n >= 2 else 0) for n in range(NMAX + 1))
    # (C) 奇异向量维数 = dim ker( all L_m : V_n -> V_{n-m}, m=1..n )
    sing = {}
    for n in range(NMAX + 1):
        if n == 0:
            sing[n] = 0
            continue
        rows = []
        for mm in range(1, n + 1):
            # 形状 (len(B[n-mm]), len(B[n]))：列 = V_n 的基
            M = [[Fraction(0)] * len(B[n]) for _ in range(len(B[n - mm]))]
            for j, wj in enumerate(B[n]):
                d = mul_left(cval, hval, mm, {wj: Fraction(1)})
                for w, v in d.items():
                    if w in B[n - mm]:
                        M[B[n - mm].index(w)][j] = v
            rows.extend(M)
        rr = frac_rank(rows, len(B[n])) if rows else 0
        sing[n] = len(B[n]) - rr
    res["singular_dims"] = sing
    # (D) chi
    if Fraction(hval) == Fraction(5, 8):
        w2 = B[2]
        i2, i11 = w2.index((-2,)), w2.index((-1, -1))
        chi = [Fraction(-3), Fraction(2)]
        q1 = sum(chi[i] * G[2][i][j] * chi[j] for i in range(2) for j in range(2))
        res["chi_int"] = [str(x) for x in chi]
        res["chi_norm_gram"] = str(q1)
        acc = {}
        for w, v in reduce_seq(cval, hval, (-2,)).items():
            acc[w] = acc.get(w, Fraction(0)) + Fraction(-3) * v
        for w, v in reduce_seq(cval, hval, (-1, -1)).items():
            acc[w] = acc.get(w, Fraction(0)) + Fraction(2) * v
        keys = sorted(acc)
        q2 = Fraction(0)
        for u in keys:
            for v in keys:
                q2 += acc[u] * acc[v] * gram_entry(cval, hval, u, v)
        res["chi_expansion"] = {str(k): str(v) for k, v in acc.items()}
        res["chi_norm_direct"] = str(q2)
        res["chi_agrees"] = (q1 == q2)
        for mm in (1, 2, 3):
            d = mul_left(cval, hval, mm, acc)
            res["L%d_chi" % mm] = {str(k): str(v) for k, v in d.items()}
            res["L%d_chi_zero" % mm] = (len(d) == 0)
        res["Lm2_norm"] = str(G[2][i2][i2])
        res["Lm1sq_norm"] = str(G[2][i11][i11])
        res["detG2_c0"] = None
    # det G_2 的通式检查
    if cval == "0":
        w2 = B[2]
        i2, i11 = w2.index((-2,)), w2.index((-1, -1))
        a, b, d = G[2][i2][i2], G[2][i2][i11], G[2][i11][i11]
        num = Fraction(a * d - b * b)
        res["detG2_c0_factor"] = str(num)
        H = Fraction(hval)
        res["detG2_expected_4h2_8h_minus_5"] = str(4 * H * H * (8 * H - 5))
    # (E) mod L_{-1}
    def modLm1_irred(n):
        if n == 0:
            return rank[0]
        Bn, Bp = B[n], B[n - 1]
        Rn = frac_nullspace(G[n])
        Rp = frac_nullspace(G[n - 1])
        span = [list(v) for v in Rp]
        cur = frac_rank(span, len(Bp)) if span else 0
        comp = []
        for i in range(len(Bp)):
            e = [Fraction(0)] * len(Bp)
            e[i] = Fraction(1)
            if frac_rank(span + [e], len(Bp)) > cur:
                span.append(e)
                comp.append(i)
                cur += 1
                if cur == len(Bp):   # 张满 V_{n-1}（不是 cur==rank[n-1]！）
                    break
        assert cur == len(Bp), (n, cur, len(Bp))
        cols = []
        for i in comp:
            d = mul_left(cval, hval, -1, {Bp[i]: Fraction(1)})
            vec = [Fraction(0)] * len(Bn)
            for kw, kv in d.items():
                if kw in Bn:
                    vec[Bn.index(kw)] = kv
            cols.append(vec)
        Rb = [list(v) for v in Rn]
        rr = frac_rank(cols + Rb, len(Bn)) if (cols or Rb) else 0
        e_n = rr - len(Rn)
        return rank[n] - e_n

    res["modLm1_irred"] = {n: modLm1_irred(n) for n in range(0, NMAX + 1)}
    res["modLm1_irred_seq_1_4"] = [res["modLm1_irred"][n] for n in (1, 2, 3, 4)]
    vm = {}
    for n in range(NMAX + 1):
        if n == 0:
            vm[0] = 1
            continue
        cols = []
        for w in B[n - 1]:
            d = mul_left(cval, hval, -1, {w: Fraction(1)})
            vec = [Fraction(0)] * len(B[n])
            for kw, kv in d.items():
                vec[B[n].index(kw)] = kv
            cols.append(vec)
        r = frac_rank(cols, len(B[n])) if cols else 0
        vm[n] = len(B[n]) - r
    res["modLm1_verma"] = vm
    res["modLm1_verma_seq_1_4"] = [vm[n] for n in (1, 2, 3, 4)]
    d = rank
    q_from_d = {n: (d[n] - d[n - 1] if n >= 1 else d[0]) for n in range(NMAX + 1)}
    res["q_from_rank_minus_rank_prev"] = q_from_d
    res["q_from_d_seq_1_4"] = [q_from_d[n] for n in (1, 2, 3, 4)]
    res["is_0011_arithmetic_restatement"] = (
        [q_from_d[n] for n in (1, 2, 3, 4)] == [res["modLm1_irred"][n] for n in (1, 2, 3, 4)])
    res["Lm1_injective_on_quotient"] = all(
        res["modLm1_irred"][n] == rank[n] - rank[n - 1] for n in range(1, NMAX + 1))
    return res


def main():
    out = {"marker": MARK, "nmax": NMAX,
           "route": "independent greedy insert/push normal-ordering recursion"}
    runs = {}
    runs["c0_h58"] = run("0", "5/8", "c=0,h=5/8 (thermal)")
    runs["generic_c12_h37"] = run("1/2", "3/7", "generic control (1/2,3/7)")
    runs["c0_h0"] = run("0", "0", "c=0 vacuum h=0")
    out["runs"] = runs
    with open("v768_v1_result.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


if __name__ == "__main__":
    main()
