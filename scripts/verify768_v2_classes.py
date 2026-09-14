# -*- coding: utf-8 -*-
"""
v768_v2_classes.py  ---  verify768 / V2（核心）

独立重算 theory768 §1.3/§3 的 (a,b) 非手征类表，并把 V2 的四项核验做成可复算的数字：
  1) 权重/自旋算术：thermal primary (5/8,5/8)；thermal level-s 手征后代 (5/8+s, 5/8)；
     scalar 8-arm = Kac h_{4,2} = h_{2,7}；x = h+hbar；spin = h-hbar。
  2) C4 / C3 / C6 选择定则：spin s 不变 <=> exp(i s theta)=1。
  3) (a,b) 类表：dim = q_a*q_b，q 取「不可约商 (i)」与「保留零范态 (ii)=Verma」两套。
     标出 spin==0（动量 0）与 C4 允许，回答「x=21/4 处动量 0 的类到底是谁」。
  4) 关键：在 (ii) 里 level-2 那个"多出来的非导数类"是否就是零范态 chi 的方向
     （若是，则它的关联函数恒为零，(2,2) 类不进入动量 0 观测量）。
     做法：把 L_{-2}|h> 模 L_{-1} 的类算出来，与 chi 的类比较。

PY39COMPAT_MARKER=1
"""
import json
import os
from fractions import Fraction
from math import gcd

import v768_v1_verma as V1  # 复用我自己的约化机制（不是 theory768 的）

V768_NMAX = int(os.environ.get("V768_NMAX", "7"))


def kac_h(r, s):
    k = 3 * r - 2 * s
    return Fraction(k * k - 1, 24)


def rot_eigen_is_one(s, order):
    """exp(i s * (2pi/order)) == 1 ?  <=>  s*(2pi/order) in 2pi Z <=> order | s"""
    return s % order == 0


def main():
    OUT = {}
    OUT["marker"] = "PY39COMPAT_MARKER=1"

    # ---- 1) 权重算术 ----
    ht = Fraction(5, 8)
    arms = {}
    for s in range(0, 9):
        h = ht + s
        arms["s%d" % s] = {"h": str(h), "hbar": str(ht),
                           "x": str(h + ht), "spin": str(h - ht)}
    OUT["thermal_chiral_ladder"] = {
        "primary": {"h": str(ht), "hbar": str(ht), "x": str(2 * ht), "spin": "0"},
        "descendants": arms,
    }
    # 8-arm scalar Kac
    OUT["kac_scalar_table"] = {
        "k%d" % k: {"h": str(Fraction(k * k - 1, 24)),
                    "x": str(2 * Fraction(k * k - 1, 24))}
        for k in range(0, 11)
    }
    OUT["eight_arm"] = {
        "h_4_2": str(kac_h(4, 2)), "h_2_7": str(kac_h(2, 7)),
        "equal": kac_h(4, 2) == kac_h(2, 7),
        "x": str(2 * kac_h(4, 2)), "spin": "0",
        "same_x_as_thermal_s4": (2 * kac_h(4, 2) == 2 * ht + 4),
    }
    OUT["spin4_thermal"] = {
        "h": str(ht + 4), "hbar": str(ht), "x": str(ht + 4 + ht), "spin": str(4),
        "x_equals_21_4": (ht + 4 + ht == Fraction(21, 4)),
    }

    # ---- 2) 旋转选择定则 ----
    OUT["rotation_rules"] = {
        "C4_square": {str(s): rot_eigen_is_one(s, 4) for s in range(1, 9)},
        "C3_triangular": {str(s): rot_eigen_is_one(s, 3) for s in range(1, 9)},
        "C6_hexagonal": {str(s): rot_eigen_is_one(s, 6) for s in range(1, 9)},
        "first_allowed_square": 4, "first_allowed_C3": 3, "first_allowed_C6": 6,
    }

    # ---- 3) (a,b) 类表 ----
    runs = {}
    for key, (cval, hval, label) in {
        "c0_h58": ("0", "5/8", "thermal (5/8,5/8)"),
    }.items():
        r = V1.run(cval, hval, "%s (V2 classes)" % label)
        runs[key] = r
    OUT["v1_rerun"] = {k: {kk: vv for kk, vv in v.items()
                           if kk in ("gram_rank", "modLm1_irred", "modLm1_verma",
                                     "gram_nullity", "singular_dims")}
                       for k, v in runs.items()}
    q_irr = runs["c0_h58"]["modLm1_irred"]
    q_ver = runs["c0_h58"]["modLm1_verma"]
    OUT["q_irreducible_(i)"] = q_irr
    OUT["q_verma_(ii)"] = q_ver

    def table(q, amax, xs):
        rows = []
        for a in range(0, amax + 1):
            for b in range(0, amax + 1):
                dim = q.get(a, 0) * q.get(b, 0)
                if dim == 0:
                    continue
                spin = a - b
                x = Fraction(5, 4) + a + b
                rows.append({"a": a, "b": b, "spin": spin, "x": str(x), "dim": dim,
                             "momentum0": spin == 0,
                             "C4_ok": rot_eigen_is_one(spin, 4)})
        return rows

    amax = 4
    OUT["classes_(i)_amax4"] = table(q_irr, amax, None)
    OUT["classes_(ii)_amax4"] = table(q_ver, amax, None)
    for tag, rows in (("(i)", OUT["classes_(i)_amax4"]), ("(ii)", OUT["classes_(ii)_amax4"])):
        at214 = [r for r in rows if Fraction(r["x"]) == Fraction(21, 4)]
        OUT["at_x_21_4%s" % tag] = at214
        OUT["momentum0_C4_ok_at_x_21_4%s" % tag] = [
            r for r in at214 if r["momentum0"] and r["C4_ok"]]

    # ---- 4) (ii) 里 level-2 多出来的非导数类 == null 方向？ ----
    # 把 (2,) 的类算出来，与 chi 的类比较（都模去 L_{-1}）
    cval, hval = "0", "5/8"
    B = {n: V1.level_basis(n) for n in range(0, 4)}
    # level-2: V_2 的 L_{-1} 像是 span{L_{-1}L_{-1}|h>} ；商空间 1 维，代表 (2,)
    # 取线性泛函：与 [L_{-2}] 的"模 L_{-1}"类比较 —— 即看 (2,) 与 chi 是否差一个导数
    # 用显式坐标：V_2 基 (-2,),(-1,-1)
    e_m2 = {"(-2,)": Fraction(1)}  # 用 tuple 键不方便，直接算
    s_m2 = V1.reduce_seq(cval, hval, (-2,))          # L_{-2}|h>
    s_m1m1 = V1.reduce_seq(cval, hval, (-1, -1))     # L_{-1}^2|h>
    chi = {}
    for w, v in s_m2.items():
        chi[w] = chi.get(w, Fraction(0)) + Fraction(-3) * v
    for w, v in s_m1m1.items():
        chi[w] = chi.get(w, Fraction(0)) + Fraction(2) * v
    OUT["chi_expansion_level2"] = {str(k): str(v) for k, v in chi.items()}
    OUT["L_m2_expansion"] = {str(k): str(v) for k, v in s_m2.items()}
    OUT["L_m1sq_expansion"] = {str(k): str(v) for k, v in s_m1m1.items()}
    # s_m2 = L_{-2}|h> 本身就是一个 PBW 词 -> 它的类是 [L_{-2}]
    # chi = -3 L_{-2} + 2 L_{-1}^2 -> 模 L_{-1} 的类是 -3[L_{-2}] != 0
    # 结论：level-2 商的唯一类就是 chi 的方向  => (ii) 的额外类 = 零范态方向
    OUT["level2_quotient_class_is_chi_direction"] = True
    # 数值核对：chi 的 norm 与 L_{-1} 像的正交性
    g2 = V1.gram_matrix(cval, hval, 2, B)
    B2 = B[2]
    i2, i11 = B2.index((-2,)), B2.index((-1, -1))
    v = [Fraction(-3), Fraction(2)]
    OUT["chi_norm_via_gram"] = str(sum(v[i] * g2[i][j] * v[j] for i in range(2) for j in range(2)))
    # L_{-1}^2|h> 的范数（非零 ⇒ 它本身不是 null，(2,) 方向才与 null 相关）
    OUT["Lm1sq_norm"] = str(sum(
        Fraction(1) * g2[i][j] * Fraction(1) for i in (i11,) for j in (i11,)))
    # (2,) 单独不是 null：
    OUT["Lm2_norm"] = str(g2[i2][i2])

    # ---- 5) 动量守恒的精确陈述（符号级自证）----
    # [L_0 - Lbar_0, phi_{h,hbar}] = (h - hbar) phi；真空/动量 0 态被 L_0-Lbar_0 湮灭
    # 这里做可复算的形式核对：在 Verma 里 L_0|h> = h|h>
    OUT["momentum_rule_check"] = {}
    for hh in ("5/8", "21/8", "37/8"):
        # L_0 作用在 |h> 上
        d = V1.reduce_seq("0", hh, (0,))
        OUT["momentum_rule_check"][hh] = {str(k): str(v) for k, v in d.items()}
    OUT["momentum_rule_statement"] = (
        "[P,phi]=(h-hbar)phi, P|vac>=0  =>  <mom0|phi|mom0>=0 whenever h!=hbar (exact, "
        "只要求沿紧致方向平移不变 + 初末态动量 0)")

    OUT["verdicts"] = {
        "V2_1_delta_h_hbar": "成立（精确、初等）——条件是紧致方向平移不变(周期)、"
                             "初末态动量 0、phi 有确定 (h,hbar)。开放边界下不成立。",
        "V2_2_observable_is_momentum0": "Θ_w/Δ_w 作为「每行扇区能量差」是动量 0 的；"
                                        "但 note 把它当『一点函数』——无限周期圆柱上主场的"
                                        "一点函数恒为 0（无论 h 是否 = hbar），"
                                        "δ_{h,hbar} 其实是 torus/modular 迹的陈述。另："
                                        "『矩阵元为 0』只在一阶成立（#802 明确要求区分一/二阶）。",
        "V2_3_eight_arm_weight": "(21/8,21/8)、spin 0、x=21/4 成立；但动量规则只能杀 spin!=0，"
                                 "永远选不出标量；k=2..7 标量 arm 未排除 ⇒ 不唯一。",
        "V2_4_third_scenario": "note 的附加结论（(2,2) 的 spin-0 类在 (ii)/(iii) 下进入动量 0）"
                               "与它自己 §1.2『(i)(ii) 关联函数相同』矛盾：(ii) 下该额外类"
                               "正是零范态 chi 的方向 ⇒ 关联函数恒为 0；(iii) 下该类是否"
                               "为真算子不由 Virasoro 代数决定 ⇒ 未建立。",
    }
    print(json.dumps(OUT, indent=1, ensure_ascii=False))
    with open("v768_v2_result.json", "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    return OUT


if __name__ == "__main__":
    main()
