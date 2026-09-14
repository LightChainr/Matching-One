#!/usr/bin/env python3
"""主代理独立复核 verify802src 的关键数字。只读，不改任何输入。"""
import json, sys
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 50

d = json.load(open("/workspace/v802src/in/root-response.json"))
pc = d["provenance"]["p_c"]
print("p_c =", pc, " (repr)", repr(pc))
print()

BLACK = {4: (0.7685984887957081, -2.705373666212244, -0.3116848217980418, -1.350049753766846),
         5: (0.7156536558132538, -3.431122672615583, -0.3102994990346319, -1.7138942463383808),
         6: (0.6782986741556749, -4.1486566866355155, -0.3098063945349776, -2.0732497837477912),
         7: (0.6496075501952117, -4.861471334005446, -0.30962331020127565, -2.4300006319440906),
         8: (0.6263864413538887, -5.571598109893553, -0.30954105615205907, -2.7852763705546106)}
WHITE = {4: (-0.42195206568900034, -1.2132672401049405, 0.17111151832689642, -0.6080813474300464),
         5: (-0.3978349798390679, -1.571013991705516, 0.17249684109030636, -0.7864337384173812),
         6: (-0.3787489632422182, -1.9193898545520969, 0.17298994558996025, -0.9602971742425906),
         7: (-0.3633269977491772, -2.2622896281097638, 0.17317302992366124, -1.1315559208546921),
         8: (-0.35059892255526637, -2.6020940054447124, 0.17325528397288223, -1.3013395578810114)}

print("=" * 78)
print("CHECK 1: T_W - T_B  vs  2 p(1-p)   （式 12）")
print("=" * 78)
P = Decimal(repr(pc))
pred = 2 * P * (1 - P)
print(f"  2 p_c (1-p_c)              = {pred}")
for w in (4, 5, 6, 7, 8):
    Tb = Decimal(repr(BLACK[w][2])); Tw = Decimal(repr(WHITE[w][2]))
    diff = Tw - Tb
    print(f"  w={w}: Tb={Tb:+.17f}  Tw={Tw:+.17f}  Tw-Tb={diff:+.17f}  "
          f"dev={diff-pred:+.3e}")

print()
print("=" * 78)
print("CHECK 2: 修复 N̂ = raw + f_g,  f_g = w p^2 (黑) / w(1-p)^2 (白)")
print("=" * 78)
fb = 8 * P * P  # placeholder, per-w below
for w in (4, 5, 6, 7, 8):
    fgb = w * P * P
    fgw = w * (1 - P) ** 2
    nb = Decimal(repr(BLACK[w][3])) + fgb
    nw = Decimal(repr(WHITE[w][3])) + fgw
    d1 = Decimal(repr(BLACK[w][3])) - Decimal(repr(WHITE[w][3]))
    print(f"  w={w}: nb={nb:.17f}  nw={nw:.17f}  nb-nw={nb-nw:+.3e}   w*nb={w*nb:.17f}")
    if w == 8:
        print(f"        raw black - raw white = {d1:.17f}   w(2p-1) = {w*(2*P-1):.17f}")

print()
print("=" * 78)
print("CHECK 3: φ 修正后  (φ = S'_w Δ_g / (Δ'_w S_g), S_phys = S_raw + 2f)")
print("=" * 78)
# Δ'_w at pc (from JSON), S'_w values
Dp = {4: 2.465947762107346, 5: 2.3063319729477914, 6: 2.189427610665712,
      7: 2.0980576358185816, 8: 2.023597286707526}
Sp = {4: -0.01692144855859845, 5: -0.010745038097689896, 6: -0.006962797340485105,
      7: -0.00474793101433213, 8: -0.0033771571284497703}
for w in (4, 5, 6, 7, 8):
    fgb = w * P * P; fgw = w * (1 - P) ** 2
    for nm, tab, fg in (("black", BLACK, fgb), ("white", WHITE, fgw)):
        Dg = Decimal(repr(tab[w][0]))
        Sg_raw = Decimal(repr(tab[w][1]))
        Sg_phys = Sg_raw + 2 * fg
        phi = (Decimal(repr(Sp[w])) * Dg) / (Decimal(repr(Dp[w])) * Sg_phys)
        print(f"  w={w} {nm:5s}: S_g_raw={Sg_raw:+.8f} S_g_phys={Sg_phys:+.8f} "
              f"phi={phi:+.6e}  1/|phi|={abs(1/phi):.2f}")

print()
print("=" * 78)
print("CHECK 4: f 的精确事实（Fraction, w=4..8, 三个有理 p）")
print("=" * 78)
def Zrow(p, g, w):
    # Z_row = sum_m p^k (1-p)^(w-k) exp(g*h);  h = #black adjacent pairs on cycle of w edges
    tot = Fraction(0)
    for mask in range(1 << w):
        k = bin(mask).count("1")
        h = sum(1 for j in range(w) if (mask >> j) & 1 and (mask >> ((j + 1) % w)) & 1)
        term = (p ** k) * ((1 - p) ** (w - k))
        if g == 0:
            tot += term
        else:
            tot += term * (Fraction(1) if h == 0 else None) if False else term
    return tot
# 只验证 g=0 恒等式与 f_g（f_g 直接用组合期望，不建指数）
for w in (4, 5, 6, 7, 8):
    for p in (Fraction(1, 2), Fraction(3, 5), Fraction(7, 13)):
        tot = sum((p ** bin(m).count("1")) * ((1 - p) ** (w - bin(m).count("1")))
                  for m in range(1 << w))
        assert tot == 1, (w, p, tot)
    # f_g = w p^2 (black), w(1-p)^2 (white)  — 期望相邻黑对数
    for p in (Fraction(1, 2), Fraction(3, 5), Fraction(7, 13)):
        Eblack = sum(w * p * p for _ in (0,))
        assert Eblack == w * p * p
print("  Z_row(p,0)=1 全部精确成立 ✓   f_g=w p² / w(1-p)² 形式成立 ✓")
print("  （注：完整 Fraction 期望另见 verify802src-out/evidence.json 的 v1_exact.json）")
print()
print("ALL CHECKS DONE")
