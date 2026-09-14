# lit-note-C.md —— 长程前沿文献定标（条目 1–8，完整版）

**作者**：lit-frontier-longb（长程检索专员）　**日期**：2026-09-14
**范围**：只做检索与写作。**未做任何计算、未做任何仓库写操作**（GitHub 只读）。不需要云机。
**继承**：`lit-note-B.md`（R/J/U）、`lit-note-B2.md`（1a/2/3/J(d)/5）、`round45-out/arms-audit.md`。
本文只补**新条目**与本轮**新发现**；冲突处**显式标出**，不重做已覆盖部分。

**格式纪律**：每条给可检索标识（arXiv id / DOI / 期刊卷页）；找不到写「未找到，已检索关键词：…」。
**加引号 = 直接引用；否则 = 转述**。「本项目内部笔记」**不作为文献证据**。

---

## 0. 总表

| 条目 | 已有结果 | 出处（可检索标识） | 关系 | 证据强度 |
|---|---|---|---|---|
| **1** square-site matching function `M_L` 的有限尺寸修正指数 | **有两份已发表测量，结论相反**。①根位移指数 `Δ₁ = 4.000 1(2)`（= 本项目「根指数 4」）②`M_L(p_c)~L^{2−x}` → `2−x = −3.42`；`p_L^⋆−p_c` 斜率 `−4.07`（推得 `w = −4.17`） | ① **Jacobsen 2015**, *J. Phys. A* **48** 454003, arXiv:`1507.03027`（§7.1，式 (36)(40)）② **Mertens & Ziff 2016**, *Phys. Rev. E* **94** 062152, arXiv:`1603.07289`（式 (38)(39)，图 6/7） | **①支持 `x=21/4 ⇒ h=21/8`；②冲突** | ① 直接数值（n ≤ 21，TM 特征值法）＋CFT 佐证；② 直接数值（L=3–7 精确 + 16/24/32/48 MC，作者自认系统太小） |
| **1′** `alpha_j=(j²−1)/12` 的出处 | **已发表**：`x_ℓ^P = (ℓ²−1)/12`（ℓ = 非重叠 path/臂数），显式「extends to half-integers the Saleur–Duplantier exponents」 | **Aizenman, Duplantier, Aharony 1999**, *Phys. Rev. Lett.* **83** 1359–1362, arXiv:`cond-mat/9901018`（式 (1)(2)） | **已蕴含**（arms-audit §1.3 的裁定获原文确认） | 已证（PRL） |
| **2** 判别臂数的 `c_slope` 渐近先例 | **未见完全对应量**。最接近的方法学先例 = **对数导数 / 导数比 FSS 估计子**（其极值随 `L^{1/ν}` 标度）＋「有效比值函数」逼近极限 | 对数导数 FSS：U(1) gauge–Higgs 双模拟论文（*Nucl. Phys. B* 通篇；「logarithmic derivatives of moments … their maxima also scale with」）；Ferrero et al. arXiv:`0803.3339`（`s^{1/ν} = 1 + x ∂_x F_ξ`）；有效比值函数：Shchur–Berche–Butera, *Nucl. Phys. B* **811** (2009) 491–518, DOI `10.1016/j.nuclphysb.2008.10.024` | **未见同口径量；方法学相邻** | 方法学先例充分，**对象不同** |
| **3** 环面按**同调秩**分辨的占据计数 `C[L,j,k]` | 已发表的拓扑分类是**绕行方向型**（`Z_0/Z_1/Z_2`），**不是**占据子图的**同调秩**；同调分辨**存在**但只在**高维** | `Z_·`：**Scullard–Jacobsen / arXiv:`2010.02887`**（式 (6) 上下文）；同调分辨：**Duncan–Kahle–Schweinhart**, arXiv:`2011.11903` (math.PR)；环面精确多项式到 L≤12：**Akhunzhanov–Eserkepov–Tarasevich**, arXiv:`2204.01517`, J. Phys. A **55** 204004, DOI `10.1088/1751-8121/ac61b8` | **口径相邻（绕行型 vs 同调秩）、对象不同** | 已发表，但**不是**本项目那张双下标整数表 |
| **4** `√(cE)·Z` ⇒ Laplace（四阶矩 6） | **命名恒等式**：**Gauss–Laplace transmutation**；标准表述「Gaussian × exponential-variance mixture = Laplace」 | **Ding & Blitzstein**, arXiv:`1510.08765`（§1）；更早 **Andrews & Mallows 1974**、**West 1987**；应用 **Park & Casella 2008**（Bayesian Lasso） | **已蕴含（恒等式，无条件）** | 已证（命名恒等式） |
| **5** 稀有大簇的稳定 CLT / 随机区间加性泛函 | **框架齐备**（marked point process CLT/stable、stabilizing functionals、Gibbs 条件原理、Palm 渐近正态），但**「随机（指数长）区间上的加性泛函 + 完整 component-Palm 锚点」的具体陈述未见** | Penrose, arXiv:`math/0410021`（stabilizing 泛函 CLT，**应用含「critical percolation 簇计数」的白噪极限**）；Basrak–Wintenberger–Zugeč, arXiv:`1903.09387`（marked Poisson cluster：CLT 或**无穷方差稳定律**）；Onaran–Bobrowski–Robert, *Electron. J. Probab.*（动态点过程**局部加性泛函 CLT**）；Kipnis–Varadhan 1986；Ferré–Hervé–Ledoux 2012（arXiv:`1201.4579`）；Prokešová–Jensen（Palm 似然渐近正态，DOI `10.1007/s10463-012-0376-7`） | **技术空白（非原理空白）** | 框架 = 已证；具体陈述 = 未见 |
| **6** 孔洞场 / cavity field / `9/8` 极点 | (a) 「纵向步长至多 1」长步机构**未见**；相邻=**长程渗流**（power-law 长边，非「步长 ≤1」）。(b) **孔洞/内边界**在渗流文献**有**，且**孔洞 ↔ 互补（白）簇的对偶**是标准的。(c) **fugacity-倾斜的活动矩阵 + 谱端点 → 极点**有样板：**R-矩阵 / Yang–Lee 零点**，`ρ(z)=Σ_λ ψ_1(λ)²/(z^{−1}+λ)`，奇异点 = 谱端点。(d) `1/9` / `9/8` **未见** | (b) **Isichenko**, *Rev. Mod. Phys.* **64** (1992) 961–1043（内 hull ↔ 互补簇）；Hu 等的 hole 幂律 + 「largest hole」分布（*Physica A* 2021, DOI `10.1016/j.physa.2021.125847`）(c) **arXiv:`0907.4037`**（R-矩阵谱倒数 = Yang–Lee 零点）；Collatz–Wielandt / Perron–Frobenius（Meyer, *Matrix Analysis*, ch. 8） | (b)(c) **方法学已蕴含**；(a)(d) **未见** | (b)(c) 已发表；(a)(d) 空白 |
| **7** first-exit 盒与向量指数矩域 | **未见**「向量指数矩域给方向范数内外界」的现成结果。可引用的样板：**Collatz–Wielandt 显式向量界**（非负矩阵 Perron 根的双侧夹逼）＋「inclusion interval」的严格理论；first-exit 尾律的标准工具 = Freidlin–Wentzell LDP | **「Upper bounds on the growth rates of hard squares and related models」**, *DMTCS*（Lemma 2 = Collatz–Wielandt：`min_i (Ax)_i/x_i ≤ λ ≤ max_i (Ax)_i/x_i`）；**Oepomo**, *Electron. J. Linear Algebra* **10** (2003) 31–45（Zbl 1022.15018，Collatz 特征值 inclusion interval）；**Lifshits–Shi**, *Bernoulli* **8** (2002) 745–765（Zbl 1018.60084，first exit + LDP 尾律）；Freidlin–Wentzell | **构件齐备，组合未见** | 构件已证；组合 = 未见 |
| **8** `#757` 投稿先例矩阵（低优先） | 几何/渗流**精确判据**文献充分：临界多项式（任意周期 basis、任意 2D 周期格）、`P_B(p)=R_2−R_0`、`P(A,B,C)=P(Ā,B̄,C̄)`、site 渗流经 covering lattice 归约；且**「full scaling law」被该领域明确列为 OPEN** | **Scullard–Jacobsen**, arXiv:`1209.1451`（J. Phys. A **45** 494004）；arXiv:`1207.3340`（deletion–contraction）；**is a graph invariant**：PRE 2012, PubMed `23214553`；**Jacobsen–Scullard 2019/2020**, arXiv:`1910.12376`（摘要明列 open question「the full scaling law」） | 「任意周期 / 精确判据」**已有**；「两次出生 / balance root / matching root」**未找到** | 部分已发表；术语未见 |

---

## 1. 条目 1（最高优先）：`h=21/8` 的独立数值证据**存在，但是「一支持一冲突」**

### 1.1 被引文献自己的数（**已逐式核对，直接引用**）

**Mertens & Ziff 2016**（*Percolation in finite matching lattices*, **Phys. Rev. E 94, 062152**, DOI `10.1103/PhysRevE.94.062152`, arXiv:`1603.07289`；下称 **MZ16**）。

MZ16 把 Sykes–Essam 关系做**有限尺寸推广**，定义 **matching function**
```
M_L(p) = N_L(p) − N̂_L(1−p) − L²χ(p)                      （MZ16 式 (15)）
```
并证明它**精确等于**若干**绕行概率之差**（式 (20)）：
```
M_L(p) = R_L^x(p) − R̂_L^x(1−p),   x ∈ {c(交叉绕行), b(两向), e(任一方向), h(横向)}
```
（式 (17) 特别给出 `M_L(p) = R_L^b(p) − R̂_L^b(1−p)`。）

**关键直接引用（MZ16 式 (38)(39) 及其后正文）**：
> `M_L(p) = A₂ L^{2−x} + 2B₁ b L^{1/ν}(p−p_c) + C₂ L^{2−y}(p−p_c)² + 2D₁b³L^{3/ν}(p−p_c)³ + …`　(38)
> 「That is, `M_L(p_c)=A₂L^{2−x}`, … In Fig. 6, using exact and Monte-Carlo data, we plot these quantities vs. `L` on a log-log plot. **These plots give `2−x = −3.42`** and `2−y = 0.705`.」
> `p_L^⋆ − p_c ∼ L^{2−x−1/ν}`　(39)
> 「The numerical value for `x` implies that this exponent has the value **`w = 2−x−1/ν = −3.42−3/4 = −4.17`**, somewhat larger than the value `4` suggested by Jacobsen [Jacobsen 2015]. In Fig. 7 we show the results for `p_L^⋆−p_c` … and find a **slope of `−4.07`** for this criterion.」
> 「If we assume that `w` is exactly `−4`, then **`x−2 = 3.25` exactly**.」

**对照 arms-audit 的链条**（`x = 2 − y_t + |根指数|`，`y_t = 3/4`；`omega = x − 2`）：

| 输入 | `x` | `h = x/2` | `24h+1` | 完全平方? |
|---|---|---|---|---|
| MZ16 实测 `2−x = −3.42` | 5.42 | 2.71 | 66.04 | ❌ |
| MZ16 实测斜率 `−4.07` | 5.32 | 2.66 | 64.84 | ❌ |
| **`21/4`（目标）** | 5.25 | **21/8** | **64** | ✅ |

→ **arms-audit §3 的数逐字复核无误**。MZ16 自己写「Presumably, larger systems are needed to find the true behavior」——**他们没有宣称 `−3.42/−4.07` 是渐近值**。

### 1.2 ⚠️ 本轮新发现（arms-audit 未覆盖）：**存在一份同模型、支持 `4` 的独立数值测量**

**Jacobsen 2015**（*Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley–Lieb algebras*, **J. Phys. A: Math. Theor. 48 (2015) 454003**, arXiv:`1507.03027`；下称 **J15**）。

J15 **§7.1 标题即 "Site percolation on the square lattice"**（= 本项目**完全相同的模型**；同为 §7.2 的才是 kagome bond）。**直接引用**：
> `p_c(n) = p_c + Σ_{k≥1} A_k n^{−Δ_k},  with 0 < Δ_1 < Δ_2 < ⋯`　(34)
> 有效指数由逐差对数导数给出：`Δ_1(n) = [log δp_c(n) − log δp_c(n−1)] / [log n − log(n−1)]`　(35)
> **`Δ_1 = 4.000 1(2).`**　(36)　「This agrees well with the value `w = 4.03 ± 0.01` … **It appears inevitable to admit that `Δ_1 = 4` exactly.**」
> `Δ_2 = 6.00(1)`　(39)　「we henceforth admit that `Δ_2 = 6` exactly.」；`Δ_3 = 8.0(5)`，「we conjecture that `Δ_3 = 8`」。
> `p_c(n) = p_c + Σ_{k≥1} A_k / n^{2(k+1)}`　(40)（即 **`Δ_k = 2(k+1)`**）
> 规模 **n ≤ 21**；`p_c = 0.592 746 050 792 10(2)`。

**CFT 依据（直接引用，§8）**：
> 「In the continuum limit there is no difference between whether the propagating cluster is an FK cluster or a dual FK cluster. Therefore `f_open(n)` and `f_closed(n)` both determine the **same critical exponent, namely `x_m`**, and they both scale like `f_1(n)` in (47).」
> 「The values `c = 0` and **`x_m = 5/48`** are of course known」
> `p_c(n) − p_c = O(n^{−4})`　(50)　「and moreover the corrections appear to be `O(n^{−6})`, `O(n^{−8})`, and so on.」

**为什么关键**：J15 的 `Δ_1` 与 MZ16 的 `w` **是同一个物理指数**（临界多项式/匹配函数的根位移指数）。而
```
w = x − 2 + 1/ν
J15:  w = 4     ⟹ x = 4 + 2 − 3/4 = 5.25 = 21/4 ⟹ h = 21/8  ✅
MZ16: w = −4.17 ⟹ x = 5.42 ≠ 21/4               ⟹ 21/8 不成立 ❌
```
MZ16 自己承认这个冲突：「somewhat larger than the value 4 suggested by Jacobsen」。

### 1.3 判定（条目 1）

> **问题**：有没有任何**已发表的、独立的**对 square-site matching polynomial（或 `M_L=P2−P0` / matching function）**有限尺寸修正指数**的测量或证明？

**答：有，且有两份，结论相反。**

| # | 被测量 | 数值 | 出处 | 对 `21/8` 的含义 |
|---|---|---|---|---|
| A | 根位移指数 `w`（= 项目「根指数」） | **4.000 1(2)**（n ≤ 21，TM 特征值法，**方形格点 site 渗流**） | **J15** §7.1 式 (36)(40)，arXiv:`1507.03027` | **支持** `x=21/4 ⇒ h=21/8` |
| B | 直接测 `M_L(p_c)~L^{2−x}` 的 `2−x` | **−3.42**（L=3–7 精确 + 16/24/32/48 MC） | **MZ16** 式 (38)，图 6 | **冲突**（`21/4` 要求 `−3.25`） |
| C | `p_L^⋆−p_c` 拟合斜率 | **−4.07**（同 B 数据）；由 B 推 `w=−4.17` | **MZ16** 式 (39)，图 7 | **冲突**（`21/4` 要求 `−4`） |

**三条并存事实**：
1. **A 与 C 是同一个指数**（都等于根位移/阈值估计子收敛指数），**A 说 4，C 说 −4.07**。二者不可能同时渐近成立。
2. **没有任何第二份**对 `M_L(p_c)` 本身指数（`2−x`）的独立测量 —— 只有 MZ16 的 **B**。所以 `−3.25`（`21/4` 的直接后果）**从未被直接测过**。
3. **两测法内核不同**：J15 用**转移矩阵特征值相等**（半无限圆柱周长 n ≤ 21，**非 Monte Carlo**）；MZ16 用**环面**精确枚举 + MC（L ≤ 48）。**几何不同**（圆柱 vs 环面），按普遍性该指数应相同。

**诚实结论**：`21/8` 这条链**不是**「无独立证据」，而是「**有一份支持（J15 `Δ_1=4`，同模型、精度高，n ≤ 21）、一份冲突（MZ16 `2−x=−3.42` / 斜率 `−4.07`，系统小、作者自认不够）**」。**判定生死需第三方**：
- 用 J15 的特征值法把 `p_c(n)−p_c` 的 `Δ_1` 推到 n > 21，看是否稳定在 4；
- 或用 MZ16 的框架把 `M_L(p_c)` 推到 L ≳ 100，看 `2−x` 趋向 `−3.25` 还是 `−3.42`。

**⚠️ 对 arms-audit §3 措辞的最小修正**：不要写「被引文献自己的数字让 `21/8` 消失」，应写「**MZ16 的小尺寸数据与 `21/4` 不符，但同一指数的独立测量 J15（同模型、n≤21）给出 `Δ_1=4.0001(2)`，与 `21/4` 一致；两者冲突，`21/8` 的生死未被这两篇判定**」。

### 1.4 附带确认：`alpha_j = (j²−1)/12` 是**已发表公式**

**Aizenman, Duplantier, Aharony 1999**（**Phys. Rev. Lett. 83 (1999) 1359–1362**, arXiv:`cond-mat/9901018`）**直接引用**摘要与式 (1)(2)：
> 「2D Percolation path exponents `x_ℓ^P` describe probabilities for traversals of annuli by **`ℓ` non-overlapping paths** … whose exponents, believed to be exact, yield **`x_ℓ^P = (ℓ²−1)/12`**. **This extends to half-integers the Saleur–Duplantier exponents for `k = ℓ/2` clusters**, yields the exact fractal dimension of the external cluster perimeter, `D_EP = 2 − x_3^P = 4/3`…」
> 式 (1)：`x_k^C = x_{ℓ=2k}^{O(N=1)} = (4k²−1)/12`（k = 簇数，ℓ = 2k = 线/臂数）；式 (2)：`D_H = 2 − x_1^C = 7/4`。

**含义**：`alpha_j=(j²−1)/12` **不是本项目自造**，是 **ADA 1999 / Saleur–Duplantier** 的已发表 path/arm 指数公式（`j` = 臂/线数）。与 arms-audit §1.3 的裁定完全一致：它与 `c=0` spinless Kac `x_bulk(k)=(k²−1)/12` **是同一个二次族**，故「8-arm ⇒ `x=21/4` ⇒ `h=21/8`」是**恒等式而非独立证据**。**arms-audit 该裁定成立，已由原始文献逐字确认。**

---

## 2. 条目 2：判别臂数的 `c_slope` 渐近先例

**待定标量**：`c_slope = h1 = (log c)_z`，其中 `c = C_L(b)` 是自归一化曲线，`b = (1/2) log(P0/P2)`。

### 2.1 最接近的方法学先例（**三族**）

**(i) 对数导数 / 导数比的 FSS 估计子（最标准的一族）**
FSS 里常用「可观测量对控制参数的对数导数」作 `1/ν` 的估计子；其**极值随 `L^{1/ν}` 标度**。直接引用一例（U(1) gauge–Higgs 双模拟，θ = π 临界端点）：
> 「we also determine ν by studying observables that have the same scaling behavior as U, for instance the **logarithmic derivatives of moments** of the topological charge. In particular we study the derivatives (44) **which have maxima that also scale with** `L^{1/ν}`. Thus these derivatives can again be fit as described in (43) and allow for an **independent determination of ν**.」
同族的显式 FSS 函数导数式（Ferrero et al., 4D Ising spin glass, arXiv:`0803.3339`）：
> `s^{1/ν} = 1 + x ∂_x F_ξ(x, s)|_{x=x_ξ(L^{-ω})}`；以及 `s^{1/ν} = 1 + g ∂_g F_g(g, s)|_{g=g(L^{-ω})}`

**(ii) 「有效比值函数」（effective ratio functions）逼近极限**
**Shchur–Berche–Butera**, *Numerical revision of the universal amplitude ratios for the two-dimensional 4-state Potts model*, **Nucl. Phys. B 811 (2009) 491–518**, DOI `10.1016/j.nuclphysb.2008.10.024`。直接引用摘要：
> 「we estimate ratios of critical amplitudes, constructing **effective ratio functions**, and computing their **limiting values at the critical point**.」
→ 即「构造一个 L-依赖的有效比值，再研究它随 L 的**逼近速率**」，与 `c_slope` 的「随尺寸的渐近」在**方法学上同型**（但对象是临界振幅比，非自归一化曲线）。

**(iii) 有限尺寸修正系数的普适比（Izmailian–Hu 型）**
**Okabe–Kawashima**, *Universal relations in the finite-size correction terms of two-dimensional Ising models*, arXiv:`cond-mat/0107514`（直接引用，转述其引用的 Izmailian–Hu, PRL **86** (2001) 5160）：
> `N(f_N − f_∞) = Σ_{k≥1} a_k/N^{2k−1}`，`ξ_N^{−1} = Σ_{k≥1} b_k/N^{2k−1}`，且 **`b_k/a_k = (2^{2k}−1)/(2^{2k−1}−1)`** 普适；`a_1 = cπ/6`、`b_1 = 2π x_H`。
→ 「两个级数系数之比是普适的」，是「比值型观测量」的**精确**样板。

### 2.2 判定（条目 2）

**未找到完全对应的已发表量**：「**自归一化曲线对平衡变量** `b=(1/2)log(P0/P2)` **的对数导数 `(log c)_z` 随尺寸的渐近**」——
**已检索关键词：** `self-normalized observable`, `log-derivative of balance observable`, `finite-size asymptotics of derivative ratios`, `ratio observable finite-size scaling`, `effective ratio function critical amplitude`, `universal amplitude ratio finite-size correction coefficients`, `logarithmic derivative finite-size scaling estimator`, `Binder cumulant derivative pseudo-critical scaling`, `finite-size scaling derivative of FSS function correlation length exponent`, `modulus/plateau approach exponent`.
**最接近的先例 = (i) 对数导数/导数比 FSS 估计子**（其极值随 `L^{1/ν}` 标度；出处见 §2.1(i)），以及 **(ii) 有效比值函数**。**投稿时建议**：把 `c_slope` 明确定位为「(i) 的一例」，并说明与 (i) 的差别在于**对数导数取自平衡变量 `b` 而非控制参数**。

---

## 3. 条目 3：环面按**同调秩**分辨的占据计数是否已有人发表

### 3.1 已发表的是**「绕行方向型」三分**，不是**「同调秩」**

**直接引用**（**arXiv:`2010.02887`**, *Critical polynomials in the nonplanar and continuum percolation models*）：
> 「All the configurations `{C}` on the torus are classified into three types as `{Z_0}`, `{Z_1}`, and `{Z_2}` according to their topological properties. … a configuration `C` belongs to `{Z_2}` if it **wraps along two different directions**, to `{Z_1}` if it **wraps along one and only one direction**, and to `{Z_0}` if it **does not wrap**. `R_2, R_1, R_0` … generally the critical polynomial is defined as `P_B ≡ R_2 − R_0`.」

同一分类亦见 **MZ16**（`R_L^e / R_L^h / R_L^s(螺旋) / R_L^b / R_L^1 / R_L^c(交叉绕行)`，arXiv:`1603.07289` §II 条目列表）。

**口径差（关键）**：`Z_0/Z_1/Z_2` 判的是**配置中是否存在某个簇的绕行**（winding 非零），即**存在性/方向**；而 `#775` 的 `j = r_black(ω)` 是**占据子图的环境同调秩** `r = dim H_1 ∈ {0,1,2}`。二者**不同**：
- 「横向与纵向由**不同簇**分别绕行」⇒ `Z_2` 成立，但**无单簇 cross-wrap**；
- `r = 2` 需两个**独立**非可缩圈，可由**两个不同簇**提供。
- 故 **`R_j`（绕行型概率）≠ `Pr[r_black = j]`（同调秩分布）**，不能把 `Z_·` 计数当 `C[L,j,k]`。

### 3.2 环面**同调分辨**研究**存在**，但对象是**高维**

- **Duncan, Kahle, Schweinhart**, *Homological percolation on a torus: plaquettes and permutohedra*, **arXiv:`2011.11903`** (math.PR, v1 2020-11-24, v4 2023-09-29)。**直接引用摘要**：
  > 「We study higher-dimensional homological analogues of bond percolation on a square lattice and site percolation on a triangular lattice. … finite cell complexes … with the topology of the torus `T^d`. When random subcomplexes induce nontrivial `i`-dimensional cycles in the homology of the ambient torus, we call such cycles **giant**. We show that for every `i` and `d` there is a sharp transition from nonexistence of giant cycles to giant cycles spanning the homology of the torus. … we prove that `p_c = 1/2` in the case of middle dimension `i = d/2` for both models. This gives finite-volume high-dimensional analogues of Kesten's theorems…」
  → **同调分辨**（giant cycle = 非平凡同调类）**是**已发表概念，且就在**环面**上；但为**高维 plaquette / permutohedral**，2D 退化为经典 Kesten。**未给**「按 `r_black ∈{0,1,2}` 与 `|ω|=k` 双分辨的精确整数表**」。

### 3.3 判定（条目 3）

> **问题**：是否有人发表过环面上按 cycle rank / homology rank 分辨的占据子集计数？

**答**：
- 「环面 + 同调分辨」研究**存在**（Duncan–Kahle–Schweinhart, arXiv:`2011.11903`），但为**高维同调渗流的相变/阈值**，**不是** 2D square-site 的 `C[L,j,k]`。
- 标准**「拓扑三分」**（`Z_0/Z_1/Z_2`）是**绕行方向型**，**口径不同于同调秩**；`#775` 要的按 `r_black`（`0/1/2`）双下标整数表——**未找到已发表版本**。
  **已检索关键词：** `rank generating polynomial`, `Tutte polynomial torus`, `cycle rank distribution`, `homology resolved counting`, `Aizenman Duplantier Aharony`, `Akhunzhanov Eserkepov Tarasevich`, `wrapping polynomial`, `topological sector counting`, `Betti number distribution subgraph`, `random subgraph homology`, `critical polynomial Z_0 Z_1 Z_2 torus`, `homological percolation torus giant cycles`, `exact percolation probabilities torus cylinder plane polynomial`, `simplicial homology random configuration torus`, `cycle rank distribution random subgraph`.
- **`#752` 的 wrapping polynomial 覆盖哪一部分**：`R_2 − R_0`（= 临界多项式 `P_B`）与 `R_j`（`j=0,1,2`，**绕行型**）已被 **Scullard–Jacobsen / arXiv:`2010.02887`** 与 **Akhunzhanov–Eserkepov–Tarasevich**（*Exact percolation probabilities for a square lattice: Site percolation on a plane, cylinder, and torus*, **arXiv:`2204.01517`**, J. Phys. A **55** (2022) 204004, DOI `10.1088/1751-8121/ac61b8`；**torus L ≤ 12 的精确多项式**）覆盖。**但它们分的是「绕行方向型」，不是「同调秩」**。若 `#752` 指这套，则**不覆盖** `C[L,j,k]`，**口径差 = 绕行型 vs 同调秩**。

---

## 4. 条目 4：`normal × exponential` 混合 = Laplace —— **命名恒等式**

### 4.1 标准出处（**直接引用**）

**Ding & Blitzstein**, *Representation for the Gauss–Laplace Transmutation*, **arXiv:`1510.08765`**（§1）：
> 「**The Gauss–Laplace transmutation** states that  `V ∼ 2Exp(1),  L|V ∼ N(0,V)  ⟹  L ∼ Laplace`,  or equivalently, if `Exp ∼ Exp(1)` is independent of `Z ∼ N(0,1)`, then **`L = √(2Exp)·Z ∼ Laplace`**.」
> 「Some proofs of the Gauss–Laplace transmutation exist in the literature (**Andrews and Mallows, 1974; West, 1987**)…」
> 「…crucial to efficiently simulate posterior distribution of the **Bayesian Lasso (Park and Casella, 2008)**, which imposes **Laplace priors**…」

同义表述（转述，佐证非唯一出处）：「compounding a Gaussian with **exponential** variance (or Rayleigh standard deviation) yields a **Laplace**」；属 **variance-gamma 族** gamma 形状参数 = 1 的特例。

### 4.2 与本项目写法的对应（**含 `c` 的版本**）

`√(cE)·Z = √(c/2)·√(2E)Z = √(c/2)·L`，其中 `L` 为标准 Laplace。取 Laplace scale `λ = √(c/2)`，密度 `(1/(2λ))e^{−|y|/λ}` 即 **`(1/√(2c))e^{−√(2/c)|y|}`** —— 与本项目写法**逐字一致**。`c = 2` 即标准 Laplace（方差 2，四阶矩 24 ⇒ **标准化四阶矩 24/2² = 6** ✅）。

### 4.3 判定

**是已命名恒等式**：**Gauss–Laplace transmutation**（normal–exponential scale mixture）。**无需额外条件**，只要 `E ⊥ Z` 且 `E ∼ Exp`；换 `Exp(λ)` 只改 `c` 的标定。
可检索标识：**arXiv:`1510.08765`**；渊源 **Andrews & Mallows 1974**、**West 1987**；应用 **Park & Casella 2008**。
→ **此条不必作为新结果证明，引用即可。**

---

## 5. 条目 5：稀有大簇的稳定 CLT / 随机区间上的加性泛函

**待定标物**（`#772` 路线 A）：supercritical site 的有限半径局部近似、**随机区间上的加性泛函**、rare barrier 的 marked point process、**稳定 CLT / 随机信息 LAN**；特别地「**完整 component-Palm 的锚点选择不能被无条件块 CLT 自动覆盖**」。

### 5.1 已有的现成框架（**逐条给标识**）

| 需要的构件 | 已有结果 | 出处 |
|---|---|---|
| 可加泛函 CLT / 函数 CLT | 可逆 Markov 链上可加泛函 CLT + FCLT | **Kipnis–Varadhan 1986**, *Comm. Math. Phys.* **104** 1–19；**Ferré–Hervé–Ledoux 2012**, *Ann. I.H.P. B* **48**(2) 396–423, arXiv:`1201.4579` |
| **随机几何上的加性泛函 + CLT（含渗流应用）** | 「stabilizing 泛函」的 LLN/CLT，**应用明列「critical percolation 簇计数」的白噪极限** | **Penrose**, *Multivariate spatial central limit theorems with applications to percolation and spatial graphs*, arXiv:`math/0410021` |
| **marked point process 的 CLT 与稳定 CLT** | 「we find sufficient conditions under which the total claim amount satisfies the **central limit theorem** or alternatively tends in distribution to an **infinite variance stable random variable**」 | **Basrak–Wintenberger–Zugeč**, *On total claim amount for marked Poisson cluster models*, arXiv:`1903.09387` |
| **（动态）点过程局部加性泛函的有限维 CLT** | 「finite-dimensional central limit theorems for **local, additive, interaction functions** of temporally evolving point processes … via a distributionally equivalent **marked point process**」 | **Onaran–Bobrowski–Robert**, *CLTs for Local Functionals of Dynamic Point Processes*, *Electron. J. Probab.*（ISSN 1083-6489） |
| **随机区间上的加性泛函（古典）** | 随机区间装箱：`N(x)` 的 LLN + **Dvoretzky–Robbins CLT**（RSA / random interval packing） | 见 **arXiv:`1311.4967`** 引言对 Dvoretzky–Robbins 的转述引用 |
| **条件（rare-event）下的极限** | Gibbs 条件原理 → 指数倾斜的 Markov 过程（driven process） | **Csiszár 1984 / van Campenhout–Cover 1981 / Dembo–Zeitouni**；**Chetrite–Touchette**, arXiv:`1405.5157`, *Ann. Henri Poincaré* **16** (2015)（亦见 `lit-note-B2.md` §4） |
| **Palm 条件下的渐近正态 / 非退化性** | Palm 似然估计子的**强相合 + 渐近正态**（Neyman–Scott、log-Gaussian Cox）；Slivnyak 区分 reduced/non-reduced Palm | **Prokešová–Jensen**, *Asymptotic Palm likelihood theory for stationary point processes*, *Ann. Inst. Statist. Math.* (2013), DOI `10.1007/s10463-012-0376-7`；Palm 教程 arXiv:`1512.05871` |
| **临界点上的非高斯「分形 CLT」** | 「**Fractal Central Limit Theorem**」(long-range correlations) holds at the unstable, critical fixed point | *Stochastic renormalization group in percolation: I*, *Physica A*（ScienceDirect S0378437102012128） |

### 5.2 判定（条目 5）

- **框架齐备**：随机几何上的加性泛函 CLT 有 **Penrose arXiv:`math/0410021`**（且**直接应用在渗流的簇计数**上）；marked point process 的 CLT/**稳定律**有 **Basrak–Wintenberger–Zugeč arXiv:`1903.09387`**；rare-event 条件化有 **Gibbs 条件原理 / Chetrite–Touchette**；Palm 条件下的渐近理论与**非退化条件**（Fisher 信息 / 二阶矩）有 **Prokešová–Jensen**。
- **未见**：「**随机（指数长度）区间上的加性泛函**」＋「**完整 component-Palm 的锚点选择**」的**具体联合陈述**（即「锚点选择」这一步在文献里通常由 **Slivnyak/Campbell 公式**处理，而**无条件块 CLT 不给锚点条件分布**——这与 spec 的判断一致）。
  **已检索关键词：** `additive functional random interval`, `marked point process CLT`, `stable CLT rare event`, `random information LAN`, `exploration process Poisson cluster`, `Fisher projection nondegeneracy`, `conditional CLT under conditioning`, `Palm conditioning`, `stabilizing functional CLT percolation`, `Palm likelihood asymptotic normality`, `fractal central limit theorem percolation`, `exponential length interval additive functional`.
- **裁定**：**技术空白（非原理空白）**。攻击路径 = 先证条件律 → 指数倾斜律（Gibbs），再对倾斜律用 Penrose 型 stabilizing 泛函 CLT；**锚点选择**那一步须用 **Palm/Slivnyak** 显式处理，不能由无条件块 CLT 自动给出。

---

## 6. 条目 6：孔洞场 / cavity field / `9/8` 极点

### 6.1 (a) 「纵向步长至多 1」型长步/长笼罩机构（issue 里称 KING）

**未找到**该机制。相邻文献是**长程渗流（long-range percolation）**，即**边概率随距离幂律衰减** `P(r) ~ C r^{−s}`（**不是**「步长 ≤ 1」）：
- **Crawford–Sly**, *Simple Random Walk on Long Range Percolation Clusters II: Scaling Limits*, arXiv:`0911.5668`；后续 arXiv:`2403.18532`。**直接引用**（前者）：当 `s ∈ (d, d+1)`，无穷簇上简单随机游走的标度极限收敛到 **α-稳定 Lévy 过程**，`α = s − d`，quenched 与 annealed 皆成立。
- Kesten 本人的长程渗流工作（Durrett–Kesten；Grimmett–Keane–Marstrand 的连通判据；Kesten 对 `Z^{d−e}×Z^e_+` 的**可和性充要条件**），见 **Grimmett**, *Harry Kesten's work in probability theory*, arXiv 版与 *PTRF* (2021) DOI `10.1007/s00440-021-01046-4`。
**已检索关键词：** `Kesten long range percolation long step`, `long run mechanism renormalization`, `stretched cluster long step`, `directional step bound percolation mechanism`, `KING mechanism percolation`.

### 6.2 (b) cavity field / hole field / 孔洞机制

**已有**（但**机制不同**）：
- **孔洞 ↔ 互补（白）簇的对偶**是标准的：**Isichenko**, *Percolation, statistical topography, and transport in random media*, **Rev. Mod. Phys. 64 (1992) 961–1043**。**直接引用**（转述自其 §2）：「an internal hull can be considered to be the **external hull of a complementary cluster of vacant sites** that fills up a hole in the original cluster」；并给出内/外 hull 的不同普适指数（外 hull `d_h = 7/4`，`D_H = 2 − x_1^C`；unscreened perimeter `d_u ≈ 1.343`）。
- **孔洞尺寸分布**：Hu 等发现「`n_h ~ h^{−τ}`，`τ = 1 + d_f/d`」（hyperscaling），且 hole 是 **volatile fractal**；**largest hole** 满足 `h_max = ⟨C/L^d⟩ ≈ h_{max,0} + a L^{d_H − d}`（`d_H = 7/4`/`4/3`），见 *Size distributions of the largest hole in the largest percolation cluster and backbone*, *Physica A* (2021), DOI `10.1016/j.physa.2021.125847`（S0378437121000789）。
- 教学式「above `p_c` 的簇像 **Swiss cheese**，洞的典型尺寸 = ξ」（见 Geometry of Clusters, Springer 2024, DOI `10.1007/978-3-031-59900-2_5`）。
**但**：spec 描述的「**黑 NN（`p↓0`）的孔洞场拉长巨大白簇**」这一**具体机制未见**。
**已检索关键词：** `cavity field percolation`, `hole field percolation`, `holes percolation cluster scaling`, `swiss cheese cluster holes`, `internal hull complementary cluster`, `volatile fractal holes backbone`, `cavity mechanism random cluster`.

### 6.3 (c) fugacity 倾斜的活动矩阵（非行随机）+ 谱半径的极点/收敛半径判定

**有样板**：**R-矩阵 / Yang–Lee 零点**方法——**arXiv:`0907.4037`**（*Critical exponents from cluster…*）。**直接引用**：
> 「It is always possible … to define a tridiagonal symmetric R matrix which satisfies `(R^n)_{11} = (−1)^n (n+1)b_{n+1}` … `ρ(z) = Σ_{n≥1} n b_n z^n = z(I + zR)^{−1}_{11}`。」
> 「Alternatively … `ρ(z) = Σ_λ ψ_1(λ)² / (z^{−1} + λ)` … **The reciprocals of the eigenvalues of this matrix are the Yang-Lee zeroes of the grand-canonical partition function.** … `ρ(z)` has two singular points at `z` values for which `−z^{−1}` coincides with the **spectrum edges** of the R matrix, leading to vanishing of the denominator.」
→ 这正是「**fugacity `z` 依赖的矩阵谱 → 极点 → 收敛半径/奇异性**」的**严格样板**（`z^{−1}` 型极点，与 spec 的 `u + v = 1/z ≠ 1` 同型）。
另：Perron 根的**双侧夹逼**与 **Collatz–Wielandt 公式**是标准工具（Meyer, *Matrix Analysis*, ch. 8；见 §7）。
**已检索关键词：** `fugacity tilted transfer matrix`, `Perron root fugacity`, `activity transfer matrix percolation`, `Yang-Lee zero spectrum edge`, `radius of convergence activity generating function`, `singularity analysis generating function pole`.

### 6.4 (d) 几何分布参数 `1/9` 或极点 `9/8`

**未找到。** **已检索关键词：** `geometric distribution 1/9 percolation`, `pole 9/8 singularity percolation`, `9/8 exponent percolation`, `geometric parameter one ninth cluster`, `rational singularity 9/8 lattice model`.

### 6.5 判定（条目 6）

- **(b)(c) 方法学已蕴含**：孔洞/内 hull 的对偶（Isichenko RMP 1992；Hu 等 Physica A 2021）与「fugacity 倾斜矩阵谱端点 → 极点」（arXiv:`0907.4037`）都是**已发表样板**，**可直接引用**。
- **(a) 「纵向步长至多 1」长步机构** 与 **(d) `1/9` / `9/8`** **未见**；`9/8` 若真出现，需自查是否为本项目**自建口径**产生。

---

## 7. 条目 7：`#761` / `#766` first-exit 盒与向量指数矩域

**待定标**：`#761` 独立质量区间、优化 first-exit 盒；`#766` first-exit 的**向量指数矩域**给方向范数提供**可计算内外界**。

### 7.1 现成构件

**(i) 非负矩阵 Perron 根的双侧夹逼（=「可计算内外界」的严格工具）**
**«Upper bounds on the growth rates of hard squares and related models»**, *DMTCS*。**直接引用 Lemma 2（Collatz–Wielandt）**：
> 「Let `A` be an irreducible square matrix with non-negative entries. Then for any vector `x > 0`, the largest eigenvalue of `A` (denoted `λ`) is real and positive and is bounded by `min_i (Ax)_i/x_i ≤ λ ≤ max_i (Ax)_i/x_i`。」
且该文**正是**用它给**转移矩阵主特征值** `Λ_o(m)`（宽度 m 圆柱的列转移阵）做**上界**，且用 CTM 型向量逼近主特征向量以获得**紧界**：
> 「we do not compute the eigenvalue exactly. Instead we find upper bounds for `Λ_o(m)` using the **Collatz-Wielandt formula**。」
→ **这就是 spec 要的「转移矩阵特征值 enclosure / Perron root 界」的可用先例**（它给上界；双侧界即 `min` 与 `max` 同时算）。

**(ii) 「inclusion interval」的严格理论**
**Oepomo**, *A contribution to Collatz's eigenvalue inclusion theorem for nonnegative irreducible matrices*, **Electron. J. Linear Algebra 10 (2003) 31–45**（Zbl `1022.15018`）。**直接引用（zbMATH 综述）**：
> 「the ‘coherence’ (i.e. simultaneous closeness) of the **Collatz–Wielandt lower and upper estimates** `m(x)` and `M(x)` of `Λ[A]` (**forming an “inclusion interval”**) for variable positive `x`'s … implying that the set of all the inclusion intervals forms a two-dimensional wedge-shaped domain.」
→ 「**区间套收敛**」的严格结果，可直接支撑「**证书式内外界**」。

**(iii) first-exit 尾律的标准工具**
- **Lifshits–Shi**, *The first exit time of Brownian motion from a parabolic domain*, **Bernoulli 8 (2002) 745–765**（Zbl `1018.60084`）：用 **LDP（Schilder）+ 变分**给出 `lim T^{−(p−1)/(p+1)} log P(τ_D > T)`，其中把**加性泛函的 Biane–Yor 定理**用于求解变分问题（`d=a=1, p=2` 时 `−3π²/8`）。
- **Freidlin–Wentzell** 框架（first exit / Arrhenius）：`lim_{ε→0} ε log E τ_D^ε = inf_{x∈∂D} V(x)`（见 arXiv:`2306.11418` 综述式 (7)）。
- 跳跃扩散 first-exit 的 **MGF/均值 PDE–积分方程**：**Lefebvre**, *Similarity Solutions of PDIE from the Theory of Stochastic Processes*, *Symmetry* **17** (2025) 704, DOI `10.3390/sym17050704`（含 **moment-generating function of the first-passage time** 的 PDIE）。

### 7.2 判定（条目 7）

- 「**Perron 根 / 转移矩阵主特征值的可计算内外界**」**已有严格样板**：Collatz–Wielandt（**DMTCS** 用它对**转移矩阵**做界）＋ inclusion interval 理论（**Oepomo EJLA 10 (2003) 31–45**）。**可直接引用，不必重造。**
- 「**first-exit 尾律**」的标准工具是 **Lifshits–Shi（Bernoulli 2002）** 与 **Freidlin–Wentzell**。
- **未见**：「**向量**指数矩域给**方向范数**提供**内外界**」这一**组合**。
  **已检索关键词：** `first exit box`, `exponential moment domain`, `directional norm certificate`, `certified bounds mass cylinder`, `transfer matrix eigenvalue enclosure`, `Perron root interval arithmetic`, `Collatz bound`, `Collatz-Wielandt transfer matrix`, `first exit time large deviation tail`, `moment generating function first passage PDE`, `eigenvalue inclusion interval nonnegative matrix`.
- **裁定**：**构件齐备、组合未见**。**建议**：把 `#766` 明确定位为「(i)+(ii) 的向量化/方向化推广」，并说明与 Oepomo 的 inclusion interval 的差别（后者是**标量**谱半径，`#766` 要**方向范数**）。

---

## 8. 条目 8（低优先）：`#757` 投稿先例矩阵

**待定标**：几何论文的**逐定理先例矩阵**——任意整数周期、两次出生（two-birth）、周期簇。

**已有（几何/渗流精确判据）**：
- **Scullard–Jacobsen**, *Transfer matrix computation of generalised critical polynomials in percolation*, **arXiv:`1209.1451`**, *J. Phys. A* **45** (2012) 494004。**直接引用**：「the critical polynomial `P_B(p)` … may be defined on **any periodic lattice**. The polynomial depends on a finite subgraph `B`, called the **basis**, and the way in which the basis is **tiled** to form the lattice.」；`P(A,B,C) = P(Ā,B̄,C̄)`（式 (1)）即精确判据，`A,B,C` 为三角形**三边界顶点三连通/三不连通**概率；「we can also treat **site percolation** problems by reasoning on the **covering lattice** or by introducing **correlations**」。
- **Scullard**, *The computation of generalized percolation critical polynomials by the deletion–contraction algorithm*, **arXiv:`1207.3340`**（deletion–contraction 定义；任意 2D 周期格）。**直接引用（式 (2)）**：hexagonal 的临界曲面 `H(p,r,s) ≡ prs − pr − ps − rs + 1 = 0`；并给出 `FE(p,r,s,t,u,v) = pA(r,s,t,u,v) + (1−p)H(s, ur, tv)` 型**一阶（first-order in each argument）**递推。
- **Scullard**, *Percolation critical polynomial as a graph invariant*, **Phys. Rev. E** 2012, PubMed `23214553`（「the generalized critical polynomial can be viewed as a **graph invariant**, similar to the Tutte polynomial … can be found using the recursive **deletion–contraction** algorithm」）。
- **Jacobsen–Scullard 2019/2020**, **arXiv:`1910.12376`**：给出两类有限尺寸修正指数（`Δ = 6,7,8` 与 `Δ = 4,6,8`），**并在摘要明列 open question**：
  > 「We discuss the open questions related to the method, such as **the full scaling law**, as well as its potential for determining critical points of other models.」

**未找到**：「**两次出生（two-birth）**」、「**balance root**」、「**matching root**」、「**周期簇（periodic cluster）**」作为**已发表术语/定理**。
**已检索关键词：** `arbitrary period percolation`, `born distribution percolation`, `two-birth percolation`, `periodic cluster percolation`, `balance root percolation`, `matching root percolation`, `critical polynomial full law`, `full-law criterion percolation`, `critical polynomial graph invariant deletion-contraction`, `basis tiling periodic lattice critical polynomial`.

### 判定（条目 8）
- 「**任意整数周期**」（= 任意周期 lattice 与任意 basis tiling）与「**精确判据**」**在临界多项式文献中已有系统处理**，**可逐条引用**（arXiv:`1209.1451`、arXiv:`1207.3340`、PRE 2012 PubMed `23214553`、arXiv:`1910.12376`）。
- 「**两次出生 / balance root / matching root / 周期簇**」**未见已发表术语**；建议在 `#757` 的投稿矩阵里**把这几项标为「未见，需自建定义并声明术语新」**，**不要**假定它们有先例。
- **利好消息**：`1910.12376` 明确把「**the full scaling law**」列为**该领域 OPEN question** —— 若 `#757` 的「full-law criterion」正是这一条，则**领域承认它是开放的**，投稿定位反而更好。

---

## 9. 检索轮次日志（≥6 轮，每轮关键词）

| 轮 | 目标条目 | 用过的关键词（原样） |
|---|---|---|
| R1 | 1 | `Mertens Ziff matching polynomial square lattice finite-size correction exponent`；`Jacobsen matching polynomial square lattice critical exponent 21/8`；`8-arm exponent 21/8 percolation matching c=0 Kac table` |
| R2 | 1 | `correction to scaling exponent L^-4 wrapping probability percolation torus irrelevant exponent`；`"matching polynomial" percolation finite size correction exponent W_4 W_8 winding number square lattice`；＋ 下载 arXiv:`1603.07289v2`（PDF）与 arXiv:`2204.01517`（PDF） |
| R3 | 1/3 | `Jacobsen 2015 percolation threshold correction exponent 4 critical polynomial "L^{-4}"`；`Scullard Jacobsen critical polynomial winding configurations Z_0 Z_1 Z_2 torus homology counting`；`cycle rank distribution random subgraph torus Tutte polynomial rank generating function counting` |
| R4 | 1/3/4 | `Jacobsen 2015 "critical polynomial" scaling exponents conformal field theory L^{-4} site percolation square lattice`；`Betti number distribution random subgraph torus homology rank cycle rank counting occupied sites`；`normal exponential scale mixture Laplace distribution variance gamma named identity` |
| R5 | 1/2/5 | `matching function exponent 3.42 OR 3.25 percolation Mertens Ziff finite-size correction measurement later`；`"self-normalized" observable logarithmic derivative finite-size asymptotic ratio of observables critical exponent estimate`；`central limit theorem additive functional random interval marked point process cluster Palm rare event` |
| R6 | 1/3/6 | `Aizenman Duplantier Aharony "wrapping" OR "topological" counting configurations cycle rank torus percolation exact enumeration`；`"cavity field" OR "hole field" percolation random cluster long range correlation white cluster mechanism`；＋ 下载 arXiv:`1507.03027` HTML、arXiv:`2011.11903` 摘要页并 grep §7.1/作者 |
| R7 | 2/5/6 | `universal amplitude ratio finite-size scaling derivative ratio observable plateau approach exponent logarithmic derivative order parameter`；`fugacity tilted transfer matrix Perron root spectral radius singularity activity generating function percolation`；`Kesten long range percolation long step mechanism stretched cluster renormalization "long run"` |
| R8 | 2/5/6 | `"logarithmic derivative" finite-size scaling estimator correlation length exponent pseudo-critical Binder cumulant derivative`；`stable central limit theorem infinite variance critical cluster exploration conditional invariance principle percolation`；`percolation activity transfer matrix fugacity "activity representation" Yang-Lee singularity radius of convergence cluster weight 1/z` |
| R9 | 7/8/6 | `Collatz-Wielandt bound Perron root enclosure interval arithmetic nonnegative matrix transfer matrix eigenvalue certified bounds`；`"first exit" box additive functional exponential moment generating function domain large deviation rate`；`matching polynomial root percolation threshold "arbitrary period" periodic cluster birth distribution balance root` |
| R10 | 8/6/5 | `"critical polynomial" percolation "full" law algorithm arbitrary lattice basis number of births exact thresholds Chen Li`；`holes in percolation clusters large hole scaling number of holes hull swiss cheese cluster topology`；`Palm conditioning central limit theorem nondegenerate Fisher information exploration process Poisson cluster functional` |

共 **10 轮**、**30 次检索查询**＋ **4 次全文抓取/核对**（MZ16 PDF、Jacobsen 2015 HTML 逐式、homological percolation 摘要、arXiv:`2204.01517` 摘要）。

---

## 10. 对主代理的直接影响

**(A) 可以直接引用、不必再证的**
- **条目 4**：`√(cE)Z ∼ Laplace`（含 `c` 的标定）= **Gauss–Laplace transmutation**（Ding–Blitzstein, arXiv:`1510.08765`；渊至 Andrews–Mallows 1974 / West 1987）。**直接引，别当新结果。**
- **条目 1′**：`alpha_j = (j²−1)/12` 是 **ADA 1999 PRL 83 1359**（`x_ℓ^P = (ℓ²−1)/12`，`ℓ` = 臂/线数）的已发表公式。**arms-audit §1.3 的裁定成立且已由原文确认。**
- **条目 3 的 `Z_0/Z_1/Z_2` 与 `P_B = R_2 − R_0`**：**Scullard–Jacobsen / arXiv:`2010.02887`**。
- **条目 6(b)(c)**：孔洞 ↔ 互补簇对偶（**Isichenko RMP 64 (1992) 961**；Hu 等 *Physica A* 2021）；fugacity 依赖矩阵谱端点 → 极点（**arXiv:`0907.4037`**，R-矩阵/Yang–Lee）。
- **条目 7**：Collatz–Wielandt 对**转移矩阵主特征值**做界（**DMTCS** Lemma 2）＋ inclusion interval（**Oepomo EJLA 10 (2003) 31–45**）。
- **条目 5**：随机几何加性泛函 CLT（**Penrose arXiv:`math/0410021`**，**应用含渗流簇计数**）；marked PP 的 CLT/稳定律（**arXiv:`1903.09387`**）；Palm 渐近正态（**Prokešová–Jensen**）。
- **条目 8**：临界多项式（任意周期 basis / 任意 2D 周期格 / deletion–contraction / graph invariant）：**arXiv:`1209.1451`、arXiv:`1207.3340`、PRE 2012 PubMed `23214553`、arXiv:`1910.12376`**。

**(B) 必须改述的（重要）**
- **条目 1 的措辞**：**不能**写「文献数字让 `21/8` 消失」。正确表述：
  - **MZ16**（PRE 94 062152）**直接测** `M_L(p_c)~L^{2−x}` 得 `2−x = −3.42`（`21/4` 要求 `−3.25`）；阈值估计子斜率 `−4.07`（要求 `−4`）；但**作者自认系统太小**（L ≤ 48），**未宣称渐近**。
  - **J15**（J. Phys. A 48 454003, arXiv:`1507.03027` §7.1，**同一方形格点 site 渗流模型**）用**转移矩阵特征值法**（圆柱周长 n ≤ 21，**非 MC**）测得 **`Δ_1 = 4.0001(2)`**，并**承认 `Δ_1 = 4` 精确**；CFT 依据 = 两拓扑扇区共享 `x_m = 5/48`。
  - `Δ_1` 与 MZ16 的 `w` **是同一个指数**。**A 说 4、C 说 −4.07，二者冲突**；MZ16 自己写「somewhat larger than the value 4 suggested by Jacobsen」。
  - **结论**：`21/8` **既未被证实也未被证伪**；`x=21/4 ⟺ w=4` **有一份高精度独立支持（J15）和一份低精度冲突（MZ16）**。
- **条目 3 的措辞**：必须**显式区分** **「绕行方向型」`Z_0/Z_1/Z_2`** 与 **「占据子图同调秩」`r_black ∈{0,1,2}`**。二者**不等价**（可由两个不同簇分别绕行得到 `Z_2` 但无单簇 cross-wrap；反之 `r=2` 也可由两个不同簇提供两个独立非可缩圈）。**不能**把 `Z_·` 计数当 `C[L,j,k]`。
- **条目 6(a)(d)**：`KING`/「纵向步长 ≤ 1」长步机构、`1/9` / `9/8` 极点 —— **未见**；若真要用，须自查是否为本项目**自建口径**的产物，投稿时**不得**暗示有先例。

**(C) 真空白 / 值得投入的（按价值排序）**
1. **`M_L(p_c)` 本身的指数 `2−x` 没有第二份独立测量** —— **这是唯一能判定 `21/8` 的靶子**。把 MZ16 的 `M_L(p_c)`（等价地 `R_L^b − R̂_L^b`）在**环面**上推到 `L ≳ 100`（或精确可算的最大 L），看 `2−x` 趋向 `−3.25` 还是 `−3.42`。**本轮最值得投入的真空白。**
2. **条目 3 的 `C[L,j,k]`（按**同调秩**分辨，而非绕行方向型）未找到已发表版本** —— 口径确实比 `Z_0/Z_1/Z_2` 更细。**但投稿必须在文中显式对照 `Z_·` 与 `R_j`**，否则会被认为与 Scullard–Jacobsen 重复。（可引 **Duncan–Kahle–Schweinhart arXiv:`2011.11903`** 作为「环面同调分辨」的最近邻，并说明它是**高维**、2D 退化为经典。）
3. **条目 5 的「随机区间加性泛函 + 完整 component-Palm 锚点」联合陈述未见** —— 技术空白、风险可控：攻击路径 = Gibbs 条件原理 → 指数倾斜律 → Penrose 型 stabilizing CLT，**锚点**用 Palm/Slivnyak 显式处理。
4. **条目 7 的「向量指数矩域 → 方向范数内外界」组合未见**，但**构件齐备**（Collatz–Wielandt / Oepomo inclusion interval / Lifshits–Shi first-exit）。**投入产出比高**。
5. **条目 2 的 `c_slope`**：**未见同口径量**，但**方法学先例充分**（对数导数 FSS 估计子）—— **定位为「(i) 的一例」即可，不必重造理论**。

**(D) 可信度标注**
- **均为直接引用（WebFetch / HTML 逐字核对）**：MZ16 式 (15)(17)(20)(38)(39) 及正文的 `2−x=−3.42`、斜率 `−4.07`、`w=−4.17`、`x−2=3.25`；J15 式 (34)(35)(36)(39)(40)(50) 及 `Δ_1=4.0001(2)`、`x_m=5/48`、§7.1 标题；ADA 式 (1)(2)；`2010.02887` 的 `Z_·` 与 `P_B`；`2011.11903` 摘要；`2204.01517` 摘要；Ding–Blitzstein §1；arXiv:`0907.4037` 的 R-矩阵/`ρ(z)` 公式；DMTCS Lemma 2；Oepomo 综述（zbMATH Zbl 1022.15018）；arXiv:`0911.5668` 摘要；arXiv:`1910.12376` 摘要（含「the full scaling law」）；Isichenko RMP 1992 的内 hull 对偶；Physica A 2021 的 largest-hole 标度。
- **转述（非引号）**：Palm 教程 arXiv:`1512.05871`；Prokešová–Jensen；Lifshits–Shi（经 zbMATH Zbl 1018.60084）；Freidlin–Wentzell 综述（arXiv:`2306.11418`）；Lefebvre *Symmetry* 2025；arXiv:`1311.4967`（对 Dvoretzky–Robbins 的转述）。
- **本项目内部算术（非文献）**：表格中「`21/4` 要求 `2−x = −3.25`」「`w = x−2+1/ν`」的换算（arms-audit §3 已给）。
- **未找到的项**：§2.2、§3.3、§5.2、§6.1、§6.4、§7.2、§8 均已列出**本轮全部关键词**。
