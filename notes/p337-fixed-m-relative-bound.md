# P337：共同 bulk 因子能否把 joint-limit 提升到固定 m？

结论：**不能据现有证明直接提升。** 若有真正的相对受限 partition 估计，bulk 确实可以消去；但当前 `exp[A_N(m)]` 是放松轮廓相容性所得的上界，并未被证明是分子、分母共有的同一因子。即使把共同 bulk 因子精确消去，并另外给予两几何对所有 h 的指数 rank1 抑制，仍不足以推出原 pooled U 消失：下文给出一个正权、单调、唯一 pooled root 的可复核反例。

本次得到三个适用于实际局部模型的有限体积相对界，其中最直接的是

\[
m^{-2W}\le Z_{\rm twist}^{\rm loc}/Z_{\rm rect}^{\rm loc}\le1,
\qquad W=25k.
\]

它把归一化损失控制在表面阶，且对 h 一致，但 `50 k log m` 的常数仍不能闭合现有 oblique Peierls 不等式。这不是固定 m 定理的反例，而是当前证明前提不足的明确边界。

## 固定来源与实际执行范围

仅阅读以下固定 Git blobs，没有追逐更新：

- `2690f665bc8029cb2370d3f1efcef5eb2853705c:notes/closed-source-poisson-double-scaling.md`，全文。
- 同提交 `notes/closed-source-oblique-twist-comparison.md`，全文。
- 同提交 `notes/closed-source-pooled-sector-odds-bound.md`，全文。
- `e17b286b29a02b7d74fa7d7489858729a8952a11:notes/closed-source-fixed-coupling-peierls.md`，全文。
- 同提交 `notes/closed-source-q-lift-and-thermal-quotient.md`，全文。

仓库：`LightChainr/Matching-One`。理论子任务只读现有证明，未访问云机、运行MC或枚举格点配置；结果随后由总任务收入研究分支。唯一数值内容是文末三个三态分布的标准库 `Fraction` 恒等式核对。

## 1. 当前 bulk 上界不能按符号名称直接约掉

Poisson note 的真实推导是

\[
Z_{\rm bad}\le 2M_h\prod_\gamma(1+u_\gamma)
                  \sum_{|\gamma|\ge R}u_\gamma,
\quad
u_\gamma=m^{2-|\gamma|},
\quad Z\ge M_h.
\]

`prod(1+u)` 来自丢掉相容性、嵌套和 cut-consistency；它不是已识别的实际小轮廓 partition。现有 rank0 下界是独立集气体，而不是同一个自由轮廓乘积。两者都具有 `exp[O(N/m²)]` 量级，并不等于它们的比值受次表面阶控制。即使某次有限阶展开把差缩到 `O(N/m⁴)`，在固定 m、N~k² 时也仍会压过 k 阶界面代价。

一个确实可以消掉 bulk 的标准相对不等式是：若实际完整模型已表示为**同一组非负活动度**的相容 polymer 气体，

\[
Z=\sum_{\Gamma\ {m compatible}}\prod_{\eta\in\Gamma}z_\eta,
\quad
Z(\gamma\in\Gamma)=z_\gamma Z_{\text{compatible with }\gamma}
                    \le z_\gamma Z,
\]

则 `P(gamma) <= z_gamma`。最后一步使用的是受限配置集包含于原配置集，且权重完全相同。

**现有 `u_gamma` 不能直接代入这个式子。** 实际源含 `h^K`、占据组件奖励及 rank 因子；翻转一个围成 a×a 占据岛的轮廓相对空配置，其精确权重是

\[
h^{a^2}m^{2-4a},
\]

不是单独的 `m^(2-4a)`。对任意 h>1，声称删除该轮廓后权重比不超过后者，已经逐配置失败。围内/围外的受限 partition 比、相变背景和嵌套奖励必须保留。本例只否定这条未经证明的删除不等式，不否定一个正确重整化后的 contour 表示可能成立。

## 2. 一个实际成立的、对 h 一致的缝合 partition 比界

令 Gt 是 `W×k` 的 shift-7k torus，Gr 是同尺寸无 twist torus；W=25k。删除 W 条竖直缝边后，两者得到同一个开缝圆柱 Gc。在同一顶点占据集合 A 上定义普通图 cycle rank

\[
\beta_G(A)=B_G(A)-K(A)+C_G(A).
\]

这里 beta 是占据图的普通循环维数，**不是** ambient image rank r。向 Gc 加回一条占据边，若连接两个组件则 beta 不变，若闭合一个环则 beta 增加1。因此

\[
\beta_{Gt}=\beta_{Gc}+j_t,\quad
\beta_{Gr}=\beta_{Gc}+j_r,\qquad 0\le j_t,j_r\le W.
\]

两种闭合都是四正则图，所以 `Bmix=4K-2B`；局部颜色模型的占据边际满足

\[
w_G^{loc}(A)=h^K m^{-Bmix+2C}
            =h^K m^{-2K+2\beta_G(A)}.
\]

于是逐配置有

\[
m^{-2W}\le w_{Gt}^{loc}(A)/w_{Gr}^{loc}(A)\le m^{2W}.
\]

对相同的全部占据集合求和，得到真正的归一化相对界

\[
m^{-2W}\le Z_t^{loc}/Z_r^{loc}\le m^{2W}.
\tag{A}
\]

这部分对所有实数 m>=1、h>0 成立。已有 Gram/PSD 证明在颜色实现的整数 m 范围另给 `Zt<=Zr`，故

\[
0\le\Delta_k=\log Z_r^{loc}-\log Z_t^{loc}\le2W\log m
             =50k\log m.
\tag{B}
\]

它不需要强制空行、不需要知道 moving root、不把几何投影塞入 unmarked transfer operator。对原 projected occupation 权重还可逐配置保留 `m^(-r)`，得到更粗的比界 `m^(-2W-2) <= w_star,t/w_star,r <= m^(2W+2)`；此处 r 可以改变，不能据此声称各 rank sector 一一对应。

**为什么仍不够：** 固定提交的斜环面证明要求

\[
\limsup \Delta_k/k < 7\tau_\infty,
\qquad \tau_\infty=-\log[3(2/m)^{1/4}].
\]

而 (B) 给出的系数是 `50 log m`，允许阈值仅为
`(7/4)log m - 7log3 - (7/4)log2`。把 (B) 代回现有 contour 界，仍不能得到 k→∞ 的衰减。这里是一个明确的失败常数，不是仅说需要更好估计。

## 3. 实际成立的 rank0/rank2 受限比界

这些界比较的是**同一 W×k 圆柱的 twisted 与 untwisted 闭合**，不是直接比较原来的 5k×5k axis 与斜伴随；不可混用几何。

令 H0 是整条接缝行全空，H2 是该行全占据。对同一占据集合：

- 在 H0 上所有跨缝占据边消失。两图占据权重相同，且 ambient rank0 判定相同。
- 在 H2 上整行已经包含水平 essential cycle。改变竖缝端点只是在这条已占据的水平行上重接，Bmix 和 C 相同；有无独立竖向 winding 也相同。因此 rank2 判定和该 sector 权重相同。

后一点也可直接检查：把 twisted seam 的一个端点换回 rectangular 端点后，用全占据行上的水平路连接原端点。这个替换保留竖向 winding；已有水平环保证 ambient rank2 与否不变。故

\[
Z_{0,t}(H_0)=Z_{0,r}(H_0),\qquad
Z_{2,t}(H_2)=Z_{2,r}(H_2).
\tag{C}
\]

给定其余占据状态，添加一站点若有 d 个占据邻居、分属 c 个组件，其局部颜色占据 odds 是

\[
h m^{-2+2(d-c)},\qquad h/m^2\le\text{odds}\le h m^4.
\]

在固定 rank0 条件下，删除占据点仍属 rank0，所以该条件只能进一步有利于空点；在固定 rank2 条件下，添加占据点仍属 rank2，所以该条件只能进一步有利于占据点。rank 投影在每个固定 sector 中为常数。逐点条件化得到

\[
P(H_0\mid r=0)\ge (1+h m^4)^{-W},
\qquad
P(H_2\mid r=2)\ge(1+m^2/h)^{-W}.
\]

结合 (C)，可得明确的相对受限 partition 不等式

\[
\left|\log\frac{Z_{0,t}}{Z_{0,r}}\right|
 \le W\log(1+h m^4),
\qquad
\left|\log\frac{Z_{2,t}}{Z_{2,r}}\right|
 \le W\log(1+m^2/h).
\tag{D}
\]

两条固定 law 都适用，前提是比较同一个 h。它们确实只付出表面阶，但相应 sector-odds 界

\[
\left|\log\frac{Z_{2,t}Z_{0,r}}{Z_{0,t}Z_{2,r}}\right|
\le W\log[(1+h m^4)(1+m^2/h)]
\tag{E}
\]

在 h~1 时系数约为 `150 k log m`，仍远不足以闭合现有 winding/denominator 竞争。甚至把右侧对 h 最小化，也只能得到 `2W log(1+m^3)`；这不是对实际 mismatch 的下界，只说明这套粗上界的最好常数仍太大。

一个额外但不能跨 root 使用的精确恒等式是

\[
Z_{j,drop}(h,m)=m^j Z_{j,star}(h,m).
\]

所以两 law 在**同一个 h** 的跨几何 sector cross-ratio 完全相同。它们各自 pooled root 的 h 不同，该恒等式不能替代在各自根处的 odds 控制。

## 4. 正性和 sector0 权重本身不消除 twist penalty

oblique note 已给 `D=25M+1`、P 固定一点并有 M 个25周期、
`T=I+D^(-2)11^T` 的正定、严格逐元素正矩阵。它与 P 对易，且

\[
R_k=\frac{\operatorname{Tr}(T^kP)}{\operatorname{Tr}(T^k)}
=\frac{(1+1/D)^k}{D-1+(1+1/D)^k}.
\]

补充一个直接的 sector0 读数：该例有

\[
w_0=\frac{M+(1+1/D)^k}{25M+(1+1/D)^k}
   =\frac1{25}+\frac{24}{25}R_k >\frac1{25}.
\]

当 D 随 k 指数增长时，Rk 仍可指数小。因此“sector0 有正权、甚至始终超过1/25”、严格 positivity、Perron 向量属于 sector0、PSD 和固定 order25 同时成立，也不能保证所需的次指数 twist ratio。给 T 乘一个正标量所产生的共同 bulk 因子，会在两 trace 中精确消去，不改变这个结论。

这沿用固定提交的反例并明确了它的 w0 值；不是实际局部颜色模型的反例。`w0>1/2+epsilon` 等更强质量界当然足够，但现有局部模型尚未从正性推出该前提。

## 5. 新反例：精确 bulk 消去 + uniform rank1 抑制仍不推出 pooled U 消失

下面反例比仅给静态 sector 分布更强：给出完整共同 thermal 参数、唯一 pooled root、全部导数和原 U 的反例。它针对逻辑推论，不声称实现了 square-lattice 源。

取偶数 N~L²、固定 m>1，令

\[
t=m^{-L},\quad a=t^4,\quad c=t^3,\quad x=h^{N/2}.
\]

每几何有三个状态，rank r=0,1,2，`q=r-1`、`E=q²`、`K=Nr/2`。给任意严格正的共同 bulk 因子 B_N(h)，例如 `exp[N f(h,m)]`，定义

\[
(Z_{0,f},Z_{1,f},Z_{2,f})=B_N(h)(a,cx,x^2),
\quad
(Z_{0,s},Z_{1,s},Z_{2,s})=B_N(h)(1,cx,ax^2).
\tag{F}
\]

所有权重都严格正；共同 bulk 精确消去。三个状态是全序，因此两 law 都关联/FKG。更强地，`K=N(q+1)/2` 给出

\[
\partial_{\log h}\langle q\rangle_g=(N/2)\operatorname{Var}_g(q)>0.
\]

两几何等权 pooled root **唯一且恰在 h=1**。此处两 total partition 完全相同，都是 `B_N(1)(1+a+c)`。对所有 h 还有

\[
\sup_{h>0}P_{1,f}(h)=\sup_{h>0}P_{1,s}(h)
=\frac{c}{2\sqrt a+c}=\frac{t}{2+t},
\]

所以 rank1 概率在每几何都对 h 一致地指数小。不是仅在一个挑选的活动度上小。

但在 pooled root，令 Z=1+a+c，有

\[
b=P_1=c/Z,\quad q_f=(1-a)/Z=-q_s,
\quad \kappa=\frac{4a+c(1+a)}{Z^2}\sim t^3.
\]

于是原来的归一化热导数比严格为

\[
\frac{U}{A_N}
=\frac{(E_{h,f}-E_{h,s})/\Delta_{angle}}
        {(q_{h,f}+q_{h,s})/2}
=\frac{2c(1-a)}{\Delta_{angle}[4a+c(1+a)]}
=\frac{2(1-t^4)}{\Delta_{angle}(1+4t+t^4)}
\longrightarrow\frac2{\Delta_{angle}}\ne0.
\tag{G}
\]

分母始终正；失败不是除零，也没有漏掉任一 normalizer 或 root motion。两个均值分别趋于±1，几何之间抵消了 pooled q，但两几何内部的方差随 b 一起消失。

此例的实际 sector cross-ratio 是 `Xi=a^(-2)=t^(-8)`，故 `|log Xi|=8L log m`。这正是独立 odds 条件缺失的地方。甚至两几何各自根的 `log h=±4L log(m)/N` 都趋于0，也仍不保证它们在 O(1/N) 共存窗口内对齐。所有 h 下的总 pressure 密度差至多 `4L log(m)/N ->0`；在 pooled root 该差更是精确为0。

因此，共同 bulk、相同热参数、总 pressure 一致、正性、FKG、唯一 pooled root 和 uniform winding 抑制都不能单独关闭原 U 的分母问题。局部模型可能通过额外结构排除 (F)，但这些结构不能被当前正性论证代替。

## 6. 可复核的小有理数核对

实际执行以下标准库脚本，直接从三态 normalized moments 算导数，再核对 (G)。没有调用格点引擎或枚举任何格点配置。

```python
from fractions import Fraction as F
N = 8
for t in [F(1, 2), F(1, 16), F(1, 256)]:
    a, c = t**4, t**3
    stats = []
    for weights in [(a, c, F(1)), (F(1), c, a)]:
        z = sum(weights)
        p = [v/z for v in weights]
        qvalues, kvals = [-1, 0, 1], [0, N//2, N]
        q = sum(v*q for v, q in zip(p, qvalues))
        e = sum(v*q*q for v, q in zip(p, qvalues))
        k = sum(v*k for v, k in zip(p, kvals))
        dq = sum(v*q*k for v, q, k in zip(p, qvalues, kvals))-q*k
        de = sum(v*q*q*k for v, q, k in zip(p, qvalues, kvals))-e*k
        stats.append((q, e, k, dq, de, p[1]))
    assert stats[0][0] + stats[1][0] == 0
    u = (stats[0][4]-stats[1][4])/((stats[0][3]+stats[1][3])/2)
    assert u == 2*(1-t**4)/(1+4*t+t**4)
    kappa = (4*a+c*(1+a))/(1+a+c)**2
    assert (stats[0][3]+stats[1][3])/2 == N*kappa/2
    print(t, stats[0][5], t/(2+t), u)
```

实际输出列为 `t, P1(root), sup_h P1, U/A (Delta=1)`：

```text
1/2    2/19               1/5    30/49
1/16   16/65553           1/33   43690/27307
1/256  256/4294967553     1/513  2863311530/1454025387
```

## 本次停止结论

固定 m 的**原斜对 pooled U**不能由已证 Poisson joint-limit 加“共有 bulk 可约去”直接宣布成立。实际可证的 (A)–(E) 只达到常数不够的表面阶；(F)–(G) 则明确排除了仅凭现有正性和小 rank1 质量补齐分母的逻辑捷径。本次理论子题到此结束，不追加耦合、模拟、扫描或后继计算清单。

独立交叉复核已通过(A)–(G)：核对了全空/全占据缝行的rank和权重对应、固定rank条件odds的方向、全域唯一root、uniform-in-h的P1最大值及三个Fraction输出。主任务另以三态协方差直接重算(G)。其中(B)的`Z_twist≤Z_rect`仍只引用整数m颜色PSD范围；三态分布始终只作逻辑反例。
