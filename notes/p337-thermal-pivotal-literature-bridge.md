# P337 thermal-pivotal 文献桥：#536 来源到 #537 当前闸门

**状态（2026-09-01）：`LITERATURE_BRIDGE_OPEN`。** 本卡只回答一个问题：二维临界渗流已有的
pivotal measure、四臂/六臂、Russo 公式、spectral sample 与 near-critical coupling，究竟能否
关闭已由仓库团队关闭为重复项的 Issue [#536](https://github.com/LightChainr/Matching-One/issues/536)
所保留的两项 midpoint 通道，以及开放 P0 [#537](https://github.com/LightChainr/Matching-One/issues/537)
的完整 original-U thermal gate。
结论是：**raw kernel 已经关闭；fused/collision 的六臂尺度输入已经具备，但尚待项目的
support reduction 才能调用；任一完整 midpoint 通道都没有关闭。** 缺口不是再测一个 raw
距离点，而是一个带 `q/E` 权重、符号和三位置 landing 信息的新不等式。

本卡没有新增枚举、Monte Carlo、数值拟合、测试或服务器任务，也不把三角格精确指数移植到
square-site/matching 模型。

## 1. 要约束的不是一个普通 pivotal count

对有限 Bernoulli 图、固定 pair kernel `G=g_xy` 和 `O in {q,E}`，沿用 #536 的定义

\[
 \Delta_zF=F(\omega^{z=1})-F(\omega^{z=0}),\qquad
 F_z^{\rm mid}=\frac{F(\omega^{z=1})+F(\omega^{z=0})}{2}.
\]

两项通道是

\[
 X^{\rm obs}_{xyz}(O)
   =(g^{\rm mid}_{xy,z}-\mathbb E_pg_{xy})\,\Delta_zO,
 \tag{1}
\]

\[
 X^{\rm ker}_{xyz}(O)
   =(O^{\rm mid}_z-\mathbb E_pO)\,\Delta_zg_{xy},
 \tag{2}
\]

且

\[
 \partial_p\operatorname{Cov}_p(O,g_{xy})
  =\sum_z\mathbb E_{\rm rest}X^{\rm obs}_{xyz}(O)
   +\sum_z\mathbb E_{\rm rest}X^{\rm ker}_{xyz}(O).
 \tag{3}
\]

第一项是 **observable/rank-pivot**：`g` 可以不变而 `q/E` 改变；第二项是
**kernel-pivot**：`q/E` 可以不变而 `g` 改变。它们共享总体、root 和代数分解，不是两张
独立统计票，也不能先取无符号 pivotal 数量再恢复符号。

## 2. 第一手文献给出的精确输入与适用边界

### 2.1 square-site 四臂：足够控制 raw kernel，不足以控制热导数

[van den Berg--Nolin, *On the four-arm exponent for 2D percolation at criticality*,
arXiv:2008.01606](https://arxiv.org/abs/2008.01606) 的 Theorem 1.1 对
**square-lattice site percolation at its critical point** 证明

\[
 \alpha_4\ge 1+\frac{\alpha_2}{2}>1.
 \tag{4}
\]

其第5节实际给出有限尺度输入

\[
 \pi_4(3n)\le \frac{C}{n}\sqrt{\pi_2(n)}.
 \tag{5}
\]

结合项目中 `g_xy` 的双四臂 support，式(5)足以推出
`E_pc |g_xy| <= C d(x,y)^(-2-eta)` 及 raw 空间绝对可和。它没有控制 Russo 求和中的
第三个位置 `z`，因而不能推出 `partial_p E[g_xy]` 或式(1)--(3)可和。特别是，(4)只给
`alpha_4>1`，不是三中心 bulk 粗界所需的 `alpha_4>4/3`。

### 2.2 三角格精确四臂/六臂：是诊断，不是 square-site 定理

[Smirnov--Werner, *Critical exponents for two-dimensional percolation*,
arXiv:math/0109120](https://arxiv.org/abs/math/0109120) 的 Theorem 4 对临界三角格
site percolation 的任意 polychromatic `j>=2` 平面臂事件给出

\[
 \alpha_j=\frac{j^2-1}{12}.
 \tag{6}
\]

所以

\[
 \alpha_4=\frac54,\qquad \alpha_6=\frac{35}{12}.
 \tag{7}
\]

这两个数有两个相反的用途：六臂指数使尺度比上的 fused/collision 和可和；但精确
`alpha_4=5/4` 仍低于 `4/3`，因此即使在三角格，朴素的三个独立四臂 bulk 账也不能关闭
绝对三位置和。它更不能被当作 Matching One 的 square-site 指数定理。

### 2.3 六臂的正增益：直接约束 fused annulus，不能替代三中心 support

[Schramm--Steif, *Quantitative noise sensitivity and exceptional times for percolation*,
arXiv:math/0504586](https://arxiv.org/abs/math/0504586) Corollary A.8 对临界
**square-bond** percolation 证明五臂概率与 `(r/R)^2` 同阶，并存在 `epsilon>0` 使六条
交替 primal/dual 臂满足

\[
 \pi_6(r,R)\le C(r/R)^{2+\epsilon}.
 \tag{8}
\]

因此

\[
 (R/r)^2\pi_6(r,R)\le C(r/R)^\epsilon,
 \tag{9}
\]

在 `R/r=2^k` 上为可和的 `C 2^(-epsilon k)`。三角格由(7)进一步给
`2^(-11k/12+o(k))`。这能接管**已经证明出现六臂 annulus** 的融合尺度；它不能从
`Delta_z g_xy != 0` 自动制造六臂，也没有直接覆盖 square-site/matching 的具体 landing
约定。

### 2.4 pivotal measure：控制正的 quad-pivotal mass，不控制 signed 联合权重

[Garban--Pete--Schramm, *Pivotal, cluster and interface measures for critical planar
percolation*, arXiv:1008.1378](https://arxiv.org/abs/1008.1378) Theorems 1.1--1.2
对临界三角格、piecewise-smooth quad 构造

\[
 \mu_\eta^Q
  =\frac{\eta^2}{\alpha_4^\eta(\eta,1)}
    \sum_{x\;Q\text{-pivotal}}\delta_x,
 \tag{10}
\]

证明它与配置联合收敛，且共形协变指数为 `3/4`。Proposition 4.4 还给出局部一、二阶矩；
若 `U` 是内面中半径 `r` 的方块而 annulus 间距为 `epsilon`，则

\[
 \mathbb E\mu^A(U)\asymp
 \frac{\operatorname{area}(U)}{\alpha_4(\epsilon,1)},
 \qquad
 \mathbb E[\mu^A(U)^2]\asymp
 \frac{r^4}{\alpha_4(\epsilon,1)\alpha_4(r,1)}.
 \tag{11}
\]

若先证明 `Delta_z q`、`Delta_z E` 可由**固定数目**的 macroscopic quad-pivotal 事件覆盖，
(11)可进入 observable-pivot 的 blockwise Cauchy--Schwarz。当前缺的正是这个 torus-homology
cover，以及 pivotal mass 与 centered signed `g_xy` 的联合矩。正的 measure 本身不给式(1)
的符号。

同一论文 Proposition 5.3 证明三角格临界普通连接概率

\[
 \mathbb P(x\leftrightarrow y)
  =(C+o(1))\alpha_1^\eta(\eta,d(x,y)/2)^2
 \tag{12}
\]

对方向一致。这是 ordinary two-point connectivity 的 raw 临界值，不是它的 `p` 导数，更不是
bilocal `g_xy`、torus rank 或 original-U 的导数定理。

### 2.5 Russo 与 near-critical arm stability：给窗口，不给加权三点和

[Nolin, *Near-critical percolation in two dimensions*,
arXiv:0711.4948](https://arxiv.org/abs/0711.4948) 提供三项直接可用的接口：

1. Theorem 1 的 Russo 公式对有限 increasing event 给
   `partial_p P_p(A)=sum_z P_p(z pivotal for A)`。取 `A={x<->y}` 就得到普通两点连接概率的
   精确导数表示，但没有自动给该 `z` 和的空间尾界。
2. Lemma 28 对 `A_plus intersect A_minus` 给“正 pivotal 概率减负 pivotal 概率”的 signed
   Russo 公式。它说明非单调臂事件必须保留符号，与 #536 的双通道语义一致；它并不等于
   式(3)，也不包含 `g/q/E` 权重。
3. Theorem 26 在 `n<=N<=L(p)` 内证明临界、近临界及相应非齐次 product measure 的固定臂
   概率一致可比；Proposition 32 给三角格

\[
 |p-p_c|\,L(p)^2\pi_4(L(p))\asymp1.
 \tag{13}
\]

论文主体为三角格。其末节说明 RSW、arm separation、characteristic-length 的方法可推广到
square lattice 的 primal/matching 配对，但精确臂指数不能随之移植。超过 `L(p)` 后，固定臂
概率的临界可比性会失败。

### 2.6 near-critical pivotal coupling：固定 lambda 可用，漂移 root 不自动可用

[Garban--Pete--Schramm, *The scaling limits of near-critical and dynamical percolation*,
arXiv:1305.5526](https://arxiv.org/abs/1305.5526) 在三角格采用

\[
 r(\eta)=\frac{\eta^2}{\alpha_4^\eta(\eta,1)},
 \qquad p=p_c+\lambda r(\eta)
 \tag{14}
\]

的 near-critical 尺度。Theorems 1.4--1.5 证明固定 `lambda` 的切片和耦合过程收敛；构造通过
对 pivotal measure 的 Poisson 标记完成。Proposition 11.6 对每个固定 `lambda` 给近临界
四臂 separation/quasi-multiplicativity，并存在与 `lambda` 无关的 `epsilon>0` 使

\[
 \mathbb P_\lambda(A_6(r,R))\le C_{R,\lambda}(r/R)^{2+\epsilon}.
 \tag{15}
\]

常数仍依赖 `lambda`。所以(15)可以把 fused-block 证明延伸到 bounded-`lambda` 窗口，不能
覆盖一个随 `L` 发散的 effective `lambda_L`。Matching One 的 pooled root 若要引用这条路线，
至少必须另证

\[
 \Lambda_L:=|p_L-p_c|L^2\pi_4^{p_c}(L)=O(1)
 \quad\text{或}\quad L\lesssim L(p_L).
 \tag{16}
\]

### 2.7 spectral sample：给平方质量和聚类，不给 cross-observable 符号

[Garban--Pete--Schramm, *The Fourier spectrum of critical percolation*,
arXiv:0803.3750](https://arxiv.org/abs/0803.3750) 对 Boolean `+/-1` crossing
observable 证明

\[
 \mathbb E|\mathscr S|=\mathbb E|\mathscr P|,
 \qquad
 \mathbb E|\mathscr S|^2=\mathbb E|\mathscr P|^2,
 \tag{17}
\]

甚至单点、双点属于 spectral sample 与 pivotal set 的概率相等；Theorem 1.1 给 crossing
spectral sample 的 sharp lower-tail control。Remark 4.6 指出“内盘有 pivotals、外环没有
pivotal”强制六臂，核心衰减正是式(9)。

这些结果针对一个 Boolean crossing 的 `hat f(S)^2`，所以保留质量而丢失符号。#536 在
`p=1/2` 的 Fourier--Walsh 基下恰好展示了缺口。若
`F=sum_S hat F(S) chi_S`，则对固定 `z`

\[
 \mathbb E[(G_z^{\rm mid}-\mathbb EG)\Delta_zO]
 =2\sum_{\substack{S\ne\varnothing\\z\notin S}}
      \widehat G(S)\widehat O(S\cup\{z\}),
 \tag{18}
\]

\[
 \mathbb E[(O_z^{\rm mid}-\mathbb EO)\Delta_zG]
 =2\sum_{\substack{S\ne\varnothing\\z\notin S}}
      \widehat O(S)\widehat G(S\cup\{z\}).
 \tag{19}
\]

所以 midpoint 两项是**两个不同 observable、相邻 Fourier 层之间的交叉重叠**，不是某个
spectral sample 的平方质量。Schramm--Steif Theorem 1.8 的 revealment bound 可控制单一
函数固定 Fourier 层的平方和；要用于(18)--(19)，还需为 dyadic pair source 与 torus
`q/E` 构造共同的低 revealment 算法并保留空间位置。square-site `p_c!=1/2` 还要先建立
相应的 biased-Fourier 版本；算术 midpoint 不是 `p`-conditional mean。

## 3. 哪些结果直接约束 #537，哪些只约束 raw kernel

| 对象 | 文献能直接给什么 | 仍缺什么 |
|---|---|---|
| raw `E|g_xy|` | square-site 双四臂加(5)给绝对可和 | 已闭合，不应再扩 raw 距离网格 |
| endpoint `z=x,y` | endpoint 恒等式降回 raw tail | 已闭合，不是 thermal gate 主体 |
| well-separated kernel-pivot | 三个互不相交圆盘内各一四臂，给 `C pi4(cR)^3` | 对 `y,z` 双重体积积分不够强；还缺共同 landing/连接增益或 signed cancellation |
| fused/collision kernel-pivot | 一旦证明尺度 `r--R` 上有六臂，(9)/(15)使融合尺度和可和 | 必须先完成 collision-tree support；六臂不能逐点假定 |
| observable/rank-pivot | 若 `Delta q/E` 有有限 quad-pivotal cover，可调用 pivotal-measure moments | 目前没有 torus homology cover，也没有与 centered `g` 的联合矩 |
| near-critical extension | `R<=L(p)` 或 bounded `lambda` 内 arm bounds 稳定 | pooled root 是否满足(16)尚未证明 |
| signed midpoint | generalized Russo 允许正负 pivotal；Fourier 形式为(18)--(19) | 无现成 cross-observable、adjacent-level、spatially localized bound |
| original-U `J2` | raw 项可由已有 tail 控制 | 仍须组合 `q/E` 两个通道、`P4` 几何投影及 `D,Y_p,Y_pp,M_pp` 尺度；文献不决定 N25 通道符号 |

## 4. 为什么“再加一个四臂”仍然不够

令 `x,y,z` 两两相距 `Theta(R)`。现有确定性 support 至多先给

\[
 \mathbb P(\Delta_zg_{xy}\ne0)\le C\pi_4(cR)^3.
 \tag{20}
\]

一个 dyadic block 中 `y`、`z` 各有 `O(R^2)` 个位置，因此

\[
 \sum_{y\sim R}\sum_{z\sim R}
 \mathbb E|\Delta_zg_{xy}|
 \lesssim R^4\pi_4(cR)^3.
 \tag{21}
\]

若形式上写 `pi4(R)=R^(-alpha4+o(1))`，式(21)只有在 `alpha4>4/3` 时衰减。
square-site 文献只证明 `alpha4>1`；更强的是，三角格精确 `alpha4=5/4` 代入后仍为

\[
 R^{1/4+o(1)},
 \tag{22}
\]

而不是可和尾。这不证明真实绝对和发散；它证明**三处局部四臂的乘积上界不足以证明收敛**。
需要从“两条 shared components 必须共同连接三个 landing 区”的几何中再提取一个全局因子，
或利用 midpoint 中心化得到 signed cancellation。

六臂结果解决的是 `r<<R` 的尺度比求和。它不解决 `r=Theta(R)` 的顶层 bulk，也不控制
`Delta_zg=0` 但 `Delta_zq/E!=0` 的 topological pivotal。

## 5. 两个可以真正判决机制的 dyadic 命题

令

\[
 A_R(x)=\{y:R<d_L(x,y)\le2R\},
\]

并只在 exact `p_c`，或在已经用(16)认证的 declared near-critical window `W_L` 中取上确界。
以下两式均要求 `L>=4R`，并分别对 `O=q,E` 成立。

### 命题 A：强的 absolute-variation 尾界

\[
 \boxed{
 \mathcal A_O(R):=
 \sup_{L,p\in W_L,x}
 \sum_{y\in A_R(x)}\sum_{z\in V_L}
 \left(
   \mathbb E_p|X^{\rm obs}_{xyz}(O)|
  +\mathbb E_p|X^{\rm ker}_{xyz}(O)|
 \right)
 \le C R^{-\delta}}
 \tag{23}
\]

对某个 `delta>0`。这是无抵消版本；成立即给两个 midpoint 通道的统一绝对 dyadic 可和，
也自动控制 `sum |E X|`。它比关闭 original-U 所必需的条件更强，因而若存在固定 carrier 的
同号下界也可能被严格否定。

一条可执行的证明路线是：

1. endpoint 项用已完成 raw tail 删除；
2. 把剩余 `(x,y,z)` 分成 well-separated、`x/z` 或 `y/z` fused、`x/y` fused、far merger；
3. fused tree 用(9)/(15)在尺度比上求和；
4. well-separated tree 必须证明一个比(20)多出至少 `R^(-1/4-epsilon)` 的共同
   landing/connection 因子（三角格指数只作为强度诊断）；
5. observable tree 必须先证明 `Delta_zq/E` 的有限 quad-pivotal/homology-arm cover。

现有文献只提供第1和第3步的外部输入，没有第4、5步。因此(23)是**明确待证命题**，不是
文献推论。

### 命题 B：逐通道 signed dyadic 尾界

若(23)过强，original-U 实际只需更弱但仍可判决的

\[
 \boxed{
 \mathcal S_O(R):=
 \sup_{L,p\in W_L,x}
 \left[
 \left|\sum_{y\in A_R(x)}\sum_z\mathbb E_pX^{\rm obs}_{xyz}(O)\right|
 +
 \left|\sum_{y\in A_R(x)}\sum_z\mathbb E_pX^{\rm ker}_{xyz}(O)\right|
 \right]
 \le C R^{-\delta}}
 \tag{24}
\]

对某个 `delta>0`。绝对值放在每个完整 channel 的 `y,z` 和之外，允许通道内部由 midpoint
中心化产生抵消，但不允许 observable 与 kernel 两列互相掩盖。它正对应 #536 需要比较的
两个机制，不会把 `abs(sum Delta g)` 冒充 `sum abs(Delta g)`。

式(18)--(19)提示一条新的证明路线：对 dyadic pair source 建立 **adjacent-level
cross-spectral/revealment inequality**，或用 multi-arm IIC/separation coupling 证明在固定
landing sigma-field 后 centered midpoint weight 的条件均值多衰减一个幂。为关闭(21)，这个
额外幂在三角格诊断下至少要超过 `1/4`。当前 GPS/Schramm--Steif 定理只控制单一 crossing
的平方谱或正 pivotal mass，没有这条 cross-observable signed 估计。

若(23)失败而(24)成立，结论应是“thermal tail 由结构性抵消关闭”，不是“没有 pivotal”。
若(24)也失败，下一结果必须在预先命名的 three-site carrier（two-bridge factorized、
`s`-transition/remote merger、或 `g` 稳定但 `q/E` 改变的 topological pivotal）上给一个固定
符号的 dyadic lower mechanism；不再添加自由 descriptor 或窗口。

## 6. 真正剩余的五个引理

1. **square-site/matching collision tree。** 把 `Delta_zg_xy!=0` 分解成有限个 landing tree，
   并证明每个 fused tree 确实含可调用(8)/(15)型六臂 annulus。
2. **torus `q/E` pivotal cover。** 把 `Delta_zq`、`Delta_zE` 覆盖成固定数目的
   primal/matching homology-arm 或 quad-pivotal 事件；4x4 的 kernel-preserving rank pivot
   说明这一步不能由 `Delta_zg` support 代替。
3. **bulk three-centre gain 或 signed cancellation。** 要么补强(20)的共同连接概率以证明
   (23)，要么证明(18)--(19)的 adjacent-level cross-spectral 衰减以证明(24)。
4. **root-window certificate。** 对每个将进入尺度判决的 pooled root 保存
   `Lambda_L=|p_L-p_c|L^2 pi4(L)` 或证明 `L<=C L(p_L)`；否则 near-critical 定理只能停在
   exact `p_c`。
5. **original-U baseline jets。** 即便(23)/(24)成立，仍要独立控制
   `D,Y_p,Y_pp,M_pp` 与面积归一化，才能把两个 `q/E` 尾界送入完整 `J2/A_N`；任何臂定理
   都不单独识别 H4、Jordan 或唯一连续场。

因此，文献给出的最强仓库判决不是“thermal gate 已关闭”，而是：**raw gate 已闭合；
fused-scale 技术已具备；唯一未决主门是 bulk/topological 的 signed 或 absolute 三位置传输。**
PR #509 `2785e3bb` / #537 新给出的 projected landing-matrix minor 是这套文献工具之前的
有限判伪层：先决定 leading three-packet 是否真是 rank-one thermal coordinate，再只为存活的
non-rank-one signed functional 或 four-packet remainder 证明相应尾界。

## 使用的第一手 arXiv 论文

- [van den Berg--Nolin, arXiv:2008.01606](https://arxiv.org/abs/2008.01606)
- [Smirnov--Werner, arXiv:math/0109120](https://arxiv.org/abs/math/0109120)
- [Schramm--Steif, arXiv:math/0504586](https://arxiv.org/abs/math/0504586)
- [Garban--Pete--Schramm, arXiv:1008.1378](https://arxiv.org/abs/1008.1378)
- [Nolin, arXiv:0711.4948](https://arxiv.org/abs/0711.4948)
- [Garban--Pete--Schramm, arXiv:1305.5526](https://arxiv.org/abs/1305.5526)
- [Garban--Pete--Schramm, arXiv:0803.3750](https://arxiv.org/abs/0803.3750)

本卡使用 `literature-search-arxiv` 工作流下载并精读上述论文；arXiv API 请求遵守限速，论文
具体复用仍应逐篇检查其许可证。
