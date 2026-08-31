# Twist 重建是否绕过 m64 稀有扇区：有界独立审查

**结论：现有正 partition 恒等式提供了准确的观测量表示，但没有给出绕过稀有 rank1 的高效估计器。** 朴素独立 partition 差分通常比直接稀有事件估计更病态；自然的共享构型抵消恰好还原 rank1 指示量。可以证明这两个结论，且可明确写出改善所需的额外输入。这里结束本次可行性审查，不启动新采样、枚举或 coupling 点。

完整阅读来源：`2690f665` 的 `closed-source-poisson-double-scaling.md`、`closed-source-oblique-twist-comparison.md`、`closed-source-pooled-sector-odds-bound.md`；`3dc47674` 的 `closed-source-hypergraph-rc-twist-projection.md`。下列 m64 数字只由已完成固定包 `p337-finite-law-window-20260831/{results/latest.json,review.json}` 的有理区间中点作算术换算，不重新评分。

## 1. 先区分“五个”与 m64

原定理的五个非负 partition 严格针对 **m=2**：四个 `Z_alpha` 加全秩 twist I。一般整数 m 使用 `m²` 个平移 twist 加 I；**m64 的 literal 构造是4097个 partition**。若能直接计算聚合量 `R=sum_{alpha≠0}Z_alpha`，可以只保留 T、R、I 三个输入，但现有恒等式没有证明这种聚合计算便宜，也没有证明足够的几何/模群对称性可把4095项压成固定少数项。改用 m2 的五个 partition 会改变固定物理源，不能用于当前 m64 结果。

令 `L_r` 为未投影 local-colour law 的 rank-r 正扇区权重，原式为

```text
T=L0+L1+L2,      R=(m²−1)L0+(m−1)L1,      I=L0,
D=T+R=m²L0+mL1+L2.
```

star 的总归一化为 D，drop 的归一化为 T；两律各自在自己的 pooled root 求值。统一记 `c=m²−1`、`J=R−cI=(m−1)L1`，则

```text
P1_star = [m/(m−1)]*J/D,
P1_drop = [1/(m−1)]*J/T.
```

## 2. 原 U 还需要两层差分

用点表示 `h*d/dh`，例如 `Tdot=sum K*weight`。对 star 令 `B=D, gamma=m/(m−1)`；对 drop 令 `B=T, gamma=1/(m−1)`。每个几何分别有

```text
Jdot=Rdot−c*Idot,
mu=Bdot/B,
C1=Cov(K,I_rank1)=gamma*(Jdot/B−J*Bdot/B²)
                  =P1*(mu_rank1−mu).
```

q 的未归一化 numerator 分别为

```text
Qstar = T−R/(m−1)−m(m−1)I,
Qdrop = T−R/(m−1)+(m−1)I.
```

因此 `q=Q/B`、`Cq=Qdot/B−Q*Bdot/B²`。先解 `mean_g q=0` 的共同 root，再用同一个原目标

```text
U/A_N = −(C1_axis−C1_tilt)/(DeltaCos4*mean_g Cq).
```

共同 `1/h` 精确取消。需要控制的差分依次是：稀有质量 `J`、热中心化 `Jdot−mu*J`、两几何的 C1 差。非负 T/R/I 不会让这些差分变为非负。若 root 也由估计值求解，另有 `delta log(h)≈−delta mean(q)/mean(Cq)`，并带入 U 的热导数；m64 当前 root/分母已有精确包围，所以即使把它们作为已知量，下面的分子困难也仍成立。

严格说，下文有条件方差界允许把真实共同 root、总 normalizer、mu 和 pooled 分母都作为 oracle 给定，只考察指定估计器类别中的分子；它们不是所有算法的下界。若实际算法没有这些 oracle，未知 root 的 U 传递还需二阶热矩 `K²I1`、`K²q` 及对应累积量来控制 `dU/d log(h)`，不能把“先求一个 root”当成没有误差或成本的操作。

## 3. 一个严格条件数与朴素独立估计的方差结论

J 对 R、I 的分量相对误差条件数为

```text
kappa_sub=(R+cI)/J = 1+2(m+1)L0/L1.
```

若 `|Rhat−R|≤eta*R`、`|Ihat−I|≤eta*I`，则严格有 `|Jhat−J|/J≤eta*kappa_sub`，且该最坏情况可达到。加上总归一化 B 的相对误差，P1 的相对误差至多 `eta*(kappa_sub+1)/(1−eta)`。在 star 的概率单位，`kappa_sub=1+2(m+1)P0/(m P1)`；drop 则为 `1+2(m+1)P0/P1`。

固定 m64 已有值给出：

| law/geometry | kappa_sub，约 | 热中心化条件数 `(mu1+mu)/abs(mu1−mu)`，约 |
|---|---:|---:|
| star/axis | 9.10067e14 | 1.65186e5 |
| star/tilted | 4.78140e19 | 21.9305 |
| drop/axis | 1.56308e16 | 2.93864 |
| drop/tilted | 2.58723e21 | 11.6821 |

star/axis 的 `mu1−mu≈0.0001513826`；因此精确得到 P1 也没有自动精确得到其热导数。两几何差需要保留共享协方差，不能把这些条件数当独立误差直接相乘。

可给出明确但**有条件**的方差下界：若 Rhat、Ihat 是独立无偏估计，各用 n 次独立样本，且其单样本相对方差均至少 v>0，则

```text
Var(Jhat)/J² = [Var(Rhat)+c² Var(Ihat)]/J²
             ≥ v*(R²+c²I²)/(n*J²)
             ≥ v*kappa_sub²/(2n).
```

即便分母 B 已知，此界也直接传给 P1。相对 RMS≤rho 必须满足 `n≥v*kappa_sub²/(2rho²)`。在四格中，`kappa_sub²/2` 从 `4.14e29` 到 `3.35e42`。**这里没有证明实际 partition 算法的 v 有正下界**，所以这些数字不是所需机器数、墙钟预算或普遍复杂度下界；它们是对“相互独立、有常数量级相对噪声的 partition 比值算法”的否定。

## 4. 自然的精确耦合恰好回到稀有事件

在同一个未投影占用构型上，非零平移 twists 的兼容数是

```text
A_nonzero=(m²−1)*1{r=0}+(m−1)*1{r=1},
A_fullrank=1{r=0}.
```

于是逐构型就有 `A_nonzero−c*A_fullrank=(m−1)*1{r=1}`；再乘 K，导数残差也精确化成 `(m−1)K*1{r=1}`。这消除了独立估计两个大数造成的 `1/P1²` 量级损失，但没有制造频繁出现的 rank1 信息。若在投影后的联合 twist ensemble 采样，占用边际正是 star，仍是其原 rank1 事件。

对此有一个不依赖细节的单样本界。任意估计变量 Y 若只在事件 A 上非零，`P(A)=p1`，且 `E Y=theta≠0`，由 Cauchy–Schwarz

```text
E Y² ≥ theta²/p1,
Var(Y)/theta² ≥ (1−p1)/p1.
```

对独立重复的平均再除以 n。令 `Y=(K−mu)1{r=1}`，即使给定精确 mu，其 C1 相对方差也受此界约束。两个几何的共同样本差仅在两者 rank1 的并集上非零时，用并集概率 `pA≤P1_axis+P1_tilt` 得到同样的界。因此，共享 twists 若只是把代数抵消做到逐样本精确，不能突破这类事件支持估计器的稀有阶。它不约束具有非事件支持输出的条件期望、确定性积分或非原测度估计器。

## 5. 哪一种新增输入才可能改善

两种具体输入足够有意义，但现有四份 note 尚未提供它们的可计算性/误差保证。

**受限扇区的桥接或重要抽样。** 需要能计算归一化权重，且对 rank1 条件律 `mu1=mu(.|r=1)` 有重叠控制的 proposal nu。若目标总归一化已知，`Y=(dmu/dnu)1{r=1}` 的相对方差精确为 `chi²(mu1||nu)`；n 次独立平均再除以 n。它可摆脱 `1/P1`，但必须提供受限/总 partition 之比的可估计桥及其重叠界。只说“在 rank1 扇区采样”不能恢复未知的 P1，也不能算出原 C1。

一个最短的信息反例：固定三个扇区内部的所有条件分布，把其未归一化总权重分别设为 `(Z0,Z1,Z2)=(1,1,1)` 或 `(1,2,1)`。任何数量的各扇区条件样本，其分布完全相同；P1 却分别为 `1/3` 与 `1/2`，且两者都有 q=0。条件热矩相同也不能修复缺失常数：一般 `C1=P1*(mu1−mu)` 随混合权重改变。这只是条件采样不能识别 sector normalizer 的信息反例，不是对固定源增加可调系数或提出新物理候选。

例如固定各扇区的 K 为0、1、4，并让共同热因子为h^K；在h=1，两套权重仍有相同的各扇区条件law和pooled零，但C1分别为−2/9与−1/4。缺少相对normalizer会直接改变需要估计的热协方差。

**Rao–Blackwell 条件积分。** 给定背景 G，需在同一带源 law 下可计算
`a(G)=E[I1|G]`、`b(G)=E[K I1|G]`，并同时得到 `E[K|G]`、`E[q|G]`、`E[Kq|G]` 以形成原 root、C1 和分母；不能用另一个无源 face law 代替。一般只知道
`0≤Var(a)≤P1(1−P1)`，两端均可达到，故“用了 RB”本身没有改善保证。若能证明 `a≤C*P1`，则 `Var(a)/P1²≤C−1`；热 numerator 还需控制 `E[(b−mu*a)²]` 相对于 `C1²`，否则 star/axis 的热消去仍可主导。当前全密度面核给出正确条件积分对象，但固定 B 行的闭式不能替代一般 B 背景下这些量的廉价、可控计算。

这是缺失假设的实质：**共享 partition 估计的协方差，或条件/桥接权重的二阶矩控制**。相同非负 partition 均值可配独立噪声，也可配完全共享噪声；甚至可以确定性精确计算、方差为零。所以不可能仅凭 partition 为正及重建恒等式推出普遍方差下界或高效性。

## 6. 其余三份 note 没有填补该算法缺口

- Poisson double-scaling 定理同时证明 rank1 消失及 `Cov(K,q)/N→1/2`，确实控制原 U 的分母和上界；它不把微小有符号差分变成可估计的大量，也不提供相对误差算法。
- pooled-sector-odds 的 Xi 界控制的是各几何内部的 q 方差，清楚保留 restricted-sector mismatch；总 partition 或压力相同不能替代该信息。在本次 m64 分母已严格正，困难仍在稀有分子。
- oblique-twist 的 order25 twist 是**空间接缝**，不是上述颜色平移 twists。其 `exp(−Delta_k)=Z_twist/Z_rect` 仍可在 PSD/可交换/有限阶前提下指数小；该 note 已给明确矩阵反例。它不能被当作免费的 partition 桥或有界方差证明。

**本轮停止结论：** 正 twist 表示值得保留为精确代数/算法入口；当前资料不足以将 m64 普通稀有抽样的停线决定改成“已有高效替代”。下一次若要启动计算，应先交付上述一种具体耦合或条件积分的二阶矩/重叠控制；本审查不据此自动建立新实验。
