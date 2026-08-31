# 用完整 q/E 生产曲线锚定 100k 簇源标记：一个可直接执行的两阶段估计

**计算已完成：这项两阶段方案没有带来稳定降噪。** [实际报告](../results/norm4-source-two-phase/REPORT.md)与[完整协方差](../results/norm4-source-two-phase/latest.json)记录一次2.52秒、零回放的计算：六N原始/两阶段SE比为0.898–1.008，补集基准的方差贡献不足总量的0.002%。旧1.9B/1B的无标记曲线已经足够精确；当前限制在簇源标记本身。下面保留执行前的推导与定义。后续直接增加旧生产端点的源标记，不把此方案再次列为待算任务。

## 1. 现有字段足够做什么

六 N 的旧生产来源固定为 `8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc:results/server-20260829/P154-norm4-production/raw/`。N65/85/130/170 使用 `n{N}_1900m.hist.csv`，N260/340 使用 `n{N}_1b.hist.csv`。每个原始 batch、方向都保留 K_minus/K_plus 边际直方图。令其累计概率为 F_minus(K)、F_plus(K)，则

```text
m(K)=E[q|K] = -1 + F_minus(K) + F_plus(K),
e(K)=E[E|K] =  1 - F_minus(K) + F_plus(K),   E=q².
```

所以完整生产的 m(p)、e(p) 及所需热导数都可从边际直方图重建，不需要原配置或 K_minus/K_plus 的新联合直方图。新数据 `results/norm4-source-thermal/raw/n{N}.csv` 已有每个 batch、方向和 K 的 `samples,sum_q,sum_e,sum_s,sum_qs,sum_es`。这里取 **bulk source S=s=CB+CW**；既有读入函数会除 N，复用时须明确单位。最后 `Udot_bulk=N*Udot_density`，原链上共同微观 fugacity 的比较用 bulk 版本。

这些字段还能精确拆出三种 q 状态的样本数和源和：

```text
n_±=(sum_e ± sum_q)/2,       n_0=samples-sum_e,
s_±=(sum_es ± sum_qs)/2,     s_0=sum_s-sum_es.
```

本方案不需要 `s²`。不建议直接在每个稀有 (K,q) 格子分别估计一个自由源参数；现有 batch 足够拟合下面的小基底，并直接估计最终读数的协方差。

## 2. 源分解与严格的 clock 零

对一个 N 的两个方向 g，固定一次系数，写成恒等分解

```text
S_g = a_g + b_N K + c_g q + d_g E + R_g.
```

这是同一个 S 的估计分解，不是新增物理源或独立证据。常数 a_g 没有归一化测度响应。**b_N 必须在该 N 的两方向相同**：`exp(t b_N K)` 只把 Bernoulli 参数变成 `logit^{-1}(logit(p)+t b_N)`。沿 pooled matching root 移动时，分子和分母的 Jacobian 相消，因此 `Udot[b_N K]=0`。逐方向分别拟合的 b_g 不具有这个共同坐标零；不能把它们也删除。

更明确地，记 `P=(first-second)/delta_cos4`、`D=bar(m')`、`B=P[e']`、`H=P[e'']`、`T=bar(m'')`，全部在未扰动 pooled root 求值。对任意固定源 X，设 `Jq=Cov(q,X)`、`JE=Cov(E,X)`，则

```text
L_Baseline[X] = N^(13/8)/2 ×
  { P[JE']/D - H bar(Jq)/D² - B bar(Jq')/D² + B T bar(Jq)/D³ }.
```

这就是 Udot，给定 baseline 后对源线性。对 X=b_N K，`Jq=b_N p(1-p)m'`、`JE=b_N p(1-p)e'`；代入上式严格为零。正确的旧同样本 scorer 原本就应满足该零，因此**仅减去 clock 不会凭空带来新的统计信息**。主要降噪机会来自更精确的 baseline 和 q-only 分量。

对 `f_g(q)=c_g q+d_g E`，三态代数 `qE=q,E²=E` 给出

```text
Jq[f] = c_g(e-m²) + d_g m(1-e),
JE[f] = c_g m(1-e) + d_g e(1-e).
```

对 p 求导只需完整生产的 m,e,m',e'；这里 c_g,d_g 在热导数中保持固定。q-only 响应因此可以直接用大档案计算；c_g,d_g 的估计误差仍来自标记子集。

## 3. 推荐估计：先去掉嵌套，再交叉拟合残差

设 A 是已标记的 n=100000 个排列，F 是 M=1.9B 或 1B 的对应完整生产块。A 恰为 F 的首段。对每个方向、K，用整数总和精确相减：

```text
sum_q(C)=sum_q(F)-sum_q(A),
sum_e(C)=sum_e(F)-sum_e(A),        C=F\A, |C|=M-n.
```

这里 C 与 A 使用不相交的原 counter，在通常 counter-PRNG 独立键假设下是独立样本块；同 seed 的四个 cyclic N 都减去同一段。几乎全部 baseline 精度得以保留，又不必假装“全样本与其子集独立”。从 C 的精确 Binomial 聚合重建 baseline、pooled root、D、B、H、T；不要直接把旧 100k scorer 的 q/E 均值替换掉而保留不匹配的协方差。

系数拟合可固定在 `p_star=0.59274605079`，采用预先固定的五折 `fold=batch_id mod 5`。在每个训练集，对两方向共同拟合上述源分解：方向专属的 a_g,c_g,d_g，加一个 pair-common b_N；用 `Binomial(N,p_star)` 给 K 加权。设计矩阵只含方向指标、K、q、E。`X'X` 中的 q²、qE、E² 可化成 E、q、E；`X'S` 只需 s、Ks、qs、Es，所以现有字段全部足够。系数只拟合一次固定 p_star 的源分解，不随评估 p 重新拟合；不额外加入系数热导数。

对每个 held-out fold h，用训练系数构造它的逐 K 残差和：

```text
sum_R  = sum_s  - a_g count - b_N K count - c_g sum_q - d_g sum_e,
sum_qR = sum_qs - a_g sum_q - b_N K sum_q - c_g sum_e - d_g sum_q,
sum_ER = sum_es - a_g sum_e - b_N K sum_e - c_g sum_q - d_g sum_e.
```

对这些残差和做精确 Binomial 聚合，得到 r,r_q,r_E 及热导数。使用补集 baseline 中心化：`Jq[R]=r_q-m_C r`、`JE[R]=r_E-e_C r`，求导也使用同一 m_C,e_C。该折的估计为

```text
Udot_h = L_C[f_h] + L_C[R_h],
Udot_two_phase = 按 held-out 样本数加权平均各折 Udot_h.
```

它等价于对 `S,qS,ES` 使用控制变量／generalized-difference 原始矩估计。例如

```text
mean(S)   ≈ a+b Np+c m_C+d e_C+mean_A(R),
mean(qS)  ≈ (a+d)m_C+b mean_C(Kq)+c e_C+mean_A(qR),
mean(ES)  ≈ (a+d)e_C+b mean_C(KE)+c m_C+mean_A(ER).
```

这些公式须逐折使用训练系数与 held-out 残差。共同 clock 项在最终 L 中直接消零，不靠四个大数数值相减。带不惩罚的 K 列、常数列的线性拟合对 `S→S+constant+bK` 等变，纯 clock 源给出零残差和零 Udot；不能随意删一项后仍声称同一源。

## 4. 系数和 baseline 的不确定性怎样进入

不能把拟合系数固定后只报告 held-out 的残差 SE。执行现有三组 source delete-one：四个 cyclic N 同删一个 1000-counter batch；N260、N340 各自一组。每次删除后**重拟合所有受影响训练折**、重算 held-out 加权和及完整链读数，补集 C 保持固定。这样保留系数误差、折之间的共享训练相关性和源分量抵消；五折不是五个独立实验。

baseline C 的误差另用原生产 batch 传播，源 A 固定，重求 pooled root 和全部系数函数。删除 A 只改变原生产第 0 批：其剩余样本数成为 18.9M（cyclic）或 9.9M（endpoint），其余批仍为 19M 或 10M。使用真实批权重的 block-influence／加权 delete-group 方差，不能仍把这一补集当 100 个等大小批次。四个 cyclic N 的原批保持共同删除；不同 seed 的两 endpoint 分组。两阶段输入已不相交，其一阶协方差贡献可相加；非线性比值的有限样本偏差仍应与误差分开描述。

若坚持直接用嵌套 F 而不减 A，对固定控制系数的标量线性化 `Y_A+beta(X_F-X_A)`，写 `R=Y-beta X,Cv=beta X`，其方差是

```text
Var = Var(R)/n + {Var(Cv)+2 Cov(R,Cv)}/M.
```

一般不能写成独立的 `Var(R)/n+Var(Cv)/M`；只有相应读数的残差与控制正交时交叉项才消失。拟合 S 的最小二乘并不保证最终 Udot 影响函数也正交。因此补集方案更容易直接、诚实地实现。

## 5. 已执行的 compute 与后续

已经完成一次**零回放的两阶段 source 分析**：读取六 N source CSV 与旧 threshold histograms，建立补集 anchors，按上述固定五折分解，输出六个 `Udot_bulk`、原 q2/Jordan source-extension 残差及 generator-drift determinant 的联合结果。报告并列保留 raw、只锚定 baseline、完整两阶段估计及共同协方差。仅替换基准打断了原同样本估计中的误差抵消，交叉拟合把误差恢复到原有量级，却没有系统改善。

这三个视图复用同一数据，不合并 p 值。若 R 仍主导 endpoint 不确定性，完整 q/E 档案已无法补上缺失的簇源信息；这一次计算会明确下一项需要的源标记，而不是继续增加准备工具。无需新模拟、服务器、s² 字段或一套新框架。
