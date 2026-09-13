# [P0 最高优先级建议] 用原 U 的完整影响函数决定采样，而不是追求等量 rank-one 事件

**日期：2026-08-31。类型：精确方法推导 + 已公布整数表上的事后可达性诊断。**
P0 标记表示本建议的审阅优先级，不升级科学证据，不创建新生产合同，也不恢复 #154/#334/F4 已停线的任务。本包不合并其他研究分支、不修改既有预测或结果；新随机样本、格点枚举、云作业均为 0。

## 判决与当前前沿

已阅读 #509 的 `f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78` 交付。原 U 的 m64 普通 iid 预算、twist 差分病态、条件扇区样本缺少 normalizer、全孔面核以及固定 m 的受限 odds 缺口均已有结果；本包不重做其首次构造。

**新增判断：star 的困难不止是 rank-one 稀有。** 把共同根、均值和热分母的估计全部纳入后，在指定的普通单 proposal、自归一化重要抽样（SNIS）plug-in 估计器类内，即使允许 oracle 最优 proposal，其原 U 的渐近方差等效负担仍约为每几何 `9.84e10`。仅用三个 rank 扇区常数做偏置已经达到 `1.066e11`，继续精调同类 proposal 的空间有限。drop 的对应 oracle 数字约为 52.6/54.4，但这不是“55 个实际样本即可”，更不允许把 drop 替换为 star 的物理目标。

这不是所有算法下界。精确条件积分、已知均值控制变量、跨几何耦合、不同估计方程或多 proposal 结构不属于上述包络。下一步应改变进入原 U 的有符号抵消方式，而不是仅把稀有扇区变成常见扇区。

## 1. 先固定同一个原始观测

每几何的有限正权族为 `pi_g(omega;xi) proportional exp(xi K) w_g(omega)`，`xi=log h`；源及其单位固定。令 `q=r-1`、`E=q^2`，每几何先独立归一化，再用固定权重 `a_g=1/2` 合并。N25 的两个角投影系数为 `c_f=1/Delta`、`c_s=-1/Delta`，`Delta=1152/625`。

在共同根处定义

```text
Q = sum_g a_g <q>_g = 0
C_O,g = Cov_g(O,K)
D = sum_g a_g C_q,g > 0
B = sum_g c_g C_E,g
R = B/D = U/A_N,    A_N=N^(13/8)/2.
```

热参数 Jacobian 在分子、分母中相消。这不是把两几何当作一个混合分布后计算协方差；原 U 的 within-geometry 语义不变。

## 2. 完整原 U 的影响函数

令 `mu_g=<K>_g`、`O_bar=<O>_g`，并定义

```text
psi_O,g = (O-O_bar)(K-mu_g)-C_O,g
D_xi = sum_g a_g cumulant_g(q,K,K)
B_xi = sum_g c_g cumulant_g(E,K,K)
R_xi = (B_xi-R D_xi)/D.
```

对每几何任意 xi 无显式依赖的 infinitesimal weight tilt `exp(s f_g)`，重新求共同根、均值及分母后，精确一阶变分为

```text
phi_g = [c_g psi_E,g - R a_g psi_q,g
         - R_xi a_g (q-<q>_g)] / D

dR/ds = sum_g <phi_g f_g>_g,    <phi_g>_g=0.
```

证明：固定 xi 时，均值变分是 `Cov(O,f)`；协方差变分的中心化核是 `psi_O`。共同根移动为 `xi_dot=-sum a_g Cov(q,f_g)/D`。对 `R=B/D` 用商法则，再加 `R_xi xi_dot`，即得上式。它是原观测的精确有限概率微分；不是新连续场或假设的物理源。

一个必须通过的控制是

```text
sum_g <phi_g (K-mu_g)>_g = 0.
```

因为共同 K 倾斜只移动热坐标，重新求根之后原 U 不变。少掉 root 项的算法通常不满足这个零。根/分母误差并非独立误差：必须保留其与分子的协方差，也不能声称重新估计根必然增加方差。

在固定 rank 扇区内，q 和 E 都是常数，所以 `phi_g` 是 K 的仿射函数。对这个明确的 U 估计器，无需为了构造 proposal 再增加 contact/clock 描述符目录。此结论不自动适用于 U_t 等更高源导数的估计器。

## 3. 三个 rank 偏置就能构成合法的重加权入口

取固定正数 `alpha_g,j`，采样

```text
nu_g(omega) proportional pi_g(omega) alpha_g,r(omega)
W_g(omega)=1/alpha_g,r(omega)
<O>_pi = <W O>_nu / <W>_nu.
```

在基准 xi0 采样时，任意邻近 xi 使用权重 `exp[(xi-xi0)K]/alpha_r`；每几何分别归一化，重建 q/E 与热矩，最后重求同一个 pooled root 和原 U。alpha 是被完全消去的采样偏置，不是新的物理耦合。

这避免了“只在各扇区内采样而不知道它们相对 normalizer”的信息缺口：需要采样的是一个跨扇区的共同正权系综，不是三个互不连通的条件链。偏置固定后，理论上可用对称单点提案和 `min(1,w_new alpha_rnew/(w_old alpha_rold))` 定义 Metropolis 核，但不可由正确平稳律推断良好混合。本包没有运行或交付高效 MCMC。

设 `p_j=P_pi(r=j)`、`M_j=<phi^2|r=j>`。若样本为 iid，plug-in SNIS 的完整一阶渐近方差系数为

```text
v(alpha) = (sum_j p_j alpha_j)(sum_j p_j M_j/alpha_j).
```

它来自 SNIS 影响核 `W phi/<W>_nu`，因此包含目标 normalizer、原 U 根和分母的估计。Cauchy 不等式给出

```text
alpha_j proportional sqrt(M_j)
v_min,sector = (sum_j p_j sqrt(M_j))^2
nu_opt(r=j) proportional p_j sqrt(M_j).
```

所以应平衡的是**原 U 影响函数的贡献**，不是扇区直方图。M_j 是 population 量；当前用完整旧整数表计算仅为 oracle 诊断，不是未来在线自适应已获保证。

若允许任意普通单 proposal、正确支持的 regular SNIS，则相同证明给下确界

```text
v_inf,SNIS = <abs(phi)>_pi^2,
nu_opt proportional pi abs(phi).
```

phi 的零点需要正支持 floor 或取极限。这个下确界只约束指定 estimator architecture；不能升级为所有无偏估计、所有 MCMC、所有条件算法的复杂度下界。

## 4. 有一个无需精确调参的 9/8 保证

若去掉一个共同倍率后，`alpha_j/sqrt(M_j)` 在 `[exp(-eta),exp(eta)]` 内，则

```text
1 <= v(alpha)/v_min,sector <= cosh(eta)^2.
```

证明：用 `omega_j=p_j sqrt(M_j)/sum p sqrt(M)` 化简为 `E_omega[r] E_omega[1/r]`；若 `a<=r<=b`，则 `r+ab/r<=a+b`，故乘积不超过 `(a+b)^2/(4ab)`。

因此把相对理想偏置取为最近的二次幂，误差在 `[1/sqrt(2),sqrt(2)]` 内，有

```text
v_dyadic <= (9/8) v_min,sector.
```

这是 iid 渐近方差保证，不约束改变 alpha 后的 MCMC 自相关时间。若实际 pilot 对 sqrt(M_j) 还有乘性误差，必须把它合并进 eta；不能直接继承 9/8。

本次表上选出的简单偏置为：

| law / geometry | rank0, rank1, rank2 bias |
|---|---|
| star / axis | `(1,2^49,1)` |
| star / tilted | `(1,2^49,1)` |
| drop / axis | `(1,2^47,1/8)` |
| drop / tilted | `(1,2^47,1)` |

这些都是对已看数据的设计诊断，不是新前瞻冻结，不授权在新数据上重新选择。

## 5. 实际旧整数表上的结果

输入逐字取自 `cae9c8997b5994c218bfe060f75656137f745755` 的
`experiments/p337-finite-law-window-20260831/inputs/{axis,tilted}.csv`；不重新枚举任何构型。权重严格保持原合同：`count*h^K*64^(-g+rho*(q+1))`，star 的 rho=0，drop 的 rho=1。这里 g 是既有 source-cost 列，不是另拟合的参数。

使用已发表有理 root bracket 的中点，不细化原生产根；120/170 位运算核对。十进制结果是高精度诊断，不是区间证书。完整结果见 [result.json](result.json)。

为便于比较，定义

```text
n_equiv = 9 (v_axis+v_tilted) / R^2,
```

即两几何各取 n 个独立样本时、对应三倍渐近标准误差的**方差等效数量**。
它不是有限样本 SNR 保证、95%覆盖率、MCMC 步数、机器预算或建议开跑的数量；尤其不能把小于100的渐近数当成有限样本承诺。

| 单 proposal / plug-in 类别 | star 的 n_equiv | drop 的 n_equiv |
|---|---:|---:|
| 原测度 iid，完整根/均值/分母影响 | 1.51803e25 | 2.52472e15 |
| 三扇区等概率 oracle bias | 1.17061e11 | 73.5279 |
| 裸 winding barrier bias：axis 64^9、tilted 64^13 | 8.80427e13 | 34171.8 |
| 最优三扇区常数 bias | **1.06600e11** | **54.4158** |
| 任意 ordinary regular single-proposal SNIS 下确界 | **9.83790e10** | **52.5631** |
| 上述二次幂 bias | 1.07259e11 | 54.9674 |

原仓库的 iid 数是 oracle-centered numerator 估计器的精确区间预算；本表第一行是包含 root/mean/denominator 的 plug-in 渐近值。两者本例数值极近，但不宣称估计器定义完全相同。

star 的理想扇区抽样概率尤其有辨识力：

```text
axis:  (0.2865897870, 0.4268244108, 0.2865858022)
tilted:(0.4999910232, 0.00001340785,0.4999955689).
```

tilted 的主要统计作用包括共同根和分母，不能从“它的 rank1 更稀有”推断应投入更大的 rank1 份额。裸 barrier bias 把 tilted 过度压向 rank1，使已需估计的根/分母信息变得低效。

这也暴露了 star 的 signed thermal cancellation：原分子很小不只是因为 P1 很小。对这个固定 estimator class，最优 sector bias 已接近无限自由 proposal 的下确界；继续同类 proposal 精调不可能带来消除 `1e11` 负担所需的数量级改善。允许别的 estimator architecture 后，这个推论不再适用。

## 6. 下一项值得最高优先级审阅的工作是什么

**不是扩大参数目录，而是检验已有全孔面核是否能降低原 U 的 conditional influence。** 对一个预先指定、按同一带源测度可精确积分的背景 G，定义

```text
phi_RB(G) = E_pi[phi(omega)|G].
```

若构造的是同一组原始矩的精确条件估计，再求同一根和比值，其一阶影响为 phi_RB。相应理想 ordinary SNIS 包络变为

```text
[E_pi abs(E_pi[phi|G])]^2 <= [E_pi abs(phi)]^2.
```

这里的 Jensen 不等式说明**条件积分有机会改变 proposal 单独不能突破的抵消负担**，但它不保证严格改善、更不保证计算便宜。现有 `(K,g,q,count)` 表不足以计算一般面背景的条件核，不能伪造此收益。

具体下一交付应报告：同一个带源全孔核的条件 q/E/K 与必要混合矩；条件 root/分母重建；`E[abs(phi_RB)]` 或完整二阶矩的实际减少；每次条件积分成本。对 MCMC，还必须分析 `W*phi/<W>` 的长程协方差，而不只报告 rank1 自相关、权重 ESS 或平稳律正确。

允许结果是“这个条件核也不改善”并停止；不把先前验证块的最显著新 descriptor 当作救场方案。drop 仅是指定的算法对照，不能替代 star 源或证明 star 的物理机制。

## 7. 核验与复算

```bash
# Python 3.9+; exact synthetic controls use standard library only.
python test_influence.py
# The numerical diagnostic was run with mpmath 1.3.0.
python -m pip install mpmath==1.3.0
python analyze.py --output /tmp/matching-one-u-influence-fresh.json
```

输出路径必须尚不存在；旧 result.json 不被覆盖。脚本不联网、不读取凭据、不生成样本。

- 两输入 SHA256 及每 K 的 `binomial(25,K)` 计数逐项检查。
- 四个有理数测试方法：逐几何中心化与共同 K 零、12 个独立 density directions 的 dual-number 原式核验、重加权归一化/二阶矩、27 个偏置误差比的 Kantorovich 界。
- 原表 120/170 位的所有 budget 字符串一致。它是精度稳定性，不是独立算法或 interval proof。
- 另对两种固定测试 weight direction、两律分别做直接重新归一化和重新求根的中心差分，四次核对误差均小于 `1e-35`。这些极小人工扰动只验证导数代码；没有新增生产耦合点或改变旧合同。
- 没有运行全仓库 CI、蒙特卡洛或新格点枚举；没有证明新尺寸表现、连续 H4 身份或所有算法下界。

## 固定来源与方法归属

1. [原 iid 可达性结果](https://github.com/LightChainr/Matching-One/blob/f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78/experiments/p337-estimator-access-20260831/RESULT.md)、[twist 审查](https://github.com/LightChainr/Matching-One/blob/f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78/notes/p337-twist-estimator-access.md)：本包的直接问题来源，不重复计作独立科学票。
2. [原 U 时钟商](https://github.com/LightChainr/Matching-One/blob/1b0ec15a/notes/p154-original-u-clock-quotient.md)：根移动与共同热坐标零是已有资产，本包把该泛函传播到完整统计影响函数与实际 proposal 包络。
3. [输入 axis](https://github.com/LightChainr/Matching-One/blob/cae9c8997b5994c218bfe060f75656137f745755/experiments/p337-finite-law-window-20260831/inputs/axis.csv)；[输入 tilted](https://github.com/LightChainr/Matching-One/blob/cae9c8997b5994c218bfe060f75656137f745755/experiments/p337-finite-law-window-20260831/inputs/tilted.csv)。SHA256 分别为 `2d23fecc98d276d9ad15ad1867199cd308f0570cb5040ef94eb6b923b4c53458`、`225031e612929ed922ba75c55e76703d59990f5283e7ac39b94f022841798da5`。
4. 最优重要抽样的一般原理不是新发明。Branchini–Elvira, [Towards Adaptive Self-Normalized Importance Samplers](https://arxiv.org/abs/2505.00372v2)，讨论 optimal SNIS proposal 及其未知量问题。
5. Owen, [Zero variance self-normalized importance sampling via estimating equations](https://arxiv.org/abs/2510.00389)，明确显示改变 estimating-equation architecture 可以越过 ordinary ratio SNIS 的非零包络。因此本包刻意不作所有算法下界。

本包新增的是：指定原 U 的完整 influence、三个扇区的目标相关最优偏置、带证明的二次幂稳健界，以及实际固定 N25/m64 旧表上的数值判别。它不另设项目状态入口。
