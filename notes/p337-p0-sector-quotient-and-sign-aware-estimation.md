# [P0 建议] 原 U 的扇区商坐标与先解析扣除、后抽样

**建议级别：P0 / 最高优先级研究评审。** 这里的 P0 标注两项建议的优先级，不把待检验的物理解释升级成定理，也不授权重新启动已停线生产。相关 #337、#275、Draft #509。

本贡献是追加的独立小包，不改主线 ledger、冻结合同、历史结果或已有 Issue 的优先级。输入截止为 `f3ecde7d` 的估计器评审及 `0dda27ba` 的 Q4 trace 交付。只使用已经完成的 N25 整数表；0 新随机样本、0 新构型枚举、0 云作业、0 官方 scorer 调用。

## 两项建议及可交付判据

**P0-A：在同一 pooled 热坐标下核对候选机制的扇区切向，而不从单个分子或 V/U 的显著性认定场。** 采用下文精确坐标，要求具名 Q1 carrier 提供每个几何的两项源切向及其热导数，连同既有归一化/根接口一起预测原 U。现已实际交付 N25 完整 Sstar 的标定：共同热坐标消去后，两个规范不变的响应项为 `+0.301509231369101` 和 `-0.175343867954932`，总和 `+0.126165363414169`。它们是响应坐标的分配，不是两个独立的物理源或场，不替代旧 frozen score。

下一份有意义的理论结果是：把一个事前指定的 carrier 放进**已经固定的 occupation continuation**，给出与竞争解释不同的这组净切向预测。若不同解释只是在同一条 q/E 曲线上重命名同一切向，则这组观测不能识别它们；不得再用新大小的同一标量测试掩盖这一点。Q4 的非零 trace 不自动提供 Q1 的完整 jet。

**P0-B：先做确定性的最小绕环层扣除，再评审一个真正的条件/桥接估计器。** 完整 m64 benchmark 证明：即使允许任意单提议分布、给出真实 root/均值/归一化，未扣除的独立普通 importance mean 仍有每几何超过 `1.3542e10` 样本的 SNR3 下界。解析扣除几何指定的最小 rank1 层后，一个理想的剩余层条件估计器的方差预算变成约 `6.8835e5` 条 iid 剩余层样本/几何。

后一个数**不是实际采样算法或生产预算**：还欠可核验的条件采样核、剩余层质量/分区函数比、root 与均值估计误差以及相关时间。P0 交付应给出这些量的可控联合误差或一个明确不可能/低收益结论；不是先开普通 MC、换一个 twist 名称或追加耦合点。保留 #154/#334/F4 停线。

## 1. 一个精确区分：rank 偏置与中间扇区活性

令每几何的三个严格正受限分区函数为 Z0,Z1,Z2，r=0,1,2，q=r-1，E=q²。定义

```
eta = (1/2) log(Z2/Z0),
rho = Z1/(2 sqrt(Z0 Z2)),   xi = log(rho),
H = cosh(eta)+rho.
```

则完整三态分布可重写为

```
P0=e^(-eta)/(2H),  P1=rho/H,  P2=e^eta/(2H),
qbar=sinh(eta)/H,  Ebar=cosh(eta)/H.
```

这是 q/E 已有信息的坐标变换，不是新状态维数或新观测。有限正温度保证本例三扇区均有正权；边界处需取极限或保留原 Zr，不能制造有定义的对数。

乘一个 rank fugacity exp(b r) 时，`eta -> eta+b`、`rho -> rho`。因此

```
I_sector=Z1²/(Z0 Z2)
```

对 rank 倾斜严格不变。现有 star/drop 两律在**同一个热参数**下正好只移动 eta；在各自 root 上不能直接把数值 rho 当相同。

相反，若插入 beta 只支持 rank1，则 `Z1 -> Z1(1+epsilon E[beta|r=1])`，Z0/Z2 不变。其固定热参数切向为

```
eta_epsilon=0,  xi_epsilon=E[beta|r=1].
```

这精确定位 `0dda27ba` 的 denominator-only trace：直接 q/E 分子为零，并不等于中间扇区分母未改变。**这不是重做那份 Q4 计算，也不是把该 trace 的 epsilon 等同于 logQ 或孔密度。**

普通 rank 投影与 rank1-only trace 是不同的有限扇区方向，但不是自动独立的连续场；共同 root 的移动还会把它们的最终读数耦合起来。

## 2. 哪些微观源信息能被这个 U 看到？

令 x=log h 为共同热坐标，微观权重 w proportional to exp(x K+t S)，且本节 S 无显式 x 依赖。对每个 rank 定义

```
k_r=E[K|r],   s_r=E[S|r].
```

直接微分给出

```
eta_x=(k2-k0)/2,          xi_x=k1-(k0+k2)/2,
eta_t=(s2-s0)/2,          xi_t=s1-(s0+s2)/2,
eta_tx=(Cov(K,S|2)-Cov(K,S|0))/2,
xi_tx=Cov(K,S|1)-(Cov(K,S|0)+Cov(K,S|2))/2.
```

eta_xx/xi_xx 分别把 k_r 替换成 Var(K|r)。对于显式依赖 x 的源，还必须加上 E[partial_x S|r] 的项，不能套用本式略去它。

所以，对 q/E 和热导数组成的原 U，一阶源作用完全由这两种条件均值对比及其热变化决定。它不识别能产生同样净切向的不同微观实现。对这些函数的重建是接口检查；只有事前约束它们的符号、幅度或跨几何关系，才是可失败的机制预测。

把 `P=1-Ebar=P1`，则

```
q_eta=Ebar-qbar²,  q_xi=-P qbar,
E_eta=P qbar,      E_xi=-P Ebar.
```

在单个几何**自身**的 q=0 根处，有特别简洁的恒等式

```
E_x/q_x = -P1 * [k1-(k0+k2)/2] / [(k2-k0)/2].
```

这是稀有质量与条件热中心偏差的乘积。原 paired root 不要求每个几何各自 q=0，因此最终计算没有使用这条简化式。

## 3. 先消去共同温度，再分配原 U 的响应

记

```
Q=mean_g(q_g),  Y=P4(E_g),  D=Q_x,
A=N^(13/8)/2,   U=A Y_x/D,   Q(x0)=0,
Jq_g=partial_t q_g,   w=mean(Jq)/D.
```

定义源对各扇区坐标在**固定 pooled Q** 下的切向

```
a_g=eta_t,g-w eta_x,g,
b_g=xi_t,g-w xi_x,g.
```

这里 w 是一个 x 的函数，必须在微分后才代入 root。于是原 U 的完整源导数精确为

```
V_eta/A = (1/D) partial_x P4[P1 q a],
V_xi/A  = -(1/D) partial_x P4[P1 E b],
V = V_eta+V_xi.
```

证明：固定 Q 的源导数为 `J_E-w E_x=P1 q a-P1 E b`；再用 `U=A partial_Q Y` 及 `partial_Q=D^-1 partial_x`。等价地展开即原有 moving-root/source/slope 四项公式。

给 S 增加任意共同温度切向 c(x)K 时，eta_t/xi_t 增加 c 倍对应 x 导数，w 增加 c；a,b 各自不变。其 x 导数中 c' 也严格抵消。因此这里两项分配不依赖该共同温度 gauge。

**不允许**先将每几何 q_g 设零、冻结 w、丢掉 w_x，或者把这两项叫作事件的因果百分比。它们不同于之前按作用量身份 `V_Sstar=2V_beta_null+V_q` 的拆分：后者拆微观源，前者拆固定 pooled Q 后的响应坐标。

### 本次实际计算

使用既有完整 `(K,g,q,count)` 表和恒等式 `Sstar=51-K-g`。在 `h=[p/(1-p)]/m`、m=exp(t) 的共同 chart 中权重是 `h^K m^-g`，所以 -g 是 score；改为原固定 p 的 Sstar 只差常数与 -K。两者均已代入核验，两项分配保持相同。

| N25, Q=1, fixed Sstar | 数值显示 |
|---|---:|
| U | +0.880466156963367675 |
| V_eta：rank 偏置的固定-Q分配 | +0.301509231369101356 |
| V_xi：中间扇区的固定-Q分配 | -0.175343867954932215 |
| V_Sstar，总和/独立直接式 | +0.126165363414169141 |

两项相反的符号由 V/A 的有理区间分别保证；正的面积因子仅用于显示小数。没有将 t=0 的分配外推到 Q4、强耦合、大 N 或一个特定 irrep。该结果是旧 exact population 的新派生读法，不增加独立证据。

## 4. 真正的单提议 importance floor：稀有性之外还有有符号抵消

沿用 `f3ecde7d` 的同一个 m64/star/N25 目标。在每几何

```
X_g=(K-mu_g) 1_(r=1),   mu_g=E_g K,
C_g=E_g X_g,   theta=C_axis-C_tilt.
```

oracle 已给真实 root、mu、归一化与 U 的分母。两几何独立，各 n 条。

对于任意一个提议概率 pi，普通无偏估计量为 `X(omega)p(omega)/pi(omega)`。由 Cauchy-Schwarz，

```
E_pi[(p X/pi)²] >= [E_p |X|]²,
Var_pi(p X/pi) >= [E_p |X|]²-[E_p X]².
```

等号在 `pi proportional to p|X|` 成立（可在 X=0 集合给零质量；或由满支持提议逼近）。这个最优提议本身需要未知的绝对矩归一化，所以这里首先是 oracle 下界。

**边界：**这是固定 integrand、单提议、普通 importance mean、独立几何的类内结果。控制变量、精确扣除、多个正负分量提议、相关耦合、条件积分、估计方程等会改变估计器，不受同一个下界约束。它不是 sign problem 的所有算法下界，也不是 SNR 到置信覆盖的转换。

完整 m64 源表给出

```
theta=+1.4356412693155306275083391892e-19,
n_IS >= 9 sum_g{(E|X_g|)²-C_g²}/theta²
      =13542251055.5519287969590743749... .
```

因此达到 SNR3 至少需 13,542,251,056 条/几何，哪怕提议分布已经类内最优。仅称“使用 importance sampling”没有回答精度问题。

## 5. 一个具体改变 integrand 的步骤：精确扣除最小绕环层

按既有绕环 barrier，不看新输出选择 cutoff：axis 最小 rank1 g=9，tilted 最小 g=13。写

```
C_g=C_min,g+E[X_g 1_(g>gmin,g)].
```

C_min 从该已知有限表精确求和。原 theta 和原 U 不变，只改变需要随机估计的余项。在这个 benchmark 中没有真实采样；也不声称大 N 的最小层计数已经计算好。

| 同一 star/m64/N25 目标，各几何等额 n | n 对应 SNR3 |
|---|---:|
| 原 iid 平均 X（复现既有预算） | 1.518025804436592e25 |
| 未扣除：单提议 importance 类内最优下界 | 1.354225105555193e10 |
| 未扣除：已知 P1 的理想 rank1 条件 iid | 1.694087689450438e10 |
| 扣除最小层：余项 importance 类内最优下界 | 5.156662505208474e5 |
| 扣除最小层：已知余项质量的理想条件 iid | 6.883500471884592e5 |

最后一行使用

```
C_hat=C_min+p_rem*(mean(K | rank1,g>gmin)-mu),
Var(C_hat)=p_rem² Var(K | rank1,g>gmin)/n.
```

所有 p_rem/C_min/mu 都被当成精确已知；root 与 pooled 分母亦同。它是**可实现估计器应对照的理想条件方差**，不是已部署的 sampling kernel 或充分的实际运行预算。

axis 上 C_min约 `+5.0665006783e-18`，剩余均值约 `-4.8975592792e-18`。大项仍然相消，但其中一项现在是确定值，不再携带 MC 方差。这就是代数扣除与单纯增加稀有扇区访问的区别。

要把这一方向转成实际算法，必须给出：

1. 余项条件分布的正确采样及真实相关时间/误差控制，不把未混合 MCMC 当 iid。
2. p_rem、normalizer、C_min、mu 与 root 的联合估计或证书；将误差传播回原 U，而不只回余项均值。
3. 更大几何上最小绕环层的计数/积分成本与未计层完整性。N25 的 oracle 表不能冒充通用算法。

外部方法可以作为工具，但“用了 MBAR/桥接/正 partition”不是本模型的方差证明。Shirts–Chodera 的 MBAR提供多平衡态资料的统计组合；Owen 的 estimating-equation/positivisation 路线也说明本节单提议下界不适用于所有设计。它们均未在本贡献中被实现或当成成功结果。

## 6. 与 Q1 carrier 和当前研究边界的关系

固定 occupation continuation `Q=exp(2t)` 的总切向可直接从 Sstar 条件矩计算。上文避免了对“整数 Q 的 twist 个数”求导，但没有因此确定总切向由哪个表示/CFT family 承担。

一个候选 carrier 若只给出 Q4 的 trace，或者只给出 R/(sqrtQ-1) 的一次零，并没有提供所需的 Q1 净切向。`0dda27ba` 已规定其完整可去 jet 和所有 baseline-preserving 条件，本贡献不重开正则 unlabelled endpoint 的选择零。

这次应优先输出的是：候选 carrier 在指定完成中的条件 rank-sector 权重及其热源 jet，去除共同温度后对上述 V_eta/V_xi 的**不同预测**。若两机制给出相同净函数，原 q/E/U 不可能单独把它们区别开；源名称与场名称不能替代预测。

弱 Q 的场身份与 m64 的算法 benchmark 是不同问题。这里没有用强耦合方差结论检验弱 Q 指数，也没有因为 P0 标注去恢复已停线 lag1/contact/F4 或新的全局生产。

## 7. 复现、证据级别与本次读取范围

本贡献的有限恒等式有完整推导；实测部分其实是**已枚举有限总体的确定性再计算**，不是新的随机实验。物理解释/未来算法仍是建议，不调整科学 claim ledger。

- 输入均来自 `cae9c8997b5994c218bfe060f75656137f745755:experiments/p337-finite-law-window-20260831/inputs/`，两 Git blob SHA 已逐字核验。
- 160 次有理 root 二分；每个区间操作作 512-bit directed dyadic rounding，最后作 outward decimal enclosure。概率按每几何分别归一化，未重复乘 binomial 系数。
- 16 项 focused tests 全通过：含整数多项式与独立 Decimal 逐行中心化计算、全部五种方差预算、温度 gauge、rank tilt、rank1-only 插入、损坏文件拒绝、root 两端异号。
- 全仓测试、CI、MC、云机、一般大 N 采样器均未在本贡献中声称已完成。
- 开发中有两次性能超时：一次临时 Decimal/高精度原型在三阶矩内重复求均值，一次初版 Fraction 区间未及时约束分母；后改为缓存均值和严格外向 dyadic rounding。另有一次 shell 命令在结果已成功写出后超时。保留这些执行边界，不把成功生产次数伪装成“从未调试”。没有新随机数据。

运行（输出路径必须不存在）：

```
python scripts/p337_sector_quotient_review.py --output /tmp/p337-sector-review.json
python -m unittest discover -s tests -p 'test_p337_sector_quotient_review.py' -v
```

机器可读区间、输入及科学边界位于 `experiments/p337-sector-quotient-review-20260831/`。程序只需 Python 标准库，兼容 Python3.9 语法；本地实际解释器见 receipt。

### 固定来源

- [当前已完成结果与实际缺口 @f3ecde7d](https://github.com/LightChainr/Matching-One/blob/f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78/docs/STATUS.md)
- [原 iid 估计器预算 @f3ecde7d](https://github.com/LightChainr/Matching-One/blob/f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78/experiments/p337-estimator-access-20260831/RESULT.md)
- [Q4 normalization-only 传递 @0dda27ba](https://github.com/LightChainr/Matching-One/blob/0dda27bab3d1b6a749a0a32b3dde666b7fe9a0dd/notes/closed-source-s4-trace-transmission-result.md)
- [Q1 complete removable jet @0dda27ba](https://github.com/LightChainr/Matching-One/blob/0dda27bab3d1b6a749a0a32b3dde666b7fe9a0dd/notes/closed-source-removable-twist-jet-interface.md)
- [正源唯一 root @85fd4923](https://github.com/LightChainr/Matching-One/blob/85fd4923/notes/closed-source-critical-root-order.md)
- [最小绕环 barrier 与强耦合层 @85fd4923](https://github.com/LightChainr/Matching-One/blob/85fd4923/notes/closed-source-size-sign-and-transmission.md)
- Shirts & Chodera, JCP129,124105(2008), [arXiv:0801.1426v3](https://arxiv.org/abs/0801.1426v3).
- Owen, *Zero variance self-normalized importance sampling via estimating equations*, [arXiv:2510.00389v1](https://arxiv.org/abs/2510.00389v1), theoretical boundary reference; not executed here.
