# [P0 建议 / exact] 真正的 Q1 激活、补全选择与局部端口型别

**建议优先级：P0，针对精确机制判别；不是新生产队列授权，也不升级连续场的证据等级。**

相关 #337、#275、Draft #267 / #509。main 基线 `580678a769b1360d701fb708e1f557a2db626f16`；有限计算的科学输入固定到 `bea717e826df5a22518774b1725ae7bcbe2cb801`。本贡献只增加文件，不覆盖冻结合同、既有结果或分支。

## 0. 已完成工作不再重新布置

对称闭合迹 `[Q-2,2]` 的 generic-Q 两连接核、N25 packing、完整归一化接口和 Q1 传递已经完成：`beta_plus(1)=-1_A-1_B`，`V_plus=-0.001904836180602413`。Draft #267 的 `5c1f9d3b7971a41d07db3c9fa4ac86529c90c199` 也已给出其真实混合导数 `+0.005036496028411871`；本次不重复该问题。

提交前又核对了执行方的新成果：`7681eedd938019d977ede41a7d74ee1b88ffbc50` 已完成对称 local four-port interaction 的 `V_old=+0.0018155512845251097`；`2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb` 已完成两插入 pole 障碍、canonical regular completion 与 `W=-0.04503611397592696`。这些都不再是首次待办。

**本次反对称有限计算是一个不同的精确对照；P0 下一步是补全选择可辨识性与具体长距离/尺度预测，不是追加一个 N25 系数来放大效果。**

## 1. 同一两个真实闭合类型的反对称系数

令 `V_Q` 是颜色置换表示，`W_Q={v:sum v_i=0}` 是标准表示，取

```
A_Q = exterior^2 W_Q,
d_minus(Q)=(Q-1)(Q-2)/2.
```

对整数 Q>=4，这是 `[Q-2,1,1]`。令 S 交换两个颜色端口，S 在本通道上作用为 -1，区别于旧 `[Q-2,2]` 的 +1。

沿用已有 N25 packing 所限定的两个连接类型：A 是两个各一次绕过指定第一 deck seam 的 essential components，闭合为 I，未滤颜色简并 Q^2；B 是一个两次绕过指定 seam 的 essential component，闭合为 S，未滤颜色简并 Q。

实际中心迹因而给出

```
beta_minus,A(Q) =  d_minus(Q)/Q^2 = (Q-1)(Q-2)/(2Q^2),
beta_minus,B(Q) = -d_minus(Q)/Q   = -(Q-1)(Q-2)/(2Q),
beta_minus,other(Q)=0.
```

A 的 `Fix(pi)^2` 在本通道的角色系数为 +1；B 的 `Fix(pi^2)` 为 -1。后者是置换幂的角色收缩，可以为负，不是负概率。rank0 的常数角色、rank2 和单次单簇的点置换角色都没有本通道。原 occupied-edge、contractible spectator 和 `Q^(-r/2)` 因子保持。

独立角色核验使用

```
chi_minus(pi)=binom(X1-1,2)-X2
            =[(X1-1)^2-(Fix(pi^2)-1)]/2.
```

这是 exterior-square character 恒等式。置换循环的 factorial moments 对加权次数<=Q 精确成立；所需最高次数为4，因而上述结果对所有 Q>=4 成立，不是对几个 Q 点插值。本包另对 Q=4,...,9 的全部共轭类作有理数求和。

**延拓约定明确：**使用上述指定的逐配置有理系数和已有正权 occupation family 延拓到 Q>0。整数值本身不能唯一排除 sin(pi Q) 一类额外解析函数；这里没有加入它们。Q1 处不声称存在一个非平凡的字面 S1 表示。

## 2. 真正基线消失的 Q1 激活

直接计算

```
beta_minus,A(1)=beta_minus,B(1)=0,
d_logQ beta_minus,A|1=-1/2,
d_logQ beta_minus,B|1=+1/2.
```

固定激活 score 是 `psi=(1_B-1_A)/2`，不是给旧对称通道手工乘 Q-1。

令 epsilon_trace 为闭合迹系数，**不是**孔密度、饱和参数或 logQ：

```
w_(Q,epsilon_trace)(omega)
 = w_Q(omega)[1+epsilon_trace beta_minus,Q(omega)].
```

在 Q=1 的邻域和足够小的 |epsilon_trace| 内，整体权重正。Q=1 时该迹对所有 p 都不改变测度，故对任意静态 O 有

```
d_logQ d_epsilon_trace <O>|_(1,0)=Cov_(Q1)(O,psi).
```

基准 Q 测度导数乘 beta_minus(1)=0；只留下固定 psi。对原 U 的非线性 root/slope 泛函，链式法则同样给出

```
K_minus := d_logQ d_epsilon_trace U|_(1,0) = V_psi.
```

**psi 所诱导的根移动和斜率变化仍全部保留。** 消失的是额外的基准 Q-变化交叉项，不是原归一化公式。若改用 `t=log m`、`Q=e^(2t)`，对应混合导数为 `2 K_minus`，没有额外 N 因子。

这是实际两参数混合导数，不是 raw-Q 分量归属。若分区函数和迹同时乘公共非零解析 prefactor c(Q)，因为迹在 Q1 为零，其一阶激活除以分区函数保持不变；旧对称通道的 raw-packet gauge 歧义不在此项发生。

## 3. 实际结果：原 global U 的激活严格为负

输入为原 N25 `(5,0)/(4,3)`、原 `DeltaCos4=1152/625`、原 pooled root 与 `A=25^(13/8)/2`。没有新增枚举、MC、root search 或 source fit。

| 固定 score / 导数 | 确定性数值显示 |
|---|---:|
| V_(1_A) | +0.001945570733316785 |
| V_(1_B) | -0.00004073455271437206 |
| 旧对称迹控制 `-V_A-V_B` | -0.001904836180602413 |
| **K_minus=(V_B-V_A)/2** | **-0.000993152643015579** |

零判决使用 `K_minus/A` 的 outward Fraction interval：

```
lo=-1808018764976707165760541220563119/170141183460469231731687303715884105728
hi=-3616037529953414331521082441126237/340282366920938463463374607431768211456 < 0
```

小数仅供显示；符号不取自浮点或拟合。两几何分别归一化。由 rank1 支持，若 `b_g=<psi>_g`，则 `j_q=-q b`、`j_E=-E b`；直接 q/E 加权分子为零并不等于归一化响应为零。

记 `Qbar=mean(q),Y=P4(E),D=Qbar_p`，在原 root，实际评分为

```
V/A = jY_p/D - Y_pp*jQ/D^2 - Y_p*jQ_p/D^2
      + Y_p*Qbar_pp*jQ/D^3.
```

输入已含 count multiplicity，不重复乘 binomial coefficient。继承的有理根区间被原整数多项式两端异号和正斜率重新验证；没有重新搜索根。这些结果和旧 symmetric control 共享同一精确总体，不是独立证据票。

## 4. P0 建议一：除去未固定的补全自由度，再提出场级预测

最新 `2ba8863f` 已明确：在 `K2+c(Q)K0` 中，正则性只迫使 `c(1)=1`。设 `c(Q)=1+alpha(Q-1)+...`，同一声明的混合响应满足

```
W_i(alpha)=W_i(0)-alpha V_old,i.
```

alpha 是不同的微观 Q 延拓选择，**不是已证明的物理规范等价**。不能拟合 alpha 救结果；同样不能把某个 W 的符号称为不依赖补全的场指纹。

一个无需拟合 alpha、也不使用病态比值的必要不变量是

```
I_ij = V_old,i W_j - V_old,j W_i.
```

代入上式，alpha 严格消失。更一般应研究 W 在 span(V_old) 的商类。若有多个未固定 counterterms，则对它们响应列的张成空间取商，保留共同协方差和确切零空间。

**下一份有价值的交付：**给出指定相互作用在匹配尺寸、同一个原 U 上的 `I_(N,cN)` 预测，或对原 root-normalized entry/completion 对比的预测，并检验其是否区分候选。两个读出必须使用同一个微观 alpha；逐尺寸或逐读出重新定义 alpha 会破坏此检验。

I 非零只排除“整条响应都是一个共同补全方向”，不能数 CFT 场；I 为零也不证明没有物理效应。若为 canonical `c(Q)≡1` 提供额外物理原则，则可直接预测其 W；正则性本身不是该原则。本建议不修改已有 canonical 评分，也不把商类误称为不同微观模型完全等价。

## 5. P0 建议二：用无 Q1 极点的对照审查局部插入型别

对 Q>=4，令

```
J_Q=I-11^T/Q,
P_minus=(I-S)(J_Q tensor J_Q)/2.
```

它满足 `P_minus^2=P_minus`、`S P_minus=P_minus S=-P_minus`。只含 Q 的逆幂，在 Q1 没有极点；任意有限、非奇异 Kronecker-delta 闭合网络都保持正则。Q1 唯一颜色时 J_Q=0，含插入的网络闭合为零，其 Q 导数却可非零，正如 K_minus 所示。

**所以“Q1 零值＋非零 Q 导数”不是 Jordan/对数混合的充分诊断；它可以来自普通表示维数的零，无需 singlet pole cancellation。** 本次只作有限模型反例，不否定可能存在的连续对数混合。

还有一个精确局部选择律。若一个**局部端口向量** l_Q 交换偶，即 `S l_Q=l_Q`，则 `P_minus l_Q=0`。更一般，若局部到 cut 的映射满足 `S B_Q=B_Q E_local`，它把局部交换偶向量送到零反对称投影。对于固定正则有理延拓，这个恒等零不能通过 Q 求导激活。

**不能把这个向量零律扩大为“任何交换不变的标量作用量都不能作用于本通道”。** 作为 endomorphism，P_minus 本身满足 `S P_minus S=P_minus`，而 `Tr(P_minus)=d_minus(Q)` 非零。交换不变算符可以作用于已有奇子空间，不等于由交换偶端点向量发射奇态。本包同时测试这两种情况。

因此，对这个反对称对照，局部识别应得到以下一种有限交付，而非重复已完成的对称四端口计算：

1. 对一个实际交换偶的局部发射及交换协变映射，证明其在本通道恒零，退休该具体识别；不扩展为全部 scalar-source/H4 机制的否定。
2. 给出明确的有序或带奇几何权的局部发射/接收，证明颜色置换和端口交换协变性，并在真实 multiplicity space 上使同一映射通过 I、S 两种外部闭合，得到 +d_minus、-d_minus 以及 rank0/rank2/单端口零控制。有限 witness 通过只是必要条件；一般 intertwiner 仍须覆盖声明的完成上下文。

闭合迹的中心投影是 endomorphism；pair insertions 常是向量、余向量或更大的 B/C contraction。应写出型别，不能用一个最终标量相符代替重数空间等式。更改几何顺序或加入奇权是一项新定义，不得追认为原无标记场。

完成这些之后才讨论空间旋转表示和尺度预测。`[Q-2,1,1]` 与 `[Q-2,2]` 不同；本次 K_minus 不提供 sqrt(N)、17/4、21/4 或 spin4 的身份。几何 P4 差商非零不是局部 spin4 证明。

本建议不自动启动 Q-grid、seam-grid、尺寸链、枚举或云生产。N25 packing 不能直接搬到更大 N；后者可能有 w|u|>2 的额外角色，旧 mod6 汇总也未必充分。P154/P334/F4 停线保持。

## 6. 核验、来源和复算

本包依赖标准库。16 项 focused tests 已通过：Q4至Q9精确角色求和，Q4/Q5投影幂等与 I/S closure，交换偶/奇端点，交换不变算符反例，输入篡改和错误根区间拒绝，全部系数重建，独立 Bernstein 形式的 moving-root 中央差商。

第一次测试全通过但有 CSV handle 的 ResourceWarning；改为 context manager 后再次全通过，科学输出未改变。完整证书复算成功。随后一个组合命令在 compileall 阶段超时，因此不声称该次完整 compileall 成功；另行检查了 Python3.9 grammar，未在本地运行3.9 runtime。未运行全仓测试、原总体枚举或完整 source-histogram Git-object 再提取。Compact rows 来自不可变 connector 内容，support totals 与已发布报告相符；提供可选完整再提取接口。

```
python scripts/p337_antisymmetric_trace_review.py --verify experiments/p337-antisymmetric-trace-20260831/result.json
python -m unittest discover -s tests -p 'test_p337_antisymmetric_trace_review.py' -v
# checkout 需包含 pinned science objects：
python scripts/p337_antisymmetric_trace_review.py --check-source-checkout /path/to/Matching-One
```

完整输入和 blob 来源见 [SOURCES.json](../experiments/p337-antisymmetric-trace-20260831/SOURCES.json)，有理证书见 [result.json](../experiments/p337-antisymmetric-trace-20260831/result.json)。

固定科学来源：

- [两种真实对称 closure 核](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/closed-source-two-trace-kernels-q1.md)
- [N25 packing 与 histogram 充分性](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/n25-winding-packing-and-pair-continuation.md)
- [稳定角色与有限 Q 边界](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/closed-source-stable-colour-character-continuation.md)
- [既有 Q1 对称迹评分](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/p337-q1-closed-trace-transmission-result.md)
- [canonical regular interaction 与补全自由度](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)

文献接口：Couvreur–Jacobsen–Vasseur, *Non-scalar operators for the Potts model in arbitrary dimension*, [arXiv:1704.02186](https://arxiv.org/abs/1704.02186)。它提供颜色与位置置换分类背景；本次不从该文移植这个正权晶格家族的普适性或具体指数。

**范围：**精确有限代数与有理评分，不是独立统计验证。P0 是建议优先级。没有修改 claim ledger、main、冻结历史、Issue 生命周期、其它 PR 优先级或服务器。
