# P337 thermal-pivotal 双通道：#536 来源卡与 #537 当前闸门

**状态（2026-09-01）。** 仓库团队已把
[Issue #536](https://github.com/LightChainr/Matching-One/issues/536) 关闭为重复项，
其正文与对话附件仍保留 finite-N25 计算的来源；当前唯一 canonical P0 是开放的
[Issue #537](https://github.com/LightChainr/Matching-One/issues/537)。本卡没有执行
这些生命周期变更。#536 的数值与反例语义仍是 **issue-only**：本轮没有独立导入
附件脚本/整数表，没有重跑、验证或提升为 `main_integrated` 事实。它们复用已完成
N25 总体，亦不是新的独立证据。开放的
[Issue #539](https://github.com/LightChainr/Matching-One/issues/539) 只保留为 P2
exact-N25 复现支撑，不与 #537 并列为第二个 P0。

上游解析边界见
[canonical pair thermal-pivotal gate](regular-pair-thermal-pivotal-gate.md)：
raw canonical kernel 在 exact `p_c` 的绝对可和性已经闭合，但 thermal 导数、
`q/E` 权重、moving root 与 original U 仍需要三位置 pivotal 控制。本卡给出该
闸门的有限图双通道接口和一次严格有界的语义复放合同；它不宣称关闭渐近闸门。

## 1. 对称 midpoint 恒等式

在有限 Bernoulli occupation 图上，令 `O`、`G` 对 `p` 无显式依赖。对站点
`z` 定义

\[
 \Delta_zF=F(\omega^{z=1})-F(\omega^{z=0}),\qquad
 F^{\rm mid}_z=\frac{F(\omega^{z=1})+F(\omega^{z=0})}{2}.
\]

`E_rest` 表示只对除 `z` 外的 occupation 取期望。Issue #536 给出的有限图
恒等式为

\[
\partial_p\operatorname{Cov}_p(O,G)
=\sum_z\mathbb E_{\rm rest}
 \big[(G^{\rm mid}_z-\mathbb E_pG)\Delta_zO\big]
+\sum_z\mathbb E_{\rm rest}
 \big[(O^{\rm mid}_z-\mathbb E_pO)\Delta_zG\big].                 \tag{1}
\]

本卡固定把第一和第二个和分别称为：

- **observable/rank-pivot channel**：kernel midpoint 乘 `Delta_z O`；
- **kernel-pivot channel**：observable midpoint 乘 `Delta_z G`。

这里的 midpoint 是对 forced-off/forced-on 两端的算术平均，**不是**按 `p`
加权的 conditional mean；若改用 `p` 加权平均，必须补交叉项。两个 channel
是同一个 covariance derivative 的代数分解，不是因果百分比、独立来源或两张
统计票。

## 2. 两个 4x4 反例固定双通道语义

两个控制均在 `4x4` nearest-neighbour torus 上，顶点编号
`i=x+4y`，marks 为 `0,2`，forced site 为 `z=5`；`z` 不与任一 mark 相邻。

### 2.1 rank 改变而 kernel 不变

在

```text
occupied-before = {3,4,6,8,10,11}
```

中加入 `5` 后，torus rank 从 `0` 变为 `1`，但 canonical pair kernel 始终
为 `g=1/4`，八端口 partition 在两端均为

```text
[0,1,2,3,0,3,4,5].
```

因此这里 `Delta_z g=0`，而全局 rank/`q/E` observer 可以 pivotal。它严格
说明：只测 kernel pivot 会漏掉 **kernel-preserving topological pivot**；局部
八端口 partition 不变不能推出 original-U 权重不变。

### 2.2 kernel 改变而 rank 不变

在

```text
occupied-before = {3,4,6}
```

中加入 `5` 后，`g:0 -> 1/4`，而 torus rank 始终为 `0`。因此 kernel 可以
pivotal，而全局 rank 不必 pivotal。

两个控制共同排除“双通道具有相同 support”或“一列可由另一列代替”的说法。
它们只固定有限图语义；不证明任一通道在临界尺度上占优。

### 2.3 shared-count 不变也可以远程改写 kernel

#537 还报告一个对每个 `L>=4` 都存在的显式有限环面族：取
`x=(0,0)`、`y=(1,1)` 和 `z=(floor(L/2),0)`，强制翻转远处 `z` 前后 torus
rank 都是 `0`、shared-component count 都是 `2`，但 canonical
`g_xy:1/4 -> 1/2`。该构造目前仍按 issue-only 边界引用；它不是概率下界。
它严格排除“只要 rank 和 shared count 不变，远处 carrier 就不能改写 kernel”
的局域化证明，并把交付 A 的 carrier 词汇从单纯的 `s`-transition 扩充为
**same-count signed rerouting**。因此至少需要同时记录 mark separation 与
pivot-to-mark distance 两个尺度。

## 3. N25 exact midpoint decomposition：issue-only 数值

Issue 正文报告对 canonical N25 axis `HNF(5,0,5)` 与 tilted
`HNF(25,18,1)` 各遍历完整 `2^25` population，并读取每个几何
`25*2^24` 个 rest-flip states。translation-anchor observable 为

\[
 G_0=\frac1N\sum_y g_{0y}.
\]

固定 pooled root、`A_N=N^(13/8)/2` 和
`DeltaCos4=1152/625` 后，正文给出：

| pair scope | observable/rank pivot | kernel pivot | combined root | total `J2` |
|---|---:|---:|---:|---:|
| all pairs | −0.005935841948386444 | +0.0004141042799437205 | +0.0000023062436033220 | −0.005519431424839401 |
| NN | −0.0008524574521935375 | −0.0008997868197966343 | +0.0000011698175873727 | −0.001751074454402799 |
| nonNN | −0.005083384496192906 | +0.001313891099740355 | +0.0000011364260159493 | −0.003768356970436602 |

正文称各列符号有有理外包络检查，并记录 Bernstein--Russo、midpoint
product、旧完整点对矩、反射及 tiny-configuration 检查。但这些检查的脚本和
整数表只存在对话附件 `matching_one_thermal_review_20260901`；本卡没有独立
导入或核验它们。因此当前可作出的最强仓库表述是：

1. 在 #536 的固定 midpoint 约定下，有限 N25 总 J2 的负值主要落在
   observable/rank-pivot 列；
2. kernel-pivot 在 NN 为负、nonNN 为正，不能用一个无类型的“kernel 占比”
   概括；
3. combined-root 列很小，但仍属于精确总和，不能静默删除；
4. 三列共享同一 N25 population、root 和 algebraic split，不是独立证据；
5. 这些 finite-N25 符号不是 `p_c` 渐近符号、CFT loading 或连续场分解。

`G_0` 的 translation anchor 只足以构造线性矩与 site-summed midpoint
channels。anchor variance 不是 full-source variance；同样，
`abs(sum_y Delta g)` 不是 `sum_y abs(Delta g)`。

## 4. endpoint 项已经从开放瓶颈中移除

因为任一 mark occupied 时 canonical `g_xy=0`，Issue #536 报告 endpoint
项精确化为

\[
 \sum_{z\in\{x,y\}}\mathbb E_p\Delta_zg_{xy}
 =-\frac{2\,\mathbb E_p g_{xy}}{1-p}.                       \tag{2}
\]

在 exact `p_c`，式(2)的 `y` 空间和由 raw summability theorem 控制。因此
endpoint raw derivative 不再是 thermal gate 的未知部分。真正剩余的是：

- `z` 位于端点邻域、连接区域或任意远 merger 区域时的 kernel pivot；
- `Delta_z g=0` 但 `Delta_z q` 或 `Delta_z E` 非零的 topology pivot；
- near-critical/pooled-root 误差与 original-U 面积归一化。

这与上游 [thermal-pivotal gate](regular-pair-thermal-pivotal-gate.md) 的
remote-merger `Delta_z g=-1/16` 反例和 `q/E`-weighted influence 对象直接
衔接。

## 5. 完整 original-U `J2/A_N` 接口

沿用现有 original-U 记号：

\[
 M=\operatorname{mean}q,\quad Y=P_4(E),\quad D=M_p,\quad
 R=\frac{Y_p}{D},
\]

并对固定 pair source `G` 定义

\[
 jM=\operatorname{mean}\operatorname{Cov}(q,G),\qquad
 jY=P_4\operatorname{Cov}(E,G).
\]

需要控制的完整归一化响应是

\[
 \boxed{
 \frac{J2}{A_N}
 =\frac{jY_p-RjM_p}{D}
 -\frac{(Y_{pp}-RM_{pp})jM}{D^2}
 }                                                             \tag{3}
\]

而不是裸的 `partial_p E[g_xy]`。式(1)必须分别用于 `O=q,E`，再通过
`P4`、geometry pairing、root 和 slope 权重进入式(3)。表中的
observable/rank pivot 与 kernel pivot 来自 `jY_p-RjM_p` 的两类 midpoint
项；combined root 保留式(3)中其余 root/centering contribution。

要证明 canonical pair 对 original U 渐近消失，必须证明带 `A_N` 的式(3)
趋零。raw kernel 可和、thermal derivative 可和、甚至 `D` 远离零，任何一项
单独都不自动给出这个结论。面积单位、near-critical root window 以及
`D,Y_p,Y_pp,M_pp` 的尺度必须明示。

对任一 centered channel `H_xyz`，还必须区分

\[
 \sum_y\left|\sum_z\mathbb E H_{xyz}\right|
 \quad\text{与}\quad
 \sum_{y,z}\mathbb E|H_{xyz}|.                              \tag{4}
\]

后者更强；其发散不能证明前一个 signed 对象发散，也不能据此排除全部抵消
机制。

## 6. 交付 A：signed/absolute three-site support

交付 A 固定 canonical `g`，不更换 completion。对 `O=q,E` 分别记录

\[
 H^{\rm obs}_{xyz}(O)
 =(g^{\rm mid}_{xy,z}-\mathbb Eg_{xy})\Delta_zO,
\]

\[
 H^{\rm ker}_{xyz}(O)
 =(O^{\rm mid}_z-\mathbb EO)\Delta_zg_{xy}.                 \tag{5}
\]

至少按 endpoint、endpoint neighbourhood、连接区域、外部 merger/rank
pivotal，以及 shared-component `s=2/3/4` 分层给出：

- pair/site 级 signed sum；
- 在 pair/site 级先取绝对值的 absolute sum；
- 各层加回总量和同一 population 的联合协方差；
- `s=2` two-bridge factorized carrier、`s`-transition/remote-merger carrier、
  `g` 不变但 `q/E` 改变的 topological carrier 三类的明确归属。

交付可以是适用条件清楚的定理、一个具体几何反例，或把唯一未控 joint
pivotal event 精确命名；不要求先拟合指数。

**A 的停止边界：** 若 signed 与 absolute 支持不能区分三类 carrier，结果
记为 `unresolved`。不新增第四个 descriptor、距离 window 或自由指数救场；
也不把 absolute bound 的失败改写成 signed divergence。

### 6.1 先做有限维判伪器，再决定哪一个 tail 值得证明

PR #509 当前 head `2785e3bb` 把完整 original-U 组合压缩为 root-conditioned
Hessian。若 `M(p(u,epsilon),epsilon)=u` 且
`Yhat(u,epsilon)=Y(p(u,epsilon),epsilon)`，则

\[
 \left.\partial_\epsilon\widehat Y\right|_u=jY-RjM,
 \qquad
 \partial_u\partial_\epsilon\widehat Y=\frac{T_N}{D}=\frac{J_N}{A_N}. \tag{6}
\]

其 signed covariance 形式把 kernel、readout、root 与 slope 放进同一个
Schur 投影：`T_t=<H,(a-Ea)S-beta B>_pool`。这意味着“raw kernel 非零”或
“某个 midpoint 列很大”都可能在完整 Hessian 中结构性抵消。

#537 因此把交付 A 的第一步收紧为一个有界精确判伪器：在相同 `C4` 与 Schur
投影后，构造 ordinary four-arm landing source/thermal transfer matrix，并检查
其所有 `2x2` minors。

- 任一 projected minor 非零，直接否定“leading three-packet 只是一个 pure
  thermal coordinate”；随后应转向该非秩一 landing functional 的 signed
  dyadic bound。
- 所有 minors 为零，只使 four-packet remainder 成为最简存活机制；它仍需
  把 `R^4*pi4(R)^4` 包络、collision tree 与 near-critical root window 连成
  定理，不能把有限秩一观察写成证明。

这个有限判伪器比直接从三个 carrier 各自开始做完整渐近证明更先验地缩小机制
空间；它不新增 production、自由 descriptor 或统计票。

## 7. 交付 B：L32/L64 各前 32 counter 的有界复放（已完成）

交付 B 只复用既有 L32/r8、L64/r16 producer 与 RNG：

1. 每个尺寸读取原 block 最前 `32` 个 configuration counters，共 `64` 个；
2. 复用每个 configuration 原有的 `32` 个 anchor/direction pairs；
3. 对每个站点强制 off/on，保留 pair-level
   `g0,g1,q0,q1,E0,E1`、component/support layer、位移与几何类；
4. absolute value 必须在 pair/site 级先取，不能对 anchor 总和后取绝对值；
5. 保存 signed/absolute moments、同配置完整联合矩、wall/CPU/RSS、线程数及
   输入/输出哈希。

语义验收边界沿用 Issue 正文：原 `g` 必须复现原 producer；midpoint
product 逐项一致；两个 4x4 控制、式(2) endpoint 因子、translation-anchor
线性矩及式(3) root/slope 全式必须一致；tiny 图可与 brute force 对照。

**B 的停止边界：** 达到 L32 `32` + L64 `32` = `64` 个原配置后立即停止。
无事件、符号混合或精度不足均不是增加 seed/counter 的理由。这是确定性语义
与成本预检，不作物理显著性、临界符号或生产预算声明。只有交付 A 已给竞争
预测、且另立冻结合同后，才讨论扩大旧数据复放或新采集；本卡不授权 GPU、
新 Monte Carlo、raw 距离网格或宏观窗口 top-up。

冻结实现、raw、描述性score与receipt现已交付于
[`results/p337-thermal-pivotal-preflight/`](../results/p337-thermal-pivotal-preflight/REPORT.md)。
64个原配置、5,242,880个pair/site callback全部通过forced/original复现、完整
carrier-mask partition与q/E midpoint恒等式。有限回放只出现18个kernel-changing
callback：L32为10个kernel-only加1个joint，L64为7个kernel-only；L64的
observable/rank原语为0，非零kernel原语落在external shell 2/3。L32的非零原语
落在shell 1，其中一个two-bridge事件属于square-NN，其余为external
shared-transition。没有观察到kernel-preserving topological callback。

这批callback数不是独立样本数；零出现也不是零概率界。它只判决了实现可行性和
有限carrier形态：在该冻结子集上，变化是稀疏、局域且以kernel-change为主，不能用
大量直接`q/E`翻转来解释Issue-only N25 observable/rank总体项。复放按合同停止，
不增加counter或seed。

## 8. 当前判决

- N25 首次 midpoint 通道拆分由已关闭 #536 保留为来源，不再作为“下一首次任务”。
- 仓库当前可安全采用的是恒等式、反例语义、#537 的理论目标与已完成 B 的有限
  输出；issue-only 数值必须继续携带未独立导入/验证标签。
- 下一高信息增益对象是 #537 的 projected landing-matrix minor 判伪器。它决定
  随后证明 non-rank-one signed landing tail，还是证明 rank-one 后的 four-packet
  absolute remainder；不是更多 counter、再拟合 raw `C(r)` 或扩大 macro ratio。
- #539 只做 P2 exact-N25 reproducibility support，且应先导入/消费 #536 已有附件，
  不把完整 `2^25` rest-flip 枚举重复一遍。
- 本卡不把 midpoint channel 称为场、因果份额或独立 evidence，也不把 remote
  merger 宣称为已证实的渐近主导机制。

上游闸门与未证明的 dyadic/near-critical 条件继续由
[regular-pair-thermal-pivotal-gate.md](regular-pair-thermal-pivotal-gate.md)
维护；[文献桥](p337-thermal-pivotal-literature-bridge.md)进一步说明四臂只关闭raw、
六臂只在已证明的collision annulus后关闭fused尺度，而bulk/topological需要逐通道
absolute命题或signed adjacent-Fourier层抵消。near-critical/pooled root必须另证
`|p_L-p_c|L^2 pi4(L)=O(1)`或`L`不超过相关长度。本卡只固定 #536 的 finite
two-channel 来源、#537 的当前判伪顺序和停止规则。
