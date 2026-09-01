# 当前成果：从取向信号到微观机制

**整理日期：2026-09-01，已读到 #275 P0 重置、#537 P1 重置、held-out N65 contact-stage 结果与两轮独立审计。** 当前唯一 P0 是原始 observable、normalizer 与候选可识别性；#537 保留完整 thermal/contact 渐近问题但不再自动扩尺寸。两个独立决策实验、齐次N50、canonical pair 有限空间块、finite rank-one/one-defect 判决和 N65 范围受限的前瞻验证均已完成并按 stop rule 收缩。[C3 相位可识别性](../experiments/p275-c3-phase-contract-20260901/REPORT.md)与[#275 判决树](../notes/p275-observable-identifiability-gate-20260901.md)记录新 P0；[下一步](NEXT-TARGETS.md)只列现在允许执行的工作。本任务继续交付 Draft #509，不合并。

## 压缩后仍须保留的结论账本

以下是已经完成的结论和语义控制，不是重新开放的计算任务：

- P43 even-channel 纠错的精确映射是 `DeltaS_cross = -DeltaS_either`；不重拟合的纠正分数为 `0.5700315436/2`。
- 独立 P43+P57 primary 合成拒绝 **global zero**；norm-5 冻结比较优先 H4 于 H12/H8，但 **child block alone** 仍与零相容。
- **N145->290** full curve 已完成；冻结`1e-10`分数仍为`9.35200/2`，但[#543零空间回溯](../experiments/p543-covariance-nullspace-audit-20260901/REPORT.md)证明解释依赖截断，故保留“该冻结分数拒绝共同倍率”，不再把被丢方向称作无害numerical null。它不再是 active compute。
- 单一 **scalar width** 和其他已测试的低阶标量捷径不足以解释较高 thermal jet；不再用自由指数追加同读出拟合。
- 有限体 **Russo** / chain-rule 的 **pivotal** 语义仍是精确控制；它本身不是新的机制证据票。
- N=26 的冻结有限族 `Beta(5,5)` 与 `Beta(7,7)` 均已被精确枚举否定；这不外推为一般代数或连续场结论。

## 当前唯一 P0：#275 原始观测量与可识别性

site-cluster colour-lift 的有限控制已经说明：未归一化 `q/E` numerator 的
pair-character 零规则不能直接传给 normalized expectation；trace、source 与
pooled moving-root `U` 也不能因相同 spin 标签而混称同一 observable。

[本分支 exact design](../experiments/p275-c3-phase-contract-20260901/REPORT.md)
完成 66 项符号检查。对同一个 phase-calibrated real C3 读数，H4/H8 在共享复振幅
或 signed-real gain 合同下可由合适第二旋转区分；任意复增益下所有两旋转设计均
不可识别。`7.5°` 可处理非零 signed-real gain，`15°` 则有精确符号别名反例。
这是一项已完成的条件设计结果，不是场身份或采集授权。P0 剩余验收是把两个
实际物理候选映射到同一 normalizer 下至少两个 raw 坐标，并用已有 covariance
做一次冻结 profile-rank 判决；若允许自由度使列空间相同，就登记不可识别并降级。

随后冻结的 2M paired N65 primitive-C3 读数在 signed-real 两候选合同内拒绝
pure H4（nominal `chi2=73.64/1, p=9.37e-18`），未拒绝 H8 alias
（`chi2=1.112/1, p=.292`）；100 个 leave-one-batch 判决同向。这个 finite observer
子门已停止、无 top-up。它没有测试混合谐波或任意复增益，也没有把 original
square-site `U` 映射到 H8，因此不完成 #275。

## 当前机制取舍：有限传递已闭合，raw 临界空间尾绝对可和

执行分支 [2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md) 已完成 canonical `Kreg=K2+K0`：所有有限网络的直接 Q1 响应为零，而固定 N25 的 `∂logQ∂epsilon U=-0.04503611398` 严格非零。旧 pure-K2 的一阶传递和固定四路径的13/8也已完成，不再等待“首次局部传入原U”。该数值依赖指定补全，不能称 completion-independent 或连续场确认。

**[空间选择规则](../experiments/p337-regular-spatial-support-20260901/RESULT.md)：两个非相邻标记若至多共享一个外部占据组件，首阶Q激活严格为零。** 对 canonical 核，占据总和后的 `Cxy=E[a_xy]` 满足 `|Cxy|≤(43/16)Pr{两点vacant且至少两个不同占据组件接触两处}`。两种独立精确算法的4140项全部一致；实际两组件构型给a=1/16，证明阈值可达到。原U的正确混合接口是完整 `W[a_xy]`，不是 `Cov(a_x,a_y)`。

执行 [`a237968f`](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md) 已用两个新依赖块拒绝 L64/距离16的有限零传递；L32/距离8到L64/距离16只保留冻结比值，不拟合指数。执行 [`410015f5`](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md) 又完成 `J2=d_logQ d_epsilon²U=-0.0055194314248394015`，相邻项为−0.0017510744544027990，其余非相邻项为−0.0037683569704366022；additive-linear与NN-contact-only两种global closure均已排除。非相邻类包含短对角，不等于宏观长程尾。

**[本次解析推进](../notes/p337-critical-spatial-summability.md)进一步排除 raw canonical interaction 的红外发散解释。** 两个不同 occupied components 同时到达两个远标记，会在两端各强制一个黑/白交替四臂事件。两个不交环带给四臂概率平方；van den Berg--Nolin 对临界方格 site 的严格 `alpha4>1` 界于是给 `E_pc|g_xy|≤C d^(-2-eta)`（某个`eta>0`）。每行绝对和一致有界，距离R之外的尾为`O(R^-eta)`。不需要方格5/4猜想，不新增距离点；停止距离网格、指数拟合、alpha扫描和“靠远距离累计放大”的解释。

Draft [PR532](https://github.com/LightChainr/Matching-One/pull/532) 的已提交增量是两共享组件端点因式分解和指定补全族的四路径3/2下界；它不给占据平均符号。其后续评论建议的uniform J2生产已被410015f5完成，不能重复启动；评论中的八通道压缩和`alpha=3/2` pair positivity尚无同一PR的新提交，保留为候选代数，不改canonical补全或当前优先级。

Draft [PR533](https://github.com/LightChainr/Matching-One/pull/533) 的 `5aa929a6` 保留了若干可信 C1 子结果：relaxed all-west Catalan/Toeplitz 单粒子恒等式、立即 `EW/WE` 回边禁制、一个合法 nonlocal one-west 反例，以及“先有指数双gap核则sine transform保留`pq`双零”的条件引理。它仍是 **P2/C0 overall**：真实beta cloud没有构造物理双gap核，rank1与rank0/2共用坐标只是断言，`w>=2`、second thermal/root/original-U和uniform remainder未闭合；`410015f5`、`2690f665`仍不在祖先链，旧`a*s`、CR与`O(Na^4)`/`O(Na^3)`矛盾仍在。Issue #542保持开放，不进入current claim。

**[thermal/pivotal 双通道审计](../notes/p337-thermal-pivotal-gate-audit.md)已经排除一个过度简化。** `d_p Cov(O,g_xy)` 精确分成 kernel reconnection 与原 rank/readout pivotal 两项；只控制 `partial_p E[g_xy]` 会漏掉第二项，并在 N9/N10/N13 的 `O=E` 精确控制中把总符号判反。Bell-8 的64,954个join和N13的1,198,080个state-edge-pair检查全部闭合；共享组件数不变时核也可换号。

#537 当前为 P1；有限 pure-thermal rank-one 路线已经按研究停止规则退休。[最初 preflight](../experiments/p537-landing-matrix-preflight-20260901/REPORT.md)在 six-block clean-two-bridge 合同下得到全部非零 minors。其后 PR #544 的 N25 radius-one collar、axis/tilted pooled-root Schur 首项 minor 为严格正的 `+2.6904188461441777e-14`；完整 fibres、聚合和 hashes 已独立复现。PR #545 给出任意 `R>=1,L>=2R+5` 的 broad four-arm pointwise family，以及 axis-L4 的 `Psi=-533831111/140737488355328` 和 matching-root 区间不变号证书。#544 最新 `e8e9c7cf` 又精确复现 axis-L4 并完成 axis-L5：`Psi4=-4.0685187141747587e-7`，而 `L^4 chi_perp` 为 `-0.66238/-0.67340`，与局部 `L^-4=N^-2` 衰减相容。

[独立范围审计](../notes/p537-finite-rank-one-decision-20260901.md)限制了这些结果的含义。#544 没有保存 `x+y+z` 全局 joint component identity/global no-extra flag，最终又把 Bell transition 汇总成 source absent/present；它证明显式 collar coarsening rank two，不证明 canonical ordinary/no-extra。#545 没有 formal Bell/no-extra 字段，并在 landing 子块内重新估计 `beta` 与 `S`；其 reduced-block 恒等式不是完整两几何总体先定 root counterterm 后的 original-U Schur summand。“JSON byte-identical”声明也未成立。L4→L5 的两个尺寸不能升级成指数或 full-U 尺度律。两条 Draft 因此关闭为 unmerged assets，不合并，也不再通过补 schema、L6、N、几何或 minor 救回 rank one。

完整渐近目标仍是把 `T_N=jY_p-R*jM_p-R_p*jM` 控制到 `o(D/A_N)`，其中 `partial_u partial_epsilon Yhat=T_N/D=J_N/A_N`，等价的 logit 形式为 `T_t=<H,(a-Ea)S-(jM/M_t)B>_pool`。固定 unmerged asset [`df4a64f6`](https://github.com/LightChainr/Matching-One/blob/df4a64f68232eec5aa5b8c8a5d920062aaa7808e/results/p537-one-defect-diagonal-edge/REPORT.md) 已在 N25 axis/tilted、`x=West(z)` 的 alternating radius-one selected sector 中保留 6,846 个 kernel-changing row classes、740,950 条 physical pair fibres：总 Schur signed mass 为 `-4.948839916450813e-6`，12 个 rank-stage×source-orbit cells、两项 row sums 与六项 column sums均严格非零；`0→1` 贡献 117.63%，`1→2` 抵消 17.63%。其冻结 contact mask 又给出 mask 0 的 **0 classes / 0 fibres**，mask 1/2 的负质量由 mask 3 的正质量部分抵消。因此这个 selected finite signed mass 由 local contact/OPE channel 承载，并由 first birth 主导；它不是完整 graph 或完整 `T_N`。

同一固定资产在上述 selected sector 内保留 contact mask、birth stage 与 source orbit 的联合 tensor。one-arm masks 1/2 只支持 NN source；所有 non-NN source 均落在 double-contact mask 3。把两个不完全相等的 one-arm masks 仅作为 single-contact aggregate 合并后，`0→1/1→2 × single/double contact` 的 `2×2` globally frozen pooled-root Schur signed table 为 `[-2.88380e-6,-2.93729e-6; -5.32257e-6,+6.19482e-6]`，determinant 有严格负区间 `[-3.3498535471290615e-11,-3.3498535471290614e-11]`。因此该有限 N25 contact×birth table 是 rank two；它仍是同一数据块的派生分解，不是独立 evidence vote 或完整 `T_N`。

[本分支的 physical one-defect gate](../notes/p537-one-defect-diagonal-edge-20260901.md)补上旧 fibres 缺失的构造性信息。第一条 axis N25 literal `z` flip 保存共同 `x+y+z` map 与单一背景，并令 rank `0→1`、Bell `9240712→6848576`、`g16 4→0`；完整总体先冻结 root/counterterm 后，source midpoint 为 `-1.0888815582478189e-11`，合计 `-8.298623728474635e-12`。第二条固定构型把 row-major `x=0,y=6,z=2` 的三点 NN 距离全部固定为 2，仍有 rank `0→1`、Bell `01203010→00102000`、`g16 8→0`、joint terminal incidence `2→1`；C4 pooled source midpoint 为 `-1.0121115955209059e-10`，full 为 `-9.586976893140449e-11`。两条的 beta-free source part 与 full weight 都有严格非零有理证书。它们说明 metric 非邻接仍可经同一 global carrier 接触；`d_NN≤1/≥2` 不能定义 contact/separated。总体资产与 literal witnesses 共同触发停止规则：blanket full-graph two-independent-defect / 自动 six-arm 路线已经否定，不再枚举完整 graph。

这些都是有限 N25 结论；df4a64f6 的 mask 0 空集只属于固定 `x=West(z)` 的 alternating selected sector。随后 frozen N65 20M block 在 canonical selected-carrier 分账中复现 `[-,-;-,+]` 与严格负 `Delta`；[完整审计](../results/p537-contact-stage-n65/REPORT.md)保留 6×6 covariance 与 positive exposure，selected total 是同一依赖块的派生诊断。它拒绝该分账内的 scalar/separable law，但 selected cells 在共同 thermal gauge 下通常移动，`theta=-1` 也只是符号象限恒等，不能升级为坐标无关算符或六票独立证据。同块 post-hoc [full-T secondary](../results/p537-full-t-transport/REPORT.md)给 `J65=-0.00162251±0.00018553`、`J65/J25=0.29396±0.03361` 的有限 original-`U` 收缩；它没有新独立 block，两点 power 仅为描述。最新 N145 200M 的冻结 full-T 区间跨过决策边界，结论是 `UNRESOLVED`、无 top-up；fixed-power 比较与 horizontalized remainder 仍只作 P2 两点 fingerprint。#537 后续只保留 remainder transport 的 proof/counterexample；两点指数、CFT 标签和 N145 外推均不改其 P1 生命周期。#539只保留P2复现支持；当前P0随机生产和云任务为空。

[固定m审查的新增结论](../notes/p337-fixed-m-relative-bound.md#7-2026-09-01实际组件气体的进一步取舍)则排除了裸组件气体在h=1使用标准非负KP判据的路线；任意非负控制函数都不能统一成立。rank2投影精确固定唯一绕行组件颜色，但实际两相内外partition比仍未控制，固定m原U定理没有被宣布完成。

## 已完成：完整齐次 N50 与有限零传递判决

| 问题 | 已得到的具体答案 | 尚未解决 |
|---|---|---|
| 原父图epsilon=1的无条件U/S响应能否直接算出？ | **[完整N50结果](../experiments/p337-homogeneous-n50-20260831/RESULT.md)：U=1.0615603877、V_S=+0.0543457827，严格有理界排除V=0。** 每图精确覆盖2^50配置，共约49.85 CPU秒、峰1.63 GiB；独立p导数核查12项通过。 | 端点正号预测存活，幅度不能原样延续。N50有限传递已完成，不自动扩N100/t/epsilon，也不确认连续H4机制。 |
| 0/1孔信息加尾概率能否控制齐次U？ | [严格信息不足见证](../experiments/p337-continuation-feasibility-20260831/THEORY.md)：两套保留真实0/1孔表的摘要补全，有相同完整q曲线、唯一root和正斜率，但U约为+10.10358/−10.07432。连符号都不能从所列约束确定。 | 补全不是原图的物理多孔律；需要利用真实连通规则约束未知层的取向热score矩，不能把该见证扩大为一切延续方法无效。 |
| 能否摆脱小孔概率展开？ | [全epsilon面核](../experiments/p337-face-kernel-20260831/REPORT.md)精确保留端口、位移、平行边及源修正；指定两孔构型的权重交叉比为e^t；固定B整行的全孔密度绕环概率已闭式求和。 | 条件核本身不替代总体；当前N50无条件值已由充分状态frontier精确求出，一般N/epsilon及有限t控制仍未完成。已有3dc47674 hypergraph/twist表示不是本次新发现。 |
| 原U的有限耦合反号 | [固定m=64复核](../experiments/p337-finite-law-window-20260831/RESULT.md)：N25原Ustar≈−5.82495e−19、Udrop≈+1.07107e−13，各在自己的共同root。执行[2690f665](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)并行交付了更强的全m>=64证书，符号不计作本次新增机制排除。 | 本次新增的是下行的采样可达性判断；固定N25半直线不等于N50饱和到齐次的延续。 |
| 十台机器能否直接测上述反号？ | 同一计数给rank1概率：轴向约1e−15，斜向约2e−20。95%机会至少见一次rank1的必要样本下界为2.28e14至4.47e19，**不启动该点普通无条件采样**。 | 该界不是估U的充分预算；不约束条件、importance或twist partition估计器，也不是墙钟时间预测。 |
| 看见rank1后，热协方差估计是否就容易？ | [直接计算原U分子的方差](../experiments/p337-estimator-access-20260831/RESULT.md)：即使给真实root、均值和分母，独立iid平均`(K−mu)I1`达到SNR3，star每几何仍需≥1.5180e25次，drop≥2.5247e15次。star/axis的条件K均值差只有约0.000151383。 | 只针对指定估计器，不是95%区间或所有算法下界。该事后可行性分析未添加原冻结合同的判决。 |
| 非负twist表示是否已解决成本？ | [条件数审查](../notes/p337-twist-estimator-access.md)：独立partition差分病态；自然逐样本抵消恰还原rank1事件。m2才是五个literal partitions，m64为4097个，除非另给聚合算法。 | 改进需要相关估计、条件积分或桥接权重的二阶矩/重叠控制，当前没有可自动运行的高效替代。 |

本次N50采用精确状态合并覆盖完整有限总体，未逐叶访问2^50个配置；无新MC或云作业。科学合同10c666b6、最终producer4ae4e710均先冻结后完整运行，合同至此结束。[连通消元和比界](../notes/p337-connectivity-reduction.md)保留可严格积分的开关及不能忽略的外部同调。m64比较先冻结`375a6f0c`再评分，原四点实验保持不变；后续预算是同一点的事后精确可行性分析，其余是明确数学构造和有界配置核验。它们不提供新的独立统计证据。

最新既有结果另已接入：[e1b96895](https://github.com/LightChainr/Matching-One/blob/e1b968959634b9b3999c727b83ed38d0b730cb20/results/defect-reweight/REPORT.md)给Xi_reweight=+4.550327123237、Xi_jump=−15.306045530801，排除jump-only；[f4057192](https://github.com/LightChainr/Matching-One/blob/f405719264c896aa873dd4aae7292795f544ba99/notes/topological-projection-reverses-global-u-tail.md)已给Sstar/Sdrop相反渐近尾。N100/N225仍是理论预测，不是已执行生产。

更新的[2690f665联合极限](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md)给N/m²→zeta<∞时的Poisson/full共存、原pooled分母控制和U超多项式衰减；它不是固定m的大N定理。该分支剩余的是固定m的oblique twist代价和受限扇区odds，不再等待首次有限反号或指定联合极限的分母界。

[固定m相对界复核](../notes/p337-fixed-m-relative-bound.md)给实际局部模型`Delta_k≤50k log m`，改进到对h一致的表面阶，但常数仍大于已有轮廓界允许的阈值。另给完整三态热族反例：即便共同bulk精确消去、rank1对所有h指数小且pooled根唯一，若几何间受限sector odds未控制，仍可有`U/A→2/Delta`。它只否定这些摘要前提的充分性，不是原格点模型的反例。

执行分支另于[0dda27ba](https://github.com/LightChainr/Matching-One/blob/0dda27ba/notes/closed-source-s4-trace-transmission-result.md)完成固定m2/Q4的S4 `[2,2]` seam插入：直接q/E分子恒零，原U仍经归一化与热导数得到`V_beta=+5.440121494634842e−6`。这个有限传递已完成，不再追加Q4点或seam；其epsilon是trace系数，不是孔密度或logQ。后续[bea717e8](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/p337-q1-closed-trace-transmission-result.md)已完成指定Q1闭合迹：固定beta1=−1_A−1_B给`V_beta1=−0.001904836180602413`，完整finite landing及有理证书均已交付。这些闭合迹之后，局部pair及其正则补全到原U的有限传递也已由上面的2ba8863f完成；尚未完成的是占据平均空间/尺度归属，不能再按旧的local-to-finite-U待办重算。raw-Q分项依赖归一化规范，不能按其符号认定连续场；上述结果也不提供m64高效估计器。

## 最新有限机制：端点闭合不能靠两个共同耦合延续

| 已完成的固定计算/推导 | 结果和真正排除的解释 | 来源 |
|---|---|---|
| checkerboard端点的闭合源 | `S=C+F4+Bvac=C+F4+T_NN-4K+2N`在N50→N25端点逐配置保持；其原U响应N25为+0.126165363414169。 | [源字典b8d043fc](https://github.com/LightChainr/Matching-One/blob/b8d043fc/notes/decimation-closed-source-and-global-u.md)、[精确值ec01768f](https://github.com/LightChainr/Matching-One/blob/ec01768f/results/p337-closed-source-n25/REPORT.md) |
| 被bare-C遗漏的F4 | 固定N25两商全枚举给`V_F4=+0.194414686460907`，对应N50端点遗漏+0.599656868156603；裸簇源端点U传输失败。 | [结果b8d043fc](https://github.com/LightChainr/Matching-One/blob/b8d043fc/results/decimation-plaquette-u/score/REPORT.md) |
| 一孔缺陷的原U混合响应 | `Xi=U_(t,epsilon)=-10.7557184075641`；`R=U U_st-U_s U_t=+27.7665635812302`有理界严格非零。纯温度内部延伸和源无关单一增益均失败，端点恒等式仍成立。 | [实际计算f5c4a74a](https://github.com/LightChainr/Matching-One/blob/f5c4a74a/results/p337-endpoint-defect/score/REPORT.md)、[解释bc17b81d](https://github.com/LightChainr/Matching-One/blob/bc17b81d/notes/checkerboard-single-defect-global-u-result.md) |
| **本次：共同温度＋同一源耦合的profile闭合** | 固定四维`f=(q₁,E₁,q₂,E₂)`、`T=∂pf,C=Cov(f,S),H=∂epsilon f`，全部四个三阶minor平方和**D3=0.000439154238009660…>1/10000**，有理界严格为正。不存在共同的`H=bT+cC`，允许b/c任意依赖p/t也不能修复。 | [结果与解释](../experiments/p337-two-coupling-closure-20260831/RESULT.md)、[全部有理界](../experiments/p337-two-coupling-closure-20260831/results/latest.json) |

这些是同一个有限N50/N25链上的相关精确计算，不是独立统计票、渐近H4或连续场身份。N25两商的Smith类不同。Xi已完成，不能从较旧的PR267 NEXT再次启动；共同clock商空间与01/02/12接口也已分别由[1b0ec15a](https://github.com/LightChainr/Matching-One/blob/1b0ec15a/notes/p154-original-u-clock-quotient.md)、[c2828e34](https://github.com/LightChainr/Matching-One/blob/c2828e34/notes/p154-lag1-current-commutator.md)交付。

**停止继续用共同温度和同一S有效耦合解释这两个几何的全部q/E响应。** 这包括任意平滑共同坐标变化，不只是一种选定拟合。共同root运动只给C/H加上T的倍数，不能消除该障碍。允许每个几何各有坐标、只匹配标量U，或引入独立新耦合，均属于此次未检验的其他问题；不能从D3直接声称识别了第三个CFT场。

合同、输入和代码先冻结于`76a070d4`，再作首次评分；只消费已公布的完整整数系数，0新样本、0新枚举、0云作业。已有混合U增益拒绝并不推出此结论：源坐标重标度可令旧R非零而新D3仍零。[计算合同](../experiments/p337-two-coupling-closure-20260831/CONTRACT.md)明确这个区别及失败边界。该结果与此前有限结果共享数据，不追认为前瞻独立确认。

## 当前决策：独立实验开始淘汰冻结预测

用户在本轮交付期间要求收缩研究自由度。本页是当前科学判断的唯一入口；PR267与Issue评论保存沿革，PR509交付资产，不分别维护竞争的“当前真相”。上面的时间记录和下面既有结论按其证据状态保留。

| 主实验 | 要消除的不确定性 | 当前状态 | 失败如何改变研究 |
|---|---|---|---|
| #154 temporal transmission | 早期隐藏结构的影响是否进入原global U，由哪个birth通道承担？ | `0820b8d2`冻结；N85=5M、N340=160M全部完成。净U导数为**0.04347±0.04363 / 0.06068±0.08266**；同时区间**[-0.07164,0.15858] / [-0.15739,0.27875]**均在±0.50内。[完整报告](../experiments/p154-prospective-transmission-20260831/REPORT.md) | **停止该lag=1源作为当前主要H4解释的优先投入**。强entry/强completion两预测被排除，双通道弱数值限制相容；不换lag、不补样、不称精确零 |
| #334 independent intervention | 既有contact机制能否预测新的coarse-state保持干预，20%残余是否有可迁移意义？ | `4b3c21b7`冻结；每N300k新prefix、120/120分片完成。R_new/R_old为N325 **0.4989 [0.4361,0.5617]**、N425 **0.5169 [0.4507,0.5831]**；两个预定范围±0.25及0.75–1.25均被排除。[固定得分](../experiments/p334-prospective-intervention-20260831/results/latest.json) | **停止这两个残余投影预测的优先投入**；新块不重拟合、不加descriptor、不把观测到的约1/2注册为救场模型 |

#334上表区间是每N97.5%渐近Student-t区间，两N按Bonferroni构成95%家族；结论条件于冻结的旧系数、均值和R_old点预测。它拒绝所声明的固定预测，不拒绝未知旧总体参数；单个signed-loading投影也不能证明或否定四feature完整充分性。旧数据只训练，未进入新得分。

#154的±为一个新批次SE，区间采用预定六坐标Bonferroni共同95%渐近覆盖。四个entry/completion区间全部落在±0.30内，W/B/C是三个互斥但不穷尽的有限尺寸数值限制，不是三套完整物理理论。结果触发已冻结的主线降级规则，未证明精确零、未否定所有滞后源，也未改变已有H4基线证据。两个生产实验都已结束；#154/#334保留一般问题并列P1，当前没有自动续跑的P0。

执行队另一个[已完成的新群体干预](https://github.com/LightChainr/Matching-One/issues/334#issuecomment-5477517774)使用不同的score-normal源，1M新prefix给T=(3.0852±0.3919)×10⁻⁸，触发其冻结的“停止完整两score标签均值闭合”规则。这个结果与本表的四feature残余检验各有源和目标，不能合并成一次重复验证；它也没有识别未扰动global anomaly的生成机制。两项结果共同要求下一步回到明确global读出的传递关系，不继续扩展contact目录。

执行队的[一次激活事件核4daae57e](https://github.com/LightChainr/Matching-One/commit/4daae57eef5c945aa050a95cd3d5d5d77582161b)也已完成；当前前瞻#154使用这个固定lag=1接口，不再重复sqrtN-lag或做lag扫描。规范化U的entry/completion读出分解共享完整源的根/分母导数，不能改称源的因果事件归因。

刚完成的组外prefix预测、条件形状和时序补观测全部归入**探索C2**。即使局部检验正确、训练测试按批分开，问题和模型仍经过同一archive生成；这些结果不充当最终独立确认。contact研究最新`323de7d5`显示原00约80%的signed loading落在四feature span内、约20%残留；它不是response R²或精确充分状态。[来源](https://github.com/LightChainr/Matching-One/blob/323de7d5ee4a980b3c77e1a972cb6c812a9f88e5/notes/p334-new64-feature-loading.md)

#275 已调为唯一 P0；#419/#370/#398及#1保留 support 资产。#275 只执行具名 observable/normalizer 的候选预测与可识别性，不恢复旧 C3/E_top 首次检测或任意新坐标。只有能说明“哪种结果停止哪个候选”的具体任务才进入P0；跨N旧数据回归、进一步 projection/Hessian/descriptor 及 generic certificate 仍属 exploratory/support。下一项实际分析按 [NEXT](NEXT-TARGETS.md) 的冻结判决包执行；本轮验证块不再用于选择新模型。

## 取向与物理响应

| 已有结果 | 证据与适用范围 | 结果位置 |
|---|---|---|
| **方位H4已有独立证据** | P43+P57独立primary合成对零为31.1857/4、p=2.81e−6；固定H4为3.4623/4、p=.484。norm-5区分所测试H8/H12；prism新增独立支持。norm-5子块本身对零仍相容。 | main：[综合报告](../notes/issue212-matching-odd-synthesis.md)；[prism原讨论](https://github.com/LightChainr/Matching-One/issues/205#issuecomment-5462845639) |
| **单一标量尺度修正已不足** | N145→290全曲线冻结单倍率9.3520/2、p=.0093；其零空间回溯显示cutoff-sensitive但默认拒绝不反转。S′纯幂律与单一rank-gap宽度均失败。norm-4 q2为20.897/2、p=2.90e−5；Jordan标量p=.067、全jet p=.054。第四代中Jordan加一个even-mode家族存活，但λ=0、1/2、1难区分，未选出次级模身份。 | [#543回溯](../experiments/p543-covariance-nullspace-audit-20260901/REPORT.md)；[已完成的开放PR273](https://github.com/LightChainr/Matching-One/pull/273)；[四代结果PR277](https://github.com/LightChainr/Matching-One/pull/277) |
| **微观簇源的偶响应已明确测得** | S=(黑NN簇数+白matching簇数)/N。P40百万N65/N85已补齐E×S及E×controls；q-fugacity补偿下C四几何均为正（110–138SE），共同raw源的H4联合p=.21156仍未分辨。full辅助p=.04911是相关视图，不能认定能量场。此C保持matching均值，但未做真实Bernoulli温度补偿，也不是原norm-4的根/斜率归一化U。 | Draft267：[最新百万偶响应](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/results/p40-even-given-odd/REPORT.md)、[原norm-4接口](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/notes/p40-even-response-norm4-interface.md)；[20k来源](https://github.com/LightChainr/Matching-One/blob/eb7ef8c9f13a88d96f32c2da62ba7ef2145cb33e/results/p154-absolute-cluster/REPORT.md) |
| **正的源/读出证据与M载荷分开** | 外部Euler源在fixed-K分解后仍有稳定约32%余项；F5已分辨两条源读出。原始matching M的同流载荷零仍相容（p=.58155），连接尚未分辨。 | [外部源](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5468416605)、[开放PR451](https://github.com/LightChainr/Matching-One/pull/451) |
| **primitive square-bond另有多character结果** | N112独立生产支持r0+r1；E_top同流方向行列式p=2.49e−5，说明两个拓扑观察方向。纯E4/E6/E4²被拒；这套square-bond观察量与square-site thermal主线分开保存。 | [N112结果](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5469692921) |

P40实际引擎在随机键中包含N：N65/N85是不同N-domain，通常PRNG独立假设下可作nominal联合统计；同N两方向共享随机数。早期仅凭seed标签判断跨N共流的说法已经更正。

P154的源全链、两阶段及百万端点已经完成。条件line响应进一步拒绝具名E-plus-clock模型，固定K/rank1空间关联也强；这两项在[Draft267固定提交](https://github.com/LightChainr/Matching-One/blob/764595ea5c838c110e416382a3a90e2ecf7297bb/results/norm4-source-line-fixed-k/REPORT.md)。本轮[角权桥](../experiments/p154-spatial-localization-20260831/REPORT.md)用同100k/1M源子集完成六N U±/v±及全部共同协方差。U−中心值均负，原总source链仍未分辨；固定K/rank1内中心化的空间源对原U严格为零，两个角权分配相反。这个精确抵消约束说明，仅加强同层O4关联不能识别原全局H4机制；后续已转入下述rank人口/进入退出及独立传递实验。源子集误差不替代原高精度普通生产。

**新的时序机制读数已经完成。** [单lag路径源](../experiments/p154-temporal-source-20260831/REPORT.md)在原2.4M排列上补出L=max(0,K−ceil√N)的早期簇数与rank。按早期rank中心化后，固定p的一阶直接响应保持早rank分布，但后期进入/退出均明确为负；早rank1不影响首次进入，只影响退出。沿共同matching根，N260的早rank0/1人口贡献为−0.06240±0.00048与+0.04733±0.00027，N340同样部分抵消。这支持rank以外的早期结构改变后期拓扑，但新的原U导数在N260/N340仍为0.843±4.882、12.249±9.922，尚未定位H4来源。这是K依赖的正路径测度源，中心化采用经验条件均值；不是旧同层源的分解，也不是逐路径hazard或有限强度结论。原三组100次共同删批与670维协方差完整保留。

## 完整分布：N100、N400、N900均已完成

三个尺度分别为2M、8M、32M共享counter块，尺度之间独立，块内形状配对。

| 读数 | 已有答案 | 现在可用的信息 |
|---|---|---|
| 共同density-map必要条件 | N100 53.914/6拒绝；N400 3.901/6、p=.690 | 精度与尺度都改变，不能把未拒绝写成模型恢复。 |
| N900冻结宽度预测 | 实测Vz=2.339461729±.120385；quarter-width p=.134967，fixed-z p=.084182 | 两者均存活；比较共享N400锚点和同一N900目标。 |
| 正三中心低矩表示 | 早侧权重依次.1806、.0654、.0320；N900未用于构造的m7/m8残差3.566/2、p=.168 | 肩部变薄、向中部重分配。三个经验中心不等于三个物理场。 |
| 两个平移的共同对称正核 | N100/N400所需第六矩为负；N900也要求−2.142812±.481800 | 所声明候选类在三个有限尺度均有障碍；继续只调权重/间距不能修复。这些探索性矩估计不等于已校准的边界证书。[N900来源](https://github.com/LightChainr/Matching-One/commit/b6db7ba57c3c5bcb6e25558b5274f08aeef1ce63) |

结果在开放[PR484](https://github.com/LightChainr/Matching-One/pull/484)：[N900完整报告](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/REPORT.md)、[第三尺度形状](https://github.com/LightChainr/Matching-One/blob/54430ea7/notes/p267-n900-three-center-shape.md)。N900无需重复启动，也没有已选出的渐近宽度指数。

## 过程与传播

| 已完成的数学/计算成果 | 直接意义 | 来源 |
|---|---|---|
| **digital Alexander与整数饱和** | M=P₂−P₀；K_minus/K_plus是两个essential births，rank-one方向固定，iota=1。有限matching根是阈值rank经Bernstein/Beta变换所得连续激活分布的等权混合中位点。 | 规则cell证明已main；一般有理/积分证明稿在[73d4960、c1a72e5讨论](https://github.com/LightChainr/Matching-One/issues/269#issuecomment-5466825850) |
| **完整birth机制与反事实干预** | 147个固定真实prefix完整law已求解。指定两例的1个/6个middle sites命中所有order≥3最小trigger，屏蔽后完整law一致。 | [147 clocks](https://github.com/LightChainr/Matching-One/blob/87b6ca5b39084c06143f31cafdaba53f90012e27/notes/p334-all147-real-prefix-clocks.md)、[middle干预](https://github.com/LightChainr/Matching-One/blob/0143632db59d867cfb658a6ad4465e5036684fff/notes/p334-middle-bridge-physical-interventions.md)，独立分支 |
| **均值clock与空间波动可区分** | uniform blockade平均响应由完整clock决定；位点影响浓度和same-mask replica包含额外空间信息。两真实prefix等生存率比较中，较平clock仍有高30.6%的E1。 | [精确噪声桥](https://github.com/LightChainr/Matching-One/blob/614eedb2429d74d6b4de7ebf15d6c8f918b54e3c/notes/p334-isoclock-positive-noise-spectrum.md)、[真实prefix结果](https://github.com/LightChainr/Matching-One/blob/795908fbc9a781a0cda704864c237deaf0327f37/notes/p334-real-prefix-iso-survival-noise-energy.md) |
| **正权传播已有具名路径** | width4/5/8实际传播已做。T4首个自相关增量为第四阶；后续慢极点/权重分解也已完成，尾部修复主要来自慢极点移动，权重变化部分抵消该作用。删除current仍保留快慢反转。 | [T4传播桥](https://github.com/LightChainr/Matching-One/blob/074a5f537caecac9cbd663dcc76ebd05ff54f302/notes/p398-width8-T4-schur-bridge.md)、[极点分解](https://github.com/LightChainr/Matching-One/commit/1f19fc1a2d9fc59dce650e95268c716762725985)，独立分支 |
| **局部高阶拓扑已有测量** | P437固定五键20k新背景给14.97SE；同块分解约99.8485%局部能量在degree≥6。 | 独立分支：[固定五键实测](https://github.com/LightChainr/Matching-One/blob/386db0a74a44be37403c666b27e1c023b81ea459/results/local-20260831/P437-N112-fixed-S5-20k/REPORT.md)、[同块分解](https://github.com/LightChainr/Matching-One/blob/888af29d58c72f113cf7cb5f80247a81a91b9273/results/p437-fixed-support-coherent-decomposition/REPORT.md)；PR437已合并的是较早的精确filter工具。 |

**P334总体配对分析也已完成。** N325/N425各20k配对counter，共40批。所有counter进入分母，目标为checkpoint rank1分层对F2及其积分的贡献，尚未覆盖rank0/rank2和完整F2。对该贡献的H4方向差，可移除后缀噪声占原个体观察方差的估计比例，canonical为49.15%/50.03%，integrated为0.816%/0.681%；比例定义为mean[(X−Y)²]/sampleVar(X)，并非标准误下降率或运行加速比。两个N的这一H4贡献均未分辨。N325/N425分别47和164个困难配对保留双向原观察，不丢样或单边替换。[完整报告](https://github.com/LightChainr/Matching-One/blob/c3bb43f1b078c5f9f76f71b25cdb3e2e331eb115/results/p334-paired-clock-loading/REPORT.md)，独立分支；后续已完成的分解见下。

**R1之后的完整P334分析也已推进。** 全A/E九层、128万辅助续接、SS/mixed/BB、接触坐标、共同Euler不可见源和全热曲线均已完成，旧mask总体rank2不等于同prefix rank2。本轮[有限q_t](../experiments/p334-finite-source-20260831/REPORT.md)在t=±1仍有可测future S(A)/D(A)响应，即时两rank与Euler增量的联合分布保持不变；使用原数据的精确importance估计，没有按新策略重采样。完整census恢复了局部检验的抽样支持：全部1502/1551双R0 prefix有两个独立允许源，固定781568条定向续接进一步给出A局部行列式均值约5.4/7.6SE的正证据，积分A对应的两个出生中心也支持局部二维响应。E/间隔和四阶平方量仍弱。[最终报告](../experiments/p334-mechanism-response-20260831/REPORT.md)明确区分E[det J(Z)]和det(E[J(Z)])，并保留原20批。

**P334的预测与条件形状已继续完成。** 同一cell00原prefix、旧8+新64续接上，以完整census源Gram G(Z)预测局部J(Z)=BG(Z)，与同为四参数的常数矩阵比较。原20批分为五折，全部删一重拟合；组外A平方误差降低38.17%±9.33个百分点、39.37%±10.59个百分点，出生中心误差降低54.42%±15.47、58.65%±15.71个百分点。这支持微观特征携带可预测信息，尚未证明跨N运输或完全闭合。[预测报告](../experiments/p334-prefix-prediction-20260831/REPORT.md)

同一批续接的无偏条件协方差检验中，minus→D的出生中心方差响应为−2.43469e−8±5.32765e−9、−1.27127e−8±3.30315e−9，两N全部20次删一保持负号。固定prefix的两出生位置若仅作确定性平移，此量应为零；该纯平移类已不足。总体形状能量仍未分辨，且不能断言每个prefix同时具有rank2和形状变化。所有读数仍是cell00对原20000分母的贡献，未把未测其他cell置零。[条件形状报告](../experiments/p334-conditional-shape-20260831/README.md)

147-prefix噪声预算仅适用于固定经验mixture。P398的93维Krylov空间和Boolean谱degree按各自生成过程解释；现有width4/8的i^j权重保持波长4，并非固定模式序号的尺度外推。

**P398固定干预与解析线性响应均已完成。** 固定η=0,±1/4干预出现cross传播，旧16维几何字典整体近似改善但未闭合；随后[精确η0导数与零频分析](../experiments/p398-linear-response-20260831/README.md)给出U′+−在旧lag括号内的反号时刻t≈1.04798965。负平稳重加权与正动态项竞争，零频两cross仍为负；16维模型积分误差0.467%/0.551%，两源模型不反号。投影借用完整π及π′，不声称盲预测、精确闭合或square-site映射；保持P2。

## 后续纠错已经改变的解释

| 旧说法或候选 | 现在应采用的结论 |
|---|---|
| P418巨大共同谱惩罚支持radius flow | batch sum统一为per-sample后四共同谱均相容；radius5单壳数值不可靠。正确归一化的P250秩与exact CRT不受此错影响。[修正报告](https://github.com/LightChainr/Matching-One/blob/e2b57aa7c5ec5c7db8cbb4f03872435f20966407/results/p418-normalized-archive/REPORT.md) |
| #43 even通道原冻结检验直接通过 | 原cross/either错配是协议失败；精确符号运输后的.57003/2属于事后确定性纠正，历史仍保留。[纠错](https://github.com/LightChainr/Matching-One/issues/133) |
| 非零q/标记耦合即可识别新场 | q仅有三值；q-only单根接触耦合已有全阶闭合。真正独立物理识别需要更丰富的微观读出。[no-go](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5463228151) |
| 一个低阶Hankel秩或形式Jordan足够 | rank随生成元与观察窗改变；旧R2幸存者被更高阶/联合map分析否定。Q=1端点也不唯一指定generic-Q切向。[完整链](https://github.com/LightChainr/Matching-One/issues/250)、[Q-lift](https://github.com/LightChainr/Matching-One/issues/333) |

## 支持资产与统一解释边界

阈值rank引擎、协方差、exact certificate、Q4表示、W5周期配对及15态terminal serial代数均已可用。#498–502新增的是有限子群/作用/理想/闭合集分类。W5已实现的两图并非自对偶，当前没有新的严格p_c界；这些资产接具体概率比较时才成为阈值推进。

#1三次四区间有限排除已经main，发布前复查时四次Jacobsen PR524也已由其他执行者合并，Mertens p-med PR525已打开。该路线没有推出次数或高度上界的理论，区间也不是已证明包含真实阈值的严格界。保留P2候选验证能力，暂停自动扩搜，具体事实、最近见证与恢复条件见[专项审查](../notes/cubic-search-review-20260831.md)。

目前尚未给出square-site p_c闭式、唯一连续场/Jordan身份、总体方向响应的几何归因。相容性不等于模型确认；有限图证明、有限尺度数据、表示论选择规则和晶格到连续场重叠各按其实际范围使用。重分析共用原块，不增加独立证据。

旧状态全文保存在[整理前固定commit](https://github.com/LightChainr/Matching-One/blob/8a68cca866d7fbca7463e2167c3ff06128d5851f/docs/STATUS.md)。本次没有改动旧结果、冻结文件或其他分支；当前分析顺序见[下一步分析](NEXT-TARGETS.md)。
