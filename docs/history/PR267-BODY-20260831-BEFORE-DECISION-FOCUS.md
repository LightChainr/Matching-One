# PR #267 body before the two-experiment focus

Captured 2026-08-31 from PR #267 while its observed head was `5db502fe4031abda153058d378033cdecd454fd0`.
The text below is preserved verbatim as historical PR-body content. Its old
next-step statements are not the current queue. Use [Next Targets](../NEXT-TARGETS.md)
and [Decision Experiments](../DECISION-EXPERIMENTS.md) for current work.

---

## Conditional winding source and fixed-K localization — 2026-08-31 i

本轮已经完成新的空间源分析，科学commit d2a3f445abcb99f9a21bc7b5cd70b5a22e6a77c7；仍交付于Draft #267。

1. **共同E＋热时钟模型不足以解释rank1内部绕环响应。** [条件绕环源报告](https://github.com/LightChainr/Matching-One/blob/d2a3f445abcb99f9a21bc7b5cd70b5a22e6a77c7/results/norm4-source-line/REPORT.md)在六个原norm4尺寸测得四实分量残差，逐N的nominal chi-square/4为3843.85、2596.77、1049.13、1378.11、8479.76、4629.15。所有根位移、rank人口校准和同批协方差均已传播；不是只看一个非零方向分量。完整396维协方差与三组100个paired delete-one向量已归档。
2. **大四重方向响应需要按几何解释。** 将每个复数响应精确乘以(a−ib)^4/N²，环面坐标系里的实幅度约.089–.110，虚部绝对值不超过.000117。这个确定性事后读法表明强实验室坐标P4含有已知环面旋转，不能直接充当原global U的格点H4算符身份。它给出了更具体的研究对象：rank1内部primitive-line组成的物理源响应。
3. **强源作用已经定位到固定K的空间排列。** [固定占据数分解](https://github.com/LightChainr/Matching-One/blob/d2a3f445abcb99f9a21bc7b5cd70b5a22e6a77c7/results/norm4-source-line-fixed-k/REPORT.md)使用同一新档案，分离同K项、跨K组成项及根移动项；同K四实分量在六N均明确非零（联合nominal chi-square=22534.76/24，最大单分量56.27SE）。固定(K,rank1)后，任意每几何各自的f_g(K,q,E)都恒定，因而这类粗粒度解释不足以覆盖实际绕环源作用。Bessel补充最大变化仅2.62e−6；288维联合协方差、完整K域及低计数覆盖已保存，没有再次回放。
4. **这轮计算已经实际完成。** 112.14秒本地CPU给2.4M旧排列补出此前缺失的物理line/source乘积，primary聚合1.31秒；N65/85/130/170各100k，N260/N340各1M，端点仍是原1000+增量9000的同编号union batches。新counter=0、服务器/GPU操作=0、测试套件=0。新旧视图和几何旋转不是新的独立样本。

总览和下一步已更新：首次条件line和fixed-K测量不再是待办。下一项直接复用现档案的W±=R(1±Re[e^-4iθ O4])/2软角权，将原global U及源切向拆成严格可加回的两项，寻找强空间源作用在全局读出中的传递或抵消。W±不是离散primitive-line类别；真实分类和birth/lifetime归因另需相应联合标记。该有限空间源结果不自动给出连续场身份或norm4缩放律。未合并PR、未关闭/锁定/降级Issue、未改写历史。

<details>
<summary>Previous PR body preserved as historical context</summary>

## Source precision and root-comoving topology — 2026-08-31 h

本轮继续实际源分析，科学 commit b77f3e2a56df2b033561e6fe9164d09ab4698949，交付仍在 Draft #267。

1. **一个新的明确机制结果：恢复matching平衡后，rank-1人口仍然改变。** [端点1M与根随动组成报告](https://github.com/LightChainr/Matching-One/blob/b77f3e2a56df2b033561e6fe9164d09ab4698949/results/norm4-source-endpoint-1m/REPORT.md)给出共同微观簇源s=CB+CW的dP1/dlogQ：六N为−.09250到−.09782，各约73–160SE；N260为−.097823±.000612，N340为−.097110±.000683。纯两方向共同affine-K热时钟预测它严格为零，因此该有限解释已经不够。pooled根仍移动约+.029，两端rank0/rank2概率等量增加；这不是能量算符身份的结论。
2. **实际增加缺失的源信息。** N260/N340各补标后续900k旧生产排列，从100k增加到1M；前100k未重读，新增随机counter=0。固定增量本地两CPU任务约98.78秒完成，聚合1.44秒；没有服务器或GPU。原H4源响应现在为−5.9723±4.2736、+11.8594±9.1981，SE由19.98/25.87下降约4.68/2.81倍。四个cyclic尺寸继续复用原100k。
3. **模型区分已经实际算过，结果保持开放。** q2/Jordan源刚性延拓的两链nominal p=.48573/.58735；独立补集锚定的共同生成元漂移为−3.9349±3.5467、−3.9476±4.1631。新增共同q、共同E的2×2响应行列式与共同q+E+clock的3×3行列式均未分辨；近暗q列使必要零条件尤其不能被当成机制识别。未扰动全生产补集q2仍有p=2.282e−5的张力，源切向弱结果不能覆盖它。
4. **已消耗完一个准备性猜想。** [两阶段基准/交叉拟合的实际结果](https://github.com/LightChainr/Matching-One/blob/b77f3e2a56df2b033561e6fe9164d09ab4698949/results/norm4-source-two-phase/REPORT.md)在零回放2.52秒分析中没有稳定降噪，raw/two-phase SE比.898–1.008；几乎全部方差来自源标记。下一步转向[具名共同源映射的可证伪预测](https://github.com/LightChainr/Matching-One/blob/b77f3e2a56df2b033561e6fe9164d09ab4698949/notes/norm4-common-source-response-determinants.md)和rank1内部的条件原始绕环方向响应，不再把普通q/E增精度列为首要行动。

完整348维联合协方差、所有paired delete-one向量、旧100k与新1M的嵌套差、原始CSV和哈希均已保存。端点基准从full1B排除完整已标记1M，避免把新900k与旧补集误算成独立；r1/r2是静态累积激活的根位移分解，不能直接称出生局部机制。没有测试套件、Issue生命周期变更、合并、历史重写或其他服务器操作。

<details>
<summary>Previous PR body preserved as historical context</summary>

## Actual source analysis — 2026-08-31 g

本轮直接完成了新的源分析，继续交付在 Draft #267，科学 commit 29f339faa71064436b4350027f258dcccc43603d。

1. [原 norm-4 两条完整链的源响应](https://github.com/LightChainr/Matching-One/blob/29f339faa71064436b4350027f258dcccc43603d/results/norm4-source-thermal/REPORT.md)：65→130→260、85→170→340，各重观察100k条旧生产排列，逐K恢复簇数与拓扑联合量。原U的根移动、热导数和斜率响应均已计入。六N共同paired-cluster fugacity根响应约dp0/dlogQ=.029（106–182SE）；这不是把E_top直接认成能量场或把有限根认成无限pc。
2. 源刚性延拓与一生成元漂移已有实际数值：q2/Jordan两链源残差nominal p=.56165/.58952，漂移行列式66.64±102.15 /65.18±104.71，当前100k子集精度尚未区分它们。共同微观源导数是v_N=N*Udot_density；不能给密度源套错跨N坐标，也不能让低精度子集覆盖旧全生产结论。
3. [热时钟商空间与源拆解](https://github.com/LightChainr/Matching-One/blob/29f339faa71064436b4350027f258dcccc43603d/notes/p40-thermal-clock-source-quotient.md)已进入实际计算：共同a+bK对原U的响应严格为零；同一个簇源分成三状态再加权和扇区内热信息。NZ N85的W=+2.088±.882，在[独立P40百万旧流](https://github.com/LightChainr/Matching-One/blob/29f339faa71064436b4350027f258dcccc43603d/results/p40-source-thermal/REPORT.md)中为−.973±.760，没有稳定同向支持。两个估计器及组分不算独立证据票。
4. [下一项可直接执行的两阶段估计](https://github.com/LightChainr/Matching-One/blob/29f339faa71064436b4350027f258dcccc43603d/notes/norm4-source-two-phase-estimator.md)：用旧1.9B/1B q/E档案减去已标记100k，取得独立补集基准，再交叉拟合源残差。尚未运行；它无需新配置，可直接回答当前终点大误差有多少能由既有数据消除。

完整JSON、方向/跨N分组协方差、所有delete-one向量、CSV、代码和运行回执已提交。四个cyclic N共享原counter域，两个HNF端点分别保留独立seed；P40使用其自身N-key分组。新增随机counter=0；本轮没有测试套件、GPU、服务器、隧道操作或任务生命周期变更。

<details>
<summary>Previous PR body preserved as historical context</summary>

## Research delivery — 2026-08-31 f / measured even tangent and completed N900

本轮交付了明确的物理源正结果，并纠正了原norm-4桥接与N900状态。

1. [P40百万配置 E_top/source 实测](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/results/p40-even-given-odd/REPORT.md)：固定matching时，绝对簇源的偶响应C在四几何均为正（110–138SE），排除该点纯q一阶拓扑切向；raw H4 p=.21156仍未解析，full辅助p=.04911保留相关视图边界，不认定能量场。缺失矩已补，49×49协方差/100delete-one/CSV/单次运行回执都在本Draft。
2. [原norm-4源接口](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/notes/p40-even-response-norm4-interface.md)：正确目标是原pooled-root、slope-normalized thermal U的源导数 L4[U_dot]，不是裸 L4[C]；链为65→130→260及85→170→340，显式保留根/斜率/三重矩响应及机制自身的source-deformation预测。
3. [N900完成](https://github.com/LightChainr/Matching-One/blob/679f96cdbc1c293203ed9652c43386a55eab8bfc/notes/frontier-increment-20260831f.md)：32M/800，Vz2.33946±.12039，两个条件预测p=.13497/.08418无赢家；旧running标签已过时。#484保留原非Draft状态，本PR267继续Draft。
4. [下一步注意力](https://github.com/LightChainr/Matching-One/blob/679f96cdbc1c293203ed9652c43386a55eab8bfc/docs/NEXT-TARGETS.md)与README/STATUS/MAP/ROADMAP/HYPOTHESIS-BOARD、ledger v29均据这些结果更新；不再把首次C测量或首次N900列作待办。

两个旧百万counter块只重观测一次，共约42秒本地CPU；没有新样本、重复科学测试、GPU、服务器或隧道操作。只做一次交付元数据/空白检查。没有合并、rebase、force-push、Issue关闭/锁定/降级/改名或labels/assignees更动。历史正文完整保留，注意力不构成许可。

<details>
<summary>Previous body preserved — historical context</summary>

## Research delivery — 2026-08-31 e / completed million-sample response and physical next readout

本轮已经完成上一版列出的P40百万样本生产分析，并据真实结果推进下一步；没有停在工具准备。

1. [P40实际报告](https://github.com/LightChainr/Matching-One/blob/a4cbf02a48c3f78ee8fb3a1e4141bd985c0bf845/results/p40-absolute-cluster/REPORT.md)：各方向绝对簇源的global matching耦合很强，但H4方向差未解析（full-control joint nominal p=.38964）；保存完整27×27协方差及对齐delete-one。逐几何补偿源与raw源分开解释，配对抵消分项不命名为物理机制。
2. 实际引擎的N随机域分离更正了旧shared-seed推断。N65/N85可在通常PRNG域独立假设下做nominal联合统计；同N方向与所有派生读数保留原依赖。
3. [下一步](https://github.com/LightChainr/Matching-One/blob/18dfe6c5bec59577638001ebf66b4b8ed37066fc/docs/NEXT-TARGETS.md)已改成具体的even-given-odd物理响应：保持matching均值，测绝对簇fugacity是否移动rank-1概率；首报共同raw源，补已有counter缺失的E*S/E*controls。不是再做首个q/source Gram、另写scorer框架或宣称能量场已识别。
4. [新前沿](https://github.com/LightChainr/Matching-One/blob/18dfe6c5bec59577638001ebf66b4b8ed37066fc/notes/frontier-increment-20260831e.md)吸收uniform-blockade均值clock闭合（d53db2f，open #484）及完整Je观察仍缺hidden reversible transport（33c6028f，branch_only），两条路线均移除已完成的首次待办。

约4MB现有矩阵、约.11秒本地CPU计算：当前无需GPU。没有新随机样本、重复科学测试或任何服务器/隧道操作。仅做一次交付元数据/空白检查。继续同一Draft PR，不合并、不重写历史；没有Issue关闭、锁定、降级、改名或labels/assignees变动。原文保留如下，注意力不构成许可。

<details>
<summary>Previous body preserved — historical context</summary>

## Research delivery — 2026-08-31 d / completed physical source and updated handoffs

在补读仓库、Issue/PR后，本轮完成了一项真正缺失的物理源分析，并把新研究结果接回总览。

1. [absolute cluster source实际结果](https://github.com/LightChainr/Matching-One/blob/eb7ef8c9f13a88d96f32c2da62ba7ef2145cb33e/results/p154-absolute-cluster/REPORT.md)：旧N65/N130各20k配置，E_top响应未解析；源保留65–70% clock/Euler之外方差。N130 matching2.47SE仅是辅助线索。CSV、JSON、39×39协方差、单次本地运行回执均已提交。
2. [前沿增量](https://github.com/LightChainr/Matching-One/blob/5a307ad1b374c173891ef71fe1ff26386a3360f6/notes/frontier-increment-20260831d.md)：全部147 clocks、实际noise budget、canonical roots、P398 current deletion及共同对称two-lobe障碍均已读实际报告并整理，不再列为首次待做。
3. [下一步注意力](https://github.com/LightChainr/Matching-One/blob/5a307ad1b374c173891ef71fe1ff26386a3360f6/docs/NEXT-TARGETS.md)：原norm-4物理身份仍第一位；直接可用的下一读数是P40百万样本q/cluster Gram，不把q当q²、N85当N130或同seed尺寸当独立证据。
4. N900作者报告running/noresult，未独立检查进程、未重复生产。

仍为同一Draft PR：没有合并、rebase、force-push、Issue关闭/锁定/降级/改名；未改labels或assignees。没有新随机样本、重复科学测试或服务器/隧道操作。原文完整保留在下方，优先级不是许可。

<details>
<summary>Previous body preserved — historical context</summary>

## Context reconciliation — 2026-08-31 c / completed work before next assignments

本轮按用户要求先补读并整理仓库、Issue、PR，再决定后续；没有开始新的科学计算。

[已完成成果与真正缺口](https://github.com/LightChainr/Matching-One/blob/7c5e36902bc4882d885c8e6ae9ef31ede9d58be6/notes/context-reconciliation-20260831c.md)、[按用途分组的23个开放PR（带时间戳）](https://github.com/LightChainr/Matching-One/blob/7c5e36902bc4882d885c8e6ae9ef31ede9d58be6/docs/ISSUE-PR-INDEX.md)和[下一步注意力](https://github.com/LightChainr/Matching-One/blob/7c5e36902bc4882d885c8e6ae9ef31ede9d58be6/docs/NEXT-TARGETS.md)已同步：

1. N400已完成；rank-clock去平滑、普通scalar及density transport分别记录。N900有既存freeze/runner，运行状态未核实，不重复启动。
2. P398已到memory/triplet geometry；不再派首次hidden-force分析。
3. P334已有完整物理/标记/条件K2时钟、冻结12前缀与47通道race、一般parallel-two-port定理；未完成的是总体方向加载等，不是首次longer-horizon。
4. tiny VJS完整三项Q导数与#370真实生产证书已完成；#493–497已合并，属于可复用支持。
5. 原norm-4次级物理响应仍为首要开放问题；所有候选和并行路线保持开放，成果引用保留branch_only/open_pr/main_integrated区别。

本PR仍Draft，不合并main、不新建重复PR、不关闭/锁定/改名或降级Issue，不改labels/assignees。保留原正文与来源历史；没有新MC、重复测试、服务器或隧道操作。

<details>
<summary>Previous body preserved — dated historical context</summary>

## Scientific output update — 2026-08-31 second continuation

本轮先吸收 #469–490 全部新正文/讨论和PR head/path，再单独读 PR #491 的物理cut结果；旧完整捕获保持原日期，不伪装成live总数。

三个新科学输出已经进入本Draft：

- N100中心/两尾分析 d973a39：中央A面积与E dipole明确，中央A dipole未解析不代表无响应。输入是PR484同2M counters，无新采样。
- P398固定读出 4846adf：C0-only landing创新显示近距离快模，不能把.90%整矩阵误差泛化成所有观察者都准确。
- [P418统一单位生产重分析](https://github.com/LightChainr/Matching-One/blob/e2b57aa7c5ec5c7db8cbb4f03872435f20966407/results/p418-normalized-archive/REPORT.md)：继承约定下原共同masked大惩罚消失；独立旧P250/CRT结果不受牵连，优化器数值边界另列。

[当前总览与下一步](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/docs/NEXT-TARGETS.md)同步PR484/485、width8 552c45d与cut-network PR491，不继续布置已完成的firstproduction/firstfilter/firstcut。六项新#370工具和有限serialmonoid工作归入可用支持，研究注意力仍面向机制和实际读数。

保留Draft/open，不合并、不关闭/锁定/改名Issue，不变更labels/assignees。没有新格点MC、测试套件或服务器/隧道操作。

<details>
<summary>Previous body preserved — historical context</summary>

## Latest scientific continuation — real analyses after context recovery

Draft #267 remains the only delivery PR, open and unmerged. The earlier context recovery is preserved below as history.

### Norm-4: a new mixed source on the archived samples

- 新的局部相互作用分析已完成：Draft #267 / [`c0880c2`](https://github.com/LightChainr/Matching-One/commit/c0880c297b40699563e8be537e777ac8cd4084c8)，[完整结果](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/results/p154-fixed-k-interaction/REPORT.md)。按原 seed/counter 重放 N65/N130 各20k配置，只新增 E_top × edge 的缺失混合矩；400个batch/方向行的原计数和rank sums完全一致，没有增加随机样本。
- 将局部 pair 源 Q 拆为 fixed-K 几何 R 与 occupation-only H：`J_Q=J_R+J_H`，H 的响应是明确的二阶热导数。P4[J_R] 两size z=.152/.891，joint chi²=.8173/2、nominal p=.66455；目前未解析，**不能**说已选中第二场，也不能说由 H 完全解释或所有局部源无效。
- P40 的 q×motif 二阶Gram不包含 q²×motif；这里 q=I2−I0，E_top=q²，不是平均matching函数的平方。因此本次不是重跑旧降噪score。Q是局部canonical势，R有全局K反项；二者的物理含义与matching图对变换在[源/读出说明](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/notes/p154-local-source-versus-connectivity-emission.md)中明确。

### P398: exact state dimension versus effective slow modes

- **正权物理两点矩阵已经完成。** branch_only [`e38fe76`](https://github.com/LightChainr/Matching-One/blob/e38fe7634354b0cb2201fa55fd9b4d37ccedeef2/notes/p398-positive-cylinder-propagation.md) 给出width4真实AP/landing C(d)及两个非零普通本征值；[`b35e100`](https://github.com/LightChainr/Matching-One/blob/b35e100a3903c706dceba57c4667386eb4510ac3/notes/p398-anisotropic-cylinder-spectrum.md)完成完整正权h/v族；[`dbd4081`](https://github.com/LightChainr/Matching-One/blob/dbd408154b4215ca41fbf26c0fd962997074f05d/notes/p398-continuous-two-channel-kernel.md)给出同宽连续距离fingerprint。旧正文“物理矩阵尚未产出”不再是当前状态。
- 新 Draft #267 / [`8f7a587`](https://github.com/LightChainr/Matching-One/commit/8f7a5875157265e32e9db08c5f3991a9b9ddb86e) 已把同一adjacent-pair/singleton微观指标与first cyclic character延拓至width5，正权42态、h=v=1/2。实际正间距两点Hankel精确秩8，传播多项式无重根；这是有限圆柱传播，不是8个CFT场。[精确结果](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/results/p398-physical-two-point/REPORT.md)。
- 消费该矩阵得到一个更有用的区别：两最慢模不拟合residue，整矩阵相对误差在d1/d2/d4为 .8996%/.1523%/.004050%；但 U=C0^-1C(d) 的 U2−U1² 相对缺口7.8838%。精确状态、有效慢模与原两读出的自主闭合不是同一个概念；whole42state过程仍是Markov。[数值分解](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/results/p398-physical-two-point/mode-visibility.md)。

当前[注意力顺序](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/docs/NEXT-TARGETS.md)与[本轮科学进展](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/notes/physical-analysis-progress-20260831.md)消费这些完成结果，不另设准备阶梯。未运行测试套件、未用服务器、未增MC样本；优先级不是许可，不关闭/锁定任务、不自动合并PR。

<details>
<summary>Earlier proposal and context review — historical text preserved</summary>

## Current handoff: full context recovery before further science

**Draft #267 remains open and unmerged.** This updates the existing scientific-information recovery and production-analysis PR; it creates no replacement PR. Priorities allocate attention, not permission.

### Start here

- [Repository Context](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/docs/REPOSITORY-CONTEXT.md): completed results versus genuinely remaining questions.
- [Next Targets](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/docs/NEXT-TARGETS.md): the single current attention board.
- [Issue/PR index](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/docs/ISSUE-PR-INDEX.md): all 464 numbered items from the initial snapshot, plus a separately dated #465–468 increment.

The initial capture contains 146 Issues, 318 PRs, 1,354 discussion comments and seven formal reviews; all bodies/discussions were read by three parallel readers. At that capture, 109 Issues and 19 PRs were open. Full captured text is preserved in the [compressed archive](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/analysis/github-context-20260831.tar.gz), with hashes and exact PR heads in the [inventory](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/analysis/github_context_inventory.json). Later #465–468 and new #13/#14/#370 discussion are bounded addenda, not a silently regenerated snapshot. Key source reports were inspected; this does not claim every code file or cited paper was reread.

### The central scientific correction

Canonical E_top Phase D and the first #370 real-production confidence application are **done**. E_top=P0+P2 is a topological coordinate, not an identified energy operator. The integrated proxy aliases its derivative, P4[S']=P4[E_top']/2; the measured mixed radius-one row was not selected. **The original norm-4 secondary physical response remains unidentified.**

The review also restores completed F5 source separation, N112 C3 E_top response, P334 trigger-graph/capacity work, P437 localized high-order topology, the ordered-filtration proxy, motif projection, coalescence, R8 thermal-null and W5 periodic gluing. Their observer, measure, dependency and integration boundaries remain explicit. Inconclusive measurements are not unrun tasks.

### Current attention

1. Identify an explicit physically distinct microscopic singlet/scale/readout against the original norm-4 residual and the surviving topology/clock plane.
2. Apply the completed rooted/landing module to a positive-weight connected two-point response matrix (#398), not more automatic jet closure.
3. Explain the P267 joint response map after amplitude-only profiling; then the P334 trigger graph's topological cut/landing sides.
4. Transport the already resolved F5/C3 observer responses and address genuine matching loading, complete-spectrum, ordered-context and concrete W5 comparison questions in parallel.

The newer #466/#468 strategy proposals are retained and cross-referenced, with stale first-step requests corrected. Newly claimed #370/#13 tool/algebra tasks may proceed, but they do not become physical-analysis prerequisites. This is not an order to stop those tasks.

### Work completed before the context-review pause, retained without rerunning

- [P267 response-ray](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/results/p267-response-ray/REPORT.md), `a09758e`: full A/E/C/W amplitude-only nominal 19.87177/3, while A/E alone remains compatible. Dependent post-reveal finite response comparison with estimated covariance, not exact finite-sample coverage or a field identification.
- [P334 capacity allocation](https://github.com/LightChainr/Matching-One/blob/research/navigation-priority-refresh-20260829/results/p334-trigger-capacity-allocation/REPORT.md), `2e32fd0` (original `52b68b3`): selected Delta W2=540 decomposes into support/block −195.657143, side constraint +649.285885 and residual organization +86.371257. Same saved graphs, not independent replication or a causal H4 share.

Earlier P334 directional allocation and P439 direct/plateau transport remain in this PR. All source statuses retain main_integrated / open_pr / branch_only / hypothesis distinctions.

### Delivery and operational scope

Navigation recovery commit: `5b8ac329ea980c362fa2cd517e152be016532536`, followed by the editorial receipt/addendum commit on this same branch. Main was brought into the Draft by ordinary merge through `e300609`; later main PR #467 at `e46b00f` is recorded separately, not claimed as already incorporated here. Ledger v25 has 247 sources, 101 scientific nodes, 15 observer sectors, 106 dependency groups and 29 attention entries. The historical file comparison remains pinned to `29f9716`.

- Selected parent Issues and key PRs receive **body-only** result cards, preserving their prior text and lifecycle. The receipt records targets and hashes.
- One delivery-only serialization/ID/reference/rank/new-source-path check and whitespace check; no scientific test suite or parent analysis rerun for this context recovery.
- No new Monte Carlo, server/tunnel operation, Issue lock/closure/rename, PR merge, rebase or force-push.
- Work stays in the isolated PR worktree; the conflicted main worktree is untouched.

The previous handoff below is historical, not a second current priority queue.

<details>
<summary>Earlier handoff — preserved historical text</summary>

## Draft scope

Scientific information recovery plus actual production-data analysis. **Keep Draft PR #267 open and unmerged.** Priority allocates attention only: no Issue locks, closures, demotions, bulk renaming or sequential approval system.

Navigation head: `087b07cb69f2481cbcd55fe2194150d7620835e5`. Incorporated main snapshot: `05af691cc7d438605c138a2f237a0523fec8ae11`, by ordinary merge. No rebase or force-push.

## Latest priority correction: results before more preparation

The overview said production-first while several machine-readable nodes still requested completed work. This update fixes both the pages and their underlying ledger:

- The first real #370 production-confidence adapter, E_top risk/hazard decomposition, P267 factorial, P250 adaptive pilot and P333/P398 coupled closure are no longer missing milestones.
- Saved P334 exact b2 and reconstructible successor H2 are no longer a proposed new acquisition. Given successor H2, the one-step clone expectation is an identity, not another mechanism experiment.
- Generic fixtures, compilers, validators, repeated test suites and full-repository lifecycle sweeps are on-demand support. Each work block starts with a named scientific result, exact object or numerical prediction; affected provenance is recorded at handoff.
- The ledger distinguishes four **runnable-now outputs**, future measurement options, on-demand support and completed reference entries. A ranked reference catalog is not a run queue.

### Current attention

1. **E_top/model survivors on real norm-4/norm-5 and P267 factorial covariance.** Compare amplitude transport with response-direction change, report survivors with uncertainty, and give numerical predictions for a separating third tau/map. Use #370 for a concrete bounded survivor question, not another synthetic ladder.
2. **P334 deterministic checkpoint replay.** Explain safe-insertion degree/2-star structure and its geometric organization on saved microscopic configurations, without new random paths.
3. **Matching/source coupling.** Obtain a same-stream M/K_A loading in a matching-sensitive geometry or a concrete line-relative response. Completed P439 direct/plateau scoring is not a new task; more K_A-only scale points do not answer the missing overlap.
4. **P398 physical response.** Consume the completed nine-mark rooted/landing accumulator in an explicit positive-weight finite-width connected two-point calculation. Report the charge-neutral 2x2 matrix, cross term and determinant.

P250 complete-spectrum/common-k36, a distinct microscopic Phase-E singlet, complex rho-C3, ordered AU/UA and concrete W5 periodic comparison remain parallel options. This update launches no new Monte Carlo production.

## Two newly completed zero-new-sample analyses

| Output | Result | Scientific consequence |
|---|---|---|
| [P334 directional allocation](https://github.com/LightChainr/Matching-One/blob/2b0844a/results/p334-fork-directional-allocation/REPORT.md), Draft `2b0844a` | Projects the parent full 22x22 covariance into total/common-gate/between/within directional responses. Within-checkpoint contrasts are `(3.61607 +/- 3.25069)e-6` at N325 and `(-1.53489 +/- 2.08635)e-6` at N425. | Exact scalar-state nonclosure is established in source `6147e22`, but this additional structure's direction response is unresolved. No unstable component/total fractions or pure-H4 identification. |
| [P439 direct/plateau transport](https://github.com/LightChainr/Matching-One/blob/8498d62/results/p439-direct-plateau-transport/REPORT.md), Draft `8498d62` | Direct max abs(z)=1.11, plateau=1.28; joint zero `chi2=4.69005/8, p=.79013`. Both loading intervals include zero. | No resolved pair of large cancelling terms was exposed. This is unresolved overlap, not proof of zero, shared dynamics or two independent sources. |

P334 inputs remain **branch_only** at `6147e22`; P439 parent `bfbceb2` is **open_pr #451**, not branch-only. Both new derived outputs are **open_pr #267**. Manifests pin inputs and hashes, paired covariance is preserved, and reused streams are not additional independent evidence. The scripts consume immutable Git source objects; those objects must be available when replaying.

## Existing results retained, not rerun

- Canonical E_top Phase D covers ten cross-geometry and six P154 norm-4 blocks. E_top is the Alexander-even topological coordinate, not an identified energy operator.
- P267 `ce01e4d` completes the 100k factorial, resolving tau-by-map interaction at `chi2=236.756/4, p=4.63e-50`; E_top alone remains unresolved.
- The completed radius-1 Phase-E mixed row `e526b9b` has unresolved/sign-changing J_bulk; A/E/J_bulk improves on A/E/C by only `1.981<4`. That row is not selected; other microscopic singlets remain open.
- #370 source `f5779b9` excludes three fixed lines on eight real production rows inside its declared Gaussian outer set, while a free ray survives. These fixed-line results do not automatically certify arbitrary physical M1/M2d/M2j classes.
- P250 Draft `eb29446` finds a 30.97% held-out state increment beyond its declared residue baseline; complete endpoint periodograms are absent and both residuals retain k=36.
- P429 Draft `2d47d72` finds H2 predictive allocations of 70.4%/97.9%; these are model-dependent, not causal or full sampling intervals. New P334 exact checkpoints now supersede its proposed H2/b2 acquisition.
- P337 risk/hazard and ambient-line pilots are complete; new state does not yet mean identified observer coupling.
- P398 branch `5389200` completes the declared 23-dimensional/nine-mark intersection; second-jet `afc619c` gives X2=0, J2=0 and an automatic higher fixed-radical Gram condition. Physical all-Gram/all-Q transfer and continuum identity remain open.

## Navigation and handoff

[NEXT-TARGETS](https://github.com/LightChainr/Matching-One/blob/087b07c/docs/NEXT-TARGETS.md) is the ranked execution entry point. README, STATUS, RESEARCH-MAP, ROADMAP, HYPOTHESIS-BOARD, the [team handoff](https://github.com/LightChainr/Matching-One/blob/087b07c/notes/production-priority-refresh-after-factorial-20260830.md) and both registries agree with it. Ledger v24 records 234 sources, 95 scientific nodes, 15 observer sectors, 100 dependency groups and 29 decision entries. Older scorer inventory counts remain explicitly historical.

## This update's checks and operational boundary

- Each new reanalysis ran once with its necessary arithmetic/covariance identity check; no parent analysis campaign or scientific test suite was rerun.
- One delivery-only YAML/ID/rank/new-source-path pass (24 new immutable paths) and `git diff --check` passed.
- No Huawei environment or tunnel was contacted, and no new Monte Carlo was started.
- Work remained in the isolated PR worktree; the conflicted main worktree was untouched.
- No Issue state changed, no PR merged, no history rewritten.

</details>


</details>


</details>


</details>


</details>


</details>


</details>


</details>


</details>


</details>

