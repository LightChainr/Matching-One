# 当前成果：从取向信号到微观机制

**整理日期：2026-08-31，22时三份意见后复核。** 两个独立决策实验已完成并执行停线；首次Xi、jump/reweight分解和两套固定源的渐近反号也已完成，不再作为待办。[三份意见的核对与实际推进](../notes/reviews-response-20260831.md)记录本轮增量；[下一步](NEXT-TARGETS.md)只保留尚未回答的问题。本任务新增结果交付Draft #509，不合并。

## 本轮：延续需要连通信息，有限反号不等于可直接采样

| 问题 | 已得到的具体答案 | 尚未解决 |
|---|---|---|
| 0/1孔信息加尾概率能否控制齐次U？ | [严格信息不足见证](../experiments/p337-continuation-feasibility-20260831/THEORY.md)：两套保留真实0/1孔表的摘要补全，有相同完整q曲线、唯一root和正斜率，但U约为+10.10358/−10.07432。连符号都不能从所列约束确定。 | 补全不是原图的物理多孔律；需要利用真实连通规则约束未知层的取向热score矩，不能把该见证扩大为一切延续方法无效。 |
| 能否摆脱小孔概率展开？ | [全epsilon面核](../experiments/p337-face-kernel-20260831/REPORT.md)精确保留端口、位移、平行边及源修正；指定两孔构型的权重交叉比为e^t；固定B整行的全孔密度绕环概率已闭式求和。 | 这是条件连通核，完整B平均、root及U误差仍未求出。已有3dc47674 hypergraph/twist表示不是本次新发现。 |
| 原U的有限耦合反号 | [固定m=64复核](../experiments/p337-finite-law-window-20260831/RESULT.md)：N25原Ustar≈−5.82495e−19、Udrop≈+1.07107e−13，各在自己的共同root。执行[2690f665](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)并行交付了更强的全m>=64证书，符号不计作本次新增机制排除。 | 本次新增的是下行的采样可达性判断；固定N25半直线不等于N50饱和到齐次的延续。 |
| 十台机器能否直接测上述反号？ | 同一计数给rank1概率：轴向约1e−15，斜向约2e−20。95%机会至少见一次rank1的必要样本下界为2.28e14至4.47e19，**不启动该点普通无条件采样**。 | 该界不是估U的充分预算；不约束条件、importance或twist partition估计器，也不是墙钟时间预测。 |

三项均无新MC、孔型全枚举或云作业。m64比较先冻结`375a6f0c`再评分，原四点实验保持不变；其余是明确数学构造和有界配置核验。它们不提供新的独立统计证据。

最新既有结果另已接入：[e1b96895](https://github.com/LightChainr/Matching-One/blob/e1b968959634b9b3999c727b83ed38d0b730cb20/results/defect-reweight/REPORT.md)给Xi_reweight=+4.550327123237、Xi_jump=−15.306045530801，排除jump-only；[f4057192](https://github.com/LightChainr/Matching-One/blob/f405719264c896aa873dd4aae7292795f544ba99/notes/topological-projection-reverses-global-u-tail.md)已给Sstar/Sdrop相反渐近尾。N100/N225仍是理论预测，不是已执行生产。

更新的[2690f665联合极限](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md)给N/m²→zeta<∞时的Poisson/full共存、原pooled分母控制和U超多项式衰减；它不是固定m的大N定理。该分支剩余的是固定m的oblique twist代价和受限扇区odds，不再等待首次有限反号或指定联合极限的分母界。

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

支持线#275/#419/#370/#398及#1保留已有产物；#275已由P1调整为P2，保持开放。只有能说明“哪种结果停止哪个候选”的具体任务才进入P0；跨N旧数据回归、进一步projection/Hessian/descriptor及generic certificate目前均属exploratory/support。下一项实际分析从[NEXT中的精确U传递接口](NEXT-TARGETS.md)出发；本轮验证块不再用于选择新模型。

## 取向与物理响应

| 已有结果 | 证据与适用范围 | 结果位置 |
|---|---|---|
| **方位H4已有独立证据** | P43+P57独立primary合成对零为31.1857/4、p=2.81e−6；固定H4为3.4623/4、p=.484。norm-5区分所测试H8/H12；prism新增独立支持。norm-5子块本身对零仍相容。 | main：[综合报告](../notes/issue212-matching-odd-synthesis.md)；[prism原讨论](https://github.com/LightChainr/Matching-One/issues/205#issuecomment-5462845639) |
| **单一标量尺度修正已不足** | N145→290全曲线单倍率9.3520/2、p=.0093；S′纯幂律与单一rank-gap宽度均失败。norm-4 q2为20.897/2、p=2.90e−5；Jordan标量p=.067、全jet p=.054。第四代中Jordan加一个even-mode家族存活，但λ=0、1/2、1难区分，未选出次级模身份。 | [已完成的开放PR273](https://github.com/LightChainr/Matching-One/pull/273)；[四代结果PR277](https://github.com/LightChainr/Matching-One/pull/277) |
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
