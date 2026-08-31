# 当前成果：从取向信号到微观机制

**整理日期：2026-08-31。** 17:20复核已更新P154全链/两阶段/百万端点和P334九层/共享标签/SS分解；[最新核验](../notes/progress-and-compute-20260831.md)注明本地提交与GitHub同步差异。[下一步](NEXT-TARGETS.md)给出团队分工，[清理记录](REPOSITORY-TRIAGE-20260831.md)记录实际Issue操作。科学完成度与main/open PR/独立分支位置分别标明，之后变化不自动包含。

## 取向与物理响应

| 已有结果 | 证据与适用范围 | 结果位置 |
|---|---|---|
| **方位H4已有独立证据** | P43+P57独立primary合成对零为31.1857/4、p=2.81e−6；固定H4为3.4623/4、p=.484。norm-5区分所测试H8/H12；prism新增独立支持。norm-5子块本身对零仍相容。 | main：[综合报告](../notes/issue212-matching-odd-synthesis.md)；[prism原讨论](https://github.com/LightChainr/Matching-One/issues/205#issuecomment-5462845639) |
| **单一标量尺度修正已不足** | N145→290全曲线单倍率9.3520/2、p=.0093；S′纯幂律与单一rank-gap宽度均失败。norm-4 q2为20.897/2、p=2.90e−5；Jordan标量p=.067、全jet p=.054。第四代中Jordan加一个even-mode家族存活，但λ=0、1/2、1难区分，未选出次级模身份。 | [已完成的开放PR273](https://github.com/LightChainr/Matching-One/pull/273)；[四代结果PR277](https://github.com/LightChainr/Matching-One/pull/277) |
| **微观簇源的偶响应已明确测得** | S=(黑NN簇数+白matching簇数)/N。P40百万N65/N85已补齐E×S及E×controls；q-fugacity补偿下C四几何均为正（110–138SE），共同raw源的H4联合p=.21156仍未分辨。full辅助p=.04911是相关视图，不能认定能量场。此C保持matching均值，但未做真实Bernoulli温度补偿，也不是原norm-4的根/斜率归一化U。 | Draft267：[最新百万偶响应](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/results/p40-even-given-odd/REPORT.md)、[原norm-4接口](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/notes/p40-even-response-norm4-interface.md)；[20k来源](https://github.com/LightChainr/Matching-One/blob/eb7ef8c9f13a88d96f32c2da62ba7ef2145cb33e/results/p154-absolute-cluster/REPORT.md) |
| **正的源/读出证据与M载荷分开** | 外部Euler源在fixed-K分解后仍有稳定约32%余项；F5已分辨两条源读出。原始matching M的同流载荷零仍相容（p=.58155），连接尚未分辨。 | [外部源](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5468416605)、[开放PR451](https://github.com/LightChainr/Matching-One/pull/451) |
| **primitive square-bond另有多character结果** | N112独立生产支持r0+r1；E_top同流方向行列式p=2.49e−5，说明两个拓扑观察方向。纯E4/E6/E4²被拒；这套square-bond观察量与square-site thermal主线分开保存。 | [N112结果](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5469692921) |

P40实际引擎在随机键中包含N：N65/N85是不同N-domain，通常PRNG独立假设下可作nominal联合统计；同N两方向共享随机数。早期仅凭seed标签判断跨N共流的说法已经更正。

旧百万q/source报告a4cbf02缺E×S，56a6267已补齐；29f339又完成六N的原U源响应。后续两阶段估计没有稳定降噪，N260/N340首1M标记与评分也已在本地b77f3e2a完成。随根移动的rank-1源响应约−.093至−.098，排除有限尺寸纯共同温度换元；原U源q2/Jordan仍未分辨（p=.48573/.58735）。共同E加时钟尚为候选，下一观测应进入rank-1内部方向，不能把已算ΔE改写成新测试。[最新来源与同步状态](../notes/progress-and-compute-20260831.md)。

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

**上述R1之后的全量问题也已推进。** C/L、direct/collective及45坐标crosswalk已完成；四状态解释99.84895%/99.86380%的R1门控积分个体方差，不能当成全局A解释率。随后完整40k路径的A/E九层分析完成，128万条辅助条件续接测到正的跨birth Gamma。30c7ddb0又完成SS/mixed/BB分解：mixed占canonical Gamma约78%，SS约7%–8%，三项按原20批精确重构总量。959a7fa2已补齐接触坐标，接下来连接这些坐标与机制贡献；不再派首次九层、首次Gamma或首次SS分解。[数值、固定提交与解释边界](../notes/progress-and-compute-20260831.md)。全部rank的精确条件空间响应仍未齐全。

147-prefix噪声预算仅适用于固定经验mixture。P398的93维Krylov空间和Boolean谱degree按各自生成过程解释；现有width4/8的i^j权重保持波长4，并非固定模式序号的尺度外推。

**P398固定速率干预在本轮完成。** η=0,±1/4完整1430态/186维扇区运行1.709秒；旧U在t>0出现cross传播，不能用共同时间换元解释。固定14→16配置函数使最大C0尺度误差改善4.36倍，但弱minus长尾相对误差仍−12.64%。投影借用各η的完整平稳分布，既不是仅用旧η=0矩阵的盲预测，也不证明闭合或square-site映射。此为有限对比，精确η0导数和零频积分未计算。[实验、代码与结果](../experiments/p398-rate-intervention-20260831/README.md)。

## 后续纠错已经改变的解释

| 旧说法或候选 | 现在应采用的结论 |
|---|---|
| P418巨大共同谱惩罚支持radius flow | batch sum统一为per-sample后四共同谱均相容；radius5单壳数值不可靠。正确归一化的P250秩与exact CRT不受此错影响。[修正报告](https://github.com/LightChainr/Matching-One/blob/e2b57aa7c5ec5c7db8cbb4f03872435f20966407/results/p418-normalized-archive/REPORT.md) |
| #43 even通道原冻结检验直接通过 | 原cross/either错配是协议失败；精确符号运输后的.57003/2属于事后确定性纠正，历史仍保留。[纠错](https://github.com/LightChainr/Matching-One/issues/133) |
| 非零q/标记耦合即可识别新场 | q仅有三值；q-only单根接触耦合已有全阶闭合。真正独立物理识别需要更丰富的微观读出。[no-go](https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5463228151) |
| 一个低阶Hankel秩或形式Jordan足够 | rank随生成元与观察窗改变；旧R2幸存者被更高阶/联合map分析否定。Q=1端点也不唯一指定generic-Q切向。[完整链](https://github.com/LightChainr/Matching-One/issues/250)、[Q-lift](https://github.com/LightChainr/Matching-One/issues/333) |

## 支持资产与统一解释边界

阈值rank引擎、协方差、exact certificate、Q4表示、W5周期配对及15态terminal serial代数均已可用。#498–502新增的是有限子群/作用/理想/闭合集分类。W5已实现的两图并非自对偶，当前没有新的严格p_c界；这些资产接具体概率比较时才成为阈值推进。

目前尚未给出square-site p_c闭式、唯一连续场/Jordan身份、总体方向响应的几何归因。相容性不等于模型确认；有限图证明、有限尺度数据、表示论选择规则和晶格到连续场重叠各按其实际范围使用。重分析共用原块，不增加独立证据。

旧状态全文保存在[整理前固定commit](https://github.com/LightChainr/Matching-One/blob/8a68cca866d7fbca7463e2167c3ff06128d5851f/docs/STATUS.md)。本次没有改动旧结果、冻结文件或其他分支；当前分析顺序见[三个问题](NEXT-TARGETS.md)。
