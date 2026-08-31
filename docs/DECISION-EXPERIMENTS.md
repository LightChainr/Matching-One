# Two decision experiments / 从档案解释转向新数据预测

2026-08-31 · 当前默认主线只有 **#154 temporal transmission** 和
**#334 independent intervention**。其他研究保留为并行探索或按需支持。
优先级表示注意力；停止规则针对写明的模型预测，不锁任务，不关闭 Issue。

## 当前判断

对“研究过拟合”的担忧有根据：完整协方差能处理读数之间的依赖，却不能
撤销由同一数据反复启发模型、选择特征和挑选观察者的过程。旧 norm4 的
2.4M 排列、P334 的原20批 prefixes 均保留为 **discovery**。新 suffix、
cross-fit 和更精确的条件均值有价值，但不变成新的 prefix 总体。

同时，项目已经开始转向：执行团队的
[独立干预冻结稿，`bc0a18c2`](https://github.com/LightChainr/Matching-One/blob/bc0a18c207e3b09f49ea6b6af6601471114d654a/notes/p334-independent-intervention-freeze.md)
已有固定模型、效应容忍度、独立随机域和生产预算。本页直接承接该冻结稿，
不另建第三套 #334 计划。后续
[dispatch，`bde1a51ca95c74448265b670ba0d9a0d87915479:notes/p334-independent-intervention-dispatch.json`](https://github.com/LightChainr/Matching-One/blob/bde1a51ca95c74448265b670ba0d9a0d87915479/notes/p334-independent-intervention-dispatch.json)
已记录封存并集成的 producer；状态仅为 **frozen before formal generation**，
不表示正式生成已经完成，也没有该独立块的结果。

| 主线 | 已完成、只作探索输入 | 当前需要的下一项交付 |
|---|---|---|
| #154 | [两条固定时钟律](../results/p154-fixed-clock-models/REPORT.md)，`3847a5cf`：纯 cos4 相对位移不符旧方向读数；标量位移暂存但 U 预测仅约.01/.018 | 给出可负担且能区分原 U 机制的预报；不要为已不符的角律追加生产，或声称8M可分辨微小标量预报 |
| #334 | hierarchy、contact loading、local A-rank、source-normal、旧prefix held-out/shape均已交付 | 按 `bde1a51c` 已封存producer执行原 `bc0a18c2` 干预；发布新100万 prefixes 的固定主统计量，不重做设计或旧档案跟进 |

## 1. #154：显式传输律已有一项淘汰，独立 U 比较仍需可分辨预报

**本轮直接分析，`open_pr #267`，`3847a5cf`：**
[固定时钟比较](../results/p154-fixed-clock-models/REPORT.md) 只比较
`e_g=−(m−delta w_g/2)F1'_g`、`x_g=−(m+delta w_g/2)F2'_g`，
其中 `w=1` 或精确 `cos4theta`，位移在根附近不随 p 变。
每 N 只用 pooled rank1 与 root 响应校准 m、delta；不拟合方向增益。
纯 cos4 律对 entry/completion 方向差的偏离为31–35个配对 SE，
其 U 预报 +2.777/−4.081 也与观测有张力。停止推广这个具体角律；
不等于排除所有 H4 机制，更不关闭 #154。

标量律方向差残差为0.19–1.21 SE，但 U 预报仅
`+.00980±.02328`、`+.01753±.02817`（N260/N340）。
[完整配对预算](../results/p154-clock-transmission-budget/REPORT.md) 已补齐：
8M新排列的边际3SE分辨率约.766/1.254；只用点预报外推就需约
489亿/410亿样本，固定校准不确定性更使两预报的3SE区间覆盖零。
这撤回了“有公式就可投产”的推断。保留标量律作为基线，下一独立
生产要有可负担的机制差异；不追加新角系数来挽救已失败模型。
以上都是同一旧块上的回顾性模型比较，非独立确认或 exact certificate。

执行团队的[单系数推导](https://github.com/LightChainr/Matching-One/blob/bde1a51ca95c74448265b670ba0d9a0d87915479/notes/p154-prospective-birth-clock-transmission-decision.md)
已完成，不能再写“没有任何数值预测”。#509 的
[较长 lag](https://github.com/LightChainr/Matching-One/blob/04743caf1450d5f88cae2747e0dbee36d7cd8ca1/experiments/p154-temporal-source-20260831/REPORT.md)
也已完成：`K−ceil(sqrt(N))` 早期源给出强 entry/exit 响应而 U 未分辨。
它不是下面的 lag=1 或平衡源，也不是新的独立排列。当前不扫描更多 lag。

### 保留的通道判别框架

固定原 square-site 方向对、lag=1、bulk 源单位和原 U；建议第一比较仅用
N130→260，作为有限规模比较，不拟合自由指数，也不把一对规模叫作完整
norm4 三代残差检验。若目标必须是 q2/Jordan 三代消除，需另给该源的三代
预测，不能继承原平衡源的缩放律。

在每个最终 K 的前一步 j=K−1，令
`epsilon=s−E[s|j,rank,g]`。将已有源分为两个可独立开启的物理开关：
`epsilon0=1_(rank=0)epsilon`、`epsilon1=1_(rank=1)epsilon`。
两者同强度相加恢复已测源的全部 q/E/U 一阶响应：早期 rank2 已吸收，
其居中源对这些未来拓扑读出为零；这不是对所有空间读出的源等价。
保持各早期层质量，随后均匀添加一个点。
使用已有三个互斥事件，不增加新的描述量：

```text
entry-source:      X0=T01+2T02, Y0=−T01
completion-source: X1=T12,     Y1= T12
Jq_i = sum_K Bin(N,p)[K] X_i(K)
JE_i = sum_K Bin(N,p)[K] Y_i(K)
```

直接02属于早期 rank0 开关；其 E 增量为零，但通过 q 改变根与归一化。
这与报告中把02分给两次 activation 的记账不同，不能混用。第二开关称
**completion**；一步协议没有测定完整 lifetime 机制。

设 bar 为两方向均值，P4 为原精确角差归一化，所有量在未扰动公共根计算。
定义 `D=bar(q')`、`B=P4[E']`、`H=P4[E'']`、`T=bar(q'')`，则

```text
pdot_i = −bar(Jq_i)/D
V_i = N^(13/8)/2 * [(P4[JE_i'] + pdot_i H)/D
                   − B*(bar(Jq_i') + pdot_i T)/D²]
V = dU/dsource = V0+V1
Rdot_i = −bar(JE_i) − pdot_i bar(E')
```

这给出了 `source → entry/completion → original U` 的显式转换器。
机制必须给出 signed T01/T12 或相应约束，再通过这个 map 预测 V。
首次 activation 对 pooled-root 贡献更大，**不推出它主导 U**；方向共同
的事件核也不自动令 V=0，因为根和分母仍会改变。

### 三个候选判别区，而非伪装成已建立的三个物理模型

事前固定有意义的 `delta_U>0`、population 最小效应 `delta_R>0`、主导
裕量 `eta>0` 与各 N 的预测符号 s_N。一个联合置信域用于两个 N 及全部
固定通道，不能逐坐标挑显著值。

- A：`|V|≤delta_U` 且 `Rdot<−delta_R`。称 population-only **总读出**，
  允许通道抵消，不能称为“完全没有方向传输”。
- B：`s_N V>delta_U` 且 `s_N V0−|V1|>eta`，两个 N 均成立。
- C：`s_N V>delta_U` 且 `s_N V1−|V0|>eta`，两个 N 均成立。

A/B/C互斥，但不穷尽所有机制。只有置信域完整进入相应区域才判定；
混合、跨尺度不一致、全模型失败或精度不足均保留原名，不补第四个模型
解释同一 confirmation 结果。零附近的不显著不能替代等价性。

**未完成项写明：** 这个 A/B/C 框架的 s_N、科学容忍度与可分辨数值
向量尚未冻结；它与上面已算过的两条具体时钟律不同。不能把候选区
写成已投产或已识别机制。下一交付是可负担的具体预报，
不是另一个 feature/lag 扫描或工具框架。

预算示例只说明成本：若选 `delta_U=0.5`，单 N、真值零、正态近似下，要
有约80%概率让95%区间整体进入该带，需 SE约0.154。由旧 SE按
`M_new≈M_old*(SE_old/SE_target)^2` 外推，N130约2.6M、N260约22M新排列。
这不是已冻结样本量；两尺度及通道联合判别需按完整协方差重新定预算。
新数据重新估计条件均值、根和误差，旧批次不并入 confirmation 分数。

## 2. #334：直接执行已有冻结的 source-normal 干预

`bc0a18c2` 已将原“继续解释 residual”收敛为更强、更窄的两项预测：

- **M0：完整两-score 标签闭合**，
  `m_C(Z,u)=c_a(Z)+b_f(Z)s_f(Z,u)+b_s(Z)s_s(Z,u)`，
  同一 prefix 的 safe classes 共用斜率。它预测 source-normal center
  响应严格零。这比第一 Jacobian 的 `J=B G` 强，不能用其失败否定后者。
- **M1：source-normal response 可迁移到新 prefixes**。四个 own-source
  center 预测为 `(4.116,3.233,3.300,3.977)*1e−8`；固定等权主量
  `T_forecast=3.6565e−8`。这是由旧数据得到的有限规模预报，不是第一性
  原理振幅或 continuum field 预测。

在全新 prefix 上，只由完整标签 census 构造

```text
G=E[ssᵀ], T_oo=s_o²−pi_a² Var_a(L_o)
phi_oo=T_oo−E[T_oo s] G⁺ s,    B(Z)=max_u |phi_oo(u)|
q±(u|Z)=[1±phi_oo(u)/B(Z)]/d
```

G⁺按确切源空间秩投影；奇异 prefix 不丢弃，B=0贡献零。
这些是合法的**非负**有限概率律，某些标签概率可以为零，不能误称严格正。
它们保持各 joint-safe class 质量、即时 rank/Euler 分布和两原 source-score
均值，改变第一 score 空间外的标签结构。各臂后接普通均匀 suffix。
`B(Z)*(F_plus−F_minus)/2` 精确无偏于 `E_uniform[phi F|Z]`，因此不是
小参数 Taylor 外推；也不是未加权的两总体均值差。

冻结生产：N325/N425，k0=193/252，每 N **20×25000 fresh prefixes**；
每个 active prefix/own-source8对± draws，non00保留为零贡献而不改分母。
新 prefix/label/suffix 随机域按 N、batch、prefix、source、rep 分离，
base seed `202608311920334`，旧数据不参与新块生成或加权。
dispatch 已明确 NePnUn（`1e313eebd1a947b8b38714aeea1404d5`）负责N325，
551oUR（`9b3eb563ddc64f2680be17501ec4fe6a`）负责N425，各14 workers；
记录的每机实际限额为14.5 CPU cores、25GiB，目录为
`/workspace/p334-independent-normal-20260831`。这些是该提交的调度记录，不是服务器锁或实时探测。

producer 提交为 `513552c77f035526efb99075b54032d288b2f4bb`，已由
`6928b3d861d2f4ce1ee93446c02e6a44e56832a6`集成进执行分支；
上述 `bde1a51c:notes/p334-independent-intervention-dispatch.json`封存其引用。
精确整数、rank-aware源投影不筛掉奇异prefix；有限两臂采用精确有理共同逆CDF耦合。
预算、seed、四个预测、主T及delta均沿用原冻结稿。**producer准备已完成；
dispatch状态仍仅是正式生成之前，不是已完成实验或 `main` 集成声明。**

固定判读为新块主量 T±3SE、`delta=1e−8`：下端>delta，停止推广 M0 的
完整闭合；上端<delta，停止推广具有实质正效应的 M1 预报；其余记未分辨，
固定预算结束。另报告相对 `T_forecast` 的差异，不重新拟合振幅。
3SE是声明的诊断区间，不是 anytime/exact certificate。
**M1失败不等于M0被证实**；例如清楚的负响应会同时反驳两项预报。

用户提到的约20%确有对应：最新
[original00 卡片](https://github.com/LightChainr/Matching-One/blob/bc0a18c207e3b09f49ea6b6af6601471114d654a/notes/p334-prefix-response-projection-scientific-card.md)
给四-contact shares约78.20–80.36%，使用旧8 clock与新64 response。
更早 receiver-R0汇总的80–99%是另一总体，不能覆盖此结果。
两者都是 signed loading，均不等于未解释方差或独立机制比例。

### 旧prefix预测与shape已经交付，仍不能替新normal实验判决

`open_pr #509`固定 `04743caf1450d5f88cae2747e0dbee36d7cd8ca1` 的
[`experiments/p334-prefix-prediction-20260831/REPORT.md`](https://github.com/LightChainr/Matching-One/blob/04743caf1450d5f88cae2747e0dbee36d7cd8ca1/experiments/p334-prefix-prediction-20260831/REPORT.md)
给出既定5折、原批次留出的 BG 相对常数模型风险改善：A为38.17%/39.37%，
出生中心为54.42%/58.65%（N325/N425）。两个模型每输出均4参数；
这是原cell00的1502/1551个prefix上的预测增益，两个N分别训练，尚无跨N零重拟合结论。

同提交的
[`experiments/p334-conditional-shape-20260831/results/REPORT.md`](https://github.com/LightChainr/Matching-One/blob/04743caf1450d5f88cae2747e0dbee36d7cd8ca1/experiments/p334-conditional-shape-20260831/results/REPORT.md)
也已交付：combined72的 `VarC.D.minus` 为
`−2.43469e−8±5.32765e−9` / `−1.27127e−8±3.30315e−9`，
保留原20000-prefix分母和20批误差；总体shape energy仍弱且未截零。
old8/new64/combined72共享原prefix，局部均值det只是追加到同factor，没有重跑。

这两份完成结果可以更新 discovery 机制解释，不应再登记为待交付。
BG只约束原两-score的第一Jacobian近似；上述shape也不是新 source-normal 主T。
它们既不直接裁决更强的完整标签闭合M0，也不把新M1预报变成独立总体确认。
固定独立块的两模型、主量和判读规则保持不变。

## 冻结、收口与仓库职责

冻结包含 source、readout、模型预测、统计单位、随机域、预算和失败/未分辨
规则，随后固定 producer SHA，再生成 confirmation 数据。探索自由保留；
结果出现后修改的模型只进入新的探索/后续预测版本，不回写本次成绩。
不要求等完整理论才做有限机制比较，也不把精确证书当作投产前置条件。

执行团队承接已有 #334生产；总览维护这两个问题与结果入口；其他团队的
旧档案交付按发现材料吸收。沟通留在仓库，不等待团队许可或逐条确认。
更新后的 Huawei Skill 已读；此前十机列表均Ready，当前未重新探测。
本任务未开机、未运行新 production，也未触碰其他作业。PR #267保持Draft。
