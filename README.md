# Matching One

研究正方晶格点渗流：有限系统为何出现稳定的四重取向修正，这种修正由什么微观结构产生，以及它怎样影响阈值逼近。

精确起点是 `p_c(square site) + p_c(NN+NNN site) = 1`。现有研究已经建立同调解释、独立H4数值证据和具名几何传播机制；正方点渗流阈值闭式与连续场身份仍开放。

## 直接开始

| 想知道什么 | 入口 |
|---|---|
| 已经做成什么，哪些解释已被后续结果改变 | [当前成果](docs/STATUS.md) |
| 下一项分析具体算什么、用哪份现有数据 | [下一步分析](docs/NEXT-TARGETS.md) |

三次穷举的范围与优先级问题见[专项审查](notes/cubic-search-review-20260831.md)。

README负责入口，STATUS负责成果，NEXT-TARGETS负责下一步分析。直接进入对应问题即可。

**两个独立决策实验都已完成并触发降级。** #154的165M新路径将该lag=1源的净U响应限制在预定弱效应范围，停止其“主要H4解释”优先投入；#334的600k新prefix同时排除两个固定残余预测。一般问题保留为P1，不自动续跑；没有追加样本或事后救场模型。当前科学判断以[STATUS](docs/STATUS.md)为唯一入口，[NEXT-TARGETS](docs/NEXT-TARGETS.md)只列执行顺序。[清理记录](docs/REPOSITORY-TRIAGE-20260831.md)保存原13项结题及优先级操作；所有已有数据、分支和PR保留。

**随后有限机制又得到一次明确排除。** checkerboard闭合源与一孔Xi都已在其他分支完成；本次直接使用其精确系数，证明离开端点时，两几何的原q/E响应不能由共同温度和同一源耦合吸收（D3严格大于1/10000）。该有限候选退休，无新枚举或云任务；不外推连续场身份。[新结果](experiments/p337-two-coupling-closure-20260831/RESULT.md)

**已补齐完整齐次N50传递。** 两几何各精确覆盖2^50配置，合计约49.85 CPU秒、峰1.63 GiB；原U=1.06156039、固定源响应V=+0.05434578，严格排除该有限点零传递。正号预测存活，不等于机制确认，不追加尺寸/源参数。[本轮结果与复核](experiments/p337-homogeneous-n50-20260831/RESULT.md)；旧m64直接采样停线及P154/P334降级保持。

**#537 的 one-defect 判决已经触发。** 在完整总体先冻结 `p/R/mu/beta` 的 N25 physical graph 上，一条真实 `z` flip 同时使 rank `0→1`、source Bell `9240712→6848576`、`g16 4→0`；恢复 geometry pool 与 C4 orbit 后，Eq.(10) 的 source midpoint 为 `-1.0888815582478189e-11`，完整 Schur 权为 `-8.298623728474635e-12`，均有严格不含零的有理区间。因此 full graph 上“两个慢变量必须由两个独立 defects 承载”以及自动领取 six-arm gain 的路线已经否定。[精确证书与判决](notes/p537-one-defect-diagonal-edge-20260901.md)同时冻结边界：该 witness 属于 contact/collision、`arm_mask=3`，不代表 ordinary separated sector，也不提供渐近下界。#537 继续作为唯一 P0，但对象已收缩为 contact/collision 对 surviving leading four-arm signed functional 的贡献；按固定 `d_NN≤1` / `d_NN≥2` 分解处理，不启动 L6、MC、距离网格或服务器生产。

**协方差零空间QA已经闭合。** [#543回溯](experiments/p543-covariance-nullspace-audit-20260901/REPORT.md)集中修复三个广义卡方实现，并仅用既存 residual/covariance 重算16个归档向量分数。15项不变；P50的默认`9.35200/2`也不变，但其丢弃方向与残差不相容，必须标为cutoff-sensitive，不能再称无害numerical null。该修复不增加证据票、样本或云任务，support线到此停止。

## 已有成果的主线

**同调结构。** 有限matching量为 `M=P₂−P₀=E[r_black]−1`；两次阈值rank分别是两个本质同调方向的出生。一般有限商的有理同调与整数饱和证明稿已经交付，完整出生时间与空间响应已有可计算对象。[证明与推论](https://github.com/LightChainr/Matching-One/issues/269)

**独立取向证据。** 同N取向对照、独立随机流和Gaussian lineages支持所测范围内的 `DeltaM ~ DeltaCos4 * N^(-13/8)`。norm-5及prism区分了所测试的H4/H8/H12候选。全曲线与norm-4结果同时要求超出一个标量振幅的结构。[证据汇总](docs/STATUS.md#取向与物理响应)

**真实机制对象。** 已求解147个选定真实prefix的完整birth clocks；空间影响可区分完整均值时钟看不到的几何。正权width8传播已有簇大小、接触和T4的明确路径解释。[过程与传播](docs/STATUS.md#过程与传播)

## 当前已完成，直接使用

- **N900完成**：32M共享counter、800批次，两个冻结宽度预测均存活；不要再按旧导航等待首次结果。[报告](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/REPORT.md)
- **P40百万样本偶响应完成**：缺失混合矩已补齐，固定matching均值的q源补偿下，四几何偶响应均明确为正；共同raw源的H4方向差仍未分辨。[报告](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/results/p40-even-given-odd/REPORT.md)
- **P418共同谱已纠正**：统一每样本单位后四个共同谱相容，旧巨大惩罚不再支持radius flow。[修正](https://github.com/LightChainr/Matching-One/blob/e2b57aa7c5ec5c7db8cbb4f03872435f20966407/results/p418-normalized-archive/REPORT.md)
- **P154角权桥已完成**：全链、百万端点、条件line/fixed-K之后，六N原U及源导数的软角分配也已交付。固定K/rank1内中心化源对原U严格为零；该约束已用于本次独立传递检验；不能继续把更强O4检测当成全局机制。[报告](experiments/p154-spatial-localization-20260831/REPORT.md)
- **P334有限空间源已完成**：九层、Gamma、接触及共同源全曲线之后，t=±1的既定q_t仍改变未来响应，同时保持即时Euler/rank联合分布。新完整census与定向续接进一步支持同prefix局部A二维成分；结果见[新交付](notes/analysis-delivery-20260831.md)。[有限源报告](experiments/p334-finite-source-20260831/REPORT.md)
- **P398精确响应完成**：继固定η干预后，η=0一阶反号时刻和零频积分已计算。平稳重加权与动态贡献竞争，16维模型积分误差约0.5%，最低两极不足；本有限模型保留P2。[实验报告](experiments/p398-linear-response-20260831/README.md)

- **P334开始给出微观预测**：固定源Gram模型对组外A响应的误差降低38–39%，出生中心响应降低54–59%；只用原数据、每模型四个参数。新增条件形状检验也支持固定prefix内出生中心方差发生变化，单纯确定性平移不足。[预测](experiments/p334-prefix-prediction-20260831/REPORT.md)、[形状](experiments/p334-conditional-shape-20260831/README.md)
- **P154早期结构影响后续rank已测得**：一个预定lag、原2.4M排列补观测给出明确的进入/退出响应，并分出早rank0/1对后期rank1人口的竞争；原H4源导数仍弱。这是新的路径源结果，不替代旧同层分解。[时序报告](experiments/p154-temporal-source-20260831/REPORT.md)

#1代数路线保留P2候选验证能力，本任务暂停新增次数/高度/区间穷举。#520–#523及随后四次区间排除由其他执行者继续交付，21时复核时main已到Yang–Zhou四次结果`580678a7`；本任务未执行这些合并。有限区间排除不等于真实阈值的无条件代数排除，亦不提供超越性证据。

## 当前研究方式

用户在2026-08-31要求结束持续扩张的探索阶段。一个新分析若不能在结果出来前说明“什么结果使哪条机制线停止”，默认列为exploratory/support，不列P0。暂不增加prefix descriptor、描述性observable、通用证书或无理论约束的全族穷举。

两项主实验各使用冻结的少数预测和全新随机块，一次采集区分量。archive可以训练模型和估算预算，不参与独立验证得分；模型、评分、样本上限和失败后的降级动作先提交，再运行。若结果否定候选，先将该候选降级，不用第四个事后模型挽救本轮得分。能力不足以形成判别预测时，明确写“未形成可执行决策实验”。

原样本、共同随机块和准确量纲继续保留。这些是正确计算的基础，也不能替代对自适应假设选择的控制。支持工具仅服务当前的具体候选；机器数量不决定任务数量。

## 背景与沿革

三队地址和分工共用[仓库协调入口](https://github.com/LightChainr/Matching-One/blob/0e9d684e88c26b904da342b4c33cdc04057a3d07/docs/TEAM-COORDINATION.md)，主要通过仓库结果交接，避免重复计算与频繁消息。

- [Research Map](docs/RESEARCH-MAP.md)、[旧Roadmap](docs/ROADMAP.md)保留路线沿革，当前顺序由README、STATUS和NEXT-TARGETS承担。
- [Draft #267](https://github.com/LightChainr/Matching-One/pull/267)保存完整研究交接并仍在更新；本次只提供简洁入口，不改动该分支。阅读具体报告时使用固定commit。
- [8月31日本轮整理前的入口](https://github.com/LightChainr/Matching-One/blob/8a68cca866d7fbca7463e2167c3ff06128d5851f/README.md)保留旧判断的时间顺序。
- 数值阈值来源仍按方法保存于[data/literature_threshold_sources.json](data/literature_threshold_sources.json)，不把参考小数当精确常数。

文档更新在独立Draft分支交付；Issue已按本轮授权结题或调整优先级。**不合并、不删除研究分支**，代码、数据和历史报告原路径保持可用。

MIT，见[LICENSE](LICENSE)。
