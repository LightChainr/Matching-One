# Matching One

研究正方晶格点渗流：有限系统为何出现稳定的四重取向修正，这种修正由什么微观结构产生，以及它怎样影响阈值逼近。

精确起点是 `p_c(square site) + p_c(NN+NNN site) = 1`。现有研究已经建立同调解释、独立H4数值证据和具名几何传播机制；正方点渗流阈值闭式与连续场身份仍开放。

## 直接开始

| 想知道什么 | 入口 |
|---|---|
| 已经做成什么，哪些解释已被后续结果改变 | [当前成果](docs/STATUS.md) |
| 下一项分析具体算什么、用哪份现有数据 | [三个分析问题](docs/NEXT-TARGETS.md) |

README负责入口，STATUS负责成果，NEXT-TARGETS负责下一步分析。直接进入对应问题即可。

**团队当前只保留两个P0：[#154微观簇源的非热响应](https://github.com/LightChainr/Matching-One/issues/154)、[#334共同位置与birth几何](https://github.com/LightChainr/Matching-One/issues/334)。** #398本轮固定干预已完成，降为P2；其余P1为储备/主线接口。已关闭13个完成或重复Issue，原24个独立成果PR保留未合并；[清理记录与优先级](docs/REPOSITORY-TRIAGE-20260831.md)。[新增分析交付与五机安排](notes/analysis-delivery-20260831.md)区分已同步成果、本地新增结果和下一步。

## 已有成果的主线

**同调结构。** 有限matching量为 `M=P₂−P₀=E[r_black]−1`；两次阈值rank分别是两个本质同调方向的出生。一般有限商的有理同调与整数饱和证明稿已经交付，完整出生时间与空间响应已有可计算对象。[证明与推论](https://github.com/LightChainr/Matching-One/issues/269)

**独立取向证据。** 同N取向对照、独立随机流和Gaussian lineages支持所测范围内的 `DeltaM ~ DeltaCos4 * N^(-13/8)`。norm-5及prism区分了所测试的H4/H8/H12候选。全曲线与norm-4结果同时要求超出一个标量振幅的结构。[证据汇总](docs/STATUS.md#取向与物理响应)

**真实机制对象。** 已求解147个选定真实prefix的完整birth clocks；空间影响可区分完整均值时钟看不到的几何。正权width8传播已有簇大小、接触和T4的明确路径解释。[过程与传播](docs/STATUS.md#过程与传播)

## 当前已完成，直接使用

- **N900完成**：32M共享counter、800批次，两个冻结宽度预测均存活；不要再按旧导航等待首次结果。[报告](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/REPORT.md)
- **P40百万样本偶响应完成**：缺失混合矩已补齐，固定matching均值的q源补偿下，四几何偶响应均明确为正；共同raw源的H4方向差仍未分辨。[报告](https://github.com/LightChainr/Matching-One/blob/56a6267d6a6826a165f93ed3a64a670ca7088180/results/p40-even-given-odd/REPORT.md)
- **P418共同谱已纠正**：统一每样本单位后四个共同谱相容，旧巨大惩罚不再支持radius flow。[修正](https://github.com/LightChainr/Matching-One/blob/e2b57aa7c5ec5c7db8cbb4f03872435f20966407/results/p418-normalized-archive/REPORT.md)
- **P154角权桥已完成**：全链、百万端点、条件line/fixed-K之后，六N原U及源导数的软角分配也已交付。固定K/rank1内中心化源对原U严格为零；下一步必须解释rank人口/进入退出，不能继续把更强O4检测当成全局机制。[报告](experiments/p154-spatial-localization-20260831/REPORT.md)
- **P334有限空间源已完成**：九层、Gamma、接触及共同源全曲线之后，t=±1的既定q_t仍改变未来响应，同时保持即时Euler/rank联合分布。新完整census使同prefix局部检验可有效使用旧labels；定向增量状态见[新交付](notes/analysis-delivery-20260831.md)。[有限源报告](experiments/p334-finite-source-20260831/REPORT.md)
- **P398精确响应完成**：继固定η干预后，η=0一阶反号时刻和零频积分已计算。平稳重加权与动态贡献竞争，16维模型积分误差约0.5%，最低两极不足；本有限模型保留P2。[实验报告](experiments/p398-linear-response-20260831/README.md)

## 分析方式

以一个科学问题和一个可交付读数组织工作：先读对应最新结果，能用现有矩就直接计算；确实缺字段时一次补齐同一问题需要的混合矩。代码、结果和短解释一起交付。

精确枚举、证书、scorer、terminal algebra都是现成支持。需要它们解决具体问题时使用；它们不构成开始物理分析的前置队列。廉价并行理论仍可继续。

新分析保留源数据、量的定义和原随机块关系。一个块的多种视图共同传播协方差；新的事后分析明确记录为事后分析。其余细节按实际科学风险处理，见[Governance](GOVERNANCE.md)。

## 背景与沿革

三队地址和分工共用[仓库协调入口](https://github.com/LightChainr/Matching-One/blob/0e9d684e88c26b904da342b4c33cdc04057a3d07/docs/TEAM-COORDINATION.md)，主要通过仓库结果交接，避免重复计算与频繁消息。

- [Research Map](docs/RESEARCH-MAP.md)、[旧Roadmap](docs/ROADMAP.md)保留路线沿革，当前顺序由README、STATUS和NEXT-TARGETS承担。
- [Draft #267](https://github.com/LightChainr/Matching-One/pull/267)保存完整研究交接并仍在更新；本次只提供简洁入口，不改动该分支。阅读具体报告时使用固定commit。
- [8月31日本轮整理前的入口](https://github.com/LightChainr/Matching-One/blob/8a68cca866d7fbca7463e2167c3ff06128d5851f/README.md)保留旧判断的时间顺序。
- 数值阈值来源仍按方法保存于[data/literature_threshold_sources.json](data/literature_threshold_sources.json)，不把参考小数当精确常数。

文档更新在独立Draft分支交付；Issue已按本轮授权结题或调整优先级。**不合并、不删除研究分支**，代码、数据和历史报告原路径保持可用。

MIT，见[LICENSE](LICENSE)。
