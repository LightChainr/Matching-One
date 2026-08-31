# 团队收敛与仓库任务清理：2026-08-31

本页记录本轮已经执行的GitHub操作；科学分工见[NEXT-TARGETS](NEXT-TARGETS.md)。关闭与优先级调整已生效，文档继续在[Draft #509](https://github.com/LightChainr/Matching-One/pull/509)审阅，未合并。

## 清理结果

| 项目 | 清理前 | 清理后 |
|---|---:|---:|
| 开放Issue | 113 | **100** |
| 标为P0的开放Issue | 57（旧标题） | **2** |
| 当前P1 | 35个旧P1标题 | **7** |
| 当前P2 | 15个旧P2标题 | **91** |
| 期初开放PR | 24 | **原24个保留，不关闭独立成果** |

本轮关闭11个已完成交付或历史综合任务、2个被其他入口承接的任务。给100个仍开放Issue设置唯一priority标签，并同步82个原标题中的优先级；优先级操作保留正文与assignees。#154/#334/#398另外增加具体团队交接，#13增加P2支持范围说明，原文均保留在下方。原113项中6项没有P0/P1/P2标题，现存100项统一按标签分级。

## 已关闭的13项

| Issue | 结题方式 | 依据与后续归属 |
|---|---|---|
| [#40](https://github.com/LightChainr/Matching-One/issues/40) | completed | 百万same-N motif评估完成，收益有限，按原规则结束pilot；后续源问题归#154。 |
| [#55](https://github.com/LightChainr/Matching-One/issues/55) | completed | 600M/design实验与评分完成，H12未分辨；按停止规则不加第三别名行。 |
| [#111](https://github.com/LightChainr/Matching-One/issues/111) | completed | Euler/Betti恒等式与控制量约化由PR196交付；方差改进支线结束，PR仍开放。 |
| [#216](https://github.com/LightChainr/Matching-One/issues/216) | completed | Q4 Jordan继承/Gram/log-slope代数完成；完整torus形状归#220，实测overlap归#154。 |
| [#269](https://github.com/LightChainr/Matching-One/issues/269) | completed | 一般有理同调与整数饱和证明已交付；分支未合并，后续birth/source归#334/#337。 |
| [#403](https://github.com/LightChainr/Matching-One/issues/403) | duplicate / not_planned | 精确H2/trigger成果保留；后续总体与空间问题收敛到#334，branching问题留#429。 |
| [#406](https://github.com/LightChainr/Matching-One/issues/406) | completed | Bochner、≥100-mode witness和档案分析完成；采用归一化修正版，后续预测归#250。 |
| [#418](https://github.com/LightChainr/Matching-One/issues/418) | completed | CRT/单位核对与档案重分析完成；旧共同谱大惩罚已纠正，后续归#250。 |
| [#486](https://github.com/LightChainr/Matching-One/issues/486) | duplicate / not_planned | cut/update目标由#487/PR491交付；后续统一到#334。 |
| [#487](https://github.com/LightChainr/Matching-One/issues/487) | completed | 原cut theorem、W2和三阶机制由PR491交付；总体/空间续题归#334。 |
| [#265](https://github.com/LightChainr/Matching-One/issues/265) | completed | 历史战略审计已吸收，具体科学问题留在相应主Issue。 |
| [#466](https://github.com/LightChainr/Matching-One/issues/466) | completed | 综合意见已由#154/#334/#398承接，不再保留平行策略队列。 |
| [#468](https://github.com/LightChainr/Matching-One/issues/468) | completed | 机制优先建议已吸收；具体问题留在#275/#334/#398/#439等入口。 |

每个关闭Issue顶部均保留完成依据、来源链接和后续归属。completed表示原任务交付结束，包含阴性或未分辨结果；不表示对应物理机制已证明，也不表示结果PR已经合并。

## 当前优先级

- **P0：#154、#334。** 本轮主要分析交付，建议合计90%的注意力。
- **P1：#12、#14、#275、#337、#398、#408、#439。** #398做一次有界参数干预；其余是储备或主线接口，不同时启动七套生产。
- **P2：其余91个开放Issue。** 保留问题与已有成果，按需支持；不默认追加大样本、设备租赁、通用证明器或同类证书。

可直接筛选：[P0](https://github.com/LightChainr/Matching-One/issues?q=is%3Aissue+is%3Aopen+label%3A%22priority%3AP0%22) · [P1](https://github.com/LightChainr/Matching-One/issues?q=is%3Aissue+is%3Aopen+label%3A%22priority%3AP1%22) · [P2](https://github.com/LightChainr/Matching-One/issues?q=is%3Aissue+is%3Aopen+label%3A%22priority%3AP2%22)。

几个具体降级理由：#158的现有效率需要数百亿量级样本，应等待新估计量；#30/#31尚无对应GPU/宽度实测，工具准备不构成新增设备理由；#146/#370已具备大量可用语义/认证支持，下一次使用由实际分析调用。#144/#244/#400/#401/#419/#429仍有实质未解问题，因此保留，没有为压低数量而关闭。

## 为什么24个PR仍保留

它们包含独立的代码、原数据、负结果、精确证明或协议，不能把“已完成分析”和“可以删除审阅入口”混为一谈。例如#267保存独有百万源响应及协方差，#273/#277保存原norm-4生产，#484/#485保存全分布与形状结果。

全部添加`status:unmerged-asset`，明确这是未合并交付，不是待重新计算的工作。没有关闭这些PR，没有改变其Draft状态、base或分支，没有开启自动合并；也没有删除worktree、代码、数据或冻结预测。

**执行期间的新增事项单列。** 后续快照截至#519，又出现8个开放代数目录PR。逐项复核时#510/#514已关闭，本次未修改；仍开放的#511/#515/#516/#517/#518/#519均加P2与未合并成果标签。它们自己声明未涉及可靠性、临界方程或严格界，因此不进入主线生产。#13顶部已明确：暂不默认扩相邻代数目录，下一项集中投入应连接具体比较不等式、局部变换或有限类反例。上述期初24项与新增项分开计数，不把动态总量写成固定24。

本页是一次已执行操作的记录，不增加新的登记或审批要求。以后结果变化时更新对应主Issue的短交接和STATUS即可，不重新建立一套战略索引。
