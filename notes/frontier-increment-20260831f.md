# N900 已完成：第三尺度约束了宽度，尚未唯一识别机制

截至本次定点捕获 **2026-08-31 07:29:32Z**，[PR #484](https://github.com/LightChainr/Matching-One/pull/484) 为 `open`、非 Draft；本笔记固定来源分支 `analysis/etop-modulus-survivors-20260831` 的完整 commit **`5f30397c5ba277fb0799fb2f7491c823de07a13d`**，按 `open_pr` 记录，不将其写成 `main` 已集成事实。此前“运行中、尚无 N900 目标结果”的正文已滞后：报告、score 和两份成功回执已经交付。本次仅阅读这些产物，未检查进程、操作服务器或另启生产。

## 已完成的物理读数

N900 使用 **32M distinct shared counters、800 aligned batches、two modulus pairs**：`tau_2i` 与 `tau_4i` 各自包含原来的两个方向，沿用精确方向归一化 `±1152/625`。这是 64M geometry-pair evaluations，独立样本块仍只有一个 32M counter block；它与 N100/N400 独立，两 modulus 之间共享随机流。几何为 N100 周期矩阵逐项乘三，保留 moduli 与有理旋转；没有新增 shear，因此不承担三形状共同 transport 的检验。

主读数是带符号 rank-step profile
`D_A = P4[A_top](4i) − P4[A_top](2i)` 在 `z=N^(3/8)(p−p_ref)` 下的中心二阶矩，`p_ref=0.59274605079`。测得 **`V_z=2.339461729 ± 0.120385154`**，保留面积归一化和同 counter batch 协方差。这是拓扑观察者的 signed-profile 宽度，不把它解释为正概率分布或局部能量场的算符识别。

两条有限尺度预测均由 N100/N400 提出，并在 N900 采集前冻结：宽度按 `N^(−1/4)` 缩放对应 `V_z(900)=V_z(400)(900/400)^(1/4)`；固定 critical-width profile 对应 `V_z(900)=V_z(400)`。

| 条件预测 | 预测值 ± anchor SE | 观测减预测 ± total SE | z | nominal p |
|---|---:|---:|---:|---:|
| quarter-power width | 2.565535388 ± 0.091547708 | −0.226073660 ± 0.151240101 | −1.494800 | 0.134967 |
| fixed critical-width profile | 2.094750873 ± 0.074748391 | +0.244710856 ± 0.141703589 | +1.726921 | 0.084182 |

## 为什么现在不能强行给出 winner

两个残差共享同一个 N900 target，也共享完全相关的 N400 anchor。令预测为 `a_i V_400`，则 `Cov(r_i,r_j)=Var(V_900)+a_i a_j Var(V_400)`；score 保存的完整比较协方差为 `[[0.0228735680, 0.0213356290], [0.0213356290, 0.0200799071]]`。两个 p 值不是两次独立支持，不能相乘、计票或当成模型后验概率。

score 未另给 prediction interval。沿用其双侧正态比较尺度，N900 观测离两预测均不足 1.96 个 **total SE**，因此两模型包含 target 与 anchor 误差的 95% 观测相容范围确实有交集，当前观测就在交集中。这一判断不能用仅含 anchor 的 `prediction_se` 区间替代，也不等价于两种机制相同。N400→N900 的有效宽度 **`0.306876850 ± 0.038610154`** 只是相邻两个面积之间的有限尺度描述，不是新渐近指数。

## 下一实质目标：让宽度机制给出额外物理预测

第三尺度的首次采集与宽度比较已经完成。下一步应让候选机制明确预测：同一 `D_A` 的宽度改变是一个内禀时钟的整体缩放，还是多个拓扑响应分量权重随尺度变化；它们对**面积、中心位置及去宽度后的形状**必须给出不同的、带共同协方差的联合响应。现有 score 已保存面积／位置／宽度的联合协方差和 800 个 batch 原始矩，现有 threshold histograms 也可复用；先消费已经交付的形状分析，再补模型真正不同的预测，不把“首次做形状”重新列成待办。

这是一项机制区分目标，而非要求再增加一个自由指数、再按同一宽度重排模型或自动追加 N900 样本。若两个候选只规定一个二阶矩，却没有 source、几何或形状响应的区别，当前产物支持的是两条条件缩放曲线，而不是已完成的微观机制模型。

## 精确来源与完成回执

- [REPORT.md](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/REPORT.md)；[score.json](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/score.json)。原始输入为同目录 `raw/tau_{2i,4i}.{hist.csv,moments.csv,metadata.json}`，SHA256 保存在 score。
- [tau_2i.receipt.json](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/logs/tau_2i.receipt.json) 与 [tau_4i.receipt.json](https://github.com/LightChainr/Matching-One/blob/5f30397c5ba277fb0799fb2f7491c823de07a13d/results/etop-n900-rank-width/logs/tau_4i.receipt.json)：均 `exit_code=0`，分别耗时 `3339.989715`、`3676.290414` 秒；同 seed `20260831141001`、counter offset `267900000000`，各自为单线程既有引擎运行。
- 预测来源 `fb1a944e1ef34e9b9dfcf32c59af25f44ce43d9a:results/p267-rank-clock-width/score.json`；采集前冻结／runner commit `ecde7c9132ed35ec1575bd82f11e816722912e6f`。本笔记未重算任何统计量。
