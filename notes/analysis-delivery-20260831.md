# 8月31日：新增分析交付与团队下一步

本轮由“数学研究俯瞰”完成具体分析，沿用“数学研究执行”和“数学研究总览”刚交付的输入。三队共用[仓库协调入口](https://github.com/LightChainr/Matching-One/blob/0e9d684e88c26b904da342b4c33cdc04057a3d07/docs/TEAM-COORDINATION.md)，后续以仓库交接为主，不另建消息循环。结果继续进入[Draft #509](https://github.com/LightChainr/Matching-One/pull/509)，不合并。

## 已完成的四项新判断

| 实际问题 | 新结果 | 产物 |
|---|---|---|
| P154：强条件空间响应是否进入原全局U？ | 六N角权分配已完成。固定K/rank1内部中心化源对原U严格为零，软两支相反；原总源导数仍未分辨。强条件O4本身不能直接充当global H4来源。 | [报告](../experiments/p154-spatial-localization-20260831/REPORT.md)、[完整同批结果](../experiments/p154-spatial-localization-20260831/results/latest.json) |
| P334：即时Euler/rank不变的空间改变，是否具有有限未来效应？ | 既定q_t在t=±1仍有可测的future S(A)/D(A)响应。完整标签census给出精确归一化；有限奇部接近一阶预测，微小非线性已量化。 | [报告](../experiments/p334-finite-source-20260831/REPORT.md)、[20批结果](../experiments/p334-finite-source-20260831/output/latest.json) |
| P334：同一prefix是否有两个未来响应方向？ | 完整census加固定定向续接后，A的局部行列式均值在两N达到5.4/7.6SE；积分A也支持两个出生中心的局部二维响应。E/间隔和四阶平方量仍弱。 | [报告](../experiments/p334-mechanism-response-20260831/REPORT.md)、[最终20批结果](../experiments/p334-mechanism-response-20260831/results-extension/score.json) |
| P398：原传播反号是否已在一阶出现，谁控制它？ | η=0解析一阶响应在t≈1.04799反号，负平稳重加权与正动态贡献竞争；零频积分仍负，16维几何模型积分误差约0.5%，最低两极不足。 | [报告](../experiments/p398-linear-response-20260831/README.md)、[完整结果](../experiments/p398-linear-response-20260831/results/latest.json) |

P154角权桥、P334有限源和P398解析响应三包已在`fb01c44aa45e4f8d37d52144e2ad7c4adfe6ce40`提交并推送；P334局部响应包随本页此次更新交付。P154使用原100k/1M源子集与三组100次共同删批，未替代高精度普通生产。有限P334使用原40k prefix及旧续接，属于配对importance估计，未按新策略直接采样；P398为同一有限模型上的解析导数方程，数值采用float64。它们都没有增加独立prefix证据。

## P334同构形局部问题：先解决估计支持，再定向增加信息

总体2×2响应矩阵rank2不能证明同一个prefix的矩阵也rank2。第一轮跨quartet无偏det/det²发现，同类标签mask只有46/36个prefix存在非共线源对，所有prefix的四quartet支持均为零；`0±0`在此表示无信息。

本轮完整census允许使用精确score `H_o(u)=pi_a[L_o(u)−mean_a L_o]`，从同类mask改为全label半差。该改进已在HZ实际完成，恢复1496/1549个prefix的四quartet支持，20批全部覆盖。精确源Gram在全部1502/1551个双R0 prefix均rank2；仅使用旧8组时，未来响应的局部det/det²仍未分辨。这区分了“允许两个独立源”和“未来读出具有两个独立响应”。

据此已完成一次固定预算增量：只对这3053个原双R0 prefix各补64个独立quartet，合计781568条条件续接。原prefix、几何、源和20批不变；新quartet用bit31随机键域，避免简单延长旧quartet索引造成键碰撞。生成用时4.164秒，分析2.220秒，预算完成后停止，没有按显著性继续追加。最终结果已保存到[局部分析包](../experiments/p334-mechanism-response-20260831/REPORT.md)。

| 最终局部目标，原20000prefix分母 | N325 | N425 |
|---|---:|---:|
| E[det J_A(p_ref)(Z)] | 1.15937×10⁻⁸ ± 2.14282×10⁻⁹ | 1.49767×10⁻⁸ ± 1.97301×10⁻⁹ |
| E[det J_integral A(Z)] | 5.32779×10⁻¹⁰ ± 7.33706×10⁻¹¹ | 5.26807×10⁻¹⁰ ± 6.77299×10⁻¹¹ |

若每个prefix的A条件均值响应均为rank≤1，这两个局部行列式均值都应为零。结果提供局部二维成分存在的数值证据，不能再完全归因于不同prefix的混合；也不说明每个prefix都rank2。积分A对应两个几何的出生中心，因此该坐标有直接机制解释。E/出生间隔的局部行列式仍未分辨；四阶平方量误差仍大，保留其负估计。新64与原8共用原prefix，不能叫独立总体复现。

## 已立即改变的分析顺序

1. P154不再反复检测同层O4关联。下一份机制预测要触及rank人口及进入/退出；必须对原U或其源导数给出可核对的改变。
2. P334的局部A二维成分已出现。执行队在本轮交付期间也已推送[协方差层级分解dc4bb041](https://github.com/LightChainr/Matching-One/blob/dc4bb041522176c7d992c56b018c9e91a96d58c4/notes/p334-birth-covariance-hierarchy.md)：总体协方差响应主要由prefix均值移动承载，固定rank-cell内还有可分辨的prefix均值项，固定prefix内部协方差变化较弱。下一步直接用已有接触/闭环及基线出生中心特征预测层内均值响应，再对照本轮局部二维结果；不重复已完成的层级分解。先消费完整census与新增续接；普通全域续接、首轮全曲线、有限t=±1与首次局部检验均移出待办。
3. P398这一有限模型的动态解释已交付，保持P2。没有明确晶格源/读出映射前，不自动扩width、扫参数或继续堆模态。

四台实际使用机器的产物均已取回，任务及各自隧道结束；2026-08-31 18:26后的舰队查询确认五台均为Ready，XP本轮始终未启动。P334的独立有理数公式复核及120项生产实现比对通过，代码、结果、输入hash与20批共同协方差均已交付。五台按需使用，核数空闲本身不构成新增样本的理由。连接技能已补入同一次会话起停后的密钥失效与恢复记录、实际运行环境及安全的进程检查方法。

仓库Issue交接：[P154角权结果](https://github.com/LightChainr/Matching-One/issues/154#issuecomment-5476888766)、[P334有限源与局部响应](https://github.com/LightChainr/Matching-One/issues/334)、[P398动态结果](https://github.com/LightChainr/Matching-One/issues/398#issuecomment-5476888961)。具体可复用输入与结论以上述四个结果包为准。
