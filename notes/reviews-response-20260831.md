# 三份意见的核对与实际推进（2026-08-31）

三份原文已全文阅读。此页是这次阅读和结果的记录，不增设当前状态入口；当前判断仍从[STATUS](../docs/STATUS.md)进入。原文涉及的两个sandbox下载包未取得，其代码/验证结果不被视为本地已复现；下列验证从原文明确给出的定义独立实现。

| 输入 | SHA256 |
|---|---|
| 意见1 | `1cda9a3b61ac7dcf32f1df533f245b19b22c9cc05f8812838efa8c0cdd9d2b42` |
| 意见2 | `a370cd9f8fb10fe7d70b6be5e210889308647a7ca55e3978fb1d3124698879fb` |
| 意见3 | `31517bb3fe3bc4ac75b9a172d1a49558bf27320e987e54a3837abaf7139c6ec2` |

## 先去掉已过时的待办

- 首次Xi、整体增益拒绝、共同热/S四profile闭合拒绝均已完成。
- 意见3建议的weighted-jump/reweight分解也已完成：总览[e1b96895](https://github.com/LightChainr/Matching-One/blob/e1b968959634b9b3999c727b83ed38d0b730cb20/results/defect-reweight/REPORT.md)给`Xi_reweight=+4.550327123237`、`Xi_jump=-15.306045530801`，合计`-10.755718407564`。执行9057325d是同一有限总体的另一实现，不重新计票。
- 齐次N25有限t四点、强源负尾、一般L首项、Sdrop正尾已完成；[f4057192](https://github.com/LightChainr/Matching-One/blob/f405719264c896aa873dd4aae7292795f544ba99/docs/NEXT-TARGETS.md)要求的是有限耦合定量判别/可达性，不再等待尾部符号。
- [3dc47674](https://github.com/LightChainr/Matching-One/blob/3dc47674899a9ca93dc7d667f03c537fdb954bbf/notes/closed-source-hypergraph-rc-twist-projection.md)已经实现hypergraph RC与twist读出。再次声称“存在四端口表示”不能算新成果。

**汇报前再次核对的进展：** 执行[2690f665](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)已把反号证明扩展到整个实数m>=64；[同提交的Poisson双极限](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md)在N/m²有界的指定联合极限闭合了oblique/pooled分母并证明U超多项式衰减。本次m64符号计算与它并行，作为复核记录；不能再把“首次有限区间反号”列为下一步。固定m的oblique twist代价与受限扇区odds仍未由这些结果解决。

上述结论按固定SHA读取；不把branch-only写成main或已合并。P154/P334/F4的既定停线保留。

## 已经实际完成的三个增量

| 评审指出的问题 | 本次交付与改变的判断 |
|---|---|
| 孔尾概率并不是U误差界 | [信息不足见证](../experiments/p337-continuation-feasibility-20260831/THEORY.md)保留真实0/1孔表，构造同q曲线、同唯一root、同正斜率的两个摘要补全，原U却分别约+10.10358与−10.07432。连分母已知也不能只靠支持/尾界定出U符号；需额外的多孔连通score矩控制。这是摘要松弛类证明，不是原图两种物理模型。 |
| 权重与连通读出必须共同消元 | [全孔密度面核](../experiments/p337-face-kernel-20260831/REPORT.md)保留面/端口/位移/平行边，给出精确cycle校正及源权重。复现原文两孔反例，并给该固定B整行任意epsilon、有限t的闭式绕环概率，免去孔阶数截断。它是条件连通律，不是总体U已求出。 |
| 原U反号是否能用普通采样分辨？ | [唯一固定m=64双律计算](../experiments/p337-finite-law-window-20260831/RESULT.md)先冻结375a6f0c后计算；符号复核并行半直线结论。新增rank1概率使普通抽样见一次稀有事件的必要预算达到10^14–10^19级，故不启动该点普通生产；不等于拒绝条件/twist等估计器。 |

没有新增孔型全枚举、随机样本或云作业，也没有把64加入旧四点实验。新固定点只给N25有限结论，既不是有限宽t区间的证书，也不是增长N下的统一尾界。

## 仓库取舍

意见1/3指出PR528会恢复旧P0计划，核对属实。本次在隔离维护工作树修正其现有导航/账本，保留历史条目和数据，不合并研究Draft。大Draft的审查入口是上表各个具体包、冻结SHA及验证脚本，不要求从完整diff重建历史。

不再追加一个总规划。下一分析从已经明确的缺口出发：把全孔面核的连接约束用于未知多孔层的取向热score矩，或为具名双律给出避免稀有扇区瓶颈的原U估计器及误差证明。它们都未因本轮有限结果自动完成；在形成可区分预测与可行预算前，P0生产队列保持空。支持线不回升为默认主线，也不以新特征补救已失败模型。
