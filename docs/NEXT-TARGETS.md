# 下一步：有限传递已补齐，只推进具体机制缺口

**2026-08-31，23时后续结果复核。** 当前结论统一读[STATUS](STATUS.md)。Xi、jump/reweight分解、强源相反尾和共同热/源profile闭合拒绝均已完成。三份意见的滞后部分已去掉，不能再按旧导航启动它们。

## 当前缺口（#337，P1；N50齐次计算已完成）

**仅靠0/1孔表、支持范围和尾概率的方案已被具体反例排除。** [本轮证明](../experiments/p337-continuation-feasibility-20260831/THEORY.md)甚至固定同一q曲线、唯一root和正斜率后，仍允许相反U。不要继续细化同类尾概率来代替连通信息。反例属于摘要松弛类，未否定原图的所有延续路线。

- **删除“先求齐次N50无条件U”的待办。** [固定完整计算](../experiments/p337-homogeneous-n50-20260831/RESULT.md)已给U=1.0615603877、V_S=+0.0543457827。有限零传递排除，正号延续存活；合同已经结束，不自动加N100、t/epsilon点或MC。
- 本次固定分解显示直接项+2.026626与斜率源项−2.082389强抵消。若理论模型只能预测局部非零而不能约束同一个root-normalized净V，就不足以成为下一生产的依据。不得把现成N50结果再拟合成新source或跨尺寸律，然后追认为前瞻预测。
- [全epsilon面核](../experiments/p337-face-kernel-20260831/REPORT.md)及[循环支撑积分](../notes/p337-connectivity-reduction.md)已完成。非循环A可全部积分，初始潜在块可保留全局span作卷积；当前有限N50也已有完整值。对尚未覆盖的一般N/epsilon，真正剩余是受控连通比较/score矩与pooled斜率，不能再细化支持尾概率或把条件B结果当总体。

共同温度＋同一S耦合的profile模型已退休，不再通过拟合b/c、改变密度项或补一个观测方向救回。额外独立耦合和几何各自的坐标尚未被否定，也不因此自动成为主线。

**有限源双律比较的全m>=64反号已由[2690f665](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)完成，暂不启动普通采样。** 本次[固定m64复核与必要预算](../experiments/p337-finite-law-window-20260831/RESULT.md)给斜向rank1概率约2e−20，仅见一次该事件的必要样本量已超过10^19。下一项若以这对作用量为目标，必须先给出保持原U语义的条件/twist等估计器及受控误差，不能直接把十台机器投入粗采样。不延伸旧N25四点峰值网格。

同一最新分支已解决N/m²有界联合极限的Poisson/full共存和pooled分母。**固定m的剩余理论问题**是控制oblique twist代价及受限扇区odds，并消去粗轮廓上界中的exp[O(N/m²)]体积因子；再引用联合极限或仅证明PSD/有限twist阶数不能完成这一步。上述进展不替代固定S、epsilon=1、多孔连通score矩的另一条延续问题。

[本次相对缝界](../notes/p337-fixed-m-relative-bound.md)已有`Delta_k≤50k log m`，常数不足，不要重复证明“只是表面阶”。同文的完整热族反例说明：即便rank1一致指数小，仍需受限sector odds控制；共同bulk消去或相同pressure不能代替它。执行分支[0dda27ba](https://github.com/LightChainr/Matching-One/blob/0dda27ba/notes/closed-source-s4-trace-transmission-result.md)的固定Q4归一化通道到原U传递也已完成；更新的[bea717e8 Q1闭合迹](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/p337-q1-closed-trace-transmission-result.md)已完成两个连接核、完整finite landing及原U严格负响应。删除“补Q1有限传递”的待办；剩余是local pair-to-cut intertwiner及尺度内容，不追加Q/seam点。

**估计器可行性已经进一步收缩。** [精确方差预算](../experiments/p337-estimator-access-20260831/RESULT.md)把原root/均值/分母给定后，m64的star普通独立热协方差平均仍需每几何至少约1.52e25样本才达SNR3。[twist审查](../notes/p337-twist-estimator-access.md)也表明：非负partition恒等式没有自动消除稀有质量、热中心化和取向差的抵消。下一算法只接受具体的相关估计/桥接重叠或条件积分二阶矩控制；仅提出“twist”“Rao–Blackwell”或“在rank1中采样”不再构成启动理由。这里没有否定全部替代算法，也不改开m2/增加参数点来回避当前预算。

## 已完成并停线的生产

1. **#154：固定新传递实验已完成，执行停线。** [165M新路径](../experiments/p154-prospective-transmission-20260831/REPORT.md)触发净U在两N均处于±0.50内的冻结规则；entry/completion各在±0.30内。停止把这个lag=1簇源列为当前主要H4解释，不接着扫lag、换源单位或补样。本轮没有选出连续场身份。
2. **#334：固定新群体实验已完成，执行降级。** [固定得分](../experiments/p334-prospective-intervention-20260831/results/latest.json)同时排除“残余投影接近0”与“旧残余点预测按±25%迁移”。封存两个失败预测，不追加prefix/尾部、不在验证块训练新模型，不把约1/2变成新的机制主张。完整contact解释仍未建立，已测局部响应保留为事实。

**现在不执行：** 新archive坐标、更多prefix特征、同块cross-size回归、泛化证书、新次数/高度/区间全族筛选。支持工作只服务一个具名机制预测的实际缺口；便宜不自动意味着P0。

结果出来后先按冻结规则否定/降级候选，保留失败和未分辨结果，不在同一验证块上加入第四个模型。需要新模型时另开未来探索阶段，其分析不能追认成本轮预测。

两个实验已完成，#154/#334一般问题保留P1，当前P0生产队列为空。#337端点研究的下一交付仍须约束**同一个global读出**。以下接口已经推导完成，直接使用，不再立项重做：令`q=-1+F1+F2`、`E=1-F1+F2`，源响应为`Jq,JE`，在共同根`mean_g q=0`处记`D=mean_g q_p`、`a=N^(13/8)/2`，则

```text
p0dot = -mean_g(Jq) / D
Ddot  = mean_g(Jq_p + q_pp*p0dot)
U     = a*P4(E_p) / D
Udot  = a*P4(JE_p + E_pp*p0dot)/D - U*Ddot/D
```

`P4`为两取向差除以冻结的DeltaCos4。对#154，`Jq=T01+2*T02+T12`、`JE=-T01+T12`，再按Bin(N,p)混合；直接0→2不能删掉。对#334，出生中心/间隔的积分响应并不决定这个含p导数、共同根和斜率归一化的泛函，因此“局部响应强”不足以推出`Udot`强。下一计算只代入一个事前指定的机制，给出它与竞争机制在本式上的符号或幅度差，以及失败后要停止的解释；无法产生区分预测时不启动下一块生产。不要从本轮验证块再选一个最显著的kernel、lag或contact特征来充当预测。

用户现已授权十台独立服务器按需直接调用，无需跨团队消息协调；每台仍核实当前任务/进程，不覆盖他人文件或停止不明进程。两个固定预算均已完成，本次精确判别没有云作业。N50的有限精确合同也已执行完毕；近期没有尚未执行的P0随机生产合同，不为凑满服务器而开新生产。理论执行交付一个余项结论，总览把结论接回STATUS，审查只核查这一个结论及停线动作；通过仓库交接，不增设会议或平行总账。

原13项Issue结题见[清理记录](REPOSITORY-TRIAGE-20260831.md)。结果通过Draft #509交付，**不合并、不删除历史数据和分支**。
