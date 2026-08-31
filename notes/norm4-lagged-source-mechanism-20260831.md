# Norm-4：一步滞后空间源与后续拓扑的最小机制读数

当前状态：**已完成探索性计算**，结果见 [REPORT](../results/norm4-lagged-source/REPORT.md)，提交 `4daae57eef5c945aa050a95cd3d5d5d77582161b`。旧2.4M排列的一步空间源能改变未来激活与moving-root rank1人口，原U方向响应仍未分辨；这批数据已经看过，不是后续模型比较的独立验证。

以下保留计算前的数学定义与数据设计；当时尚未运行计算、测试、replay或服务器操作。后续前瞻性比较见 [决策实验](../docs/DECISION-EXPERIMENTS.md)。

## 结论：下一项应测“一步前同 rank 空间源如何改变下一次激活”

推荐固定 **一个 activation 的滞后**，不扫描 lag，也不重复 same-time centered-source 的零响应。对最终占据数 K，先在 K−1 的配置中施加同 rank 内居中的 cluster 源，再按原均匀未占据点规则添加最后一个点。保留三个互斥跃迁 `0→1`、`0→2`、`1→2` 的源加权计数，足以恢复后续 matching/root、rank-1 population 和原 U 的一阶响应。

现有 seed/counter 能精确重建所需旧排列；**现有已聚合 CSV 不能直接给出这三个联合标记**。所需工作是给原排列补新时间交叉量，不是增加 Monte Carlo 随机样本，也不是重新检验静态投影恒等式。

依据 `8799dfe18782ca85ece4b69d508a27d7173bd28d:results/norm4-global-source-projection/REPORT.md`，原 bulk 源的 occupancy/rank 分量已完成；同一时刻的 spatial residual 对全局 q/E 与 U 严格零。以下对象改变的是源的作用时刻，保持读出定义及原基准曲线。

## 1. 明确的前向协议，而非把两个边际配成 joint

对固定几何 g 和原均匀排列 π，记第 j 个前缀的配置为 A_j，

\[
r_j=q_j+1\in\{0,1,2\},\qquad E_j=q_j^2,\qquad s_j=C_B(A_j)+C_W(A_j).
\]

在基准微正则分布下定义

\[
\mu_{g,r}(j)=\mathbb E[s_j\mid r_j=r,g],\qquad
\epsilon_j=s_j-\mu_{g,r_j}(j).
\]

因而 \(\mathbb E[\epsilon_j\mid r_j,g]=0\)。对每个最终 K≥1，令 j=K−1，并采用以下固定协议：

1. 保持 j 与早期 rank r 的基准概率；在该 `(j,r,g)` 层内按 \(\exp(t\epsilon_j)\) 倾斜空间配置，并在层内归一化。
2. 从剩余 N−j 个点中均匀添加一个点，使用原拓扑演化规则，读取 \(q_K,E_K\)。
3. 最终 K 的外部权重仍取 \(w_K(p)=\mathrm{Bin}(N,p)\)；K=0 定义源响应为零。

层内归一化使早期 rank 分布保持不变，均匀 continuation 核也没有改变。t=0 时，均匀 j 配置加一个均匀未占据点恰好给出原均匀 K 配置，所以整条基准 q/E 曲线、p₀ 与原 U 完全相同。其一阶路径源就是 \(\epsilon_{K-1}\)。

这是“在所选观测步之前一个 activation 注入源”的明确协议。它不是一个固定早期 Bernoulli 参数下、跨所有后续时刻共同使用的单次注入，也不是原平衡 bulk fugacity \(\exp(t s_K)\)。如果以后改为固定早期 p_source 或更长 lag，需要另写相应的联合时间权重；不能沿用本式冒充那个实验。

## 2. 真正的新约束来自 future rank，不来自重新证明 same-time 零

same-time 条件期望立即给出

\[
\mathbb E[\epsilon_j q_j]=\mathbb E[\epsilon_j E_j]=0.
\]

未来读出则为

\[
H_q(K)=\mathbb E[\epsilon_{K-1}q_K]
       =\mathbb E[\epsilon_{K-1}(q_K-q_{K-1})],
\]

\[
H_E(K)=\mathbb E[\epsilon_{K-1}E_K]
       =\mathbb E[\epsilon_{K-1}(E_K-E_{K-1})].
\]

右侧一般不受 same-time 投影恒等式约束。若“给定占据数与 rank 就已足够预测该源倾斜后的 continuation”，则这两个核为零。非零说明同早期 rank 内的空间配置通过本源与下一步拓扑跃迁发生耦合，排除这一源不敏感的 rank-only 闭合。

此结论针对指定源与均匀加点协议；不能仅凭非零就宣称未扰动 rank 序列在所有意义下非 Markov，或已识别连续极限场。未扰动的平均 rank 转移矩阵仍可能描述其自身边际；这里问的是它能否在源改变层内配置后继续预测。

## 3. 最小数据核只有三个跃迁类别

沿用精确的 \(K_1=K_{minus}\)、\(K_2=K_{plus}\)。每条排列最多需要记录两个事件：

- `01`：K=K1<K2，早期 rank 为 0。
- `02`：K=K1=K2，早期 rank 为 0，直接完成两次激活。
- `12`：K=K2>K1，早期 rank 为 1。

对每个 `(N,g,batch,K,type)`，新增两个整数总和即可：`event_count` 与 `sum_s_previous`，其中源总是同一排列上的 \(s_{K-1}\)，单位为 bulk cluster count。不要存成另一个排列的 s，也不要使用 \(s_K\)。

设保留样本数为 M，三个居中事件核为

\[
T_{01}(K)=\frac{S_{01}(K)-\mu_{g,0}(K-1)n_{01}(K)}M,
\]
\[
T_{02}(K)=\frac{S_{02}(K)-\mu_{g,0}(K-1)n_{02}(K)}M,
\]
\[
T_{12}(K)=\frac{S_{12}(K)-\mu_{g,1}(K-1)n_{12}(K)}M.
\]

则

\[
H_q=T_{01}+2T_{02}+T_{12},\qquad H_E=-T_{01}+T_{12}.
\]

首次激活与第二次完成的源核分别为 \(T_{01}+T_{02}\)、\(T_{12}+T_{02}\)。**直接 `0→2` 对 E 的增量为零，但对 q 的增量为 2；不能因为它不进入 H_E 就丢掉它，它仍影响根位移与 U 的分母。**

这里不用保存所有 \((j,K)\) 的二维矩阵，也不用估计一组自由 lag。原 full-filtration 代码在每个排列内已经同时持有 K1/K2 与全部黑/白 component 数组；在数组被压缩成 batch 总和之前，把上述事件加入新总和即可。

作为解释边界，若只选旧报告中的 rank-1 居中源 \(R_j(s_j-\mu_{g,1}(j))\)，这个一步实验只剩 T12：它不能影响发生于早期 rank-0 的首次激活。全 spatial residual 的三类核避免把该支持限制误读成“没有 birth 机制”。

## 4. 原 U 与 moving root 使用同一套基准，但必须重算新源 rootdot

记 \(Q_g(p)=\mathbb E[q_K]\)、\(\mathcal E_g(p)=\mathbb E[E_K]\) 为原基准曲线。新源均值在每个 K 都为零，故

\[
J_q^{lag}(p)=\sum_Kw_K(p)H_q(K),\qquad
J_E^{lag}(p)=\sum_Kw_K(p)H_E(K).
\]

计算 p 导数时，所有早期微正则居中均值与 lag=1 协议固定，只微分 w_K。令

\[
D=\overline{Q'},\quad B=\mathcal P_4[\mathcal E'],\quad
H=\mathcal P_4[\mathcal E''],\quad T=\overline{Q''},\quad A=N^{13/8}/2.
\]

在原 pooled p₀ 上，新源的三个联动读出为

\[
\dot p_{lag}=-\overline{J_q^{lag}}/D,
\qquad
\dot R_{1,lag}=-\overline{J_E^{lag}}-\dot p_{lag}\overline{\mathcal E'},
\]

\[
v_{lag}=\partial_tU
=A\left[
\frac{\mathcal P_4[(J_E^{lag})']+\dot p_{lag}H}{D}
-\frac{B\{\overline{(J_q^{lag})'}+\dot p_{lag}T\}}{D^2}
\right].
\]

可以复用原 central/LOO 的 p₀，并从现有 q/E profile 求该根上的基准热 jets；**不能复制原同一时刻 bulk s 的 rootdot**，因为这里的 Jq 已换成早期居中源。上述源已是 bulk 单位，不额外乘 N。保持原精确方向归一化和 N^(13/8)，不重拟合指数。

最小首项科学读出建议为 **lag=1 的 root-comoving rank-1 响应**，同时给出由同三个事件核确定的 rootdot 与 v_lag。这样，若人口响应已分辨而 v_lag 仍弱，仍能明确报告前向空间机制已到达哪个观察者层，而不把人口变化提升为 H4 算符识别。原 q2/Jordan 对同一时刻 bulk fugacity 的源延伸假说不能自动套用到这个新协议。

## 5. 已有档案究竟够什么

本次读到的 source CSV 每个 batch/方向/K 保存 `samples,q,E,s,qs,Es` 的总和；line CSV 保存同 K 的 `R,Rs,RO4,RO4s` 总和。它们足够恢复早期各 rank 的居中均值。用固定 K 的微正则矩 x=〈q〉、e=〈E〉、m=〈s〉、tq=〈qs〉、tE=〈Es〉，有

\[
\mu_{g,0}=(tE-tq)/(e-x),\qquad
\mu_{g,1}=(m-tE)/(1-e)
\]

（仅在相应状态有支持时）。分母为零意味着该早期状态没有样本，其对应事件计数也必须为零，不作伪造插补。

这些 CSV 没有 per-counter 轨迹行，也没有 \(s_{K-1}\) 与后来 K1/K2 事件的交叉总和。两个 batch 均值相乘不能替代它们；原 threshold 边际直方图及 K1/K2 的联合低阶矩也没有这个源标记。现有跨 K batch 协方差是估计量的不确定性，不能用作所缺的逐排列两时刻矩。

因此：居中均值、基准曲线、热 jets、原根及分组可直接复用；T01/T02/T12 的源加权事件总和需要新标记。它们能从**原排列重观察**得到，但不能从已保存的同一时刻边际直接聚合出来。

## 6. 保持原随机流与 paired uncertainty 的可执行选择

设计优先沿用已标记的六 N 范围，避免另引入随机流：

- N65/85/130/170：seed `2026104501`，counter `[5100000000,5100100000)`，每个 N 为 100×1000；四 N 属于同一依赖组。
- N260：seed `2026105401`；N340：seed `2026105402`。各自使用 `[8200000000,8201000000)` 的 1M 原排列，两个端点各成一组。
- 端点 batch b 严格沿用 `[8200000000+1000b, +1000)` 与 `[8200100000+9000b, +9000)` 的 union，而不是改成连续的新 10000 分块。两方向同一排列共同删除。

每次对齐留一同时删除新事件标记与旧同 K source 批次，重新估计 \(\mu_{g,r}(j)\)、事件核和新源 rootdot，使用该原留一的 p₀。报告三事件核及三个观察者读出的完整协方差；不把同一 source 的 birth/exit、两方向或不同 K 当作独立实验。

经验居中是数据适配的分解：以上层内倾斜协议的总体定义是固定微正则均值，plug-in 及其留一传播的是该均值的估计误差。低支持状态需保留计数与权重信息；同一数据构造出的 early centered-null 仍是代数，而未来事件核才是新结果。此设计未运行，采集时间、信号大小及可分辨性均未知。

## 已读来源：精确 commit:path

- `8799dfe18782ca85ece4b69d508a27d7173bd28d:results/norm4-global-source-projection/REPORT.md`：已完成 occupancy/rank 实际分解与 same-time residual 零的边界。
- `8799dfe18782ca85ece4b69d508a27d7173bd28d:analysis/norm4_global_source_projection_contract.json`：bulk 单位、固定微正则投影、原 U 与配对留一合同。
- `6f99f901cb919f488fcd6c2b4f4c357360b17764:src/norm4_source_thermal_replay.cpp`：每排列的完整 component 数组、Kminus/Kplus，随后压缩为同 K batch 总和。
- `584559582ded1f55d06f69831b3f573cebf6e673:src/norm4_source_line_replay.cpp`：物理绕环方向、rank-1 支持与既有 endpoint union 的保存方式。
- `6f99f901cb919f488fcd6c2b4f4c357360b17764:scripts/replay_norm4_source_endpoint_increment.py`：端点旧 counter 增量与不可替换的源流。
- `6bd46ad30bc8f583c3ca1f1c8a1b95e7d90571bc:scripts/analyze_norm4_source_endpoint_1m.py`：原/增量 loader、密度转 bulk 的单位关系及基准 root-comoving 读出。
- `3e6157f237242938e1c1b12415bca256b11896b0:results/norm4-source-thermal/raw/n65.csv`、`6bd46ad30bc8f583c3ca1f1c8a1b95e7d90571bc:results/norm4-source-endpoint-1m/increment/raw/n260.csv`、`91fe30f3ebbe976d65db2538093af5e2c45b11d0:results/norm4-source-line/raw/n65.csv`：本次直接读取的代表性 CSV 表头，确认按 batch/K 而非逐排列存储。
- 原 backend 的精确地址由已读 replay driver 固定为 `bfab0330f5f56ca4d746b45d737f1607e3d229a0:src/threshold_rank_orientation_mc.cpp` 与 `bfab0330f5f56ca4d746b45d737f1607e3d229a0:src/threshold_rank_integer_period_mc.cpp`；后者保留 N260/N340 的真实 HNF 几何，不能用其 a,b 标签换成 cyclic quotient。

这是一项已具备具体读出与最小新增联合量的 hypothesis，不是已完成的两时刻生产结果；也没有将新的 lagged source 等同于原 bulk-Q 源或原 U 的唯一物理机制。
