# 下一步：thermal/pivotal 传递与固定m两相权重

**2026-09-01，已读到执行 `410015f5`、PR533 `5aa929a6` 与 #537 landing preflight。** 当前结论从[STATUS](STATUS.md)进入。两个独立生产实验、齐次N50、canonical pair 的有限空间块、完整联合原U响应和第一项rank-one有限反例均已完成。下一步不增加同类有限点、距离网格、补全扫描或旧实验救值。

## 1. 当前研究问题：绝对可和的空间核能否经热响应获得尺度增强

执行[2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)已经定义并算完canonical `Kreg=K2+K0`：直接Q1响应恒零，单点混合原U响应严格非零，固定四路径双插入系数13/8。随后 [`a237968f`](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md) 拒绝了L64/距离16的有限零传递；[`410015f5`](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md) 完成完整 `J2=-0.0055194314248394015`，并以严格负的非相邻分项拒绝NN-contact-only closure。旧“local pair还没有进入有限U”“先计算uniform J2”两个待办全部删除。

**[临界空间可求和定理](../notes/p337-critical-spatial-summability.md)已经完成下一层淘汰。** 至少两个共享occupied components在每个远端点都强制交替四臂；两个不交环带和方格site严格`alpha4>1`界给

```text
E_pc |g_xy| <= C d^(-2-eta),
sup_(L,x) sum_(y!=x) E_pc|g_xy| < infinity,
tail(d>R) = O(R^-eta)
```

其中某个`eta>0`。所以固定normalization的raw canonical susceptibility不能靠宏观远尾产生发散。这个结论不使用方格5/4猜想，也不自动适用于a237的数值`p_ref`或N25 pooled root；没有近临界控制时不作外推。

新的[精确审计](../notes/p337-thermal-pivotal-gate-audit.md)已经证明，下一项不能只研究 `partial_p E[g_xy]`。对 `O=q,E`，完整 `d_p Cov(O,g_xy)` 有 kernel reconnection 与原 rank/readout pivotal 两项；tiny torus上遗漏第二项会把符号判反。普通site flip是八端口组件join，但共享组件数本身也不决定符号。

现在只保留 #537 一个P0理论闸门：**把未定义的 ordinary four-arm/no-extra-branch 变成一次可判定的有限合同，再决定rank-one路线是否已经死亡。** [现有精确preflight](../experiments/p537-landing-matrix-preflight-20260901/REPORT.md)在six-block clean-two-bridge合同下已给出全部非零minors；总和的P4/root determinant为`-2.4843232721775393e-5`，固定M的Schur residual为`-4.217141611550048e-6`，`T_p/M_p=+5.806332966676667e-6`。因此不能继续把“构造finite matrix”列为未执行，也不能先进入五/六臂概率证明。

1. 冻结allowed Bell-8 landings、全局`no-extra-branch`布尔定义、逐`z` before/after landing与rank transition、row basis和C4 character；这是一次语义交付，不是新模型目录。
2. 若正式定义包含已有clean-two-bridge事件，立刻退休“四臂leading block只是温度重参数化”，转而记录其signed functional，并只在能给原U符号或尺度预测时继续。
3. 若正式定义排除该事件，给exact producer补上述最小逐记录字段并复测一次；不得同时改变几何、N、source或投影来救rank one。
4. 结果前不启动新MC、N、距离、动量、三插入、descriptor或五/六臂概率工作。#539只是P2复现支持；PR532的两桥因式分解与3/2下界不形成生产授权。

渐近验收仍是完整 `T_N=jY_p-R*jM_p-R_p*jM=o(D/A_N)`，等价于root-conditioned mixed Hessian `partial_u partial_epsilon Yhat=T_N/D=J_N/A_N`。朴素三位置绝对计数需要`alpha4>4/3`而严格输入只有某个`alpha4>1`；有限反例说明signed cancellation不能在冻结事件前被当作既有lemma。

## 2. 保留的理论缺口：固定m的真实两相相对权重

[固定m审查](../notes/p337-fixed-m-relative-bound.md)已经给表面界、sector-odds不足反例，并进一步证明：裸组件气体的标准非负KP判据在h=1也无法对大体积统一成立，任意非负控制函数都不能补救。停止继续优化这套裸表示的短轮廓计数常数。

PR533 `5aa929a6` 的 relaxed Catalan/Toeplitz identity、立即回边禁制、nonlocal one-west反例与条件`pq`双零可保留为C1子结果。真实beta cloud没有构造物理双gap核，rank1共同坐标、`w>=2`、second thermal/root/original-U和uniform remainder仍未闭合；关键依赖也不在祖先链。因此该Draft保持P2/C0 overall，Issue #542保持开放，不把joint-limit叙述改写成fixed-`m`进展。

rank2投影逐配置等于固定唯一绕行组件颜色；若真正的相内簇尾已可求和，等面积torus的小簇贡献可逐项相消，几何差可达 `O(N exp(-c ell))`。但实际共存窗口内的内外相受限partition比和大轮廓尾尚未控制。下一理论交付只应补这个模型特有的归一化控制；仅重复共同pressure、rank1小、正性或已有Poisson联合极限不足以完成它。

该问题与canonical局部pair空间响应是不同的指定作用量问题，不互相借数值或边界。没有把固定m原U定理登记为完成。

## 已完成，不再立项

- [齐次N50](../experiments/p337-homogeneous-n50-20260831/RESULT.md)：完整父图epsilon=1/t0，U=1.0615603877、V_S=+0.0543457827，排除有限零传递，正号预测存活；合同结束，不自动加N100/t/epsilon。
- Xi、jump/reweight分解、[共同温度加同一S的四profile闭合拒绝](../experiments/p337-two-coupling-closure-20260831/RESULT.md)、全孔面核及非循环开关积分均已完成。不要重做，不用新坐标救回失败闭合。
- [全m>=64反号](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)、N/m²有界Poisson联合极限、Q4和Q1有限seam传递均已完成；它们不证明固定m大N结论。
- [普通估计器预算](../experiments/p337-estimator-access-20260831/RESULT.md)：即使给真实root/均值/分母，m64的star独立iid热协方差平均达到SNR3仍需每几何至少约1.52e25样本。[twist恒等式](../notes/p337-twist-estimator-access.md)没有自动解决抵消；当前不把十台机器投入这类采样。
- [协方差零空间QA](../experiments/p543-covariance-nullspace-audit-20260901/REPORT.md)：三个静默截断实现已统一修复，16个归档向量分数已用既存充分统计回溯；15项不变，P50仅改解释边界。#543不形成新生产或新证据线。

## 已完成并停线的生产

| 主实验 | 冻结结果 | 执行动作 |
|---|---|---|
| #154 temporal transmission | [165M新路径](../experiments/p154-prospective-transmission-20260831/REPORT.md)：两N净U在±0.50内，各entry/completion在±0.30内 | 固定lag=1簇源退出当前主要H4解释；不换lag、不补样、不称精确零 |
| #334 independent intervention | [新群体得分](../experiments/p334-prospective-intervention-20260831/results/latest.json)：同时排除残余投影接近0和旧残余点预测按±25%迁移 | 两失败预测封存；不追加prefix，不在验证块重拟合，不把约1/2注册成救场模型 |

#154/#334的一般问题保留P1；#537是唯一P0理论任务，当前P0随机生产为空。F4等原停线不变；#275/#419/#370/#398、#539和代数全族筛选保持support，不通过增加archive坐标、descriptor或generic certificate恢复优先级。

## 原U接口直接复用

这一步已经推导完成，不再立项。令 `q=-1+F1+F2`、`E=1-F1+F2`，源为指定的a、`Jq=Cov(q,a), JE=Cov(E,a)`，在共同根 `mean_g q=0` 处记 `D=mean_g q_p`、`A=N^(13/8)/2`：

```text
p0dot = -mean_g(Jq) / D
Ddot  = mean_g(Jq_p + q_pp*p0dot)
U     = A*P4(E_p) / D
Udot  = A*P4(JE_p + E_pp*p0dot)/D - U*Ddot/D
```

P4是两取向差除冻结的DeltaCos4。新双插入Q激活代入a_xy，保留两几何各自的中心化、根移动和分母。单点a_x与a_y的协方差不能替换它。两个局部lambda独立；共同epsilon/N下，无序对对 `∂logQ∂epsilon²U` 的贡献是 `2W[a_xy]/N²`，裸Q1的epsilon响应仍为零。

十台独立服务器已获按需使用授权；不为填满机器而开新生产。实用计算任务先检查所选机器当前进程，不覆盖不明任务。通过仓库交接，减少跨团队消息；本轮理论与有限代数没有云任务。

结果在Draft #509交付；维护PR528只同步导航。**不合并，不删除历史数据、冻结合同或分支。** 原Issue清理见[记录](REPOSITORY-TRIAGE-20260831.md)。
