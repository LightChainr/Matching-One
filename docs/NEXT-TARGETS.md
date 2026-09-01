# 下一步：thermal/pivotal 传递与固定m两相权重

**2026-09-01，已读到执行 `410015f5`、临界 raw-kernel 定理与 thermal/pivotal 双通道审计。** 当前结论从[STATUS](STATUS.md)进入。两个独立生产实验、齐次N50、canonical pair 的有限空间块和完整联合原U响应均已完成。下一步不增加同类有限点、距离网格、补全扫描或旧实验救值。

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

现在只保留 #537 一个P0理论闸门：**为两个通道写出真实三位置support，并把完整 `T_N=jY_p-R*jM_p-R_p*jM` 控制到 `o(D/A_N)`。** 朴素绝对三臂组计数需要`alpha4>4/3`，严格方格site输入只有某个`alpha4>1`，不能复用raw pair的同一论证。验收和停止规则为：

继续使用已经闭合的 root-conditioned 坐标。若 `M(p(u,epsilon),epsilon)=u`、`Yhat(u,epsilon)=Y(p(u,epsilon),epsilon)`，则有限体积精确恒等式为 `partial_u partial_epsilon Yhat=T_N/D=J_N/A_N`。有限landing的第一项可证伪候选已经判决：N9普通四臂、无额外分支的C4轨道在实际finite matching root上给出精确非零`2 x 2` minor；其minor多项式与matching-root多项式在`Q[p]`中互素，`E -> E-Rq`的root Schur剪切保持行列式不变。N16三纤维证书又给出参数无关的`det=-chi/2`。因此“三组四臂leading landing纯属温度坐标”的结构性lemma为假，不能把余项自动提升成`R^4*pi4(R)^4`四包络。

执行顺序随之改变：停止五/六臂替代证明与扩大N25枚举，直接把两个必需通道组合成已具名的signed ordinary-landing functional——kernel reconnection与rank/readout pivotal，经同一C4和root Schur投影后的空间和。下一交付只应给这个functional的符号抵消规律、尺度界或exact-`p_c`到pooled-root运输；若拿不到足以达到`T_N=o(D/A_N)`的速率，就交付最小未控landing类并停止。有限minor已经完成其停线任务，不再重复扫描更多tiny torus。

1. 若完整两通道和root/slope运输给 `T_N=o(D/A_N)`，停止把这一local interaction作为原anomaly的尺度放大机制，转为有限局部修正资产。
2. 若只有一个具名pivotal/landing通道逃过界，必须先给同一原U的符号或尺度预测，才冻结一次新读数；不能从结果中再选bridge irrep。
3. 若现有严格臂界不足，交付最小未控joint event或signed cancellation条件即停止；不能用未经证明的5/4替代，也不能把上界发散写成实际发散。
4. 在该闸门前不启动新MC、N、距离、动量、三插入、alpha/positivity completion或descriptor。#539只是P2可复现支持；重建旧N25 `J2`不算关闭渐近闸门。PR532已提交的两桥因式分解与3/2下界保留，但不形成生产授权。

## 2. 保留的理论缺口：固定m的真实两相相对权重

[固定m审查](../notes/p337-fixed-m-relative-bound.md)已经给表面界、sector-odds不足反例，并进一步证明：裸组件气体的标准非负KP判据在h=1也无法对大体积统一成立，任意非负控制函数都不能补救。停止继续优化这套裸表示的短轮廓计数常数。

PR533 `a7680426` 的 Bessel determinant `I0(2c)^2-I1(2c)^2>0` 可保留为受限 directed/noncrossing 的 C1 子引理；它的 full-lattice moving-root 负号仍依赖未证的 uniform connected-polymer/root-shift 界，且只覆盖 bounded `L/m`。因此该 Draft 保持P2/C0 overall，不把它的 `m→∞` joint limit改写成fixed-`m`进展，也不启动已被其新head超越的#534。

rank2投影逐配置等于固定唯一绕行组件颜色；若真正的相内簇尾已可求和，等面积torus的小簇贡献可逐项相消，几何差可达 `O(N exp(-c ell))`。但实际共存窗口内的内外相受限partition比和大轮廓尾尚未控制。下一理论交付只应补这个模型特有的归一化控制；仅重复共同pressure、rank1小、正性或已有Poisson联合极限不足以完成它。

该问题与canonical局部pair空间响应是不同的指定作用量问题，不互相借数值或边界。没有把固定m原U定理登记为完成。

## 已完成，不再立项

- [齐次N50](../experiments/p337-homogeneous-n50-20260831/RESULT.md)：完整父图epsilon=1/t0，U=1.0615603877、V_S=+0.0543457827，排除有限零传递，正号预测存活；合同结束，不自动加N100/t/epsilon。
- Xi、jump/reweight分解、[共同温度加同一S的四profile闭合拒绝](../experiments/p337-two-coupling-closure-20260831/RESULT.md)、全孔面核及非循环开关积分均已完成。不要重做，不用新坐标救回失败闭合。
- [全m>=64反号](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md)、N/m²有界Poisson联合极限、Q4和Q1有限seam传递均已完成；它们不证明固定m大N结论。
- [普通估计器预算](../experiments/p337-estimator-access-20260831/RESULT.md)：即使给真实root/均值/分母，m64的star独立iid热协方差平均达到SNR3仍需每几何至少约1.52e25样本。[twist恒等式](../notes/p337-twist-estimator-access.md)没有自动解决抵消；当前不把十台机器投入这类采样。

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
