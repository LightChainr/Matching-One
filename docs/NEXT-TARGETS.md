# 下一步：固定补全的多组件空间传递

**2026-09-01。** 当前结论从[STATUS](STATUS.md)进入。两个独立生产实验已经按冻结规则停线；齐次N50和局部pair的有限原U传递也已完成。下一步不再增加同类有限点或重开旧实验。

## 1. 当前研究问题：多组件空间结构如何进入同一个原U

执行[2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)已经定义并算完canonical `Kreg=K2+K0`：直接Q1响应恒零，单点混合原U响应严格非零，固定四路径双插入系数13/8。旧“local pair还没有进入有限U”的待办删除。具体cut/连续场归属仍不能从这些有限结果直接推出。

**[本次空间结果](../experiments/p337-regular-spatial-support-20260901/RESULT.md)进一步排除了单一共享组件传播。** 首阶Q激活必须由至少两个外部组件连接两个非相邻标记；canonical核给 `Cxy=E[a_xy]` 和 `|Cxy|≤(43/16)Pr{两点vacant且至少两不同占据组件接触两处}`。实际两组件见证表明阈值可达到，但不保证占据平均非零。4140项已全部计算并独立复核，下一次不再扩这一核目录。

下一项真正需要的结论是：**同一固定核的带符号多组件空间概率，在距离增长时保留哪一项，以及它通过W[a_xy]是否保留在原U。** 先使用已给出的概率支撑和精确系数约束符号抵消；不能把普通连通概率、Cov(a_x,a_y)或单一见证1/16当成答案。若某个近似只允许一个共享组件，它对本响应给零，已应停用；不能事后添加自由系数救回。

只有在理论给出可互相区分的空间/原U预测后，才冻结一个新的生产问题。未形成这种预测前，不启动距离网格、counterterm扫描或更大N；也不把“局部相关非零”升级为P0。允许的正则补全共享上述选择规则，但43/16和单点W只属于canonical补全，不能换系数保持原结论。

## 2. 保留的理论缺口：固定m的真实两相相对权重

[固定m审查](../notes/p337-fixed-m-relative-bound.md)已经给表面界、sector-odds不足反例，并进一步证明：裸组件气体的标准非负KP判据在h=1也无法对大体积统一成立，任意非负控制函数都不能补救。停止继续优化这套裸表示的短轮廓计数常数。

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

#154/#334的一般问题保留P1，当前P0生产为空。F4等原停线不变；#275/#419/#370/#398和代数全族筛选保持support，不通过增加archive坐标、descriptor或generic certificate恢复优先级。

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
