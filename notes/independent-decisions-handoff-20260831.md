# Independent experiments now remove specific mechanism predictions

历史记录：本页保留 #154 正式评分前的交接状态；完整生产及主/次级判决现见
[最终交接](independent-decisions-final-20260831.md)。下文“待交付”不再是当前任务。

## 两项 #334 决策已完成，#154 承接固定新块

项目已经走到独立数据上的模型淘汰：**完整两-score 条件标签均值闭合失败；
四-contact 残差的“接近零”和“旧固定幅度原样迁移”预报也均失败。**
两项 #334 实验的源与主量不同，结论不矛盾，不能合成一个“统一 residual”。
#154 已有冻结的新165M排列实验；下一交付是它的原始主判读和预先声明的
次级时钟线比较，不再增加描述量或新生产。

本页把主判决、定义、边界和下一动作放在一起。它引用未合并分支，不把
分支结论升格为 main，也不关闭或锁定 Issue。当前唯一队列仍是
[NEXT-TARGETS](../docs/NEXT-TARGETS.md)。

## Score-normal 干预排除完整两-score 均值闭合

执行队的 [score，`1164ba91`](https://github.com/LightChainr/Matching-One/blob/1164ba91035bb0ee37bcea52f700312475c257b5/results/p334-independent-normal-intervention/score.json)
给出固定四个 own-source 出生中心响应的等权平均：

`T=(3.0852005663 ± 0.3918738407)×10⁻⁸`，
`T±3SE=[1.9095790443,4.2608220883]×10⁻⁸`。

区间下端超过冻结的 `delta=10⁻⁸`，触发
`stop_complete_two_score_label_closure`。旧点预报 `3.6565×10⁻⁸`
落在该区间内；这项正响应的迁移未被排除，不等于确认唯一机制。

其源由二阶密度 score 去掉原两个 score 的分量得到 `phi`，满足
各 joint-safe class 内均值零，且 `E[phi*s_f]=E[phi*s_s]=0`。
有限两臂 `q±=(1±phi/max|phi|)/d` 直接抽样；主量中保留
`max|phi|*(C_plus−C_minus)/2` 的权重。它恰为 normal 响应，无Taylor外推。
这些概率非负，某些标签可为零。

被排除的是具体完整条件均值表示
`m_C(Z,u)=c_a(Z)+b_f(Z)s_f(Z,u)+b_s(Z)s_s(Z,u)`，同一prefix的
classes共用两斜率。第一Jacobian近似 `J=BG`、class-dependent斜率或
其他label结构仍可能成立。结果不是第三个连续场的证据，也不等于
未扰动global Matching异常的来源识别。

每 N325/N425 独立采样 `20×25000=500k` 新prefix，完整人口分母不变；
00数量为36938/38876，每 own-axis 固定8组配对。旧20批没有合并入新得分。
原始块先提交于
[`f1b36436`](https://github.com/LightChainr/Matching-One/tree/f1b36436acef4ce9935df23b22b1dc53f109bdfb/results/p334-independent-normal-intervention)，
再执行冻结评分。两台生产含编译耗时12.578/14.739秒，40批exit0。
[最终说明，`d0a9daf1`](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md)
记录取回和恢复Ready；此前dispatch中的“正式生成前”是不可变历史记录，
已被此完成交付接续，不应继续显示为当前状态。

## 四-contact 固定残差预报失败，残差并未消失

#509 的 [独立实验，`14b2c98e`](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md)
使用原一阶 loop-score `H_s`，严格正的 `q±=(1±H_s/8)/d`；它保持
class质量和即时拓扑，但没有要求同时保持原score均值。
其主量是固定四-contact 预测之后的clock-loading协方差投影：

`R=pi00*[2 Cov(mu_C,tau_C−rhat_C)−0.5 Cov(mu_W,tau_W−rhat_W)]`。

这不是上面的 `E[phi C]`。模型、旧训练均值及R_old点值已在
[`4b3c21b7`](https://github.com/LightChainr/Matching-One/blob/4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277/experiments/p334-prospective-intervention-20260831/CONTRACT.md)
冻结，新块不重拟合。

- N325：`R_new=(2.5359933±0.1388471)×10⁻⁹`；相对旧固定点值为
  `.4988857`，预定区间 `[.4360616,.5617098]`。
- N425：`R_new=(2.2061875±0.1228952)×10⁻⁹`；相对旧固定点值为
  `.5169035`，预定区间 `[.4506760,.5831311]`。

两区间均避开 C0 的 `[-.25,.25]` 和 C1 的 `[.75,1.25]`。
因此两条预报都退出本次比较。每N使用t59的97.5%渐近区间，两个主量
按Bonferroni共同95%；推断**条件于旧训练点值**，没有加入旧训练不确定性。
所以不能改写成“真实新旧总体残差证明减半”，也不事后注册一个1/2模型。

每N `60×5000=300k` 全新prefix，baseline与contrast尾部独立生成，
共120分片。它和执行队的新块分属不同随机域；四个新N-domain合计
1.6M fresh prefixes，但不同估计目标不合成一个效应或证据票数。
同N的两个receiver、所有源、A/E/K1/K2和尾部保持共同协方差。
固定预算已结束，未来传输问题不再自动触发本实验加样或第五descriptor。

## #154 的独立实验已经冻结，不是旧小效应预算的重跑

[正式冻结，`0820b8d2`](https://github.com/LightChainr/Matching-One/blob/0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f/experiments/p154-prospective-transmission-20260831/CONTRACT.json)
与后续 `14b2c98e` 的authorization明确 N85=5M、N340=160M新排列，每N200批、共九分片。
仍为 lag=1 的早期rank内中心化bulk CB+CW源；对每个最终K在K−1注入。
新块内重估条件均值和根，p微分只作用于Binomial权重，旧数据不合票。

其六个主坐标是两个N的entry、completion、net U源导数，共享完整源的
根移动和分母导数。三个有限数值模板为：双通道各在±.30；entry≥.60且
completion在±.30；反之。两N共同满足才可保持相应模板；净U等效带为
±.50。预定共同区间使用 `Phi^-1(1-.05/(2*6))`，不作实时显著性续采。
这些是非穷尽数值假说，不是三个已建立的物理理论。

上一轮 [scalar .01/.018预算](../results/p154-clock-transmission-budget/REPORT.md)
依然适用于那个小效应候选；它并不否定此处以较大响应区间为目标的新实验。
我们已淘汰的 [纯cos4相对时钟](../results/p154-fixed-clock-models/REPORT.md)
也不等同于下述两个纯contrast first-jet限制。

执行队的 [完整传输图，`c2828e34`](https://github.com/LightChainr/Matching-One/blob/c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847/notes/p154-birth-clock-transmission-map.md)
把四个有效birth clocks写为 `a00+eta*a10+tau*a01+eta*tau*a11`。
整个共同a00函数从原U中抵消，三个contrast值及其p导数可传输。
N340的a10名义gain约246.4，是a01的约915倍；**这是未标记基线的灵敏度，
不是已测源幅度或新预测效应**。因此小的a01方向不能代表全部temporal机制。

[次级冻结，`83f3eba8`](https://github.com/LightChainr/Matching-One/blob/83f3eba88d7f1290704f82610c28669dc5e12f3c/notes/p154-clock-line-secondary-freeze.md)
只比较预先给出的纯、局部平坦M10/M11增益线，在同一完整新块上采用
`(C_e*v_completion−C_c*v_entry)/hypot(C_e,C_c)`，不除以随机entry。
与源读数在同一delete-one中重估基线gain，四个残差共同区间；两条线均
含零，不排除不等于确认。它是相关的次级解释，不回写主矩形判决，
不新增第三个主实验或改变预算。

## 下一次交付的完成条件

1. #154 九个固定分片全部取回后，发布一次原合同主评分及完整六维协方差。
   不替换源、不挑N、不放宽矩形；混合、全部失败和未分辨照原名报告。
2. 承接已经预输出冻结的两条时钟线评分，不再建立重复scorer；与原六坐标
   保存共同协方差。两条都失败就记录失败，不在同块拟合替代线。
3. #334 两项实验标为完成、具体候选停止推广；进入global observer的
   传输仍是未解问题，不把它写成现成新生产任务。

## 交付与范围记录

本总览没有重采样、重跑科学测试或覆盖他人的结果；只读已提交科学交付及
指定P154目录的运行元数据。2026-08-31约11:58 UTC，Zy/XP原生产进程已
不在进程表，指定目录分别已有3/2个run回执；HZ/TV各1回执已取回。
Tg的一次检查遇到helper并发维护，未中断或重启任何进程；该情况不等于
作业失败。完整结果仍以执行方九分片正式交付为准。未阅读部分新P154源
结果、未做部分评分；各队继续通过仓库交接。

呈递采用技术报告的“判决→源/主量→不确定性→下一动作”结构。只有三个
不同单位的核心比较，采用数值段落而不画共用坐标图，以免误示同一估计量。
