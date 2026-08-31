# 两条主线已给出独立判决：先兑现失败规则，再选择新机制

2026-08-31。接收 #509 `f4999e29` 和执行分支 `612df8ec` 后更新。
[评分前交接](independent-decisions-handoff-20260831.md)保留为历史；其中
“#154待交付”已被本记录替代。以下是现有正式交付的综合解释，没有重评分。

## 总判断

用户提出的“不要把新增分析当成减少机制不确定性”已经转成实际决策。
两条主线的三项独立实验都完成了：#154排除两个强传输数值模板，#334排除
完整两-score均值闭合及两项固定contact残差预报。现在兑现这些结果：
该lag1源退出主要H4解释的默认优先投入；三个已完成实验均不自动续采。
这是具体模型的注意力调整，不关闭、锁定Issue，也不否定整条理论路线。

还不能把“剩下一个数值模板”写成“唯一物理机制”。W/B/C是非穷尽的有限尺度
响应区域；#334两个实验更是在不同源/投影上作判决。强行凑成三套互斥完整
理论，会掩盖这个区别。

## 1. #154：165M全新排列已完成，弱传输停止规则成立

正式交付：[`f4999e29612da16a3650f24d124fb59137f053d7:REPORT.md`](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md)，
同目录 `PROSPECTIVE_RESULT.json` 保存完整32×32、主6×6协方差及每N200个
配对删批向量。冻结 `0820b8d2`，N85=5M、N340=160M，九分片全部取回。

| N | entry导数及同时区间 | completion导数及同时区间 | net导数及同时区间 |
|---:|---|---|---|
|85|.044492 [−.043502,.132487]|−.001022 [−.076653,.074609]|.043470 [−.071640,.158581]|
|340|−.053571 [−.220264,.113122]|.114247 [−.048797,.277290]|.060675 [−.157394,.278745]|

六个区间使用冻结的Bonferroni渐近正态95%家族，z=2.6382572735。
四个通道区间完全进入W的±.30带；B的entry≥.60与C的completion≥.60
分别在两个N被排除。两个net区间完全进入预设±.50带，因而原停止规则成立。
W记为not_excluded，不由此识别物理理论或严格零。

源是K−1的early-rank内中心化bulk CB+CW，再加入一个均匀空位的lag1条件政策。
entry/completion是原U在同一完整源下的**读出分解**，共享rootdot和Ddot，
并非分别开关两个事件源后的因果归因。旧样本不进入点估计或协方差。
结论只约束这个源、这两个N和冻结分辨率；并非所有source、所有尺度或H4机制。

## 2. 次级时钟线也已完成，没有“存活时钟”救场

[`612df8ec1cbe3be3938ee2e1f6183a1aefc6510b:secondary result`](https://github.com/LightChainr/Matching-One/blob/612df8ec1cbe3be3938ee2e1f6183a1aefc6510b/notes/p154-clock-line-secondary-result.md)
给出M10/M11的四个线残差；全部同时区间包含零，两条限制均not_excluded。
两条线本来都经过零，弱响应不能选择其中之一。不能把它们用来覆盖主停止规则，
也不拟合一个可填满二维读出平面的任意两幅度混合。

这项次级规则在 `83f3eba8`、实现 `4129d47d` 于读取新结果前提交，是生产期间/
之后的outcome-blind次级注册，不是另一套生产前主实验。它复用正式源响应，
只按同批删除重算未标记基线gain；共同10×10协方差保留与主六坐标的依赖。
大gain说明潜在灵敏度，不说明实际源占据该方向；旧cos4相对时钟律也不是M11。

## 3. #334：两个独立结果相容，但不能合成一票

| 实验 | 已交付判决 | 仅约束的模型 |
|---|---|---|
|1M新prefix的source-normal干预，`1164ba91/d0a9daf1`|T=(3.08520±.39187)×10⁻⁸，3SE下界>10⁻⁸；正预报3.6565×10⁻⁸未排除|完整条件label均值由两个原score、跨class共用斜率表示的闭合被排除；不是所有第一Jacobian近似失败|
|600k新prefix的四-contact残差预报，`14b2c98e`|R_new/R_old区间为[.4361,.5617]/[.4507,.5831]；零附近与旧幅度附近两个冻结带均失败|条件于旧固定训练值的loading预报；不是新“精确一半”定律或全部四-feature响应闭合检验|

来源：[normal报告](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md)、
[contact报告](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md)。

新的[精确投影crosswalk `ccabada3`](https://github.com/LightChainr/Matching-One/blob/ccabada318b1eeb12ae28d53391b13ab44c116d2/notes/p334-independent-interventions-crosswalk.md)
写成 `m^c=Hᵀb+r`、`E[Hr|Z]=0`。第一响应看 `tau=Gb`，normal源只看
`E[phi*r|Z]`；四-contact的R又是跨prefix的固定预测残差loading投影。
它们没有换算恒等式。把两种结果都叫“那个20% residual”会造成错误的机制合并。

## 下一步：不给失败参数化续命，也不拿空闲服务器填满队列

当前三个主实验与次级解释都已完成，没有待补的首次生产、正式评分或新独立票。
下一项值得投入的科学工作是一个具名微观模型对**原global U**的定量传输律，
而不是对本次弱响应或contact约一半再拟合。

可用现成的 `source → birth-clock first jet → U` 精确图作接口：模型必须给出
它实际激发的contrast值/导数，再落到既定几何上的响应预报。只列gain、把
`tau=Gb`当作拟合恒等式、或证明局部response非零，都没有提供这项预报。
这是下一主实验的科学内容要求，不是新的通用工具工程或任务许可。

#275/#370的现成投影/证书优先服务这个具名模型；不扩充generic certificate目录。
complex-C3、Gaussian/annulus、三角跨微观不变量与connectivity radical仍可并行
探索，但不冒称对本次新数据的预注册确认。本轮不新增第四模型、不换lag、不续采。
将来有独立理论导出的新实验可以重新排优先级，失败路线本身没有被永久禁止。

## 交付与计算边界

九片生产墙钟约590.215秒；主评分内部0.644秒，调用封装1.365秒，正式调用一次。
JSON SHA256为 `2c02bf3214ba4c9b31e8ad7ae65addf6d1ff0a82882518d1151556e7d0ce6821`。
同提交 `DELIVERY_VERIFICATION.json` 和 `SCORE_RUN.json` 为来源；本总览没有复跑。
来源 `CLOUD_COMPLETION.json` 在20:02:29 CST记录生产团队五机已返回Ready、现场
十机均Ready。这是有时间戳的交付记录，不是本总览的新实时查询或关机操作。

本轮按更新Huawei Skill只读过指定P154目录的运行元数据；后续直接消费已提交包。
未启动新生产、改密钥、停止其他作业或操作Issue生命周期。结果在 #509/执行分支，
这里只把精确来源和判决接入Draft #267，未宣称已集成main。表格用于直接比较同一
冻结读出和不同实验投影；不另做仪表板或扩大验证任务。
