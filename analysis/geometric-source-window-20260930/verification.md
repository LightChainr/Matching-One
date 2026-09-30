# 固定窗口几何源：旁路核验

2026-09-30，独立支持worker。结论：**所审版本未发现实质公式或实现错位**。
源归一化、unordered next-h集合条件平均、原排列延续、三个条件risk均值及整批协方差实现同一固定合同。
这是源码审查加一个指定构型的精确分数对照，不是生产数据、云端构建或执行结果的认证。

审查时HEAD为`d9824b5bf478fc8acad238be2b9f31ce99787966`；下述五文件相对
`5929403a002c6e956cb74f1b6bb8d827bcc84109`无变化。只写本文件，无commit/push。

## 发现与对应

| 核对项 | 结论与依据 |
|---|---|
| 固定合同 | contract为square NN occupied-site torus，L=512、N=262144、b=155385、h=735、终点156120。run.py:51–52的实参及score.py:13的M=106759、H=735一致。输入均为safe-insertion-production旧14批，未混入独立completion-hazard复核块。 |
| kick归一化 | 完成点总质量c/m；安全点总质量(s/m)×Z/Z=s/m，合计1。θ=0每空点概率1/m；安全点P′₀(v)=(d_v−2e/s)/m，完成点导数0。因此likelihood score为d_v−2e/s，已包含正确归一化。最终均匀源期望不能再额外乘s/m。 |
| unordered集合平均 | note:15–36正确。均匀h-subset B的概率为1/C(m,h)，kick后为Σ[v∈B]Pθ(v)/C(m−1,h−1)，比值为(m/h)ΣPθ(v)。存活B全部由初始安全点构成，导数score恰为Σd_v/h−2e/s；不存活B乘survival后为0。度数均取初始A，既不沿途重算，也不把两步曲率外推至735步。 |
| e、degree及平方和 | engine.cpp:23–27从去重的无序协同边给两端各加1；故Σ_safe d=2e。degree每个前缀清零，sum_degree_squared为Σd²；乘法先转uint64_t。score.py:39先用Python整数计算sΣd²−4e²，随后才除法，没有先以浮点做这两个大数的相减。 |
| 原排列及真实延续 | engine.cpp:18–21与旧safe-insertion engine.cpp:68–76使用同一完整Fisher–Yates、mt19937_64 seed和bounded算法。旧probe_rng独立；旧probe只读取const order并修改Geometry副本，不消耗排列rng、不改order。新counter也不消耗rng、不修改g；故order[b:b+h]恢复同一完整排列的原后缀。 |
| 终点存活及score | engine.cpp:28–33累加恰好h个后缀站点的**初始**degree，并对原g继续逐点插入。达到rank2后跳过后续Geometry更新不会改变终点survival，因为rank单调；window_sum仍累加全部h项，随后由survival=0消去。score.py:41实现指定条件平均，不是非零θ有限差分。 |
| 三个条件risk均值 | score.py:30–50只排除rank_b≠1，所有rank1前缀均进入同一risk分母；未按终点存活或s>0再筛选。s=0时三项均为0，但该前缀仍计入risk，符合条件于原rank_b=1的目标。总体点值为各批分子和/各批risk和，不是等权平均各批条件均值。 |
| 两步伴随量 | −(sΣd²−4e²)/(s m(m−1))=−s Var_safe(d)/(m(m−1))。它与有限窗口响应及baseline_survival按同一批、同一rank1风险人口联合保存；不是额外独立证据。 |
| batch covariance | score.py:47–59每次删除整批后重新加总分子和risk分母，得到14个三维比率向量。以删除向量均值为中心、系数13/14构造完整3×3协方差，SE取对角平方根；没有丢掉跨指标协方差，也没有误把14批当140000个独立删除单位。 |

实现边界：通用note允许h=1；当前engine要求h≥2，这不影响固定h=735合同。
run.py和score.py的参数是匹配当前合同的硬编码，本核验不自动适用于未来修改后的参数。
run.py的children_peak_rss还可能包含编译子进程峰值，因而可能保守限制worker；不改变估计量。
这些均不构成当前科学读出的阻碍，不要求增加测试或重跑生产。

## 唯一实际计算对照

使用本地研究解释器`/Users/lc/python-envs/research-py311/bin/python -B -`，代码从stdin运行；
只用标准库`fractions.Fraction`、`itertools.combinations`、`math.comb`，不写临时文件。
通过独立的格点提升坐标/回路绕数求rank，没有调用旧derive.py、旧census或生产scorer。

固定square L4，v=x+4y，唯一输入为既定构型
`A={0,1,4,5,8,12}`，**只核对h=3**：

- m=10、c=0、s=10；协同边为(2,3)、(6,7)，e=2。
- 四个degree=1、其余六个degree=0；Σd²=4、mean_d=2/5；安全score总和精确为0。
- 只计算120=C(10,3)个终点集合，存活概率精确为5/6。
- 直接kick导数：对10个首点及各自36个后续二元素集合，按
  `P′₀(v)=(d_v−2/5)/10`与后缀概率`1/C(9,2)`累加360项，得到**−7/180**。
- 条件平均读出：对相同120个终点集合累加
  `survival×(Σd_v/3−2/5)/C(10,3)`，得到**−7/180**。
- 两结果以Fraction精确相等，且与已写出的A、h=3数值相符；没有数值差分或容差判定。

这个单构型对照核实归一化单kick与集合平均公式的连接；它没有编译/执行本次C++，
也不声称已验证L512生产行、735步结果或云运行状态。当前源码对应性由上表的逐段审查提供。

## 实际操作范围

已读指定五文件；为确认源与依赖语义，另读源kick定义、旧safe-insertion engine及被include的pair/geometry接口。
做了上述唯一精确小对照，并只读检查版本差异及下面的文件hash。
无SSH、无云端操作、无旧档案数据重放/评分、无旧L4 census、无三构型全滞后重跑、无扩样或cutoff搜索。
本核验未发现需立即中断主线的实质问题；生产执行及结果核对仍由主线程负责。

审查文件SHA256（区分主线程后续修改）：

```text
notes/geometric-source-window-identity-20260930.md
6b0c6f7e0b2b69aa106bc1256a697fc805907868a8fb3349a7cc4606bddf21be
analysis/geometric-source-window-20260930/engine.cpp
1269a3fda3d7a9dcd302fa872ad4437488aeb4af5d5ffd254f2b3f31e754d0b1
analysis/geometric-source-window-20260930/score.py
2c8fb9f5a6c31b8af74d497a04d74e78ea2e5dabe7fc7e2789cd3f7279d1769e
analysis/geometric-source-window-20260930/run.py
314d4929d8ff45c6e74464cad5a65061510fb331f4fb6521b3000e9d94f7d0f0
analysis/geometric-source-window-20260930/contract.json
8f044bc068e9edc09b8a1b56acb0d0cd399431b42d14c4d3c60cc77eab5892eb
```
