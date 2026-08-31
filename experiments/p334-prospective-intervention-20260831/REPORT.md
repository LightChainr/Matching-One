# P334：独立 prefix 干预排除了两个冻结的残差预测

**剩余 clock-loading 投影在新群体中仍为正，但只有旧固定点预测的约一半。两项预定候选都被排除。** 本实验使用全新的300,000个prefix/尺寸，先保存contact预测，再生成独立baseline与正概率干预的条件尾部。未改模型、descriptor、预算或判据；“约一半”只描述本次结果，不建立新的比例模型。

## 唯一预定主读数

每N取两个receiver的配对平均。R是完整人口加权的单一残差投影：

`R = pi00 * [2 Cov_00(mu_C,tau_C-rhat_C) - 0.5 Cov_00(mu_W,tau_W-rhat_W)]`。

这里的预测器、四个contact变量、训练均值和旧R点值全部在F0冻结。主推断条件于这些固定值，旧数据不进入验证得分。

| N | R_new（±一个新60-batch SE） | R_new/R_old | 预定97.5%区间 |
| --- | ---: | ---: | ---: |
| 325 | (2.5359933 ± 0.1388471)e−9 | 0.4988857 | [0.4360616, 0.5617098] |
| 425 | (2.2061875 ± 0.1228952)e−9 | 0.5169035 | [0.4506760, 0.5831311] |

冻结的C0容差是比值[−0.25,0.25]，C1容差是[0.75,1.25]。两N的主区间均与两条带不相交，因此按预定规则分别排除“这个残差投影消失”和“旧点残差在±25%内传递”。每N使用t59的97.5%区间，两项主读数采用Bonferroni共同95%覆盖。没有以某个辅助方向替换主读数。

四个receiver结果分别为：N325 first `(1.9991762±0.1673358)e−9`、second `(3.0728104±0.2149378)e−9`；N425 first `(1.1693240±0.1688134)e−9`、second `(3.2430510±0.1580975)e−9`。同N的相关性通过共同60批factor保留，未将四行作独立合票。

## 正干预与物理读数

新实验保持原有两取向几何、k0、完整prefix采样协议。00中使用既有exact-census源 `H_s=pi_a(L_s-mean_a L_s)`，joint-degree class内精确中心化。正政策是 `q±=(1±H_s/8)/d`，每类质量不变、概率严格为正。同class的两个arm逐对保持k、两rank、两degree及两图Euler值，改变loop/contact组成。

通过同class maximal-coupling消去共同零差分量，实际生成正负残余分量的成对尾部；这是Rao–Blackwell化正政策contrast，不把所有尾部说成来自完整q±的直接iid抽样。单位响应 `tau=Delta/(2epsilon)` 与旧H响应**严格相同**，没有有限指数源的Taylor截断。以下都是tau；实际q+对q−的有限效应为表中数值乘1/4。

预定辅助own-source读数（完整新人口分母；误差和所有cross-source方向在JSON保留）：

| N / receiver | tau A(p_ref) | tau E(p_ref) | tau K1/(N+1) | tau K2/(N+1) |
| --- | ---: | ---: | ---: | ---: |
| 325 first | −3.057579e−5 | 7.034975e−6 | 2.451343e−6 | 3.523728e−6 |
| 325 second | −3.090389e−5 | 7.021971e−6 | 2.466638e−6 | 3.527064e−6 |
| 425 first | −3.226028e−5 | 7.865001e−6 | 2.328744e−6 | 3.330296e−6 |
| 425 second | −3.245795e−5 | 7.767146e−6 | 2.323796e−6 | 3.325292e−6 |

方向是出生中心延后、固定p的A下降。粗状态相同仍可出现明确未来响应。固定四contact模型只有C/W及其线性变换的数值预测；**没有据此补造A(p_ref)/E(p_ref)预测系数**。这些辅助读数是同一批birth paths的相关函数，不增加独立实验数。

## 新样本、预算与核验

- 每N完整300,000个新prefix，60批×5,000，无按00数量补抽；00实际为N325 **22,037**、N425 **23,687**，其他八cell保留并对该干预响应置零。
- 每00 prefix固定32条独立uniform baseline尾部；每个物理source固定64个正负contrast。总计 **1,463,168条baseline + 11,704,832条contrast arm = 13,168,000条共同尾部排列**，每条同时计算两取向。N325两个source各有2个零质量情形，按冻结规则响应精确0、无需尾部，没有另抽prefix；N425无零质量source。
- 新prefix counters：N325 `[53032500000,53032800000)`；N425 `[53042500000,53042800000)`。条件流使用冻结bit63地址域。命名域分离避免复用原stream；不将PRNG种子差异宣称为数学独立性证明。
- 原冻结18个文件hash始终匹配；120个shard及全部原始文件取回并逐文件SHA256核验，预测文件均早于本shard的tail阶段。预算与同class/rank保持检查全部通过。
- 导出的旧β/均值重构旧捕获与剩余loading的最大误差为3.31e−24；旧2,048个prefix重构、126个旧00完整census与原描述符核对通过。生产前仿射政策和RB恒等式做了独立有理数检查。

## 冻结时序与资源

F0为 `4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277`，提交时间 **2026-08-31 19:36:13 CST**。首个prefix作业于19:38:29.630开始；首份预测19:38:30.575封存；最后一份预测19:39:13.645封存；全部尾部于19:39:19.839前完成。每shard的精确事件时间与预测hash保存在原始run receipt中。

HZsCM6运行0–29 shards/尺寸，55.729秒；TV2N0X运行30–59，55.033秒；两者均exit0、14workers、BLAS1。两机实际cgroup限额为14.5CPU/25GiB；采样到的容器总内存峰值约1.180GB/1.174GB（含系统服务），单child maxRSS约108,428/107,800KiB，不冒称并发进程RSS总和。编译器GCC10.3.1，云端Python3.9.9/NumPy1.26.4/SciPy1.13.1。

全部120shards回收后，冻结的最终scorer于19:43:02运行**一次**，0.737秒、exit0。机器上的本实验科学进程均已结束；按主任务随后安排，**HZ/TV保持Running并保留既有隧道，交给P154复用**，未在本包中谎报Ready。此生命周期交接不改变实验预算。P154完成后，主任务已将本轮原账号五台全部恢复Ready并关闭自有隧道，见[最终云状态](../p154-prospective-transmission-20260831/CLOUD_COMPLETION.json)。

## 可以据此说什么

独立prefix群体和正干预表明：固定四项线性预测之外的这一个clock-associated响应投影仍存在；旧残差的固定数值没有通过预设精度的前瞻检验。这限制了可传递预测，尚不确定失败来自旧训练估计误差、固定线性形式或更具体的contact动力学。R接近零也可能由抵消造成，所以本实验不将任一结果解释成完整response充分状态、因果唯一机制或连续场数。主检验没有把旧训练不确定性加入验证得分，也没有执行无条件新旧总体差的检验。没有结果后拟合1/2、增添descriptor、追加样本或重开Hessian分析。

## 文件

- [CONTRACT.md](CONTRACT.md)：冻结方案；[FREEZE.json](FREEZE.json)：参数与原文件hash。
- [results/latest.json](results/latest.json)：主/辅助点值、60批LOO与共同factor；SHA256 `1a2f569fdf9ab7c3d1a8c886828807f562fa98e458c1db06ef90faba66aa3386`。
- [execution/INDEPENDENT_REVIEW.json](execution/INDEPENDENT_REVIEW.json)：独立复核源归一化、原始抽样、共同协方差和判定；可携带统计复算JSON逐字一致。
- [execution/RECOVERY.json](execution/RECOVERY.json)、[execution/SCORE_EXECUTION.json](execution/SCORE_EXECUTION.json)：两机运行/取回、唯一评分与时序。
- [results/sufficient_statistics/](results/sufficient_statistics/)：120个分片充分统计原样入Git，可在干净checkout用冻结评分器复算完整结果，无需取回236MB路径原始数据。
- [RAW_MANIFEST_SHA256.tsv](RAW_MANIFEST_SHA256.tsv)：全部1,320个原始/分片统计文件；原始文件保留在本地`production/`（按冻结.gitignore不自动入Git）。远端完整包及逐文件checksum也均已保留。
- [DELIVERY_MANIFEST_SHA256.tsv](DELIVERY_MANIFEST_SHA256.tsv)：Git交付文件清单，另存而不改F0 manifest；不包含另由RAW_MANIFEST记录的production路径文件。

复算命令（输出目录需不存在；只重复确定性评分，不产生新随机样本）：

```bash
python score_prospective.py --input results/sufficient_statistics --output /tmp/p334-recomputed-result
```

该命令在本实验目录中运行，需要NumPy/SciPy；结果应与已保存JSON逐字一致。原始路径仍完整保留，压缩的充分统计不替代原始数据的来源清单。
