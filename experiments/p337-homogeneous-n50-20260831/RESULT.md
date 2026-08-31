# 完整齐次 N50：原 U 的固定源传递严格为正

**已补齐原父图 epsilon=1、t=0 的完整有限计算。** 固定 `(5,5)/(1,7)` 两几何、原 q/E、同一 S* 和共同根，不使用孔尾近似或条件 B 总体。结果为

| 原读数 | 精确计算的十进制显示 |
|---|---:|
| 共同 matching 根 p | 0.59275940130675926587 |
| U | **1.0615603876876551** |
| V=d_t U | **+0.0543457826695583** |
| logit 热坐标的 pooled slope | 1.8356031693516655 |

严格判决使用有理区间，完整数据在[score.json](results/score.json)。特别是

```text
V/A50 >= 1885252764556849639086318759550302603 / 10^40 > 0,
A50 = 50^(13/8)/2.
```

因此冻结的有限零传递假设 V=0 被排除；从饱和端点正响应延续而来的
V>0 预测没有被否定。**这不是机制确认或连续 H4 身份。** 原饱和端点
V约+0.3891471785，不能将其幅度原样用于这里，也不从这两个数拟合尺度律。
本合同在此结束，不追加 N100、t/epsilon 网格、新源或MC。

## 实际传递中的抵消

预先固定的四项分解全部在 `z=logit(p)` 坐标、原 U 单位下报告：

| 项 | 对 V 的贡献 |
|---|---:|
| 直接 mixed 热/源响应 Y_zt/D | +2.0266262078 |
| 共同根移动 z_t Y_zz/D | +0.0569741811 |
| 斜率的源变化 −Y_z M_zt/D² | −2.0823887351 |
| 斜率的根移动 −Y_z z_t M_zz/D² | +0.0531341289 |
| **合计** | **+0.0543457827** |

直接项与斜率项大幅抵消，必须保留原归一化、共同根及热分母。
两个根移动项分别依赖热坐标，合计 V 不依赖 p/logit 的选择。
它们是同一个精确响应的分解，不是四份独立证据，也不从中选择新描述量。

## 完整总体与实际成本

采用已证明充分的黑 NN 连通性／同调商空间状态，并将 K 收入值多项式。
每个几何精确覆盖 **2^50=1,125,899,906,842,624** 个占据配置；每个 K 的
计数均为 binomial(50,K)。保存的是完整 `(K,q,count,sum_S)`，不是采样。
实现没有逐一访问2^50个叶子，而是严格合并具有相同未来响应的状态。

| 几何 | CPU秒（至计算回执） | 峰RSS MiB | 最大边界key数 |
|---|---:|---:|---:|
| (5,5) | 14.580789 | 927.515625 | 1,009,330 |
| (1,7) | 35.266000 | 1667.593750 | 2,009,158 |

两项在本机串行完成，合计约49.85 CPU秒；每项预算固定为120 CPU秒、
4096 MiB和500万边界key。逐层质量检查、CSV、参数、源码与二进制hash见
[first](results/first/table.json)／[second](results/second/table.json)。
上述CPU不含编译和早期探针，不能当作整个研究过程的总耗时。
无新随机路径、云作业或对旧 production 的补样；十台云机查询时均Ready。

最初双颜色Python表示及K在key的C++版本触及资源门，未出完整N50表。
最后只改K存储，保留其全部系数；[有界对照记录](feasibility/polynomial/REPORT.md)
说明如何从部分可计算推进到本次完整结果。旧探针为历史可计算性回执，
不再作为待执行实验；当前唯一生产入口是 `produce-polynomial.py`。

## 冻结和独立核查

- 科学合同／评分器：`10c666b65566b25ddb8eaa02219947a9c5a261f2`。
- 初版producer：`962adad6`；最终K多项式producer及完整预算：
  `4ae4e710404e40ac31eaf962680aaceed2543ddf`。科学合同未修改。
- N9/N13直接 lifted-homology oracle全表通过；两个N25完整表逐整数
  复现a70eeff0独立枚举。旧N25 U与V也按独立p导数及主logit矩公式复现。
- 两N50表各完成全部binomial计数与2^50总质量检查。
- [独立复核](REVIEW.md)使用直接Bernstein p一/二导数，Decimal120/160
  的12项读数全部落入主Fraction包络；两种坐标四项合计相同。
  Decimal交叉核查不替代严格有理证书。输入与冻结文件hash一致。

没有改变旧数据、m64合同或N25四点实验。N25和N50是各自完整有限图的
确定性结果，不是新独立随机证据，也不证明一般 epsilon／一般 N 的延续。

## 复算入口

在仓库根目录执行，输出须选择未存在的新目录：

```sh
python3 experiments/p337-homogeneous-n50-20260831/produce-polynomial.py --freeze-commit 4ae4e710404e40ac31eaf962680aaceed2543ddf --geometry first --output-dir /tmp/n50-first-new
python3 experiments/p337-homogeneous-n50-20260831/produce-polynomial.py --freeze-commit 4ae4e710404e40ac31eaf962680aaceed2543ddf --geometry second --output-dir /tmp/n50-second-new
python3 experiments/p337-homogeneous-n50-20260831/score.py --freeze-commit 10c666b65566b25ddb8eaa02219947a9c5a261f2 --first /tmp/n50-first-new/table.json --second /tmp/n50-second-new/table.json --output /tmp/n50-score-new.json
```

纯复核可直接用已保存的两张表及 `review.py --case n50-target`；无需重跑
producer。可计算性旧脚本会写各自探针文件，若要重现其历史过程，先复制到
新的临时目录，勿覆盖本包归档回执。
