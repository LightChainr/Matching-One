# P337 thermal/pivotal preflight producer

**状态：冻结的64-counter复放与描述性评分均已完成。** 本接口只复放既有
L32/L64随机流各自最前32个配置，总计64个。seed、counter数、32个原pair和
Bell8 kernel全部固定；命令行不能提供另一个seed或样本数。可复现输出见
[`REPORT.md`](../results/p337-thermal-pivotal-preflight/REPORT.md)、
[`latest.json`](../results/p337-thermal-pivotal-preflight/latest.json)与
[`run.json`](../results/p337-thermal-pivotal-preflight/run.json)。

冻结合同在
[`analysis/p337_thermal_pivotal_preflight_contract.json`](../analysis/p337_thermal_pivotal_preflight_contract.json)，
producer在
[`scripts/p337_thermal_pivotal_preflight.cpp`](../scripts/p337_thermal_pivotal_preflight.cpp)。
它是Issue #536交付B的语义与成本preflight，不是新生产、显著性分析或
L32/L64的第三份证据。

## 1. 它实际复放什么

输入严格沿用`a237968f`的两个流：

| L | r | seed | counter |
|---:|---:|---:|---:|
| 32 | 8 | `2026083123593201` | 0–31 |
| 64 | 16 | `2026083123596401` | 0–31 |

每个配置仍以x-fast row-major顺序，每站点消耗一个`mt19937_64` word，
`word < 10934234699625173385`时occupied。原32个pair仍是4×4 translation
anchors各自的水平与竖直`L/4`位移，顺序为`j,i,H,V`。

对每个原配置、每个站点z，producer分别构造`z=0`与`z=1`：

1. 重建square-NN黑组件和square-matching白组件；
2. 用digital Alexander恒等式得到原`q=r−1`与`E=q²`；
3. 对32个pair逐一得到`g16_0/g16_1`、`shared0/shared1`与Bell8 key；
4. 检查由原z bit选中的forced state逐pair恢复独立算出的原配置；
5. 在pair/site层计算midpoint product primitives及其绝对值，再聚合。

occupied marked endpoint仍返回`g16=0, shared=−1`。vacant port使用真实物理
edge ID；在r≥8的冻结pair上八个incident edges互不重叠，因此与旧producer
的八singleton语义逐字相同，同时避免以后误把共享edge静默拆开。

## 2. midpoint单位与不可省略的恒等式

令`dg16=g16_1−g16_0`、`dq=q1−q0`。输出保存

```text
q_observable_num32 = (g16_0+g16_1) dq
q_kernel_num32     = (q0+q1) dg16
q_product_num16    = q1 g16_1-q0 g16_0
```

并在每个pair/site强制

```text
q_observable_num32+q_kernel_num32 = 2 q_product_num16.
```

`E`有完全相同的三列与恒等式。observable/kernel numerator的物理单位是
`1/32`，product-delta numerator是`1/16`，因为`g=g16/16`且midpoint另除2。
所有`abs_*`列都在单个pair/site取绝对值后才累加；代码没有先做anchor sum再
取绝对值。

这些是未中心化原语。Issue #536中的完整通道还要用总体`E[g]`、`E[q/E]`
中心化，并进入

```text
J2/A_N=(jY_p-R*jM_p)/D-(Y_pp-R*M_pp)*jM/D^2.
```

64个preflight配置不能估计这些总体量，producer也不把它们冒充完整J2 scorer。

## 3. carrier是可交叉的有限标签

CSV保存完整四bit mask，而不是强迫四类互斥：

| bit | 标签 | 有限定义 |
|---:|---|---|
| 1 | `two_bridge_persistent` | 两forced state均endpoint-vacant、`shared0=shared1=2`、至少一个g非零且至少一个primitive改变。 |
| 2 | `shared_transition_or_merger` | 两endpoint均有效、shared数改变、至少一个g非零且primitive改变。 |
| 4 | `kernel_preserving_topological` | `dg16=0`、共同g非零，但`dq`或`dE`非零。 |
| 8 | `kernel_changed` | `dg16!=0`；另以count拆成`kernel_only`与`joint`。 |

bit 1与bit 4可以重叠，这正是必须保留的信息。完整mask作为row key，所以每个
pair/site仍只进入一个聚合row；按某个bit再次求和时才得到非互斥carrier view。
`mask=0`也保留。标签是有限callback分类，不是因果比例或连续场。

## 4. 空间分层与输出

每个pair/site先算

```text
d=min(d_Linf(z,x),d_Linf(z,y)).
```

`d=0`是shell 0；`d>=1`时shell k为
`[2^(k−1), min(2^k−1,L/2)]`。另存relation mask：bit 1/2分别是x/y
endpoint，bit 4/8分别是x/y的四邻接square-NN。这里NN不是L∞的八邻域。

一次运行产生：

- `PREFIX.config.csv`：64个原配置的K、q、E、occupation fingerprint、原32pair
  g16总和与shared分层，可用于核对旧producer语义；
- `PREFIX.shell.csv`：grain为`configuration × pair × shell × relation mask ×
  carrier mask`，保存signed/absolute primitives、状态值及identity controls；
- `PREFIX.metadata.json`：seed/count、线程数、两尺寸和总wall/CPU、peak RSS、
  carrier/单位与证据边界。

实现对全部约524万pair/site callback逐项计算，但不序列化一个可能数百MB至GB
的detail CSV；它保留pair ordinal并在预定shell/relation/mask内存充分统计量。
若后续确实需要逐z取证，应先修改合同说明具体缺少哪一项，而不是在结果出来后
无界打开raw dump。

## 5. 已执行命令与结果

最终runner核对kernel SHA256
`36ae069d370b1d7a4398861c928afb41aa76885c8895c696b1bc0c97e9c314fd`后，在本地
Apple ARM64、10线程运行：

```bash
clang++ -O3 -std=c++17 -pthread \
  scripts/p337_thermal_pivotal_preflight.cpp \
  -o /private/tmp/p337-thermal-preflight.1Xhiii/p337_thermal_pivotal_preflight

/private/tmp/p337-thermal-preflight.1Xhiii/p337_thermal_pivotal_preflight \
  --kernel analysis/regular_pair_spatial_kernel.tsv \
  --output-prefix /private/tmp/p337-thermal-preflight.1Xhiii/preflight \
  --threads 10
```

线程数只改变速度并写入metadata，不改变任务划分或输出顺序。实现按站点并行、
按配置串行；RNG先单线程生成完整配置，因而不会因线程调度改变既有counter。
主成本确为每个配置的`2N`次完整拓扑重建，即每尺寸`O(32N²)`；kernel约32MB，
其余内存主要是每线程DSU scratch和一个配置的`N×32`forced pair结果。

producer本身不调用shell计算SHA256。runner必须另建receipt，记录compiler/target、
binary/kernel/output SHA256、确切命令、exit code、wall/CPU/RSS。任何失败都按原
输入报告，不更换seed、窗口或counter。

producer内部计时为3.43953 wall seconds、25.9780 CPU seconds，peak RSS
40,173,568 bytes；随后immutable raw搬入仓库并只运行一次冻结描述性scorer：

```bash
python3 scripts/analyze_p337_thermal_pivotal_preflight.py \
  --config results/p337-thermal-pivotal-preflight/raw/preflight.config.csv \
  --shell results/p337-thermal-pivotal-preflight/raw/preflight.shell.csv \
  --metadata results/p337-thermal-pivotal-preflight/raw/preflight.metadata.json \
  --output-json results/p337-thermal-pivotal-preflight/latest.json \
  --output-md results/p337-thermal-pivotal-preflight/REPORT.md
```

scorer严格核对64 configs、两个固定seed/counter集合、每个config/pair恰好N个
site、完整carrier-mask partition、pair/site absolute下界和全部midpoint identity；
只汇总L/shell/mask/bit/endpoint-NN视图，不做显著性、总体centering或full J2。
任一输出已存在时拒绝覆盖。

runner可把三份immutable raw从临时执行目录复制到上述`raw/`目录。metadata仍
保留原绝对路径；scorer要求basename不变、把旧路径写入`relocated_from`，并对
仓库中的实际输入重新计算SHA256，而不把安全搬迁误判为语义变化。

64个配置、32 pair/config以及5,242,880个pair/site callback全部通过结构检查；
q/E midpoint恒等式的残差总和与逐项最大绝对值均为0。有限样本里只有18个
`kernel_changed` callback：L32有11个（10个kernel-only、1个joint），L64有7个
（全为kernel-only）。它们全部同时属于`two_bridge_persistent`或
`shared_transition_or_merger`；没有观察到`kernel_preserving_topological` callback。
这不是该carrier的零概率结论。

空间上，L32全部非零q/E原语落在shell 1；其中一个two-bridge kernel-only事件来自
square-NN关系，余下10个shared-transition事件来自external关系。L64的7个事件全为
external shared-transition，5个位于shell 2、2个位于shell 3。endpoint原语全部为0。
未中心化总和显示L64 observable/rank列为0而kernel列非零；L32只有一个joint事件使
observable/rank列非零。**因此这次preflight最有价值的判决是：现有有限回放支持一个
稀疏、外部、kernel-change主导的接口实现，而不能把N25的负observable/rank总体项解释成
大量局部q/E翻转。** 仓库后续已把路线收紧为开放 #537：先对同一 C4/Schur 投影下的
ordinary four-arm landing transfer matrix 检查全部 `2x2` minors，再只为存活的
non-rank-one signed functional 或 four-packet remainder 做 tail；不是扩大 counter
或继续描述性分类。

## 6. 内置语义控制与剩余边界

实际preflight启动后会先执行两个Issue #536的L4控制：

- `{3,4,6,8,10,11}`、marks 0/2、z=5：ambient rank为`0→1`，所以
  `q=r−1`是`−1→0`；`g16:4→4`，key为`[0,1,2,3,0,3,4,5]`；
- `{3,4,6}`：ambient rank恒为0，所以`q:−1→−1`；`g16:0→4`。

失败信息会逐字打印实际`q0/q1`、`g16_0/g16_1`、`shared0/shared1`与两key，
便于区分rank/q混写、physical-edge alias和kernel lookup问题。还逐项检查
`s<=1 => g16=0`、forced/original复现和q/E midpoint identities。
第一次执行在第一个L4 gate、任何科学配置复放之前暴露了ambient rank与
`q=r−1`的记号混写；只修正该语义断言并保留所有冻结输入后，第二次执行通过
两个L4 gate及全部逐项控制。两次编译/运行尝试和首个失败均记录在
[`run.json`](../results/p337-thermal-pivotal-preflight/run.json)，没有以换seed、
counter或kernel掩盖失败。

当前唯一需要显式保留的语义边界是：carrier masks按有限`shared`数与kernel变化
定义，不能单凭64配置把`shared_transition`升级为渐近remote-merger主导，也不能
把`kernel_preserving_topological`升级为独立场。结果无事件或不分辨仍按64配置
停止；扩大旧数据回放或启动新生产都需要另一份carrier-specific合同。
