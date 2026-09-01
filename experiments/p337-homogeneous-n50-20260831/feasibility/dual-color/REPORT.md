# N50 齐次全孔端点的精确可计算性：有界实测

**结论：没有发现现成、已跑通的齐次 N50 数据或 transfer 引擎。双颜色 frontier + Euler 恒等式确实能精确保留所需六个矩，tiny 对照已逐整数通过；但当前 Python 表示在 N50 第19/20层已达到资源门，不能把它宣布为全程预算可控的 N50 求解器。** 没有计算或发布 N50 的 U/source-response。

全部内容只写 `/tmp/p337-n50-feasibility-20260831/`。没有修改冻结文件、启动云机、MC、增加孔阶、commit 或 push。三段探针记录的 CPU 合计 **12.400078秒**（不含少量启动/导入），最高进程 RSS **400.03125 MiB**，随后全部退出；安全限制为55秒初始 CPU硬限、较短进程CPU门和400/500 MiB RSS门。

## 1. 固定仓库快照与去重结果

- 研究工作树：`f3ecde7da04d9e01047d1a8bc7eb27d7d048fa78`。
- 新理论/程序快照：`0dda27bab3d1b6a749a0a32b3dde666b7fe9a0dd`。
- `SOURCES.json` 记录9个实际复制的既有 Python 文件及 SHA256；所有副本来自固定 f3ecde7d Git blob，不依赖移动中的工作树。

实际入口检查：

| 文件/数据 | 实际能力 | 与本目标的差别 |
|---|---|---|
| `scripts/digital_alexander_quotient_frontier.py` | `build_state_table` 遍历 `range(1 << N)`，随后在 Boolean subset lattice 上计数全部排列 | “frontier”指小阶研究前沿，不是边界状态压缩；N50仍有2^50表项 |
| `scripts/noncrossing_connectivity_codec.py` | width≤8的非交叉 partition 编解码 | 自身合同明确不含 transfer weights、torus sectors、row propagation |
| `scripts/planar_transition_table.py` | detach/join 的小宽度状态转移表 | 不含占据/空置双颜色、周期闭合、q/S* 和热权重 |
| `scripts/transfer_resource_probe.py` | 单命令时间/RSS记录 | 不是 transfer engine |
| `0dda:scripts/p337_s4_trace_exact.cpp` | 固定N25、黑DSU增加模2/3 seam flags，仍在完整二叉子集树访问2^25叶子 | `static constexpr n=25`、数组25；不是N50压缩算法 |
| `0dda:scripts/p337_closed_source_exact.cpp` / `finite_exact.cpp` | N25 `(5,0)/(4,3)` 的source矩/直方图 | 固定N25、完整2^25遍历 |
| `0dda:scripts/p337_endpoint_defect_exact.cpp` / `defect_reweighting_exact.cpp` | 实际N50图，但只放开25个B点；A饱和或固定一个A孔 | 是零/一孔端点输入，不是epsilon=1全部50点自由 |
| `results/p337-s4-trace-transmission/` | Q4、N25规范seam投影的已完成响应 | 已完成，不应重做；不是t=0齐次N50 |

对两个固定树的 `scripts/src/analysis/notes/results` 名称与正文搜索没有找到本目标的完整 N50 `(K,q,S*)` 数据。`positive-decimation-source-mainline.md` 已明确饱和端点不等于齐次 N50/N100。该判断针对以上固定快照，不是对未来提交的断言。

## 2. 为什么不需要完整 (K,q,S*) 直方图

这里固定 epsilon=1、t=0，因此50个站点具有同一个 Bernoulli p，原 q/E 不变。对于 source 的**一阶响应**和所需任意热导数，只需每个 K 的六个整数：

\[
M_{ij}(K)=\sum_{|A|=K} q(A)^i S_*(A)^j,
\quad i=0,1,2,\quad j=0,1.
\]

即 `count, sum_q, sum_q2, sum_S, sum_qS, sum_q2S`。这不是新观测量；它们正是原 q、E=q² 与同一个 source 的必要原始矩。

用 `p^K(1-p)^(50-K)` 合成它们，得到原 separately normalized q/E、S、qS、ES；source导数是对应协方差。热导数通过这些固定次数的 Bernstein 多项式求导，不需要额外状态。最后必须按原共同 root 与 slope motion 计算

\[
z_t=-M_t/M_z,
\quad
\frac{D_{S_*}U}{A_{50}}
=\frac{Y_{zt}+z_tY_{zz}}{M_z}
 -\frac{Y_z(M_{zt}+z_tM_{zz})}{M_z^2}.
\]

`z=logit(p)`，`M=(q_first+q_second)/2`，`Y=(E_first-E_second)/Delta`，父图固定顺序 `(5,5)/(1,7)` 有 `Delta=-1152/625`。本次没有执行这个 N50 root/response 评分。

完整 S 分布在有限 t 重加权时才需要；本目标 t=0一阶 source 响应无需把 S 放入状态 key。

## 3. 已实现并验证的具体 frontier 算法

实现文件：`probe.py`（q在key的初版）与 `moment_probe.py`（六矩版）。固定顶点顺序就是已有 `integer_period_torus` 的 quotient-key 排序，不做物理参数选择或孔阶扫描。

边界保留仍有未来 NN/NNN邻接的已处理站点，记录：

1. 每个边界站点的颜色，以及黑 NN / 白 matching 连接组件的规范化 partition；正/负标签分别表示黑/白。没有强行用平面 Catalan 状态删掉 torus 连接。
2. K。

在已证明适用的 honest square-cell quotient 上：

\[
q=C_B-C_W-K+B-F_4,
\qquad
S_*=C_B+C_W+F_4+2N-4K+B.
\]

这里 B 是占据NN边数，F4是全占据四角面数。探针也检查两个N50图以及tiny图的NN/matching图简单、每个面四角不同，且所有面计数在最后一个顶点出现时结算。组件离开边界后不再可能与未来连接，其历史贡献已进入值，不需要继续保留组件标签。

局部增量完全由边界状态和新颜色决定：

| 事件 | dq | dS |
|---|---:|---:|
| 新黑点 | 0 | −3 |
| 新白点 | −1 | +1 |
| 新占据NN边 | +1 | +1 |
| 合并两个黑组件 | −1 | −1 |
| 合并两个白组件 | +1 | −1 |
| 完成一个全黑面 | −1 | +1 |

初值为 `q_partial=0, S_partial=2N`。同状态合并严格加和。与当前纯黑/gain方案的初值2N+1并不矛盾；这里包含白组件，最终两种source恒等式在tiny对照中逐配置一致。

六矩版不把 partial q 放进 key。若值为 `(n,Q,Q2,S,QS,Q2S)`，转移 dq=d、dS=e，则

```text
n'   = n
Q'   = Q+d*n
Q2'  = Q2+2*d*Q+d*d*n
S'   = S+e*n
QS'  = QS+d*S+e*Q+d*e*n
Q2S' = Q2S+2*d*QS+d*d*S+e*Q2+2*d*e*Q+d*d*e*n
```

这是精确有限阶多项式平移，没有额外闭合假设。所有探针值使用 Python 整数；移植时不能仅凭2^50 fits64就把高阶中间矩用64位累加。本任务的保守界例如 `|partial q|<=5N`、`|partial S|<=10N`，故六矩和可超过64位，但远低于signed128容量。

## 4. 实际执行结果

既有入口 `planar_transition_table.build_width_table(8)` 已实际运行：

- **1430状态，22,880条序列化记录，1,977,229字节，0.605602 CPU秒**。
- 记录 hash 与原小宽度合同吻合；它只证实 planar detach/join 已有实现，不意味着1430就是本N50 torus状态数。

独立tiny核对使用固定源码 `integer_period_torus.classify_configuration`，从真实 lifted homology 直接得到 rank，不通过 frontier Euler式计算 q。

| tiny | 所有配置 | q-key版峰状态 | 六矩版峰状态 | 逐整数核对 |
|---|---:|---:|---:|---|
| N9 `(3,0)` | 512 | 256 | 256 | 完整(K,q,count,sumS)与六矩全部一致 |
| N13 `(3,2)` | 8192 | 2048 | 2048 | 完整(K,q,count,sumS)与六矩全部一致 |

移走 q-key 的实际收益发生在后段：N13第12层由 **2040→1736态**，末层 **22→14态**。早期峰值不变，因为还没有大量历史连接完成。

N50测量严格停在资源门，**没有完整N50矩，更没有N50科学评分**：

| 父图 / 实现 | 完整图最大边界宽度 | 最后完整层 | 该层状态 | 停止条件 |
|---|---:|---:|---:|---|
| `(5,5)` q-key小探针 | 20 | 16 | 65,536 | 50,000状态门 |
| `(1,7)` q-key小探针 | 16 | 16 | 65,536 | 50,000状态门 |
| `(5,5)` 六矩探针 | 20 | 19 | 524,288 | 400,000状态门，阶段峰RSS384 MiB |
| `(1,7)` 六矩探针 | 16 | 19 | 381,520 | 处理第20层期间，进程RSS400.03125 MiB停止 |
| `(1,7)` q-key比较 | 16 | 19 | 381,520 | 300,000状态门，独立进程峰RSS265.484375 MiB |

N50 `(1,7)` 第19层两种key的状态数**完全一样**；六矩版此时没有减态，增加的值字段反而更耗内存。不能把N13后段的收益直接移植成N50预算保证。

`(5,5)` 的固定顺序在前19层仍有19个已处理点全部留在边界，524,288正好是2^19。此处首先遇到的是边界颜色本身的指数增长，不是一个已测得的额外homology状态爆炸。`(1,7)` 同层已把524,288个prefix合并为381,520态，但尚未跨越整体峰值。

所有 RSS 是实际进程的 `getrusage` peak；同一进程内后续阶段继承早期高水位，不能把它误当该阶段单独RSS。q-key比较在独立进程执行，因此265.484375 MiB是其独立峰值。RAM门检查周期0.05秒，400.03125 MiB是门触发后的实际采样高水位，低于500 MiB安全上限。

## 5. 可复算产物与本次停止点

- `SOURCES.json`：固定代码来源及各文件hash。
- `probe.py`, `probe_result.json`：既有入口、独立tiny oracle、q-key和首轮N50资源层记录。
- `moment_probe.py`, `moment_probe_result.json`：六矩平移、tiny完整矩、N50完整层及内存停止记录。
- `q_key_comparison.json`：同一N50 `(1,7)` 第19层的独立q-key对照。
- `MANIFEST.json`：以上产物及固定源码的SHA256。

纯tmp可复算命令（各脚本自带资源门）：

```sh
python3 /tmp/p337-n50-feasibility-20260831/probe.py
python3 /tmp/p337-n50-feasibility-20260831/moment_probe.py
```

第一段会重跑已有tiny对照和有界资源测量，不是生产；第二段读取第一段tiny结果后核对六矩。**不应把删掉资源门后直接运行此Python原型当成已验证的全N50计划。**

可交接的具体结论是：原全孔核在epsilon=1/t=0确有一个无需2^50叶子逐一处理的精确状态合并算法，且本次已验证它保留原q/E/S与热导数所需信息；现有仓库没有已完成的N50引擎/数据，Python双颜色表示的实测资源门在第19/20层触发。全程峰值和耗时仍未知。本次按根任务要求停止，不再扩算法或增加测量。
