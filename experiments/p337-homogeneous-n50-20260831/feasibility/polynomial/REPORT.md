# 固定 N50 的 K 多项式压缩：有界对照

**结果：移去 key 中的 K 后，两幅 N50 父图都在原资源预算内完成到第40层。** 这是可计算性进展；本探针没有完成第50层、生成N50终表或做科学评分。根任务已将独立多项式 producer 冻结为 `4ae4e710404e40ac31eaf962680aaceed2543ddf`；其后完整目标由根任务负责，不属于本报告结果。

只改 K 的存储：rank0二维gain、rank1垂直投影、rank2丢弃gain，以及同一物理边位移、规范化和顶点顺序全部保留。值是动态连续 K 区间上的 `(count,sumS)`；occupied 分支把指数加1。没有固定51单元数组，也没有新孔阶、源、观测量或参数扫描。

每个系数使用 signed64，所有相加与 `sumS+dS*count` 先经 signed128 并检查范围。N≤50、保守 `|partialS|≤6N+1`，使总绝对和小于 `(6N+1)2^50 < 2^63`。整数范围门还覆盖 packed signed16 gain。

## 已完成的独立整数核对

N9/N13 使用既有直接 lifted-homology oracle；N25 使用固定 `a70eeff09f51ce2fa0fea5ae637e9191efbf2e1f` 的完整 `(K,g,q,count)`，通过 `S*=51-K-g` 得到每 `(K,q)` 的 count/sumS。四项输出均逐整数完全相等，并检查全部K边际为binomial、总count为2^N。

| 图 | 原K-key峰态 | 多项式boundarykeys峰态 | 原CPU秒 | 多项式CPU秒 | 原RSS MiB | 多项式RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| N9 (3,0) | 162 | 81 | .005448 | .006524 | 1.48 | 1.50 |
| N13 (3,2) | 1010 | 437 | .003993 | .003614 | 1.72 | 1.64 |
| N25 (5,0) | 68,305 | 14,237 | .178955 | .071762 | 25.17 | 9.47 |
| N25 (4,3) | 59,474 | 10,701 | .123493 | .036062 | 18.59 | 6.97 |

CPU与RSS均来自各独立子进程的 `wait4`；不是估算。

## N50只跑到第40层

本次每图：35秒CPU软门、40秒OS硬限、1792 MiB RSS门、5e6 boundarykeys门、max40层；串行运行，无完整目标授权参数。

| 图 | 原K-key探针停止 | 多项式停止 | 多项式CPU秒 | 多项式RSS MiB |
|---|---|---|---:|---:|
| (5,5) | 第33层后，20秒CPU门；峰4,901,320态 | 正常达到第40层门 | 8.362659 | 660.0625 |
| (1,7) | 第28层后，构建29层时触5e6态门 | 正常达到第40层门 | 20.126842 | 1287.703125 |

具体压缩规模：

- (5,5) 第39层峰 **842,882 boundarykeys / 9,329,835非零Kcells**，值区间占149,277,360字节。第40层降为546,214 keys / 6,482,388 cells / 103,718,208 valuebytes。
- (1,7) 第40层为 **1,730,742 boundarykeys / 19,952,451非零Kcells**，值区间占319,239,216字节。此时keys/cells仍在增长，不能说已经越过全程峰值。
- `valuebytes` 只计该完成层的连续系数数组；实际RSS还包含两个轮换map、key、vector及hash分配开销。逐层文件也记录capacity_valuebytes，避免把逻辑值大小当总内存。

同一物理状态的不同 K 从多个hash节点合并到一个动态值区间。新的非零Kcells仍是需要保留的准确整数信息；压缩没有通过删除物理状态或截断K支持来节省内存。

尚余10层，完整峰值与终表此前没有测得。本报告不将部分层的速度外推成全程保证，也不把根任务随后启动的完整生产预先写成完成。

## 产物与停止状态

- `frontier.cpp`：独立C++17多项式实现；旧 `/tmp/p337-black-cpp/frontier.cpp` 保持原样。
- `geometry_*.txt`、`INPUTS.json`：固定Python prepare导出的边与frontier，含physical dx/dy及来源hash。
- `expected_*.json`、`check_*.csv`：四个完整旧输入/小图对照。
- `RUNS.json`：代码/binary hash、每项argv、CPU、RSS、退出状态及整数核对。
- `probe_5_5.jsonl` / `probe_1_7.jsonl`：每层boundarykeys、nonzeroKcells、valuebytes及RSS。
- `probe_*.json`：`complete=false`、`pre_final_layer_gate`，没有N50 CSV。

编译使用 `clang++ -O3 -std=c++17 -Wall -Wextra`，无编译警告。CLI允许根任务后续明确给出至多600秒CPU、18432 MiB RSS的有限预算；**本次实际只执行上述小预算**。完整N50还需外层冻结验证driver提供真实authorization commit；C++仅记录该值，不声称自己核验Git。

全部探针子进程均由自己的runner wait4回收并正常退出。没有自己的后台计算、云机或隧道待清理；未修改研究树或冻结copy。按根任务指令，本项到此停止，不追加算法或测量。
