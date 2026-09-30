# 固定负向目标：独立既有档案复核

唯一科学任务：将前块的固定预测 `deltaY < 0` 应用于旧 completion-hazard
档案的14×10000个排列，补读精确协同边数e。该档案其它结果早已读过；
本次是此前未测过e的固定目标读出，不是全新前瞻采样或原封未看的holdout。

已完成：deltaY=-0.003199368397±0.000620090643（14批jackknife SE），
负向预测在此独立既有档案中再现。前块为-0.003449385668±0.000772252763；
本块减前块=+0.000250017271±0.000990397261，两块点差仅0.2524个比较SE，未合并。
140000行全部对应；共同支持54120/54723=98.8981%，1202个共同格。
具体数字、原始统计及结论边界见[结果](results/RESULT.md)和[完整JSON](results/result.json)。

参数不变：square site、L=512、N=262144、a=154646、b=155385。
只在rank1且有安全空站点的前缀内，以精确整数方向D和当前完成数c分格：
早组J1≤a，晚组a<J1≤b，Ybar=2sum_e/[(N-b-c)n]，w=nE*nL/(nE+nL)。
这里c是完成数，不是旧manifest中第三个时间cutoff 156120。

保存deltaY、deltaE、deltaZ2、earlyY、lateY及5×5完整协方差。
14次对齐整批删除均重新求共同支持和权重；点估计不作jackknife偏差修正。
deltaZ2=-2deltaE/[(N-b)(N-b-1)]，是同一几何公式的派生量。
不合并两独立块；仅用两块点差及sqrt(SE1²+SE2²)比较。

## 来源与行对应

- `inputs/archive/`：completion-hazard-production-20260929旧manifest及14个压缩生产批次。
- `inputs/reference-result.json`：前块原结果的未改动副本；只读取，不重新评分。
- `source/`：两个原engine.cpp的逐字节副本，保留include相对结构。
- `source-record.json`：实际源hash、准备时HEAD、各源文件最近提交、批次来源。
- `data/`：本次精确重放原始行、每批命令和退出状态、完整run.json。
- `results/result.json`：完整评分、整数格充分统计、批次依赖、协方差和删除结果。
- `results/RESULT.md`：简短数值报告。

旧rank_b=0若J1>b、=2若J2≤b、否则=1；只把旧rank0的未来J1映射为0。
其余J1和dx_b、dy_b、nu_b原样逐行核对。旧排列的Fisher–Yates消耗独立mt19937_64
流，tau使用另一time_rng；沿用旧seed重放，不改变几何或随机数实现。
无伪造四探针，无原scorer调用，无新seed或独立样本。

## 实际运行命令

本地准备与最终评分使用研究Python3.11.15；只需标准库：

```sh
/Users/lc/python-envs/research-py311/bin/python analysis/exact-pair-independent-archive-20260930/runner.py prepare
/Users/lc/python-envs/research-py311/bin/python analysis/exact-pair-independent-archive-20260930/score.py
```

云端仅TV2N0X（4a8d1d443419434889e49148ed0a7ba6），新路径
`/workspace/Matching-One-TV2N0X/exact-pair-independent-20260930`。
启动前实时Ready且隧道已停；Running后仅一次ssh-key-reset，成功且密钥mode600。
初始进程清单无其他研究作业。实测16 affinity CPU、14.5核配额、25GiB内存；使用14 workers。

容器重启后系统g++未保留。首次benchmark命令在读取编译器版本时终止，尚未重放任何行。
随后把编译工具安装到本次路径的独立installroot；未更改现有研究目录或系统包集合：

```sh
yum -y --installroot=/workspace/Matching-One-TV2N0X/exact-pair-independent-20260930/toolchain --releasever=2.0 --setopt=reposdir=/etc/yum.repos.d --setopt=logdir=/workspace/Matching-One-TV2N0X/exact-pair-independent-20260930/toolchain-logs --nodocs install gcc-c++ time
```

第一次安装命令因log目录尚未创建立即退出，创建该目录后安装成功。
GCC10.3.1在该目录内chroot编译原源，Python3.9.9驱动重放；具体编译命令及二进制hash见benchmark.json。

```sh
cd /workspace/Matching-One-TV2N0X/exact-pair-independent-20260930
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 runner.py benchmark
nohup env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 runner.py replay > run.log 2>&1 < /dev/null &
```

一次32行成本/对应检查，使用旧batch0同seed，耗时1.020989秒，峰值75960KiB，全部映射对应。
此重复前缀的科学权重为0。完整重放PID1015，14批全部成功，墙钟520.875642秒。
原始行及每批命令均已取回。
完整评分合并执行140000行对应检查，不额外重跑原枚举、旧评分或全套测试。

此次仅写本目录，无commit/push、Issue操作、现有engine/文档/ledger修改或PR合并。
PR838读取时为OPEN/Draft；资源收尾事实见execution.json。

2026-09-30 02:25:19 UTC（北京时间10:25:19）已实时核实TV2N0X为Ready、己方隧道已停。
关机前无其他研究进程，远端数据与任务内编译环境保留。本范围未完成事项：无。
