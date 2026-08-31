# P154：独立 lag1 传递检验完成，执行冻结停线规则

独立新样本排除了“进入主导”和“完成主导”两项有限尺度数值预测。两个尺度的四个通道区间全部落入双通道弱模板；净响应区间也都完全落入冻结的 ±0.50 带，因而**停止把此 lag1 源作为当前主要 H4 解释继续优先投入**。不追加样本、不更换 lag/source，不增加第四个救场模板。

冻结版本为 [`0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f`](https://github.com/LightChainr/Matching-One/commit/0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f)。[CONTRACT.json](CONTRACT.json) 的七个固定文件逐字核对了该提交；授权与哈希见 [authorization.json](authorization.json)。新生产在冻结和明确授权之后开始，完成全部九片、核验后只运行一次主评分。

| N | 读数 | 点估计 | SE | 六坐标同时区间 |
|---:|---|---:|---:|---|
| 85 | 进入贡献 | 0.044492 | 0.033353 | [-0.043502, 0.132487] |
| 85 | 完成贡献 | -0.001022 | 0.028667 | [-0.076653, 0.074609] |
| 85 | 净响应 | 0.043470 | 0.043631 | [-0.071640, 0.158581] |
| 340 | 进入贡献 | -0.053571 | 0.063183 | [-0.220264, 0.113122] |
| 340 | 完成贡献 | 0.114247 | 0.061800 | [-0.048797, 0.277290] |
| 340 | 净响应 | 0.060675 | 0.082657 | [-0.157394, 0.278745] |

六坐标采用同一个 Bonferroni 渐近正态 95% 家族，临界值 2.6382572735。W 要求两个 N 的每个通道都在 [−0.30,0.30]；本次四个区间均整体包含于该带。B 要求进入贡献 ≥0.60 且完成贡献在弱带，C 交换两通道；B、C 的主导要求各在两个 N 都被排除。W 的状态仍记作 `not_excluded`，区间整体相容不等于识别了一种物理理论。净响应判据独立预定为两 N 的区间均落入 [−0.50,0.50]，本次明确满足。

这里检验的是同一个一步条件政策：在 K−1 的早期 rank 层内以簇数 s=C_B+C_W 作指数倾斜并层内归一化，再均匀加入一个空位；一阶 score 为 s−E[s|K−1,r,g]。条件均值与 pooled matching root 在全新数据和每次删批中重估，旧发现期数据未并入新点估计或协方差。

读数按 F₁=I(K≥K₋)、F₂=I(K≥K₊) 定义，q=−1+F₁+F₂、E=1−F₁+F₂。令 A=N^(13/8)/2、D=mean(q_p)、P₄[f]=(f_first−f_second)/(2304/1445)，则 U_entry=−A·P₄[F₁,p]/D、U_completion=+A·P₄[F₂,p]/D。两部分的源导数共用完整源的 rootdot 和 Ddot，逐次回加为原 U 的导数；0→2 事件仍保留在 q 与 root/slope 响应里。这是一个读数在同源下的贡献分解，不能用来认定源对两个事件的独立因果责任。

每个 N 有 200 个等权批次；同一排列的两个方向、全部 K、全部通道始终配对。两个 N 使用显式不同随机 domain，分别删批；完整 32×32 协方差、两组全部删批向量和主六坐标 6×6 协方差都保存在 [PROSPECTIVE_RESULT.json](PROSPECTIVE_RESULT.json)。独立复核已从原始统计重算32个原定坐标和400个删批向量，协方差、回加及冻结判定一致，见[INDEPENDENT_REVIEW.json](INDEPENDENT_REVIEW.json)。新样本实际方差均低于含 1.25 保护因子的规划方差，比例为 0.626–0.890；无需也没有补样。

统计边界集中如下：阈值是预定科学决策分辨率，区间是基于 200 批删批方差和正态近似的渐近区间，不是严格有限样本证书。本结论只约束该 lag1 条件政策在 N85/N340 的响应，不声称严格为零、任意 N 的缩放、其他源或整个 H4 理论已被否定。

实际生产使用 ZyTrST、TgFr7R、XPk2PZ、HZsCM6、TV2N0X 五台 ARM 容器：每台实测 CPU 配额 14.5 核、内存 25 GiB，使用 14 worker、BLAS=1、GCC 10.3.1。N85 为 5,000,000 条，N340 为 160,000,000 条（八片各 20,000,000），共 165,000,000 个新排列；无重复分片、失败重跑、预算外样本。生产墙钟 590.215 秒。冻结前 QA 另复现了 11,000 条旧排列，未混入生产。

所有新增原始统计均已取回到 [production/](production/)：九个 `.csv.gz` 和对应 `.run.json`，包含每批、每方向、每 K 的 q,E,s,qs,Es、01/02/12 计数及前一步源之和。每片 receipt 记录精确随机区间、源 domain、冻结版本、实际编译命令、二进制及 raw/gzip 哈希。传输包哈希见 [DELIVERY_ARCHIVE_SHA256.json](DELIVERY_ARCHIVE_SHA256.json)，逐项核验及计数见 [DELIVERY_VERIFICATION.json](DELIVERY_VERIFICATION.json)；五机日志与环境在 [runtime/](runtime/)。主评分约 0.644 秒，唯一调用的版本和命令见 [SCORE_RUN.json](SCORE_RUN.json)。

2026-08-31 11:58:53 UTC 交接时，五台均已完成队列且无 P154 活进程，所有结果已取回；机器与其单隧道仍保持 Running 交主任务统一关闭。这是关闭前的历史交接。主任务随后逐台核查并关机、停止自有隧道；2026-08-31T20:02:29+0800的[最终云状态](CLOUD_COMPLETION.json)确认本次使用五台均Ready，现场十台也均Ready；第二账号只读查询，未由本任务开关。

旧档案仅用于冻结前规划和 QA。精确发现期 lag1 CSV 与摘要来自 `4daae57eef5c945aa050a95cd3d5d5d77582161b` 的 `results/norm4-lagged-source/`，已原样置于 [inputs/](inputs/)，来源/哈希见 [SOURCES.json](SOURCES.json)。`old_profiles.npz` 的每 N 数组为 (100,2,N+1,5)，末维依次是 CSV 字段 `sum_q,sum_e,sum_s,sum_qs,sum_es`；方向按 first/second、行按 batch/K 排列。原始 bulk CSV 固定于 `7da1eeb0e51cf430987dbf204d23713c2ab5a46c`：N65/85/130/170 取 `results/norm4-source-thermal/raw/nN.csv`；N260/340 按同一 batch/方向/K 将该文件与 `results/norm4-source-endpoint-1m/increment/raw/nN.csv` 的五个整数和逐项相加（每批 1000+9000）。原始 CSV 哈希另列于 [OLD_RAW_SOURCES.json](OLD_RAW_SOURCES.json)。这些旧文件及旧 root 不被 fresh scorer 读取。

两个 vendor 是 `bfab0330f5f56ca4d746b45d737f1607e3d229a0` 的 `src/threshold_rank_orientation_mc.cpp`、`src/threshold_rank_integer_period_mc.cpp` 原样副本。冻结新 producer 使用它们的几何、同调并查集和 counter RNG，不依赖其他工作树。

从包含本结果的干净 checkout 复算，只需复制本目录到临时位置、在副本移开已保存输出，再对九个已交付原始统计运行冻结 scorer；不需要云机或生成任何排列：

```bash
python3.11 -m venv /tmp/p154-score-venv
/tmp/p154-score-venv/bin/pip install numpy==2.3.5 scipy==1.16.3
```

在仓库根目录运行：

```python
from pathlib import Path
import shutil, subprocess, tempfile
source = Path("experiments/p154-prospective-transmission-20260831")
copy = Path(tempfile.mkdtemp(prefix="p154-recheck-")) / "package"
shutil.copytree(source, copy)
(copy / "PROSPECTIVE_RESULT.json").rename(copy / "PROSPECTIVE_RESULT.saved.json")
subprocess.run(["/tmp/p154-score-venv/bin/python", str(copy / "score_production.py"),
                "--input-dir", str(copy / "production")], check=True)
```

复算后的路径和计时元数据可不同，数值、协方差和冻结判定应一致。原包的唯一生产评分结果不覆盖；冻结脚本、合同不修改。
