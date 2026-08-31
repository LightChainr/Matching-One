# P398 固定总速率干预

**已完成：** 1430态完整有限模型在η=0,±1/4上的旧读出响应，包含每个η的平稳分布、完整生成元与旧14/16维几何比较。实际远端计算1.709秒。结果和科学边界见 [results/REPORT.md](results/REPORT.md)，机器结果见 [results/latest.json](results/latest.json)。

独立实验包，不修改 Matching-One 的既有工作树。数学配置见 `EXPERIMENT.json`；旧状态操作、配置函数与基准结果冻结于 `1f19fc1a2d9fc59dce650e95268c716762725985`，来源与 SHA256 见 `frozen_sources.json`。`frozen_model.py` 是从所列原文件逐函数原样抽取，不运行原模块的额外任务。

只需 Python 3、NumPy、SciPy。实际执行环境为 Python3.9.9、NumPy1.26.4、SciPy1.13.1；Huawei ZyTrST ARM64单线程。原始运行命令：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 run_rate_intervention.py --output results
```

结果包含 `results/latest.json`、`results/REPORT.md`、三个完整稀疏 forward 生成元（乘4后的整系数）和三个完整 character-i 扇区的矩阵。程序发现已有 `latest.json` 时拒绝覆盖，不启动参数扫描、MC 或 width 扩展。

本包已包含该次结果。如明确需要复现，用同一个命令指定一个新的相对目录，例如 `--output reproduced-results`；程序始终固定同样的三个η与旧lags，不覆盖本次归档。`SHA256SUMS`提供包内文件校验；`results/latest.json`保留实际执行时的原输入SHA，包括执行后已补充说明的README及已更新最终包校验清单的原哈希；执行脚本和frozen_model哈希没有变化。

本地仅做过 `--validate-width4` 小规模实现校验：生成元精确 K 共轭余量0、平稳共轭余量2.78e−17、same-even/cross-odd余量2.29e−16。未在本地运行 width8 主实验。

新干预同时改变生成元与平稳分布。原14/16维子空间按固定配置函数定义，在新平稳内积下比较；此约定不等于冻结旧η=0矩阵。两个旧93维ray在η≠0时允许耦合，完整传播使用186维不变扇区。

本实验只适用于有限宽度 join/detach 模型；不是 square-site norm4、CFT或连续极限识别。
