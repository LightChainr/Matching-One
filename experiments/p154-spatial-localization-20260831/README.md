# P154 spatial localization and global-U bridge

Completed results: [REPORT.md](REPORT.md), [results/latest.json](results/latest.json). This experiment consumes existing conditional-line/source marks; it does not repeat their replay, first line-source comparison, or first fixed-K localization.

The bundle is self-contained. [SOURCES.json](SOURCES.json) maps every compressed input to its fixed GitHub commit/path and SHA256. `inputs/*.gz` contains unchanged CSV/JSON bytes compressed for transport. `frozen_moments.py` preserves the original Binomial integration function. No git access is needed to run the bundle.

Required: Python3.9 or later, NumPy and SciPy. Executed with Python3.9.9, NumPy1.26.4, SciPy1.13.1 on Huawei ARM64, one BLAS thread:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python analyze_angular_bridge.py --output results
```

The shipped `results` is immutable to this script. An explicitly requested reproduction can use another relative output path, for example `--output reproduced-results`; it still analyzes the same six sizes and same100 aligned omissions, with no new configurations or root finding.

The frozen analysis is [EXPERIMENT.md](EXPERIMENT.md). It includes the soft U+/U− and bulk-source allocations, fixed q2/Jordan chain combinations, and a rank1/fixed-K centered-source decomposition whose global-U derivative is structurally zero. The latter is a conditional projection, not an independent physical intervention. All correlations and source-subset precision limits remain explicit in the report.
