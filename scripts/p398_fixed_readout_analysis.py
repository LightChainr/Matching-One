#!/usr/bin/env python3
"""Consume saved width-five C(d)/residues for named physical A/L readouts.

No transfer/eigenvalue recomputation, fit, parameter scan, or simulation.
The exploratory null filters follow one prescribed residue-algebra rule;
they are never selected or optimized against the resulting lag errors.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def cyclotomic(value):
    zeta = np.exp(2j * np.pi / 5)
    return sum(float(Fraction(coefficient)) * zeta ** power
               for power, coefficient in enumerate(value["coefficients"]))


def exact_matrix_embedding(rows):
    return np.array([[cyclotomic(value) for value in row] for row in rows], dtype=np.complex128)


def pair(value):
    return [float(np.real(value)), float(np.imag(value))]


def decode_pair(value):
    return complex(value[0], value[1])


def covariance(source, matrix, readout):
    return source @ matrix @ readout.conj()


def display(value):
    return f"{value.real:.12g}{value.imag:+.12g}i"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-input", type=Path, default=ROOT / "results/p398-physical-two-point/latest.json")
    parser.add_argument("--mode-input", type=Path, default=ROOT / "results/p398-physical-two-point/mode-visibility.json")
    parser.add_argument("--output-json", type=Path, default=ROOT / "results/p398-fixed-readout/latest.json")
    parser.add_argument("--output-md", type=Path, default=ROOT / "results/p398-fixed-readout/REPORT.md")
    args = parser.parse_args()
    started = time.perf_counter()
    exact_bytes, mode_bytes = args.exact_input.read_bytes(), args.mode_input.read_bytes()
    exact, modes = json.loads(exact_bytes), json.loads(mode_bytes)
    exact_sha, mode_sha = hashlib.sha256(exact_bytes).hexdigest(), hashlib.sha256(mode_bytes).hexdigest()
    correlations = [exact_matrix_embedding(value) for value in exact["scientific_result"]["correlations_by_lag_0_to_8"]]
    c0 = correlations[0]
    residues = [np.array([[decode_pair(value) for value in row] for row in mode["residue_A_L_rows_columns"]])
                for mode in modes["modes"]]
    eigenvalues = [decode_pair(mode["eigenvalue"]) for mode in modes["modes"]]

    # w^T O defines an actual linear combination of the existing local
    # Fourier observables, O=(A,L). No replacement state or extra mark.
    # X1^T G1=0 and G_r conjugate(Y_r)=0 follow from rank-one residues.
    beta = c0[1, 0] / c0[0, 0]
    raw_readouts = {
        "A": np.array([1, 0], dtype=complex),
        "L": np.array([0, 1], dtype=complex),
        "J": np.array([-beta, 1], dtype=complex),
        "X1": np.array([1, -residues[0][0, 0] / residues[0][1, 0]], dtype=complex),
        "Y1": np.array([1, -np.conj(residues[0][0, 0] / residues[0][0, 1])], dtype=complex),
        "Y2": np.array([1, -np.conj(residues[1][0, 0] / residues[1][0, 1])], dtype=complex),
    }
    definitions = {
        "A": ("A", "Original adjacent-pair Fourier readout; fixed microscopic indicator."),
        "L": ("L", "Original singleton/landing Fourier readout; fixed microscopic indicator."),
        "J": ("L-beta A; beta=C_LA(0)/C_AA(0)", "Equal-time-only landing innovation, E[J conjugate(A)]=0. No lag or spectrum used in this rule."),
        "X1": ("A-(G1_AA/G1_LA)L", "Exploratory source-side null of the slowest residue, X1^T G1=0."),
        "Y1": ("A-conjugate(G1_AA/G1_AL)L", "Exploratory future-readout null of the slowest residue, G1 conjugate(Y1)=0."),
        "Y2": ("A-conjugate(G2_AA/G2_AL)L", "Exploratory future-readout null of the second-slowest residue, G2 conjugate(Y2)=0."),
    }
    readout_info, normalized = {}, {}
    for name, raw in raw_readouts.items():
        variance = covariance(raw, c0, raw)
        if variance.real <= 0:
            raise ArithmeticError(f"Nonpositive physical readout variance: {name}")
        normalized[name] = raw / np.sqrt(variance.real)
        readout_info[name] = {
            "formula": definitions[name][0], "meaning": definitions[name][1],
            "raw_coefficients_A_L": [pair(value) for value in raw],
            "equal_time_variance_before_normalization": pair(variance),
            "unit_variance_coefficients_A_L": [pair(value) for value in normalized[name]],
            "rule_status": "model_derived_exploratory_filter" if name.startswith(("X", "Y")) else "fixed_microscopic_or_equal_time_rule",
        }
    channel_definitions = [
        ("AP_auto", "A", "A"),
        ("landing_auto", "L", "L"),
        ("AP_to_landing", "A", "L"),
        ("landing_to_AP", "L", "A"),
        ("equal_time_landing_innovation_auto", "J", "J"),
        ("slow1_source_null_auto", "X1", "X1"),
        ("slow1_biorthogonal_cross", "X1", "Y1"),
        ("split_slow12_null_cross", "X1", "Y2"),
    ]
    channels = []
    for name, source_name, readout_name in channel_definitions:
        source, readout = normalized[source_name], normalized[readout_name]
        projected_residues = [covariance(source, residue, readout) for residue in residues]
        lag_rows = []
        for lag, matrix in enumerate(correlations):
            full = covariance(source, matrix, readout)
            terms = [coefficient * eigenvalue.conjugate() ** lag
                     for coefficient, eigenvalue in zip(projected_residues, eigenvalues)]
            top_two, reconstructed = sum(terms[:2]), sum(terms)
            signal = abs(full)
            lag_rows.append({
                "lag": lag,
                "unit_equal_time_variance_covariance": pair(full),
                "absolute_signal_in_unit_variance_units": float(signal),
                "two_slowest_prediction": pair(top_two),
                "two_slowest_absolute_error_in_unit_variance_units": float(abs(full - top_two)),
                "two_slowest_relative_error": float(abs(full - top_two) / signal) if signal else None,
                "all_modes_reconstruction_absolute_error": float(abs(full - reconstructed)),
                "all_modes_reconstruction_relative_error": float(abs(full - reconstructed) / signal) if signal else None,
                "individual_mode_contributions": [pair(value) for value in terms],
            })
        channels.append({
            "name": name, "source": source_name, "readout": readout_name,
            "projected_residues_in_unit_variance_units": [pair(value) for value in projected_residues],
            "lag_rows": lag_rows,
            "interpretation": "The first two modes are suppressed on opposite ends by construction; the nonzero remaining absolute signal is the result, not the near-100-percent truncation error." if name == "split_slow12_null_cross" else
                              "The slowest mode is suppressed by a model-derived numerical filter, not by a previously proved physical symmetry." if name.startswith("slow1") else
                              "Named physical readout; no lag-dependent coefficient selection.",
        })
    indexed = {channel["name"]: channel for channel in channels}
    source_null = np.linalg.norm(normalized["X1"] @ residues[0])
    future_null1 = np.linalg.norm(residues[0] @ normalized["Y1"].conj())
    future_null2 = np.linalg.norm(residues[1] @ normalized["Y2"].conj())
    elapsed = time.perf_counter() - started
    result = {
        "schema": "matching-one.p398-fixed-readout.v1",
        "status": "completed_numerical_readout_analysis_of_saved_width5_model",
        "sources": [
            {"path": str(args.exact_input.relative_to(ROOT)), "sha256": exact_sha, "schema": exact["schema"]},
            {"path": str(args.mode_input.relative_to(ROOT)), "sha256": mode_sha, "schema": modes["schema"]},
        ],
        "definition": {
            "model": "width5, Q=1, h=v=1/2, positive discrete T=HV, first C5 character; existing AP A and singleton/landing L only",
            "physical_observable": "O_w=w_A A+w_L L; real and imaginary parts are real, bounded functions of the same physical frontier configuration",
            "source_readout": "S_xy(d)=x^T C(d) conjugate(y); each of x,y separately normalized to stationary variance one",
            "absolute_signal": "|S_xy(d)|/sqrt(Var(X)Var(Y)), already realized by the unit-variance coefficients; this is not a significance or sample-size estimate",
            "two_mode_prediction": "Reuse the saved residues and two largest-|lambda| modes with no re-fit",
            "null_selection": "One fixed algebraic rule from G1/G2, applied before lag errors; source and future nulls can differ because residues need not be Hermitian",
            "spectral_status": "Numerical complex128 projection of a saved exact physical model. No exact algebraic null or independent validation is claimed for the rounded filter coefficients.",
        },
        "readouts": readout_info,
        "channels": channels,
        "numerical_diagnostics": {
            "inherited_mode_diagnostics_not_recomputed": modes["numerical_diagnostics"],
            "source_X1_null_G1_norm": float(source_null),
            "future_Y1_null_G1_norm": float(future_null1),
            "future_Y2_null_G2_norm": float(future_null2),
            "max_projected_all_modes_absolute_reconstruction_error": max(row["all_modes_reconstruction_absolute_error"] for channel in channels for row in channel["lag_rows"]),
            "split_null_max_relative_all_modes_reconstruction_error_lags_1_to_8": max(row["all_modes_reconstruction_relative_error"] for row in indexed["split_slow12_null_cross"]["lag_rows"][1:]),
        },
        "new_remote_context_read_not_computed": {
            "commit": "552c45d7595ebcb0d04555cec03b2a5bfd8da44a",
            "ref": "theory/p398-width8-source-spectrum-20260831",
            "path": "notes/p398-width8-source-spectrum.md",
            "relation": "Already has duality-protected two rays and 93 propagation dimensions per ray, but for width8 continuous G=sum(J-I)+sum(D-I) with fixed i character. The width5 discrete T=HV numerical filters here are not those symmetry-protected rays.",
        },
        "boundaries": [
            "Exact width5 propagation rank eight, an effective two-slowest-mode whole-matrix approximation, and a readout-specific fast-mode filter are mutually compatible.",
            "Exploratory mode-null filters are algebraic consequences of this same model, not independent evidence, held-out validation, or a universal fixed microscopic symmetry.",
            "The split slow1/slow2 filter loses its leading signal under two-mode truncation by construction; its absolute normalized surviving signal quantifies the physical magnitude of that loss.",
            "Normalized covariance amplitudes do not specify experimental detectability without a sampling covariance and coefficient-calibration uncertainty.",
            "No continuum-field count, full-Q lifting result, energy-operator identification, or claim of non-Markovianity of the underlying 42-state process follows.",
        ],
        "execution": {
            "command": " ".join([sys.executable, "scripts/p398_fixed_readout_analysis.py", *sys.argv[1:]]),
            "python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine(),
            "elapsed_seconds": elapsed, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "transfer_recomputed": False, "eigenvalues_recomputed": False, "new_monte_carlo_samples": 0,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    output_sha = hashlib.sha256(args.output_json.read_bytes()).hexdigest()
    overview = "\n".join(
        f"| {channel['source']} → {channel['readout']} | {channel['lag_rows'][1]['absolute_signal_in_unit_variance_units']:.9g} | "
        f"{100 * channel['lag_rows'][1]['two_slowest_relative_error']:.7g}% | "
        f"{channel['lag_rows'][4]['absolute_signal_in_unit_variance_units']:.9g} | "
        f"{100 * channel['lag_rows'][4]['two_slowest_relative_error']:.7g}% |"
        for channel in channels
    )
    coefficient_rows = "\n".join(
        f"| {name} | `{display(raw[0])}` | `{display(raw[1])}` | "
        f"{readout_info[name]['equal_time_variance_before_normalization'][0]:.10g} |"
        for name, raw in raw_readouts.items()
    )
    double = indexed["split_slow12_null_cross"]["lag_rows"]
    double_rows = "\n".join(f"| {row['lag']} | {row['absolute_signal_in_unit_variance_units']:.12g} | "
                             f"{100 * row['two_slowest_relative_error']:.9g}% | "
                             f"{row['all_modes_reconstruction_relative_error']:.4g} |" for row in double[1:])
    ap = indexed["AP_auto"]["lag_rows"][1]
    landing = indexed["landing_auto"]["lag_rows"][1]
    innovation_channel = indexed["equal_time_landing_innovation_auto"]
    innovation = innovation_channel["lag_rows"][1]
    innovation_d4 = innovation_channel["lag_rows"][4]
    null1 = indexed["slow1_biorthogonal_cross"]["lag_rows"][1]
    report = f"""# P398：固定 AP/landing 读出与可见快模

## 直接结果

width5 的两慢模近似是**观察者相关**的，而非只由精确传播秩决定。
同一个 d=1，整矩阵误差约 0.8996%，但固定 AP 自相关误差为
**{100 * ap['two_slowest_relative_error']:.7g}%**，固定 landing 自相关为
**{100 * landing['two_slowest_relative_error']:.7g}%**。
只按等时 C(0) 去掉 AP 共线部分的 landing 创新 J，其误差为
**{100 * innovation['two_slowest_relative_error']:.7g}%**，归一绝对信号
`{innovation['absolute_signal_in_unit_variance_units']:.10g}`。
这个 J 不读取非零距离或谱来选择系数，是比全矩阵范数更具体的现成物理读出。

J 的第三模（λ≈{eigenvalues[2].real:.12g}）归一 residue 为 {decode_pair(innovation_channel['projected_residues_in_unit_variance_units'][2]).real:.12g}，
d=1 的实际贡献为 {decode_pair(innovation['individual_mode_contributions'][2]).real:.12g}，而完整信号为 {innovation['absolute_signal_in_unit_variance_units']:.12g}。
因此近距离快模可见性并非由双端 null 规则强制造成。
J 也不是精确 slowest-null：其第一模 residue 约 {decode_pair(innovation_channel['projected_residues_in_unit_variance_units'][0]).real:.9g}，足够远时会重新占优。
d=4 的 {100 * innovation_d4['two_slowest_relative_error']:.4g}% 相对误差发生在抵消后的 {innovation_d4['absolute_signal_in_unit_variance_units']:.3g} 微小信号上，不能当作大绝对效应。
这说明同一个具名读出的有效模组成还依赖距离。

再按既有谱系数构造 slowest-null 的 source/readout X1→Y1，d=1 的误差为
**{100 * null1['two_slowest_relative_error']:.7g}%**，绝对信号
`{null1['absolute_signal_in_unit_variance_units']:.10g}`。
让 X1 在源端消去 mode1、Y2 在读出端消去 mode2，则 d=1 剩余归一绝对信号为
**{double[1]['absolute_signal_in_unit_variance_units']:.10g}**，d=4 为
`{double[4]['absolute_signal_in_unit_variance_units']:.10g}`。
两模截断会丢掉该剩余信号，但这是滤波规则的构造结果，**不是新的模型否证或独立验证**。

## 固定规则，不按距离误差优化

全程只消费既有 `C(d)` 和八个 residue，不重算 transfer、本征值或平稳分布。
物理量仍为原来的 adjacent-pair A 与 singleton/landing L；所有新组合的实部、虚部
都是同一有限正权模型上的有界配置函数，没有额外 mark。

以下原始系数以 `O=w_A A+w_L L` 约定；计算时每个源与读出分别除以其等时标准差。
所以表中的绝对信号是 `|E[X(0) conjugate(Y(d))]|/sqrt(Var X Var Y)`，
而不是被放大系数制造出的任意振幅。

| 名称 | A 系数 | L 系数 | 原始等时方差 |
|---|---|---|---:|
{coefficient_rows}

- `J=L-beta A`，`beta=C_LA(0)/C_AA(0)`；因此与 A 等时正交。
- `X1=A-(G1_AA/G1_LA)L`，满足 source 侧 `X1ᵀG1≈0`。
- `Yr=A-conjugate(Gr_AA/Gr_AL)L`，满足 future 侧 `Gr conjugate(Yr)≈0`，r=1,2。
- X/Y 是**模型导出的探索性滤波器**，不是事先独立验证的观测，也没有挑选最有利的 lag。
  条件中的近似号表示保存谱系数与滤波系数为 complex128；未声称 exact algebraic null。
- source-null 与 readout-null 一般不同；把它们强行写成同一个实系数会丢失本模型的有序读出信息。

## 具名读出的实际误差与绝对信号

两慢模固定为已存 `|λ|` 最大的两项，residue 不重新拟合。

| 归一源 → 读出 | d1 绝对信号 | d1 两模误差 | d4 绝对信号 | d4 两模误差 |
|---|---:|---:|---:|---:|
{overview}

全套 d=0…8 复数 covariance、每模 residue、每项贡献和绝对/相对遗漏在
[latest.json](latest.json)。AL 与 LA 各自保留，没有假设有限过去边界过程可逆。

## 双端 slow1/slow2 消去：真正剩下多少信号？

X1→Y2 的两个最大模分别在两端被消去。故近 100% 截断误差本身是恒等性，
不是科学惊喜；下列绝对量与实际数值精度才描述剩余快模是否微小。

| d | 归一绝对信号 | 两模相对遗漏 | 八模数值重建相对误差 |
|---|---:|---:|---:|
{double_rows}

滤波后远距离信号很小；这不是采样显著性或所需样本量的估计。
真实采集还会受配置采样 covariance 和系数校准误差影响，本次没有制造这些未存信息。

## 三层对象应当分开

1. **精确有限状态传播：**原结果的正间距 Hankel 秩为 8；底层 42-state T 仍为 Markov。
2. **有效慢模近似：**两慢模可很好重建整体矩阵，但固定创新读出与消去投影的误差不同。
3. **连续极限场：**本计算没有识别场数、Jordan 身份、full-Q 物理提升或 norm-4 能量投影。

窄幅核对远端时已读到 `552c45d7595ebcb0d04555cec03b2a5bfd8da44a`
（`theory/p398-width8-source-spectrum-20260831`）的
`notes/p398-width8-source-spectrum.md`：它在 width8 连续生成元、固定 i 字符下
已经得到 Kreweras 保护的两 ray、每 ray 93 个传播方向及具名投影。
这里是 width5 离散 T=HV、ζ₅ 字符，**不是重做该结果**，也不把本次数值 null
冒充相同的对偶保护。

## 数值与来源

- X1 source-null residue 范数：`{source_null:.5g}`；Y1/Y2 future-null 范数：
  `{future_null1:.5g}` / `{future_null2:.5g}`。
- 全部投影的八模重建最大绝对误差：
  `{result['numerical_diagnostics']['max_projected_all_modes_absolute_reconstruction_error']:.5g}`。
- 原模式分解 condition number 与 residual 原样引用，不重算谱；本输出是数值解释，非新 exact certificate。
- 原精确 JSON SHA256：`{exact_sha}`。
- 原 mode-visibility JSON SHA256：`{mode_sha}`；两份输入均未改动。
- 本输出 JSON SHA256：`{output_sha}`。
- Python {platform.python_version()} / NumPy {np.__version__} / {platform.machine()}，{elapsed:.4f} 秒；
  无 MC、无测试套件、无额外 width/参数扫描。

```bash
{sys.executable} scripts/p398_fixed_readout_analysis.py
```
"""
    args.output_md.write_text(report, encoding="utf-8")
    print(json.dumps({
        "d1_summary": [{"channel": channel["name"],
                        "absolute_signal": channel["lag_rows"][1]["absolute_signal_in_unit_variance_units"],
                        "two_mode_relative_error": channel["lag_rows"][1]["two_slowest_relative_error"]}
                       for channel in channels],
        "split_null_d4_signal": double[4]["absolute_signal_in_unit_variance_units"],
        "split_null_d8_signal": double[8]["absolute_signal_in_unit_variance_units"],
        "diagnostics": result["numerical_diagnostics"],
        "output_json_sha256": output_sha,
    }, indent=2))


if __name__ == "__main__":
    main()
