#!/usr/bin/env python3
"""Numerical mode visibility from P398's existing exact width-five JSON only.

No transfer construction, stationary solve, exact-engine rerun, fitting, scan,
or simulation. The two retained modes are fixed by largest |eigenvalue|.
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


def field_value(value):
    zeta = np.exp(2j * np.pi / 5)
    return sum(float(Fraction(coefficient)) * zeta ** power
               for power, coefficient in enumerate(value["coefficients"]))


def field_matrix(rows):
    return np.array([[field_value(value) for value in row] for row in rows], dtype=np.complex128)


def complex_pair(value):
    return [float(np.real(value)), float(np.imag(value))]


def matrix_pairs(matrix):
    return [[complex_pair(value) for value in row] for row in matrix]


def scientific(value):
    return f"{value:.6g}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "results/p398-physical-two-point/latest.json")
    parser.add_argument("--output-json", type=Path, default=ROOT / "results/p398-physical-two-point/mode-visibility.json")
    parser.add_argument("--output-md", type=Path, default=ROOT / "results/p398-physical-two-point/mode-visibility.md")
    args = parser.parse_args()
    started = time.perf_counter()
    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    source_sha = hashlib.sha256(source_bytes).hexdigest()
    scientific_source = source["scientific_result"]
    block = field_matrix(scientific_source["backward_block"])
    readouts = field_matrix(scientific_source["readout_vectors_A_L"]).T
    weights = np.array([float(Fraction(value)) for value in source["physical_model"]["charged_orbit_stationary_weights"]])
    stored_correlations = [field_matrix(value) for value in scientific_source["correlations_by_lag_0_to_8"]]

    eigenvalues, right_vectors = np.linalg.eig(block)
    order = np.argsort(-np.abs(eigenvalues), kind="stable")
    eigenvalues = eigenvalues[order]
    right_vectors = right_vectors[:, order]
    modal_readouts = np.linalg.solve(right_vectors, readouts)
    left_readout_factors = readouts.T @ (weights[:, None] * right_vectors.conj())
    residues = [np.outer(left_readout_factors[:, mode], modal_readouts[mode, :].conj())
                for mode in range(len(eigenvalues))]
    # C(d)=R^T W conjugate(B^d R), hence conjugate(lambda)^d.
    mode_contributions = [[residue * eigenvalue.conjugate() ** lag
                           for residue, eigenvalue in zip(residues, eigenvalues)]
                          for lag in range(9)]
    errors = []
    for lag, (terms, exact_numeric) in enumerate(zip(mode_contributions, stored_correlations)):
        all_modes = np.sum(terms, axis=0)
        top_two = np.sum(terms[:2], axis=0)
        norm = float(np.linalg.norm(exact_numeric, "fro"))
        errors.append({
            "lag": lag,
            "stored_correlation_frobenius_norm": norm,
            "two_slowest_relative_frobenius_error": float(np.linalg.norm(top_two - exact_numeric, "fro") / norm),
            "all_modes_relative_reconstruction_error": float(np.linalg.norm(all_modes - exact_numeric, "fro") / norm),
            "two_slowest_correlation": matrix_pairs(top_two),
            "all_modes_correlation": matrix_pairs(all_modes),
        })
    block_norm = np.linalg.norm(block, 2)
    eigen_residual_matrix = block @ right_vectors - right_vectors * eigenvalues[None, :]
    relative_eigen_residuals = [float(np.linalg.norm(eigen_residual_matrix[:, i]) /
                                    ((block_norm + abs(eigenvalues[i])) * np.linalg.norm(right_vectors[:, i])))
                              for i in range(len(eigenvalues))]
    condition = float(np.linalg.cond(right_vectors, 2))
    reconstructed_block = (right_vectors * eigenvalues[None, :]) @ np.linalg.inv(right_vectors)
    block_error = float(np.linalg.norm(block - reconstructed_block, "fro") / np.linalg.norm(block, "fro"))
    modes = []
    for i, (eigenvalue, residue) in enumerate(zip(eigenvalues, residues)):
        modes.append({
            "mode": i + 1,
            "rank_by_eigenvalue_modulus": i + 1,
            "eigenvalue": complex_pair(eigenvalue),
            "eigenvalue_modulus": float(abs(eigenvalue)),
            "inverse_decay_rate_rows": float(-1 / np.log(abs(eigenvalue))),
            "residue_A_L_rows_columns": matrix_pairs(residue),
            "residue_frobenius_norm": float(np.linalg.norm(residue, "fro")),
            "relative_eigen_residual": relative_eigen_residuals[i],
            "individual_contribution_frobenius_relative_to_C_by_lag_0_to_8": [
                float(np.linalg.norm(mode_contributions[lag][i], "fro") / errors[lag]["stored_correlation_frobenius_norm"])
                for lag in range(9)
            ],
        })
    u1 = np.linalg.solve(stored_correlations[0], stored_correlations[1])
    u2 = np.linalg.solve(stored_correlations[0], stored_correlations[2])
    semigroup_defect = u2 - u1 @ u1
    relative_semigroup_defect = float(np.linalg.norm(semigroup_defect, "fro") / np.linalg.norm(u2, "fro"))
    c0_condition = float(np.linalg.cond(stored_correlations[0], 2))
    elapsed = time.perf_counter() - started
    result = {
        "schema": "matching-one.p398-width5-mode-visibility.v1",
        "status": "completed_numerical_interpretation_of_existing_exact_model",
        "source": {
            "path": str(args.input.relative_to(ROOT)) if args.input.is_relative_to(ROOT) else str(args.input),
            "sha256": source_sha,
            "schema": source["schema"],
            "exact_positive_separation_hankel_rank": scientific_source["positive_separation_block_hankel"]["rank"],
            "exact_characteristic_square_free": scientific_source["characteristic_square_free"],
            "unchanged_input": True,
        },
        "definition": {
            "scope": "width=5, Q=1, h=v=1/2, the existing positive past-connectivity A/L readouts only",
            "precision": "NumPy complex128; eigenspectrum and residues are numerical, not additional exact certificates",
            "R": "8x2 matrix whose columns are A and L on the charged orbit representatives",
            "W": "diag(charged_orbit_stationary_weights)",
            "eigen_decomposition": "B=V diag(lambda) V^-1",
            "residue": "G_r=(R^T W conjugate(V))[:,r] outer conjugate((V^-1 R)[r,:])",
            "correlation": "C(d)=sum_r G_r conjugate(lambda_r)^d",
            "two_mode_selection": "The two eigenvalues of largest modulus, fixed before reading approximation errors; no fitted amplitudes or optimized pair.",
            "relative_error": "||C_two(d)-C_stored_exact_embedded(d)||_F / ||C_stored_exact_embedded(d)||_F",
            "contribution_norm_boundary": "Individual Frobenius contribution norms are not probabilities and need not sum to one; complex matrix residues can cancel.",
        },
        "numerical_diagnostics": {
            "eigenvector_matrix_condition_number_2": condition,
            "max_relative_eigen_residual": max(relative_eigen_residuals),
            "relative_block_reconstruction_error_frobenius": block_error,
            "max_all_modes_correlation_reconstruction_error_lags_0_to_8": max(row["all_modes_relative_reconstruction_error"] for row in errors),
        },
        "modes": modes,
        "lag_errors": errors,
        "two_readout_semigroup": {
            "definition": "U_d=C(0)^-1 C(d); defect=U_2-U_1^2",
            "U1": matrix_pairs(u1), "U2": matrix_pairs(u2),
            "U1_squared": matrix_pairs(u1 @ u1),
            "defect": matrix_pairs(semigroup_defect),
            "relative_frobenius_defect_over_U2": relative_semigroup_defect,
            "C0_condition_number_2": c0_condition,
            "interpretation": "Nonzero defect rules out one time-homogeneous linear 2x2 propagator reproducing C(0), C(1), C(2). The underlying 42-state transfer remains Markov; this does not establish non-Markovianity of every nonlinear process built from the A/L pair.",
        },
        "interpretation_boundaries": [
            "Exact propagation rank eight and a good effective two-slowest-mode approximation are compatible.",
            "A small Frobenius error for the whole 2x2 matrix does not guarantee the same relative error in a cancellation-sensitive projected observable.",
            "The eigenvalue ordering is about decay in rows of this finite cylinder, not the count or dimensions of continuum fields.",
            "This is not a full-Q lift, a new width/anisotropy scan, or an energy-operator identification for Matching One.",
        ],
        "execution": {
            "command": " ".join([sys.executable, "scripts/p398_width5_mode_visibility.py", *sys.argv[1:]]),
            "python": platform.python_version(), "numpy": np.__version__,
            "machine": platform.machine(), "elapsed_seconds": elapsed,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "exact_engine_rerun": False, "new_monte_carlo_samples": 0,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    output_sha = hashlib.sha256(args.output_json.read_bytes()).hexdigest()
    mode_rows = "\n".join(
        f"| {row['mode']} | {row['eigenvalue'][0]:.12g} {row['eigenvalue'][1]:+.3g}i | "
        f"{row['residue_frobenius_norm']:.9g} | {row['individual_contribution_frobenius_relative_to_C_by_lag_0_to_8'][1]:.7g} | "
        f"{row['individual_contribution_frobenius_relative_to_C_by_lag_0_to_8'][8]:.7g} |"
        for row in modes
    )
    error_rows = "\n".join(f"| {row['lag']} | {row['two_slowest_relative_frobenius_error']:.10g} | "
                           f"{100 * row['two_slowest_relative_frobenius_error']:.8g}% | "
                           f"{row['all_modes_relative_reconstruction_error']:.3g} |" for row in errors[1:])
    report = f"""# P398 width-5：八个精确传播成分是否仍有有效两慢模？

## 直接答案

精确传播秩 8 和近似两慢模可以同时成立。
只保留模最大的两个本征值（约 `{eigenvalues[0].real:.12g}`、`{eigenvalues[1].real:.12g}`），
不用重新拟合 residue，完整 2×2 两点矩阵的 Frobenius 相对误差在
d=1 为 **{100 * errors[1]['two_slowest_relative_frobenius_error']:.8g}%**，
d=2 为 **{100 * errors[2]['two_slowest_relative_frobenius_error']:.8g}%**，
d=4 为 **{100 * errors[4]['two_slowest_relative_frobenius_error']:.8g}%**，
d=8 为 **{100 * errors[8]['two_slowest_relative_frobenius_error']:.8g}%**。
这里的误差来自实际 A/L 矩阵，而不是仅看本征值的大小。

这是对[原精确 JSON](latest.json)的 numerical 解释；没有重跑精确引擎，
没有做新 width、参数扫描、MC 或自由幅度拟合。
原结果的 `rank=8` 与 square-free 多项式原样保留。

## 每个传播模的 residue

记 `R=(A,L)`、`W=diag(5π_orbit)`、`B=V diag(λ) V⁻¹`。

`G_r=(RᵀW conjugate(V))[:,r] outer conjugate((V⁻¹R)[r,:])`，
`C(d)=Σ_r G_r conjugate(λ_r)^d`。

模按 `|λ|` 降序固定排序；取前两个，没有根据误差重新选对。
全部 2×2 复数 residue 保存在 [mode-visibility.json](mode-visibility.json)。
下表的 contribution 是单项 `||G_r conjugate(λ_r)^d||_F/||C(d)||_F`，
不是概率、不能当作百分比分账，因为不同矩阵项可能相消。

| 模 | λ（数值） | residue 的 Frobenius 范数 | d=1 contribution | d=8 contribution |
|---|---:|---:|---:|---:|
{mode_rows}

## 两慢模的实际矩阵误差

分母是精确 JSON 嵌入 complex128 后的真实 `||C(d)||_F`。

| d | 两慢模相对误差 | 百分比 | 全八模重建相对误差 |
|---|---:|---:|---:|
{error_rows}

## 两个读出是否本身构成二维自主传播？

以 `U_d=C(0)⁻¹ C(d)` 定义两点有效传播，数值结果为

`||U₂−U₁²||_F / ||U₂||_F = {relative_semigroup_defect:.12g}`
（**{100 * relative_semigroup_defect:.8g}%**）；`cond₂(C(0))={c0_condition:.9g}`。
U₁、U₂ 与完整 defect 矩阵已保存在 JSON。

这个 defect 直接说明：不存在一个固定 2×2 线性传播矩阵同时重建
C(0)、C(1)、C(2)。它不意味着底层 42-state 正权 transfer 失去 Markov 性，
也不证明所有由 A/L 构造的非线性过程必定非 Markov。
整体两慢模近似是否好与两个原始读出是否精确自主闭合，不是同一个判断。

## 数值状态与解释范围

- 本征向量矩阵 condition number（2-norm）：`{condition:.9g}`。
- 最大相对本征方程 residual：`{max(relative_eigen_residuals):.6g}`。
- B 的本征分解重建相对误差：`{block_error:.6g}`。
- d=0…8 全八模的两点重建最大相对误差：
  `{result['numerical_diagnostics']['max_all_modes_correlation_reconstruction_error_lags_0_to_8']:.6g}`。
- 所有谱、residue 与误差是 NumPy complex128 数值，不冒充新增 exact certificate。
- 整体 Frobenius 误差小，并不保证某个有精细抵消的投影读出也具有同样小的相对误差。
  本表支持的是指定两点矩阵的有效慢模近似，不是任意观察者都能只保留两场。
- 有限圆柱传播秩、慢模近似与连续极限场识别是三个不同问题。
  这组结果不识别 norm-4 的 E_top 微观能量投影，不推及全部宽度或 full-Q。

## 来源与复现

输入 SHA256：`{source_sha}`，原文件未改动。
本输出 JSON SHA256：`{output_sha}`。
Python {platform.python_version()} / NumPy {np.__version__} / {platform.machine()}，
本次轻量数值解释 {elapsed:.4f} 秒。

```bash
{sys.executable} scripts/p398_width5_mode_visibility.py
```
"""
    args.output_md.write_text(report, encoding="utf-8")
    print(json.dumps({"two_slowest_eigenvalues": [complex_pair(value) for value in eigenvalues[:2]],
                      "condition_number": condition,
                      "max_relative_eigen_residual": max(relative_eigen_residuals),
                      "relative_semigroup_defect": relative_semigroup_defect,
                      "two_mode_relative_errors_lags_1_to_8": [row["two_slowest_relative_frobenius_error"] for row in errors[1:]],
                      "input_sha256": source_sha, "output_json_sha256": output_sha}, indent=2))


if __name__ == "__main__":
    main()
