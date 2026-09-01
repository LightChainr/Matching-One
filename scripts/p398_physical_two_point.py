#!/usr/bin/env python3
"""One exact width-five extension of P398's positive-cylinder A/L readouts.

No simulation, parameter scan, new marks, or signed-module interpretation.
All transfer, stationary, closure and correlation ranks use exact arithmetic.
The field basis is (1,z,z^2,z^3), z^4+z^3+z^2+z+1=0.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import cmath
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

from noncrossing_connectivity_codec import canonical_rgs, noncrossing_states


ROOT = Path(__file__).resolve().parents[1]
WIDTH = 5
DENOMINATOR = 2 ** (2 * WIDTH)
SOURCES = [
    {
        "commit": "e38fe7634354b0cb2201fa55fd9b4d37ccedeef2",
        "ref": "origin/theory/p398-positive-cylinder-propagation-20260831",
        "path": "notes/p398-positive-cylinder-propagation.md",
        "status": "branch_only",
        "role": "Completed width-four positive measure, A/L emissions and two-point propagation.",
    },
    {
        "commit": "b35e100a3903c706dceba57c4667386eb4510ac3",
        "ref": "origin/theory/p398-anisotropic-cylinder-spectrum-20260831",
        "path": "notes/p398-anisotropic-cylinder-spectrum.md",
        "status": "branch_only",
        "role": "Completed width-four positive h/v spectrum; no repeat scan here.",
    },
]


def rational_solve(matrix, rhs):
    """Small rational solve, also used to invert cyclotomic elements."""
    rows = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix, rhs)]
    n = len(rows)
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        factor = rows[col][col]
        rows[col] = [x / factor for x in rows[col]]
        for i in range(col + 1, n):
            factor = rows[i][col]
            if factor:
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[col])]
    out = [F(0)] * n
    for i in range(n - 1, -1, -1):
        out[i] = rows[i][-1] - sum(rows[i][j] * out[j] for j in range(i + 1, n))
    return out


def integer_solve(matrix, rhs):
    """Fraction-free elimination of the single finite stationary system."""
    rows = [list(row) + [value] for row, value in zip(matrix, rhs)]
    n = len(rows)
    previous = 1
    for col in range(n - 1):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        value = rows[col][col]
        for i in range(col + 1, n):
            for j in range(col + 1, n + 1):
                numerator = value * rows[i][j] - rows[i][col] * rows[col][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Nonexact Bareiss division")
                rows[i][j] = quotient
            rows[i][col] = 0
        previous = value
    out = [F(0)] * n
    for i in range(n - 1, -1, -1):
        out[i] = (F(rows[i][-1]) - sum(rows[i][j] * out[j] for j in range(i + 1, n))) / rows[i][i]
    return out


@dataclass(frozen=True)
class K:
    c: tuple[F, F, F, F]

    @staticmethod
    def scalar(value):
        return value if isinstance(value, K) else K((F(value), F(0), F(0), F(0)))

    def __bool__(self):
        return any(self.c)

    def __add__(self, other):
        other = K.scalar(other)
        return K(tuple(x + y for x, y in zip(self.c, other.c)))

    __radd__ = __add__

    def __neg__(self):
        return K(tuple(-x for x in self.c))

    def __sub__(self, other):
        return self + (-K.scalar(other))

    def __rsub__(self, other):
        return K.scalar(other) + (-self)

    def __mul__(self, other):
        other = K.scalar(other)
        coefficients = [F(0)] * 7
        for i, x in enumerate(self.c):
            if x:
                for j, y in enumerate(other.c):
                    if y:
                        coefficients[i + j] += x * y
        for power in range(6, 3, -1):
            value = coefficients[power]
            if value:
                for j in range(4):
                    coefficients[power - 4 + j] -= value
        return K(tuple(coefficients[:4]))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = K.scalar(other)
        return self * inverse(other)

    def conjugate(self):
        return sum((value * POWERS[(-i) % WIDTH] for i, value in enumerate(self.c)), ZERO)

    def approximate(self):
        z = cmath.exp(2j * cmath.pi / WIDTH)
        return sum(float(value) * z ** i for i, value in enumerate(self.c))


ZERO = K.scalar(0)
ONE = K.scalar(1)
ZETA = K((F(0), F(1), F(0), F(0)))
POWERS = [ONE]
for _ in range(1, WIDTH):
    POWERS.append(POWERS[-1] * ZETA)


@lru_cache(maxsize=None)
def inverse(value):
    if not value:
        raise ZeroDivisionError("Zero cyclotomic inverse")
    columns = [(value * POWERS[j]).c for j in range(4)]
    matrix = [[columns[j][i] for j in range(4)] for i in range(4)]
    return K(tuple(rational_solve(matrix, [1, 0, 0, 0])))


def matrix_product(left, right):
    return [[sum((x * right[k][j] for k, x in enumerate(row) if x), ZERO)
             for j in range(len(right[0]))] for row in left]


def matrix_vector(matrix, vector):
    return [sum((x * y for x, y in zip(row, vector) if x and y), ZERO) for row in matrix]


def rank(matrix):
    rows = [list(row) for row in matrix]
    pivot_row = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_inverse = inverse(rows[pivot_row][col])
        rows[pivot_row] = [x * pivot_inverse for x in rows[pivot_row]]
        for i in range(pivot_row + 1, len(rows)):
            factor = rows[i][col]
            if factor:
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def column_rank(columns):
    return rank(list(map(list, zip(*columns))))


def characteristic_polynomial(matrix):
    """Descending coefficients by exact Faddeev--LeVerrier."""
    n = len(matrix)
    work = [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]
    coefficients = [ONE]
    for k in range(1, n + 1):
        work = matrix_product(matrix, work)
        coefficient = -sum((work[i][i] for i in range(n)), ZERO) / k
        coefficients.append(coefficient)
        for i in range(n):
            work[i][i] += coefficient
    return coefficients


def polynomial_gcd(left, right):
    """Ascending coefficient Euclidean algorithm over Q(zeta_5)."""
    def trim(poly):
        while poly and not poly[-1]:
            poly.pop()
        return poly

    left, right = trim(list(left)), trim(list(right))
    while right:
        remainder = list(left)
        leading_inverse = inverse(right[-1])
        while len(remainder) >= len(right):
            shift = len(remainder) - len(right)
            coefficient = remainder[-1] * leading_inverse
            for i, value in enumerate(right):
                remainder[i + shift] -= coefficient * value
            trim(remainder)
        left, right = right, remainder
    if not left:
        return []
    leading_inverse = inverse(left[-1])
    return [value * leading_inverse for value in left]


def rotate(state):
    return canonical_rgs((state[-1],) + state[:-1])


def detach(state, site):
    labels = list(state)
    labels[site] = max(state) + 1
    return canonical_rgs(labels)


def join(state, site):
    left, right = state[site], state[(site + 1) % WIDTH]
    return canonical_rgs(left if x == right else x for x in state)


def emissions(state):
    adjacent = sum((POWERS[j] for j in range(WIDTH) if state[j] == state[(j + 1) % WIDTH]), ZERO)
    singleton = sum((POWERS[j] for j in range(WIDTH) if state.count(state[j]) == 1), ZERO)
    return adjacent, singleton


def construct():
    states = noncrossing_states(WIDTH)
    index = {state: i for i, state in enumerate(states)}
    n = len(states)
    counts = [[0] * n for _ in states]
    for source_index, source in enumerate(states):
        distribution = Counter({source: 1})
        for operation in (detach, join):
            for site in range(WIDTH):
                following = Counter()
                for state, multiplicity in distribution.items():
                    following[state] += multiplicity
                    following[operation(state, site)] += multiplicity
                distribution = following
        for target, multiplicity in distribution.items():
            counts[index[target]][source_index] = multiplicity
    stationary_matrix = [[counts[i][j] - DENOMINATOR * (i == j) for j in range(n)] for i in range(n)]
    stationary_matrix[-1] = [1] * n
    stationary = integer_solve(stationary_matrix, [0] * (n - 1) + [1])

    seen, orbits = set(), []
    for state in states:
        if state in seen:
            continue
        orbit, member = [], state
        while member not in orbit:
            orbit.append(member)
            seen.add(member)
            member = rotate(member)
        orbits.append(orbit)
    charged_orbits = [orbit for orbit in orbits if len(orbit) == WIDTH]
    representatives = [orbit[0] for orbit in charged_orbits]
    backward = [[sum((F(counts[index[state]][index[representative]], DENOMINATOR) * POWERS[k]
                      for k, state in enumerate(orbit)), ZERO)
                 for orbit in charged_orbits] for representative in representatives]
    readouts = [list(values) for values in zip(*(emissions(state) for state in representatives))]
    orbit_weights = [WIDTH * stationary[index[state]] for state in representatives]
    exact_conditions = {
        "columns_stochastic": all(sum(column) == DENOMINATOR for column in zip(*counts)),
        "nonnegative_transfer": all(value >= 0 for row in counts for value in row),
        "stationary_normalized": sum(stationary) == 1,
        "stationary_full_support": all(value > 0 for value in stationary),
        "stationary_residual_zero": all(sum(counts[i][j] * stationary[j] for j in range(n)) == DENOMINATOR * stationary[i] for i in range(n)),
        "rotation_invariant_stationary": all(stationary[index[rotate(state)]] == stationary[index[state]] for state in states),
        "transfer_rotation_equivariant": all(counts[index[rotate(target)]][index[rotate(source)]] == counts[index[target]][index[source]] for target in states for source in states),
        "rotation_covariant_emissions": all(emissions(rotate(state)) == tuple(ZETA * x for x in emissions(state)) for state in states),
        "zero_means": all(sum((stationary[index[state]] * emissions(state)[a] for state in states), ZERO) == ZERO for a in range(2)),
    }
    # This identity certifies the orbit reduction on every physical state,
    # not merely on the representative rows used to construct its matrix.
    images = [matrix_vector(backward, readout) for readout in readouts]
    exact_conditions["backward_readout_reduction_on_all_states"] = all(
        sum((F(counts[index[target]][index[state]], DENOMINATOR) * emissions(target)[a] for target in states), ZERO)
        == POWERS[k] * images[a][i]
        for i, orbit in enumerate(charged_orbits) for k, state in enumerate(orbit) for a in range(2)
    ) and all(
        sum((F(counts[index[target]][index[orbit[0]]], DENOMINATOR) * emissions(target)[a] for target in states), ZERO) == ZERO
        for orbit in orbits if len(orbit) == 1 for a in range(2)
    )
    if not all(exact_conditions.values()):
        raise ArithmeticError(f"Physical construction identity failed: {exact_conditions}")
    return states, counts, stationary, orbits, representatives, backward, readouts, orbit_weights, exact_conditions


def serialize(value):
    if isinstance(value, K):
        approximate = value.approximate()
        return {"coefficients": [str(x) for x in value.c],
                "approximate": [approximate.real, approximate.imag]}
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


def matrix_display(matrix):
    return "\n".join("  " + " | ".join(f"{value.approximate().real:.12g}{value.approximate().imag:+.12g}i" for value in row) for row in matrix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-json", type=Path, default=ROOT / "results/p398-physical-two-point/latest.json")
    parser.add_argument("--output-md", type=Path, default=ROOT / "results/p398-physical-two-point/REPORT.md")
    args = parser.parse_args()
    started = time.perf_counter()
    states, counts, stationary, orbits, representatives, backward, readouts, weights, identities = construct()
    print(f"width=5 physical states={len(states)}, charge-one dimension={len(backward)}; stationary solved", flush=True)
    powers = [readouts]
    for _ in range(8):
        powers.append([matrix_vector(backward, vector) for vector in powers[-1]])
    ranks = [column_rank([vector for pair in powers[:k + 1] for vector in pair]) for k in range(8)]
    correlations = [
        [[sum((weights[i] * readouts[a][i] * power[b][i].conjugate() for i in range(len(backward))), ZERO)
          for b in range(2)] for a in range(2)] for power in powers
    ]
    hankel = [[correlations[i + j + 1][a][b]
               for j in range(4) for b in range(2)] for i in range(4) for a in range(2)]
    hankel_rank = rank(hankel)
    characteristic = characteristic_polynomial(backward)
    ascending = list(reversed(characteristic))
    derivative = [degree * ascending[degree] for degree in range(1, len(ascending))]
    gcd = polynomial_gcd(ascending, derivative)
    determinant_c0 = correlations[0][0][0] * correlations[0][1][1] - correlations[0][0][1] * correlations[0][1][0]
    elapsed = time.perf_counter() - started
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    codec_path = ROOT / "scripts/noncrossing_connectivity_codec.py"
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result = {
        "schema": "matching-one.p398-physical-two-point.v1",
        "status": "completed_exact_width_five_positive_cylinder_analysis",
        "sources": SOURCES,
        "definition": {
            "width": WIDTH, "Q": "1", "h": "1/2", "v": "1/2",
            "measure": "Independent square-bond cylinder, stationary connectivity through the past; read immediately after H.",
            "transfer": "T=H V, H=product_i((I+J_i)/2), V=product_i((I+D_i(1))/2); columns are source states.",
            "A": "sum_j zeta_5^j 1{j connected to j+1 in the past-frontier partition}",
            "L": "sum_j zeta_5^j 1{j is a singleton in the past-frontier partition}",
            "extension_rule": "Keep the same local indicators and first nontrivial cyclic Fourier character. At width four zeta_4=i this is exactly the existing A/L pair; at width five the cyclic group and character differ.",
            "field_basis": ["1", "zeta_5", "zeta_5^2", "zeta_5^3"],
            "field_relation": "zeta_5^4+zeta_5^3+zeta_5^2+zeta_5+1=0",
            "correlation": "C_ab(d)=E[O_a(X_0) conjugate(O_b(X_d))], connected because both exact stationary means vanish.",
            "backward_block": "(T^t f)(representative_i)=sum_j B_ij f(representative_j), f(R^k representative_j)=zeta_5^k f(representative_j).",
            "not_the_signed_module": "This positive 42-state Q=1 transfer is not the older signed 23-dimensional retained-mark module, and is not a full-Q lift of that module.",
        },
        "physical_model": {
            "state_count": len(states), "states": states, "transfer_count_denominator": DENOMINATOR,
            "transfer_counts_columns_source": counts, "stationary": stationary,
            "rotation_orbit_sizes": [len(orbit) for orbit in orbits],
            "charged_representatives": representatives, "charged_orbit_stationary_weights": weights,
            "exact_construction_identities": identities,
        },
        "scientific_result": {
            "charge_one_dimension": len(backward), "backward_block": backward,
            "readout_vectors_A_L": readouts,
            "krylov_rank_by_max_power_0_to_7": ranks,
            "two_readout_span_closed": ranks[1] == ranks[0],
            "positive_separation_block_hankel": {"block_rows": 4, "block_columns": 4, "lags": "1..7", "rank": hankel_rank, "matrix": hankel},
            "characteristic_polynomial_descending": characteristic,
            "characteristic_derivative_gcd_ascending": gcd,
            "characteristic_square_free": len(gcd) == 1,
            "correlations_by_lag_0_to_8": correlations,
            "determinant_C0": determinant_c0,
            "identification_boundary": [
                "Ranks refer to exact finite-width propagation of the two specified physical readouts, not to a number of continuum fields.",
                "The natural width change does not assume two-state closure; failure of closure rules out that finite closure extension, not the width-four result.",
                "A full positive-separation Hankel rank is evidence in the actual two-point sequence, not just an ambient state count.",
                "Square-freeness, when present, excludes a Jordan collision only in this width-five charge-one block at h=v=1/2.",
                "No inference to all widths, other charges, the Matching One norm-4 residual, or an energy-operator projection follows from this one calculation.",
            ],
        },
        "execution": {
            "python": sys.version, "executable": sys.executable, "machine": platform.machine(),
            "elapsed_seconds": elapsed, "checkout_head": head,
            "script_sha256": script_sha,
            "codec_sha256": hashlib.sha256(codec_path.read_bytes()).hexdigest(),
            "command": " ".join([sys.executable, "scripts/p398_physical_two_point.py", *sys.argv[1:]]),
            "new_monte_carlo_samples": 0, "parameter_points": 1, "widths_executed": [5],
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(serialize(result), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    json_sha = hashlib.sha256(args.output_json.read_bytes()).hexdigest()
    closed_text = "仍然闭合" if ranks[1] == ranks[0] else "不再闭合"
    simple_text = "是（与导数 gcd 为常数）" if len(gcd) == 1 else f"否（与导数 gcd 次数 {len(gcd) - 1}）"
    mode_text = (f"正间距 Hankel 满秩且特征多项式 square-free，因此两个已指定读出的整个矩阵序列\n"
                 f"共同探测到 {hankel_rank} 个不同的有限宽度传播本征成分；这不等于 {hankel_rank} 个连续极限场。"
                 if hankel_rank == len(backward) and len(gcd) == 1 else
                 "实际可观测序列的传播维数以下方 Hankel 秩为准；不能仅凭环境空间维数认定成分数。")
    report = f"""# P398：width-5 真实正权两点传播

## 科学结果

同一局部 A/L 读出在 width-5 **{closed_text}为二维传播**。
精确 Krylov 秩（最大传播次数 0…7）为 `{ranks}`；只用正间距
`C(1)…C(7)` 组成的 4×4 个 2×2 block Hankel 矩阵秩为 **{hankel_rank}**。
完整 charge-one 块的维数为 {len(backward)}，特征多项式是否 square-free：{simple_text}。

这是使用已有物理接口得到的有限宽度科学结果，不是再造 rooted closure 工具。
{mode_text}

## 固定的测度与自然延拓

- Q=1 独立 square-bond cylinder，width=5，h=v=1/2。
- 正概率的 42 个 circular-noncrossing 过去连通状态，T=H V，先 vertical 再 horizontal，
  在 H 后读取；每一步是 10 个独立 Bernoulli bond bits 的精确求和，没有 Monte Carlo。
- `A=Σ_j ζ_5^j I(j~j+1)`；
  `L=Σ_j ζ_5^j I(j 为 singleton)`，j=0…4。
  两个指示函数原样保留，仅以当前圆周群的第一个非平凡 Fourier character 加权。
  width-4 的 ζ₄=i 正好还原已有 A/L；C4 和 C5 的 charge-one 并非同一群表示。
- 8 个长度 5 的轨道贡献 charge-one；两个旋转不动状态不贡献该 charge。
  平稳分布满支撑、两个读出的均值精确为零。
- 全过程在 Q(ζ₅) 中精确计算。JSON 的每个场元素依次保存
  `(1,ζ₅,ζ₅²,ζ₅³)` 的四个有理系数；小数仅供阅读。
- 这是普通的正权 42-state 模型，不把旧 signed 23-state retained-mark module
  当作随机 transfer，也没有完成旧模块的 full-Q 物理提升。

`C_ab(d)=E[O_a(X_0) conjugate(O_b(X_d))]` 是 connected 两点量。
正间距 Hankel 秩直接来自两点序列，因而不会把 8 维环境空间自动当成可观测秩。
自然延拓本来不保证二维闭合；本计算检验的是该特例能否随宽度保持。

## 实际两点读数

`C(0)`（行列次序 A,L）：

```text
{matrix_display(correlations[0])}
```

`C(1)`：

```text
{matrix_display(correlations[1])}
```

`C(8)`：

```text
{matrix_display(correlations[8])}
```

`det C(0) ≈ {determinant_c0.approximate().real:.12g}`。
完整 d=0…8 矩阵、正间距 Hankel、精确 transfer、平稳分布和特征多项式都在
[latest.json](latest.json)，无需凭幅度拟合来认定传播秩。

## 消费的已有结果与解释边界

1. `e38fe7634354b0cb2201fa55fd9b4d37ccedeef2`（`branch_only`），
   `notes/p398-positive-cylinder-propagation.md` 已完成 width-4 的正权两点 A/L 矩阵，
   两个传播本征值为 `(3±√5)/64`。本次不重算它。
2. `b35e100a3903c706dceba57c4667386eb4510ac3`（`branch_only`），
   `notes/p398-anisotropic-cylinder-spectrum.md` 已处理 width-4 完整 h/v 正权族、
   边界和 signed Jordan 特例。本次不再做 anisotropy 扫描。
3. 本结果只回答 width-5、h=v=1/2 的上述 first-character 两点传播。
   不把非二维结果写成 width-4 失效，不把 simple spectrum 泛化为所有宽度无 Jordan，
   不据此识别 Matching One 的 norm-4 场或 E_top 的微观能量投影。
4. 真正新增的信息是：固定同样的微观读出后，宽度变化是否迫使有限传播模型
   增加成分；无需增加任意新 mark 来制造该结果。

## 复现

```bash
{sys.executable} scripts/p398_physical_two_point.py
```

本次实际运行：Python {platform.python_version()} / {platform.machine()}，
{elapsed:.3f} 秒，单一 width、单一参数点、零 MC 样本；没有全仓库测试。
checkout HEAD：`{head}`（输出为本 Draft 工作树新增结果）。
脚本 SHA256：`{script_sha}`。
JSON SHA256：`{json_sha}`。
"""
    args.output_md.write_text(report, encoding="utf-8")
    print(json.dumps({"krylov_ranks": ranks, "positive_separation_hankel_rank": hankel_rank,
                      "square_free": len(gcd) == 1, "elapsed_seconds": elapsed,
                      "json_sha256": json_sha, "output_json": str(args.output_json),
                      "output_md": str(args.output_md)}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
