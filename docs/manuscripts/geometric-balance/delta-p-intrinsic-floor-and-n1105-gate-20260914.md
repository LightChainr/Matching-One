# dpfloor：δp 的内禀一致性底噪 `F` 与 N1105 的 go/no-go

日期：2026-09-14 ｜ 机器 `DevEnvC_NePnUn`（账号2，16 vCPU / 32 GiB，Python 3.9.9）
全部算术在云机执行，本机只做读写/上传下载。**仓库写操作：无**（GitHub 只读）。
**没有跑 N1105，没有开任何大算。** 最贵的单个运算是 `axis_n9` 的 w=9 稠密特征分解（212 s）与
`slope52_n2` 的 ARPACK 复算（508 s）；全流程合计约 1 小时墙钟，远小于一次 N1105。

交付：`floor.json`（每条路径的原始数字 + `F` + `S/F`）、`scripts/`（11 个自包含脚本）、
`raw/`（全部机器可读输出与日志）。

---

## 0. 必须回答的那一句

> ### 「N1105 应该现在开跑吗？」
>
> **按现有配置（shipped settings）开跑：不应该 —— 那是 NO-GO。**
> **把两个配置改掉之后开跑：可以 —— GO，安全系数 `S/F ≈ 5.8×10³`。**
>
> 具体地：
>
> | 问题 | 数 | 判定 |
> |---|---|---|
> | `n1105mix` 报的「系统误差 2.10e−12」 | 是 **100 % 的 `p_c` 记账口径差**（可消除），不是任何实现底噪 | 归零 |
> | 紧配置下的实现底噪 `F` | **4.68e−16**（p 单位） | — |
> | 信号 `S`（`C_mix=0.0073`） | 2.71e−12 | **`S/F = 5793` ⇒ GO** |
> | 信号 `S`（`C_mix=0.0898`） | 3.33e−11 | **`S/F = 71185` ⇒ GO** |
> | shipped 配置下的复现离散度 `F_shipped` | **6.31e−12**（几何相关，最大 6.31e−12；shipped 文件自身对紧值的偏差最大 **9.29e−12**） | **`S/F_shipped = 0.43` ⇒ NO-GO** |
>
> **要改的两行（以及怎么知道改好了）见 §4。** 一句话：**判据本身可达，但必须先统一 `p_c` 并把
> 特征求解器/求根器的公差收紧；用 shipped 的那套设置去跑 N1105，先天测不出 `P0=P8`。**

---

## 1. `Ω(axis, ℓ=8)` 两份文件差 `8.6e−9` 的成因判定 —— **(c) 定义/记账口径，且具体就是 `p_c`**

### 1.1 数字（可复现：`scripts/fl_convention.py` → `raw/step0_convention.json`）

两份文件指的是：

| 记号 | 文件 | 生产者 | 用的 `p_c` |
|---|---|---|---|
| 文件 A | `oblique-spin4-controls.json` | rev769 的 `scripts/oblique_charge_transfer.py` | **0.59274605079** |
| 文件 B | `closure-amplitude-raw.json` | sector802b 的 `s802b_c1_amplitude.py`（`PC=0.59274605079210`） | **0.59274605079210** |

（生产者出处见 §2.0：rev769 的脚本在本机 `rev769-repo/scripts/oblique_charge_transfer.py`，
默认 `--reference-pc 0.59274605079`。这是**已被核实的代码级出处**，不是猜测。）

同一几何 `(1,0), n=8`（= 轴向 w=8，ℓ=8，ℓ⁴=4096）：

```
Ω := -(p_root - p_c) * ell^4
Ω_文件A = 0.3001967221684936
Ω_文件B = 0.30019673077049447
差      = 8.602000889368355e-09          （相对 2.865e-08）
```

### 1.2 严格分解（残差恰好为 0）

```
Ω_B - Ω_A = -[(p_root,B - p_root,A) - (p_c,B - p_c,A)] * ell^4
```

| 项 | 值 | 占比 |
|---|---|---|
| `p_c,B - p_c,A` = 0.59274605079210 − 0.59274605079 | **+2.0999868510784836e-12** | — |
| 来自 `p_c` 不同的贡献 `= Δp_c·ℓ⁴` | **+8.601546142017469e-09** | **99.9947 %** |
| `p_root,B - p_root,A` | **−1.1102230246251565e-16**（1 ulp） | — |
| 来自 `p_root` 不同的贡献 `= −Δp_root·ℓ⁴` | +4.547473508864641e-13 | 0.0053 % |
| **残差** | **0.0（精确）** | — |

⇒ **`8.6e−9` 全部由「两份文件用了不同的 `p_c`」解释，剩余 4.5e−13 只是根值的 1 个 ulp。**
`n1105mix` 的 note §1.2/§3.3 说「该差**不能**由两处 `p_c` 约定解释」是**错的**：
两处 `p_c` 相差 2.1e−12（15955 个 ulp），乘 ℓ⁴=4096 恰好是 8.6e−9。它把
「δp 差 = 2.1e−12」与「`p_c` 差 = 2.1e−12」当成了巧合，其实那是同一个恒等式。

### 1.3 排除另外三个候选

* **(a) 闭合幅度 `A_r` —— 排除。** `Ω`（以及 `Δ_w`、`p_root`）只由 **Perron 根**决定；
  sector802b 自己在 `s802b_c1_amplitude.py` 的 C3(a) 就写明「`Delta_w` is built from Perron
  ROOTS only, so an amplitude cannot enter it」。`A_open` 等是诊断量。本报告用
  `p_ch` 与 `p_c` 重算 `Ω_B`，与文件里印的 `w4_times_offset` **精确相符（差 0.0）**，
  完全不需要任何幅度量。
* **(b) 求解器容差/停止条件 —— 不是主因，但确有 0.005 %。** 两份文件的**特征分解代码不同**
  （文件 A 用 ARPACK `tol=1e-10`，文件 B 用稠密 `scipy.linalg.eig`），
  但它们在同一个几何上给出的 `p_root` 只差 1 ulp。真正的容差代价不在这一对文件上，
  而在 shipped 配置的**复现性**上（§3.3）。
* **(d) 数值噪声 —— 排除。** 残差为 0，即这个 8.6e−9 里没有第三样东西。

### 1.4 定义核对（别把定义差当数值差）

`n1105mix` §1.1 发现的定义错位被**独立复核确认**：

| 字段 | 实际定义 | 与重算的最大差 |
|---|---|---|
| 文件 A `A_estimate` | `−(p_root − p_c)·ℓ⁴ / cos4θ` | **0.0** |
| 文件 A `shift_times_ell4` | `(p_root − p_c)·ℓ⁴ = −Ω` | （与 `A_estimate` 差恰好一个 `cos4θ`） |
| 文件 B `w4_times_offset` | `−(p_root − p_c)·ℓ⁴ = Ω` | **0.0** |

⇒ 两份文件印的是**不同字段**，但只要都换算成 `Ω`，剩下的就是 `p_c`。

---

## 2. 独立路径清单（每条：算法是什么、在哪一步不同、是不是真独立）

### 2.0 先核代码（任务 §2.3 第 1 条）

* **`sector802-out/scripts/sector802_lib.py` 与 `sector802b-out/scripts/sector802_lib.py`
  逐字节相同**（md5 均为 `4646de953442c8cc71ae255f6cc3eec2`），而 `s802b_c1_amplitude.py`
  直接 `from sector802_lib import enumerate_states, aggregate, build_matrices, perron_fh`。
  ⇒ **产生文件 A（`root-response-raw.json` / `root-response.json`）与文件 B
  （`closure-amplitude-raw.json`）的是同一份 rank 生产**；两者的差**只反映调用参数与后处理**
  （文件 B 多算了 `solve_root` 的精确根与幅度解剖）。**它们不是两条独立路径，本报告明说。**
* 引擎：两者都用 `tagged_winding_span.py`（pinned #739 副本），`sector802` 的 note §0 与
  sector802_lib 的 docstring 都写明了。
* 本报告 reran 的 `PATH A` 就是这份库（md5 校对过，见 `raw/log_pathA.txt` 前的 md5 输出）。
* rev769 的 `oblique_charge_transfer.py` 是**另一套**自动机：SL(2,Z) Bezout 斜基、
  `(ds,dt)` 边集、**frontier 记忆 = max dt**、**没有 winding 元组**（用一个"矛盾合并"
  `dsu.bad` 当场丢弃整行）。⇒ 它与 `sector802_lib` 在**状态表示与 winding 判定**上不同。

### 2.1 路径表

| 路径 | 算法 | 与别的路径**在哪一步**不同 | 是真独立路径？ | 覆盖几何 |
|---|---|---|---|---|
| **A_lib** | `#739` 引擎 BFS 自动机（winding 元组）→ 稠密 `R` → `scipy.linalg.eig` → 二分/brentq | 基线；**就是文件 A/B 的那份代码** | **与文件 A/B 不独立**（同库） | 轴向 w=4,6,8 |
| **B_oblique** | Bezout/SL(2,Z) 斜基自动机（`dsu.bad` 当场拒环）→ 稠密或 ARPACK → brentq | **自动机不同**（状态表示 + winding 判定 + 前沿记忆），求解器/求根器也不同 | **是**（与 A 无共享代码） | 全部 13 个几何（含 4 条斜向） |
| **M_mine** | **双覆盖**前沿划分：状态 = 提升到 2w 个位置的连通划分；**winding ⟺ 提升把 `X` 与 `X+w` 认成同一块** | **第三种自动机**，winding 判定方式与前两者都不同 | **是** | 轴向 w=3,4,5,6 |
| **H_ld** | 由整数重数表在 **float128** 里装配 `R` → **幂迭代**求 Perron 根 → float128 secant 求根，**完全不用 LAPACK** | **算术级独立**（其余全用 float64/LAPACK） | **是（算术）** | 轴向 w=4,6,8 |
| **S_solvers** | 同一张 `R`（来自 A）配 5 种特征求解器（`scipy.linalg.eig` / `numpy.eigvals` / `scipy.eigvals` / ARPACK 1e−10 / ARPACK 1e−14） | **只换线性代数与求根器** | **不是独立路径**（同自动机、同矩阵），只测求解器噪声 | 轴向 w=4,6,8 |
| （旁证）rev769 `diagonal_charge_transfer.py` | 同源但"两行记忆"的 (1,1) 变体 | 与 B 同源 | 不作独立路径 | (1,1) n=2..5 |

### 2.2「同一个对象」的证书（比态数更强）

`scripts/fl_matrix_cert.py` → `raw/step_matrix_certificate.json`：在 w=3,4,5,6、两个扇区上比较
A/B/M 三条自动机的**安全块全谱**：

| w | 扇区 | 安全态数 (A/B/M) | `max\|Δλ\|` A−B | A−M | B−M |
|---|---|---|---|---|---|
| 3 | G4/G8 | 7/7/7 | **0.0** | **0.0** | **0.0** |
| 4 | G4/G8 | 19/19/19 | **0.0** | **0.0** | **0.0** |
| 5 | G4/G8 | 51/51/51 | **0.0** | **0.0** | **0.0** |
| 6 | G4/G8 | 141/141/141 | **0.0** | **0.0** | **0.0** |

⇒ 三个独立实现在 w≤6 上生成**同一个矩阵**（态数相同、谱逐位相同）。
这把"路径差"干净地还原成"求解器/求根器/算术差"。

---

## 3. 底噪 `F` 的数字

### 3.1 定义

> **`F` := 在同一几何、同一 `p_c`、同一定义下，各独立路径给出的 `p_root` 的最大两两离差（以 p 为单位）。**

被采纳的独立路径集合：`{A_lib, B_oblique_tight, M_mine}`（三条不同实现），
以及再加 `H_ld`（float128 算术）的版本。全部用**同一套紧设置**（稠密特征分解；ARPACK 仅在
矩阵过大时使用且 `tol≤1e-14`；求根 `xtol=1e-16`）。

### 3.2 结果（`floor.json`）

| 量 | 值 |
|---|---|
| `F_p_implementation_float64`（三条 float64 实现） | **2.22e−16**（= 1 ulp） |
| `F_p_including_arithmetic_path`（+ float128） | **4.68e−16** |
| **采纳的 `F`** | **4.68e−16** |
| `F_Ω` 换算到 ℓ=8（`F·8⁴`） | 1.92e−12 |
| 求解器变体离散度（同矩阵，w=4/6/8） | 3.33e−16 / 8.88e−16 / **5.55e−16** |
| ARPACK 自身 λ 误差（探针，w=8） | ≤3.50e−15 相对 ⇒ **≤1.73e−15 在 p 上** |
| 逐几何（有 ≥2 条独立路径的） | w=4 **0.0**；w=5 **0.0**；w=6 **2.22e−16**；w=8 **0.0**（A 与 B 逐位相同） |

逐几何明细（p_root，取自 `floor.json → geometry_table`）：

| 几何 | ℓ | A_lib | B_oblique_tight | M_mine | H_ld | `F`(采纳) |
|---|---|---|---|---|---|---|
| axis_n4 | 4.000 | 0.5914171708531382 | 0.5914171708531382 | 0.5914171708531382 | 0.5914171708531385 | 3.19e−16 |
| axis_n5 | 5.000 | — | 0.5922358232050272 | 0.5922358232050272 | — | 0.0 |
| axis_n6 | 6.000 | 0.5925073562056412 | 0.5925073562056414 | 0.5925073562056414 | 0.5925073562056417 | 4.68e−16 |
| axis_n8 | 8.000 | 0.5926727605746269 | 0.5926727605746269 | — | 0.5926727605746273 | 3.94e−16 |

**`axis_n8`（就是 2.1e−12 那个几何）：A 与 B 两条独立实现的 `p_root` 逐位相同（离差 0.0）；
唯一非零离差来自 float128 路径，3.94e−16（≈3.5 ulp）。**

### 3.3 它不是全部 —— **shipped 配置的代价比 `F` 大 4 个数量级**

把 B 路径用**出厂设置**（ARPACK `tol=1e-10` + `scipy.brentq` `xtol=3e-11`，即文件 A 的实际设置）
重跑一遍（`raw/step_shipped_vs_repro.json`）：

| 几何 | ℓ | shipped − tight (p) | shipped 设置内部离散 (p) |
|---|---|---|---|
| axis_n8 | 8.000 | −1.00e−16 | −5.69e−13 |
| axis_n9 | 9.000 | +2.10e−15 | −1.30e−12 |
| diag_n4 | 5.657 | +4.31e−13 | −3.76e−14 |
| diag_n5 | 7.071 | +2.80e−13 | −6.18e−13 |
| slope21_n3 | 6.708 | −1.84e−14 | −1.30e−13 |
| slope21_n4 | 8.944 | **+1.71e−12** | −4.82e−12 |
| slope31_n3 | 9.487 | −5.80e−14 | −1.74e−12 |
| slope32_n2 | 7.211 | −1.36e−13 | +1.50e−15 |
| slope52_n2 | 10.770 | **+9.29e−12** | −6.31e−12 |

* **`max |shipped − tight| = 9.29e−12`** —— shipped 文件自己的根值，最坏处偏离紧值 9.3e−12，
  是信号 `S = 2.71e−12` 的 **3.4 倍**。
* **`max |loose − tight| = 6.31e−12`** ⇒ 出厂设置的**复现性只有 ~6e−12**。
* 附带：`scipy.brentq` 的**默认** `xtol=2e-12` 就足以把 `p_root` 推偏 4.33e−13（w=8），
  7.04e−14（w=6）—— 见 `raw/pathS_solvers.json` 的 `p_root_brentq_default`。

### 3.4 它是否覆盖 `S`？

```
S = 2.71e-12（小幅度）       ⇒  S/F = 5793      ⇒ GO
S = 3.33e-11（大幅度）       ⇒  S/F = 71185     ⇒ GO
S = 2.71e-12 / F_shipped     ⇒  0.43            ⇒ NO-GO（若用出厂设置）
```

---

## 4. 裁定：GO / NO-GO / 临界

**裁定：GO —— 但有条件；退回条件只有一个（配置），不是判据本身。**

1. **「2.1e−12 的不可消除底噪」不存在。** 它是 `p_c` 记账差（§1）。判据并没有在几何之前
   撞上精度墙；它撞上的是**两份文件用了两个不同的 `p_c`**。
2. **真底噪 `F = 4.68e−16`，`S/F ≈ 5.8×10³`。** 所以 `P0=P8` 的 `~1e−12` 判据在
   算术上是**可达的**（余量近 4 个数量级）。
3. **但若沿用 shipped 的那套设置，就是 NO-GO**（`S/F_shipped = 0.43`）。差别全在：

### 修什么（两行级改动）

* **(i) 四个取向统一 `p_c`，并且只用一个值。** 这是**正确性问题**，不是精度问题：
  `Ω = −(p_root − p_c)ℓ⁴`，`p_c` 差 ε 会给所有取向加同一个纯 H0 位移 `ε·ℓ⁴`，
  在 ℓ=33.24 时 ε=1e−12 ⇒ Ω 位移 1.2e−6，直接把 `P0=P8` 打穿（而 H8 不动）。
* **(ii) 收紧求解器与求根器：** 状态数 n≤~2500 用稠密 `scipy.linalg.eig`；
  更大时 ARPACK `tol ≤ 1e-13`；求根用 `xtol ≤ 1e-15`（**不要** scipy `brentq` 默认 2e−12，
  **更不要** shipped 的 3e−11）。四个取向必须用**同一套**设置。

### 怎么知道修好了

* 对**同一几何**用 §2 里 ≥2 条独立路径复算 `p_root`，离差必须 **≤1e−15**（本报告的 F 检验）；
* 把 `(1,0)` 轴向当控制几何：三条独立自动机的安全块**全谱逐位相同**（§2.2 已证 w≤6）；
* 开机前先跑一次 N=325（ℓ=18.03，3 取向）的 H0/H4/H8 版本 + 同样的 F 检验；若它也过，
  再上 N=1105。N=325 比 N=1105 便宜得多且条件数更好（L2=0.598 vs 0.967）。

---

## 5. 不被 `F` 覆盖的系统偏差（**换路径测不出来**）

`F` 只测「实现噪声」。下面这些是**所有路径共有**的，必须单独记账：

1. **`p_c` 的取值/统一性。** `F` 的定义里 `p_c` 是**共同常数**；一旦四个取向用了不同 `p_c`
   （shipped 文件正是如此），产生的 H0 污染 `ε·ℓ⁴` 完全在 `F` 之外，且量级足以击穿判据。
2. **有限 ℓ 的模型误差（ℓ^{−8} 污染）。** `n1105mix` 在轴向上测到 `C = 3.54±0.45`（8σ）；
   这是**模型**误差，换代码路径一点都不会变。跨 ℓ 拟合 `B(θ)` 的偏差不进入 `F`。
3. **态空间完备性 / 任何截断。** 本报告的所有自动机都是**完备枚举**（`state_cap` 从未触发），
   所以这里这一项为零；但 N1105 若用 `state_cap`、记忆截断或任何近似，那是新的系统偏差，
   `F` **测不出来**。
4. **斜向 twist 的实现路径。** 斜向几何（`(a,b)` 的 Bezout 斜基 + 前沿记忆）在本报告里
   只有 `B_oblique` **一条**路径（`diagonal_charge_transfer.py` 与它同源）。
   `F` 对**轴向**是三条独立路径互证，对**斜向**实际上只有一条 ⇒ 斜向 twist 的正确性
   没有被 `F` 覆盖。
5. **谐波分解/投影后处理。** 从 `Ω(θ,ℓ)` 抽 `B0,B4,B8` 的拟合、投影权重、跨 ℓ 外推
   都是后处理模型误差，不在 `F` 内。
6. **求根器停止条件**（这一项**在** `F` 内，但极易被配置放大：默认 brentq 就够造成 4e−13）。

**举证上的红线**：本报告**不用** `F` 去担保任何会被上面 1–6 击穿的结论。§4 的 GO 只针对
「**在统一 `p_c` + 紧公差 + 完备自动机**这一前提下，四个取向的 `p_root` 能一致到 5e−16」。
第 1 项已由本报告修掉；第 2 项（ℓ^{−8}）仍然是判定 `P0=P8` 的**主要剩余风险**。

---

## 6. 「若要判 `P0=P8`，还缺什么」

1. **ℓ^{−8} 污染的定量控制。** 现在只有 4–5 个宽度、4 个局部斜率、2 种外推
   （`n1105mix` §2 已复现：换横坐标得 3.995，加 1/w² 得 4.336；包含 17/4=4.25）。
   N1105 的四个取向在**同一个 ℓ** 上，这能去掉"跨 ℓ 外推"的一部分，但
   `B(θ)` 里仍混着 `C(θ)/ℓ²`（ℓ^{−8}）——**需要 ≥3 个 ℓ 或独立的 `C(θ)` 估计**。
2. **N1105 四个取向的可行性尚未实测。** `(24,23)` 的短边 23（四方向短边 4/9/12/23），
   状态数未知；`n1105mix` §3.2 已指出"去掉 `(24,23)`"会让条件数烂 8.5×。
   这需要在**统一 `p_c` + 紧公差**下先做一次**小规模 scout**（例如先只做
   `(33,4)` 与 `(24,23)` 两个方向，看状态数与墙钟）。
3. **斜向 twist 自动机的独立复核。** §5 第 4 条：斜向只有一条路径。建议把
   `M_mine` 的双覆盖构造推广到斜基（或反之，把 Bezout 自动机写成第二份独立实现）。
4. **N=325 先导**（ℓ=18.03，3 取向，H0/H4/H8）—— 便宜、条件数更好，可作为
   「四取向 `p_root` 一致到 ≤1e−15」这一 F 检验的第一次真实彩排。

---

## 7. 复现方式

```bash
H=~/.workbuddy/skills/connect-huawei-codebuddy/scripts/huawei
export PYTHONPATH=/workspace/mo/compat:/workspace/mass761/pylibs/pylibs_local:/workspace/dpfloor/scripts
cd /workspace/dpfloor
python3 scripts/fl_convention.py                      # §1 约定分解
python3 scripts/fl_path_lib.py /workspace/sectorA/in/engine 0.59274605079210 4,6,8 out/pathA_lib.json
python3 scripts/fl_path_oblique.py 0.59274605079210 tight out/pathB_oblique_tight.json
python3 scripts/fl_path_oblique.py 0.59274605079210 tight out/pathB_oblique_tight2.json slope21_n4,slope31_n3,slope32_n2,slope52_n2
python3 scripts/fl_path_oblique.py 0.59274605079      loose out/pathB_oblique_loose.json
python3 scripts/fl_path_mine.py 0.59274605079210 3,4,5,6 out/pathM_mine.json
python3 scripts/fl_hp.py 0.59274605079210 4,6,8 out/pathH_longdouble.json
python3 scripts/fl_solvers.py 0.59274605079210 4,6,8 out/pathS_solvers.json \
        scipy_eig,numpy_eigvals,scipy_eigvals,arpack_1e10,arpack_1e14
python3 scripts/fl_matrix_cert.py                     # §2.2 同对象证书
python3 scripts/fl_arpack_probe.py                    # §3.3 ARPACK 公差代价
python3 scripts/fl_probe_hp.py                        # H_ld 收敛证书
python3 scripts/fl_shipped_vs_repro.py                # §3.3 shipped vs 紧
python3 scripts/fl_assemble.py                        # F 与 S/F
```

环境：Python 3.9.9 + numpy 1.26.4 / scipy 1.13.1 / mpmath 1.4.1
（`PYTHONPATH` 指向 `/workspace/mo/compat` 的 `int.bit_count` 垫片与
`/workspace/mass761/pylibs/pylibs_local` 的现成科学栈；**未安装、未修改任何系统包**）。
本报告**未写入** `/workspace/sectorA/`、`/workspace/v802src/`、`/workspace/n1105mix/`、
`/workspace/verify768/`；只写 `/workspace/dpfloor/`。

## 8. 交付物

| 文件 | 内容 |
|---|---|
| `dpfloor-out/note.md` | 本文件 |
| `dpfloor-out/floor.json` | 逐路径原始 `p_root`/`Ω` + 逐几何 `F` + `S/F` + 判定 |
| `dpfloor-out/scripts/*.py` | 11 个自包含脚本（云机上实际跑过的版本） |
| `dpfloor-out/raw/*.json` | `step0_convention` / `pathA_lib` / `pathB_oblique_{tight,tight2,loose}` / `pathM_mine` / `pathH_longdouble` / `pathS_solvers` / `step_arpack_probe` / `step_matrix_certificate` / `step_hp_certificate` / `step_shipped_vs_repro` |
| `dpfloor-out/raw/log_*.txt` | 全部运行日志 |
