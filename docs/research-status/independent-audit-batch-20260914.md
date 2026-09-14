# 独立核验批次：2026-09-14 外包核验产出（第三方审计，非自有推导）

**提交者身份**：这是**操作者侧的独立核验批次**，不是项目所有者的推导。
本批全部内容来自**独立子代理**在**孤立云机**上跑的验证/复算/审计，目标是
「把合作者交付的包真的跑起来 + 对承重恒等式做对抗性审计」。

> ⚠️ **口径声明（重要）**
> 1. 本批是**第三方复现与审计**，**不声称**任何新的物理结论。
> 2. 所有数字都由子代理在云机脚本中产出；本机 Mac 未做任何计算。
> 3. 严格区分 **精确恒等式 / 有限枚举事实 / 有限宽度拟合 / 条件定理 / 假设 / 猜想**。
> 4. **不与 STATUS 科学判决冲突、不改 `#275` 原始 U 合同、不改任何原件。**
> 5. 本分支**未被合并**；本批**未新增 issue、未删除数据、未操作服务器**。
> 6. 部分交付的**笔记由子代理在被中断前写出**；被中断者的**原始输出**由操作者原样归档，
>    并在文件名中标注 `raw` / `telemetry`。

---

## 一、本批包含什么（按与当前 P0 的关联度排序）

### A. `#802`（P0）—— 匹配根的切向—法向联合响应

| 文件 | 内容 | 强度 |
|---|---|---|
| `docs/manuscripts/geometric-balance/root-tangential-normal-response-20260914.md` | 控制表逐位复现 + `Δ_w`/`Δ'_w` + `(T,N)` 联合响应 + 四项裁定 V1–V4 | **控制表复现：有限枚举事实**；指数：**有限宽度拟合** |
| `results/research-dispatch/sector-root-response-20260914.json` | 同上，机器可读 | —— |
| `docs/manuscripts/geometric-balance/closure-amplitude-and-delta-exponent-20260914.md` | 闭合幅度 `A_r` 与 `Δ_w` 的**独立性**（幅度进不去 `Δ_w`）+ 外推指数 | **精确恒等式 + 有限宽度外推** |

**要点（可复核）**
- 控制表 `w·I4`/`w·I8`（w=4..8）逐位复现，最大偏差 **4.87e-11**
  （控制表只印 10 位小数 ⇒ 实为复现到其全部打印精度）。
- `Δ_w(p_c)` = `+3.275880e-03 … +1.483071e-04`（w=4..8，全非零）；
  `Δ'_w(p_c)` = `+2.4659 … +2.0236`（全正）。导数用 **Feynman–Hellmann 解析式**。
- **结构性结论**：一阶下 `N = 0 ⟺ g` 与 `∂_p` 平行 `⟺ T = −1`。
  即「法向响应为零」与「纯热重新调定」在一阶是同一件事；
  盲点的定量内容是倍数比 `φ ≈ 1.9e-4`（内禀读出比切向分量大约 5×10³ 倍）。
- **`Δ_w` 只由 Perron 根构成 ⇒ 闭合幅度无法进入它**（`identical_Delta_after_amplitude_removal = true`；
  有限 `m` 迹估计 `D(m) − Δ_w → 0` 如 `e^{−mΔ}/m`）。
- **独立验证**：Alexander duality **222720 次检查 0 违例**；反射律 `P8_p(k)=P4_p(2−k)`
  **12 例最大差 0.000e+00**；FH 实现 vs 精确代数恒等式 **4.13e-13**。
- **不能宣称**：`Δ_w`/`Δ'_w`/根偏移**同源，是一项数据生产不是三份独立证据**；
  `w≤8` 的拟合**不是** exponent measurement；未识别 8-arm / spin±4。

### B. `#768`（P0）—— 首个非公共扇区修正的通道与预测

| 文件 | 内容 |
|---|---|
| `docs/manuscripts/geometric-balance/first-noncommon-correction-channels-20260914.md` | G1 四件事分开 / G2 等级计数 / G3 预测像 |
| `results/research-dispatch/theory-768-channels-20260914.json` | 通道清单与 `(T,N)` 预测表 |
| `results/research-dispatch/theory-768-level-dimensions-20260914.json` | 模去 `L_{−1}` 的等级维数 |
| `results/research-dispatch/theory-768-virasoro-20260914.json` | 等级计数的精确有理数输出 |

**要点（可复核）**
- **`0,0,1,1` 成立**，但**只对「把 level-2 奇异向量商掉的不可约商」这一个模**：
  `p(n)=1,1,2,3,5,7,11`；`rank G_n = p(n)−p(n−2)`；奇异向量只在 level 2（维数 1）；
  `χ = −3L_{−2}+2L_{−1}²`，**范数恰为 0**。
- **null 的去向**：(i) 与 (ii) **代数上都能站住，但关联函数完全相同**（零范且被正模湮灭）；
  **只有 (iii) 改变物理**，而 `c=0` 渗流的物理模**就是 (iii)**（能量算符是零范态且是对数多重态 bottom field）。
- **G1：只证明消去两类**——identity 家族（dual-even，在差里恒等消去）与所有 **C4 禁止**通道
  （spin `1,2,3 mod 4`）。**`k=2..7` 的标量 arm 通道（`x=1/4,2/3,5/4,2,35/12,4`）未被证明消去**。
- **G3 的关键（本轮可操作的排除）**：在平移不变圆柱/环面上，权重 `(h,h̄)` 的一点函数 `∝ δ_{h,h̄}`。
  因此**热族 level-4 手征后代（spin ±4, `x=21/4`）在动量 0 观测量中矩阵元为 0**；
  而 scalar 8-arm（`(21/8,21/8)`，spin 0）不零 ⇒ **预测像不同**，构成可行的排除。
- **不能宣称**：不能宣称 level-4 类的**差矩阵元非零**；不能宣称 8-arm 或 spin4 谁胜出；
  不能宣称 `w^{−17/4}` 已被解释。**`∂⁴ε` 显式验证其类为 0，不得当作插入。**

### C. `#775`（P1，当前最小产物）—— 精确 rank 表 L=3..6

| 文件 | 内容 |
|---|---|
| `results/geometric-consistency/rank-sector-C-L5-20260914.json` | **L=5 全枚举**（和 = 2²⁵ 精确） |
| `results/geometric-consistency/rank-sector-C-L6-20260914.json` | **L=6 全枚举**（和 = 2³⁶ 精确） |
| `results/geometric-consistency/ml-exponent-20260914.json` | `M_L(p_c)` 与 `2−x` 拟合 |
| `results/geometric-consistency/ml-definition-check-20260914.json` | **定义核验**（L=3 与 MZ16 多项式逐项相同） |

**要点**
- `C[L,j,k]` 对 **L=3,4,5,6** 精确整数；L=5 和 = `2²⁵`、L=6 和 = `2³⁶`（无漏配置）；
  并查集 vs 通用覆盖 BFS **0/25 万**差异；L=6 交叉 **0/5 万**差异。
- **L=4 的 `rank_totals` `[36559, 19932, 9045]` 与独立仲裁的 `W4−W8` 分布逐位吻合。**
- **定义已澄清**：L=3 时 `P2−P0` 与 **Mertens–Ziff 2016** 的多项式**逐项相同**
  `[−1,0,0,6,0,0,0,−18,18,−4]`；交叉绕行配置数 91 = `C2_total`。
- `M_L(p_c)`（L=3..6）= `0.02496990, 0.01033324, 0.00445693, 0.00236193`；
  拟合 `rms_log`：`−3.25`(Jacobsen) `0.0428`、`−3.42`(MZ16) `0.0396`、
  `−4.0` `0.0226`、free `−4.205` `0.0193`(dof=1)。
  ⇒ **两个候选都不被支持，数据偏好更陡；`L≤6` 判不了 `21/8`。**
- **不能宣称**：`L≤6` 只是**有限宽度兼容性**，不是指数测量；L=7/8 未产出
  （`#650` 已取消该阶梯）。

### D. `21/8` 的文献定标（**更正一条仓库内判词**）

| 文件 | 内容 |
|---|---|
| `docs/manuscripts/geometric-balance/literature-calibration-frontier-20260914.md` | 长程检索（10 轮 / 30 次查询 / 4 次全文核对） |
| `docs/manuscripts/geometric-balance/correction-21-8-20260914.md` | 操作者侧更正说明 |

**要点**：文献里**有两份互相冲突的测量**，不是一致的否定。
- **支持整数 4**：**Jacobsen 2015**, *J. Phys. A* **48** 454003, arXiv:`1507.03027` §7.1 式 (36)——
  **方形格点 site 渗流（同一模型）**、周期 TL 转移矩阵特征值、**n ≤ 21**、**非 MC**，
  `Δ₁ = 4.000 1(2)`；原文「**It appears inevitable to admit that `Δ₁ = 4` exactly**」；式 (40) `Δ_k = 2(k+1)`。
  该文同时给出 `p_c = 0.59274605079210(2)`（即本项目惯用的 `p_ref`）。
- **冲突**：**Mertens & Ziff 2016**, *PRE* **94** 062152, arXiv:`1603.07289`，`2−x = −3.42`、斜率 `−4.07`，
  **作者自认「larger systems are needed」**。
- ⇒ 正确表述是「**生死未判定**」，**不是**「已被文献否掉」。
- 另：`√(cE)·Z ~ Laplace` 是**命名恒等式**（Gauss–Laplace transmutation，Ding–Blitzstein arXiv:`1510.08765`）。
- `#775` 那张按**同调秩**分辨的整数表**无人发表过**（已发表的 `Z_0/Z_1/Z_2` 是**绕行方向型**）。

### E. P1/P2 已完成或负面结果的归档

| 文件 | 内容 |
|---|---|
| `docs/manuscripts/geometric-balance/gamma-w-parameter-dependence-crosscheck-20260914.md` | **独立确认**「`γ_w−κ = A e^{−cw}` 里的 `c` 不普适」；并澄清「受限转移阵」在该参数族下**是恒等** |
| `docs/manuscripts/geometric-balance/cavity-field-98-negative-20260914.md` | `#779`：`z_c(p)` **强 w-依赖**；w=8 时 `z_c(1/8)=1.17878`（距 `9/8` 差 4.8%）；`Z_p(z)` **未**达 `1/(9−8z)` ⇒ **负面** |
| `docs/manuscripts/geometric-balance/query-hole-definition-difference-20260914.md` | `#777`：`S_D(0)` 与八点多项式**不一致但属定义差异**（真查询孔洞多出 4 个十点矩形） |
| `docs/manuscripts/geometric-balance/ward-kb-locking-20260914.md` | `#774`：四项 Ward verdict 全 compatible（w=2..8）；`residual_covariance_scaled` **命名不清**（不是 (6.6) 的 `Σ`） |
| `docs/manuscripts/geometric-balance/audit-780-monotone-fragmentation-20260914.md` | `#780` 包审计：9 条声称无错误结论；`216`/`49` 是**参数网格大小而非独立证据** |
| `docs/manuscripts/geometric-balance/reproduce-genealogy-noise-20260914.md` | 谱系/噪声包 layer1 逐字复现通过 |

---

## 二、本批**没有**完成的事（明确列出）

1. **`#801` 的七机作业映射**：本批**未**回填到 issue（操作者的红线原为禁止仓库写操作）。
   **操作者侧的映射表已单独落盘**，待授权后回填。
2. **`#780` 的共同标签并合缺口**：一项子代理任务**被中断，零产出**；未计入本批。
3. **时钟可实现性 / 有限误差**（`#780` 的下一份产物）：**未做**。
4. **`p_0 → p_c` 的有效范围**：**未做**（属待研究的有效范围边界）。
5. 本批**不含**任何 L5–8 尺寸阶梯（`#650` 已取消）。
