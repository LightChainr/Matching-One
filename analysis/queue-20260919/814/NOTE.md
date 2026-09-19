# #814 — Δ₁=4 机制闸门：结论 (c)

日期：2026-09-19。消费者：总控 #813/T1；实际 square-site matching 根的形状机制判定。只新增本票文件，不修改历史声称。

## 判决

**选择 (c)：现有材料不足以判定实际格点模型的 E₄ 形状律。** 不是 (a)，也没有足够证据选择 (b)。T6 不执行；E₄ 及其全部下游继续冻结。条件性 Ward 公式不是实际格点机制已识别。

## 1. 必须更正的两个前提

Jacobsen [J] §7.1 的 Δ₁=4.0001(2) 与 4,6,8,… 是数值外推及其猜想序列。§8 式(43)由两扇区共同磁指数只得到临界自由能密度差 o(n^-2)；作者紧接着明确指出，需要进一步成分才能得到式(44)的根漂移 O(n^-4)。所以“共同磁指数相消，因此已经解释出四次幂”过强。仓库 `fixed-width-charge-jacobsen-bridge` §7 已经保留了这个限制；本次不是首次在仓库识别这个缺口。

其次，仓库 `thermal-null-ward-root-mechanism` 的实际候选是**热场 ε 的四级后裔**，不是单位族的 T²+T̄²。两者自旋同为 ±4，但总标度维数不同。把扇区相消与 spin-4 写成互斥选项并不成立：相消可以消去共同低阶项，而另一个匹配奇后裔控制首个非共同项。也不能反向从“可能兼容”跳到 (a)。

## 2. 根漂移、自由能密度和行自由能必须分开

以下是明确假设下的标度记账，不是 square-site 普适性的独立证明。取热维数 x_t=5/4、y_t=2-x_t=3/4。圆柱上临界扇区自由能密度差若首个非共同贡献为 L^-x，而非退化热斜率为 L^-x_t，则根漂移为 L^-(x-x_t)。固定模参数环面的无量纲差 X=P₂-P₀ 同样给出：修正 L^(2-x) 除以热斜率 L^y_t。

| 对象 | x | RG 指数 2-x | 若其非共同且首项非零，根漂移指数 x-x_t |
|---|---:|---:|---:|
| 单位族 spin-4 准初级，简写 T²+T̄² | 4 | -2 | 11/4 |
| 热族四级 spin-4 后裔 | 21/4 | -13/4 | 4 |

因此自由能里 L^-4 与阈值根里 L^-4 不是同一句话。热候选在每行对数特征值差中对应 L^-17/4，在每行热斜率中对应 L^-1/4，比值才是 L^-4。这与仓库 bridge 的记账相容；不把该文件的小宽度诊断认作本轮独立证据。

共同 x_m 只取消首个共同磁性项。它本身没有证明所有较低维无关场的两扇区振幅相同，也没有把长圆柱的一个特征值结论提升为有限模参数的完整投影 trace 结论。

## 3. 已有条件性 E₄ 公式，及其尚未闭合的接口

采用仓库笔记的坐标与正规化：周期 u,v，τ=v/u，Im τ>0；ε 的 h=h̄=5/8，χ=(L₋₂−(2/3)L₋₁²)ε，且

    U₄ = L₋₄ ε + (80/87)L₋₂ χ + L₋₁ V₃.

对无接缝/接触异常、插入位置平移不变的周期 Ward 泛函 F，

    G₄(u,v) = (π⁴/45)u⁻⁴E₄(τ),
    F(L₋₄ ε) = 3hG₄ F(ε).

再要求目标泛函的 F(L₋₂χ)=0，便有

    F(U₄+Ū₄) = (π⁴/12)Re[u⁻⁴E₄(τ)]F(ε).

这组公式取自已读仓库条件性推导；本轮不声称重新完成 PBW/Zhu 推导。只有再建立实际格点的源映射、同阶排他性、共同微观耦合及根的非退化性，才推出

    p*Λ − pc = A Re[u⁻⁴E₄(τ)] + o(|u|⁻⁴).

A 必须在同一格点模型和坐标约定下跨形状相同，不能每个 τ 另拟合。一般 u=ℓe^(iθ) 时形状量是 ℓ^-4 Re[e^-4iθ E₄(τ)]；固定面积 N=|u|² Imτ 时应把 (Imτ)² 计入 N^-2 的系数。N 不是线性尺寸。

[S2] 的 Ising spin-4 与 [K] 的自由模型 Kronecker 展开均不能补上这个非自由、带 rank 投影的源接口。模权重、晶格四重旋转和一条圆柱极限不唯一指定完整环面响应。

## 4. 最小额外输入：不再要求一轮同轴尺寸拟合

**输入 A：实际投影下一个响应恒等式。** 对 square-site 的 X=P₂−P₀，在同一热正规化下定义临界插入泛函 F₋（包括概率正规化的导数）。计算/证明其应力张量沿环面运输的 seam/contact 项，并判定

    DΛ = F₋(U₄) − 3hG₄(u,v)F₋(ε)

及反手征差额是否消失。在上述 Ward 条件下 DΛ=(80/87)F₋(L₋₂χ)，但不能先假设实际 rank 投影满足这些条件。只需目标响应的 null 解耦，不必先解决整个 c=0 对数模块；奇异/零范不等于该投影下一点函数解耦。

**输入 B：格点首项与排他性。** 把实际匹配奇格点修正投影到该 U₄+Ū₄，确定一个跨几何共同的耦合；证明维数小于 21/4 的非共同项消失，并排除或分别计算同阶 scalar、其他模块和接缝贡献。共同磁指数不是这一投影的替代品。

**输入 C：根与几何极限。** 证明 F₋(ε) 对目标根非退化，误差对固定 τ 的适当紧集一致；若使用 Jacobsen 圆柱标定 A，还需单独控制 finite-aspect 到 cylinder 的极限交换。A、B 决定形状接口，C 使其可用于根。

具体阅读入口是仓库 Ward 笔记 §3、§5、§7 中未完成的接口；它援引的 Gaberdiel–Lang arXiv:0810.0106 式(2.10)及 Javerzat–Picco–Santachiara arXiv:1907.11041v2 §4.2 可用于核对普通 Ward/null 约定，但再读这两篇本身不会自动证明 actual-site 的 A/B/C。本轮两条外部文献仅记 [LIT]，未读取其正文。

## 5. T6 和分级

T6 继续冻结，本轮不对任何三个几何根求值。其代数消元不能把未经识别的 F_i 变成物理预言。解冻需满足 (a) 的实际源条件，不是仅把上节条件公式重写一次。

R1、R2 一级优先级和 R3–R5 的工作排序不变；R5 中受本闸门约束的 E₄ 子线仍冻结。这里不修改任何历史结果的声称级别。`two-stage-angular-radial-improvement` 的有限角向代数与历史数值保留，不能借本轮对照为其径向机制解冻。

## 6. 文献与读取边界

- [J] J. L. Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley–Lieb algebras*, arXiv:1507.03027v1。**PRIMARY_TEXT_READ**：§7.1 收敛/外推，式(28)–(37)；§8 式(40)–(44)，尤其印刷页22的限制。相关PDF页已截图核读。https://arxiv.org/abs/1507.03027
- [S1] J. Salas and A. D. Sokal, *Universal Amplitude Ratios in the Critical Two-Dimensional Ising Model on a Torus*, cond-mat/9904038v2。**PRIMARY_TEXT_READ**：§4.1 的修正来源和拟合边界。此篇不提供 actual-site 热后裔识别。https://arxiv.org/abs/cond-mat/9904038
- [S2] J. Salas（单独作者）, *Exact Finite-Size-Scaling Corrections … II. Triangular and hexagonal lattices*, cond-mat/0110287v2。**PRIMARY_TEXT_READ**：引言印刷页3–4，含方格 T²+T̄²、y=-2 的讨论。不能将此编号署名为 Salas–Sokal，也不能将 Ising 的结论直接迁移到 site percolation。https://arxiv.org/abs/cond-mat/0110287
- [K] E. V. Ivashkevich, N. Sh. Izmailian, C.-K. Hu, *Kronecker's double series and exact asymptotic expansions for free models of statistical mechanics on torus*, cond-mat/0102470v3。**PRIMARY_TEXT_READ**：引言及 §III.A，式(15)–(17)，Kronecker/θ 与自由 Ising、dimer、Gaussian 模型的适用范围。https://arxiv.org/abs/cond-mat/0102470

没有提出“未见先例”或新颖性判断。检索返回异常/无关条目没有用作不存在先例的证据。

## 本轮执行范围

实际执行：读取 #813、#814、RESEARCH-ATLAS（包括负面结果表）、AGENTS 及三份指定分支笔记；上述四篇论文指定正文段落和PDF截图；运行本目录 `python check.py` 的六项 Fraction 精确指数/系数记账，输出 result.json；生成 SHA-256 清单。未运行枚举、转移矩阵、指数拟合、T6、PBW/Zhu 检查、父PR的测试或全仓CI。Δ₁=4.0001(2) 为 [J] 引用数字；历史根与振幅未重算。没有修改原结果、冻结 scorer 或状态/路线图文件，没有合并PR。
