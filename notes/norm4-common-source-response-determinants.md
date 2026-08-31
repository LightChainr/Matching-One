# 共同源的有限响应行列式：三态组成与方向响应

本笔记固定本轮主分析的具名假说，不增加模拟、源拟合或自由指数。全部响应均在同一 N、同一个 pooled matching 根及同一估计器合同下计算。

## 固定读出、源和线性响应

两方向记为 g，`q=rank_black−1∈{−1,0,1}`、`E=q²`；令 `m_g=<q>`、`e_g=<E>`。上横线表示两方向平均，`P4=(first−second)/Δcos(4θ)` 使用原精确方向归一化。所有下式在 `bar m(p0)=0` 的 pooled 根求值，撇号表示 p 导数。

固定行顺序为

\[
 \left(\overline{P_1},\ W,\ U\right)
 =\left(1-\bar e,\ \mathcal P_4[e],\ A_N B/D\right),
 \quad A_N=N^{13/8}/2,
\]

其中 `D=bar m′`、`B=P4[e′]`、`H=P4[e″]`、`T=bar m″`。固定源列顺序为 **bulk `s=CB+CW`、无量纲 q、无量纲 E**，分别对应 `exp(t s)`、`exp(λq)`、`exp(ηE)`；读入的 density 源 `S=s/N` 必须先乘 N 转成 s 响应，q/E 列不乘 N。

对固定源 X，记 `Jq[X]=Cov(q,X)`、`JE[X]=Cov(E,X)`，则根移动为 `p0_dot[X]=−bar Jq[X]/D`，三个读出的 moving-root 响应是

\[
 r_X=-\overline{J_E[X]}+\frac{\bar e'}D\overline{J_q[X]},
 \qquad
 w_X=\mathcal P_4\!\left[J_E[X]-e'\frac{\overline{J_q[X]}}D\right],
\]

\[
 u_X=A_N\left\{
 \frac{\mathcal P_4[J_E'[X]]}D
 -\frac{H\overline{J_q[X]}}{D^2}
 -\frac{B\overline{J_q'[X]}}{D^2}
 +\frac{BT\overline{J_q[X]}}{D^3}\right\}.
\]

## 无需新混合矩的 q/E 对照列

三态代数 `q²=E、qE=q、E²=E` 给出每方向的精确表达式：

| 源 X | Jq[X] | JE[X] | Jq′[X] | JE′[X] |
|---|---|---|---|---|
| q | e−m² | m(1−e) | e′−2mm′ | m′(1−e)−me′ |
| E | m(1−e) | e(1−e) | m′(1−e)−me′ | e′(1−2e) |

因此 q/E 两列只需现有无源 q/E thermal jets；bulk s 列只需现有逐 K 的 `s,qs,Es` 及热导数。原 [two-phase core](../scripts/norm4_source_two_phase_core.py) 已含这些三态公式和 U 的线性源泛函。无需 `s²`、新回放或将条件源均值指数化为有限 λ 权重。

## 两个单源假说与一个共同三态假说

方向专属常数 `a_g` 在归一化测度下响应为零。两方向共同的 `b_N K` 只平移 Bernoulli log-odds；沿 pooled 根移动时，这三个读出的响应均严格为零。b、c、d 在每个 N 内必须方向共同，并在 p 导数中保持固定；这里不额外假设它们跨 N 相同。

**Hq：共同 q＋热时钟。** 若 s 在上述读出的一阶映射可由 `a_g+b_N K+c_N q` 表示，则

\[
 \Delta_q=r_su_q-r_qu_s=0.
\]

**HE：共同 E＋热时钟。** 若相应映射可由 `a_g+b_N K+d_N E` 表示，则

\[
 \Delta_E=r_su_E-r_Eu_s=0.
\]

这两个 2×2 行列式使用固定行 `(rank1root,U)`；未知 c 或 d 已被消去，不需要用响应结果拟合幅度。只测 rank1root 一个读出通常允许调整幅度，无法完成这种排除。

**H3：共同三态＋热时钟。** 任意 `c_N q+d_N E` 一般可以覆盖两个读出的响应平面；因此排除 Hq 或 HE 不等于排除全部共同三态投影。固定第三行 `Wroot=P4[E]` 后，H3 的必要条件是

\[
 \Delta_3=\det
 \begin{pmatrix}
 r_s&r_q&r_E\\
 w_s&w_q&w_E\\
 u_s&u_q&u_E
 \end{pmatrix}=0.
\]

显著非零可排除 `a_g+b_N K+c_N q+d_N E` 对这三个读出的共同一阶映射。零结果仅表示该必要条件尚未被排除；q/E 对照列若退化，行列式本身也会缺乏辨别力。

## 解释与协方差合同

保留同方向对、同 counter 组、同根合同的共同协方差；源与基准贡献均传播，多个行列式和读出不叠加为独立证据。方向专属 `a_g` 常数允许存在；自由的 `b_g`、`c_g`、`d_g` 会扩大假说空间，上述共同系数行列式不能排除该空间。

这项分析区分“有限三态组成改变”与“超出共同三态＋热时钟投影的方向响应”。W/U 仍是全局方向观察者，非零行列式不自动识别能量场、局部空间算符或连续极限生成元。它为既定六 N 主分析提供一个明确、可证伪的一阶比较，而非新增计算清单。
