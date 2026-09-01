# Canonical pair 的 thermal-pivotal 闸门仍然开放

有限 midpoint 恒等式、两个独立通道反例、已关闭 Issue #536 保留的 N25
分解、开放 Issue #537 的当前 P0 与 L32/L64 有界 counter 复放合同，见
[双通道科学卡](p337-thermal-pivotal-two-channel.md)。

**判决。** 临界点上 raw canonical kernel 的绝对可和性，不能推出它的
热导数绝对可和，也不能据此解析排除 original-U 的尺度增强。障碍不是一个
缺失的局部 descriptor，而是一个明确的三点 pivotal susceptibility：翻转一个
任意远的站点可以合并两个同时接触两标记的 occupied components，从而改变
`g_xy`。现有严格 `alpha_4>1` 只控制两个端点的 raw support；Russo 求和
再积分 pivotal 位置后损失两个空间幂，即使利用 pivotal 处额外四臂，现有
指数界仍不足以证明可和。

因此本闸门的结果是 **OPEN_WITH_EXPLICIT_REMOTE_PIVOTAL**，不是解析
排除。当前 canonical P0 是 #537；#536 只是 issue-only 来源，#539 是 P2
exact-N25 支撑。本轮没有改变这些 Issue 的生命周期或优先级标签。下一动作先用
root-conditioned Hessian 的 projected landing-matrix minors 判伪 pure-thermal
rank-one 机制，再只为存活结构证明带 `q/E` 权的三点 influence tail。未得到它
之前，不把 raw 距离曲线、另一个 completion、三插入或未冻结的宏观比值当作
机制判决。

本说明没有新增样本、枚举、测试、服务器任务或数值拟合。

## 1. 有限环面上的精确 Russo--Margulis 接口

在有限顶点集 `V_L` 上令

\[
 \nabla_z f(\omega)=f(\omega^{z=1})-f(\omega^{z=0}).
\]

对任意有界的 signed observable，Bernoulli 乘积测度给出精确恒等式

\[
 \partial_p\mathbb E_p f=\sum_{z\in V_L}\mathbb E_p\nabla_z f.       \tag{1}
\]

这里取固定 canonical `Kreg=K2+K0` 的两点核 `g_xy`，并在任一标记
occupied 时令它为零。已完成的有限核满足

\[
 -\frac94\le g_{xy}\le\frac{43}{16},\qquad
 |\nabla_z g_{xy}|\le\frac{79}{16}.                         \tag{2}
\]

于是

\[
 |\partial_p\mathbb E_p g_{xy}|
 \le I_L(x,y;p):=\sum_z\mathbb E_p|\nabla_z g_{xy}|.          \tag{3}
\]

式(3)是所需 thermal gate 的无抵消版本。证明 signed 导数较小并不足以
控制 original-U；后者还需要带 `q/E` 的版本，见第5节。

### 1.1 pivotal support 的确定性内容

设 x、y 非相邻。记 `S_xy(omega)` 为“两标记均 vacant，且至少两个不同
occupied NN components 同时接触两处端口”的事件。已有 exact support
theorem 给

\[
 g_{xy}\ne0\Longrightarrow S_{xy}.                          \tag{4}
\]

若 `z` 不是 x 或 y 且 `nabla_z g_xy !=0`，则至少成立：

1. `S_xy(omega^{z=0})` 或 `S_xy(omega^{z=1})`；否则两值都为零；
2. 两个八端口 partition 不同；否则 exact Bell8 lookup 给相同值；
3. 在 `z=0` 的状态中，z 的 occupied 邻点至少属于两个不同、且均连接到
   某个 marked port 的黑 NN components；占据 z 会合并它们；
4. 这两个黑 components 在到达 marked-port 区域前必须保持分离。平面
   matching 对偶因而在它们之间给出两条白臂。

所以远离两个标记时，z 处至少有黑、白、黑、白四条交替臂，直到首次到达
标记邻域。z=x 或y的两个特殊项只由原双端点四臂事件控制，不产生位置和。

这比“某一状态有 `g_xy !=0`”更强，但仍不把 influence 限制在 x--y 的
直径内：连接端口的 components 可以先绕到任意远处再被 z 合并。

## 2. raw 双四臂证明为何不能直接微分

令 `r=floor(d_L(x,y)/4)`。临界 raw 证明给

\[
 \mathbb E_{p_c}|g_{xy}|
 \le C\,\pi_4(2,r)^2.                                      \tag{5}
\]

[van den Berg--Nolin, arXiv:2008.01606](https://arxiv.org/abs/2008.01606)
Theorem 1.1 与第5节给出的有限尺度形式是

\[
 \pi_4(3r)\le C r^{-1}\sqrt{\pi_2(r)}.                     \tag{6}
\]

RSW 在几何级环带上只需给某个 `rho>0` 使
`pi_2(r)<=C r^(-rho)`；于是

\[
 \pi_4(r)\le C r^{-1-\rho/2},\qquad
 \mathbb E|g_{xy}|\le C r^{-2-\rho}.                       \tag{7}
\]

这足够对 y 求和。它没有控制式(3)中的 z 求和。

### 2.1 只使用原 support 会精确损失两个空间幂

对 `d(z,{x,y})=O(r)` 的 `O(r^2)` 个站点，有限能量把固定 z 状态的
`S_xy` 概率与未条件概率相差至多常数。因此仅由(4)--(7)得到

\[
 \sum_{z:\ d(z,\{x,y\})=O(r)}
       \mathbb E|\nabla_zg_{xy}|
 \le C r^2\pi_4(r)^2
 \le C r^{-\rho}.                                          \tag{8}
\]

相较 raw `r^(-2-rho)`，Russo 位置积分恰损失 `r^2`。对距离 r 的 y
还有 `O(r)` 个方向位置；(8)远不足以给 `sum_y |partial_p E g_xy|<infinity`。
这不是导数发散的下界，只是证明：原双四臂 upper bound 本身不能闭合热闸门。

### 2.2 利用 pivotal 自己的四臂，仍缺一个定量指数

考虑 x、y、z 两两相距同一量级 R 的 bulk block。在三个互不相交的固定
比例圆盘内，`nabla_zg_xy !=0` 强制三个四臂事件。乘积测度给出严格的
局部上界

\[
 \mathbb P(\nabla_zg_{xy}\ne0)
 \le C\pi_4(cR)^3.                                         \tag{9}
\]

一个 dyadic R-block 有 `O(R^2)` 个 y 和 `O(R^2)` 个 z，所以其绝对
influence 总量至多

\[
              C R^4\pi_4(cR)^3.                            \tag{10}
\]

若写成假设性的幂 `pi_4(R) approximately R^(-alpha_4)`，要靠(10)证明
bulk blocks 可和，至少要求

\[
                         \alpha_4>\frac43.                  \tag{11}
\]

不需要假定指数存在也能看见现有缺口：把(6)代入(10)只给

\[
 R^4\pi_4(R)^3\le C R^{1-3\rho/2}.                         \tag{12}
\]

要让每个 dyadic block 衰减，需要 `rho>2/3`。文献输入只给某个
`rho>0`，没有给这个强度。因此 `alpha_4>1` 不是 thermal summability
定理。这里也没有把三角格的精确指数移植到方格。

### 2.3 相撞中心还需要六臂融合界

式(9)只覆盖三点良好分离。若 `s=d(x,z)<<R=d(x,y)`，两个近邻四臂中心
在尺度 s 以上融合；标准 arm-separation 路径需要控制从 s 到 R 的至少
六臂事件。相应的充分型 dyadic 和为

\[
 \sum_R R^2\pi_4(R)
   \sum_{s\le R}s^2\pi_4(s)^2\pi_6(s,R).                  \tag{13}
\]

远处 z、而 x/y 先融合的情形是同一多尺度问题。当前输入既没有(13)所需的
六臂 landing estimate，也没有把所有 collision trees 统一成一个三点
influence bound。因而最小缺失几何信息是：

- bulk：第三个、位于 pivotal z 的四臂事件及足够强的三点积分界；
- fused blocks：六臂事件的定量界；
- far blocks：这些臂跨越 `max(d(x,z),d(y,z))` 的一致 landing/连接界。

仅增加一个 raw 距离点不会补上任何一项。

## 3. 一个任意远的 merger pivotal：局域化的明确反例机制

存在一族完全有限、可画在任意大方形或环面上的配置，说明不能把 z 限制在
`O(d(x,y))` 内：

1. 取两个非相邻 vacant marks x、y；
2. 放置两个互不相交的 occupied corridors A、B，使 A 和 B 都分别接触
   x 与 y；在两端使用已实现的 two-component `0011` landing；
3. 把两条 corridor 继续延伸到一个任意远的 vacant site z，并令 z 的两个
   邻点分别属于 A 与 B；以 vacant/matching corridor 保持 A、B 分离。

在 `z=0` 时两个 shared components 的 canonical coefficient 是已知的
`g_xy=1/16`。在 `z=1` 时 A 与 B 被合并，只剩一个 shared occupied
component，exact support theorem 给 `g_xy=0`。因此

\[
                       \nabla_zg_{xy}=-\frac1{16}.           \tag{14}
\]

z 与 x、y 的距离可以任意大。这个构造不是概率下界，也不证明导数发散；
它严格否定“g 只依赖标记间局部区域”以及任何用有限 dependency radius
删除远 pivotal 的证明。它把可能的逃逸机制唯一地命名为：**两个长程共享
components 的远端合并 pivotal**。

#537 进一步给出一个更强的 issue-only all-`L` rerouting 族：对每个 `L>=4`，
可令 `x=(0,0)`、`y=(1,1)`、`z=(floor(L/2),0)`，使翻转前后 rank 都为 `0`、
shared-component count 都为 `2`，而 `g_xy:1/4 -> 1/2`。因此 carrier 不能只
按 `s` 是否变化分类；same-count signed rerouting 同样可在远处发生，证明必须
同时保留 mark separation 与 pivot distance。

PR #509 当前 head `2785e3bb` 还对 Bell8 join 做了完整有限审计：4,140 个
partitions、64,954 个允许的 2/3/4-block joins 中有 29,970 个非零差，其中
9,952 个为正、20,018 个为负，并给出 `|Delta g|<=17/4` 的抽象有界包络。
这些计数固定了有限 carrier 的符号复杂性与 envelope；它们不提供临界出现概率
或渐近主导性。

## 4. 什么条件足以关闭 scalar thermal gate

定义尾部三点 influence

\[
 \mathcal I(R)=
 \sup_{L,x}\sum_{y:\ d_L(x,y)>R}\sum_{z\in V_L}
       \mathbb E_{p_c}|\nabla_zg_{xy}^{(L)}|.                \tag{15}
\]

若能证明存在 `delta>0`

\[
                  \mathcal I(R)\le C R^{-\delta},          \tag{16}
\]

则(1)、Tonelli 与 dominated convergence 给

\[
 \sup_{L,x}\sum_y|\partial_p\mathbb E g_{xy}^{(L)}|<\infty
\]

及一致尾界。式(10)--(13)给出证明(16)必须支付的 arm 账，而非已经完成的
证明。反过来，若某个 dyadic block 的 signed influence 有稳定非零下界，
则 scalar thermal derivative 不可和；这仍需符号信息，不能从绝对上界反推。

还需要一个 near-critical 版本。original pooled roots `p_L` 和旧生产的
`p_ref` 未被现有定理认证为字面 `p_c`；要应用到它们，(16)必须在声明的
near-critical window 内一致成立。只在 `p_c` 的点态证明不能静默替代。

## 5. original-U 四项：两项受 raw bound，两项逃逸

对任一 pair source S，现有完整线性泛函为

\[
 \frac{\mathcal L[S]}{A_N}
 =\underbrace{\frac{jY_p}{D}}_{\text{direct}}
 -\underbrace{\frac{Y_{pp}jM}{D^2}}_{\text{root motion}}
 -\underbrace{\frac{Y_pjM_p}{D^2}}_{\text{source slope}}
 +\underbrace{\frac{Y_pM_{pp}jM}{D^3}}_{\text{root slope}}, \tag{17}
\]

其中 `jM`、`jY` 是 q/E 与 S 的中心化响应，`D=M_p`。对一个固定宏观
window `d(x,y)>=cL`，raw theorem 给

\[
 \mathbb E|S_{2,\mathcal W}|, |jM|, |jY|
       \le C N^{-1}L^{-\eta}.                              \tag{18}
\]

四项的状态如下：

| original-U 项 | pair-source 因子 | raw bound 是否控制 | 尚缺什么 |
|---|---|---|---|
| direct `jY_p/D` | `p` 导数后的 E-weighted pair source | 否 | `nabla_z g` 与 `g*nabla_z E` 的带权三点 influence |
| root motion `-Y_pp jM/D²` | 未微分的 `jM` | 是，按(18) | 另需 baseline `Y_pp/D²` 的尺度界 |
| source slope `-Y_p jM_p/D²` | `p` 导数后的 q-weighted pair source | 否 | `nabla_z g` 与 `g*nabla_z q` 的带权三点 influence |
| root slope `+Y_p M_pp jM/D³` | 未微分的 `jM` | 是，按(18) | 另需 baseline `Y_pM_pp/D³` 的尺度界 |

混合根移动 `p_star^(logQ,epsilon²)=-jM/D` 同样由 raw source bound 控制
其分子。这里“控制”只指 pair-source 空间尾；baseline 分母和热 jets 自身
仍需已有有限尺寸尺度输入，不能由(5)推出。

尤其，证明 scalar `sum_y |partial_p E g_xy|` 可和仍不是充分条件。
对 `O=q,E`，

\[
 \nabla_z(Og)=O(\omega^{z=1})\nabla_zg
              +g(\omega^{z=0})\nabla_zO,                   \tag{19}
\]

而中心化还含 `partial_p E[O]` 和 `partial_p E[g]`。第二项允许 z 改变
torus rank/q/E 而不改变局部 kernel；它有自己的 topological pivotal
support。中心化项还要求 raw pair tail 乘上 O 的总 influence。关闭
original-U gate 所需的一个充分对象因此是

\[
 \mathcal I_O(R)=\sup_{L,x}\sum_{d(x,y)>R}\sum_z
 \mathbb E\left[|\nabla_zg_{xy}|+
       \max(|g_{xy}^{0}|,|g_{xy}^{1}|)|\nabla_zO|\right]
 +\sup_{L,x}\sum_{d(x,y)>R}\mathbb E|g_{xy}|
       \sum_z\mathbb E|\nabla_zO|,
 \quad O=q,E.                                               \tag{20}
\]

要求(20)对 q、E 都有一致衰减，并结合 baseline jets，才足以把 canonical
interaction 降为有限局部修正。

## 6. 改变下一动作的最终判决

PR #509 `2785e3bb` 给出与式(17)等价但更适合判伪的压缩。令
`M(p(u,epsilon),epsilon)=u`、`Yhat(u,epsilon)=Y(p(u,epsilon),epsilon)`，则

\[
 \left.\partial_\epsilon\widehat Y\right|_u=jY-RjM,
 \qquad
 \partial_u\partial_\epsilon\widehat Y=\frac{T_N}{D}=\frac{J_N}{A_N}. \tag{21}
\]

而 `T_t=<H,(a-Ea)S-beta B>_pool` 把 source、readout、root 与 slope 放入同一
Schur-projected signed covariance。故第一个 prospective discriminator 不是直接
追求最强 absolute tail，而是在同一 `C4`/Schur 投影后构造 ordinary four-arm
landing source/thermal transfer matrix并检查全部 `2x2` minors。

1. **不能解析排除。** `alpha_4>1` 完成 raw summability，但 Russo 位置
   积分与 remote merger pivotal 留下真实的 thermal escape。
2. **不启动描述性 raw 距离网格。** 它已经被 raw theorem 回答，也不能
   测量(20)。
3. **下一项先做一个有限判伪器：** 任一 projected `2x2` minor 非零即排除
   leading three-packet 是 pure thermal coordinate，并把后续定理聚焦到该
   non-rank-one signed landing functional；若全部 minors 为零，才推进
   four-packet remainder 的 `R^4*pi4(R)^4` absolute 包络。两条分支都还必须
   结算 bulk common landing、fused collision、far rerouting、`g` 不变的
   topological pivotal 与 near-critical root window，且不能让两列互相抵消后
   才报“可和”。若 tail 与 baseline jets 受控，canonical pair 的渐近放大机制
   应正式降级；若失败，必须给一个具名 carrier 的有符号 lower mechanism 与
   original-U 尺度预测。
4. **有界语义preflight已经完成并停止。** [64-counter结果](../results/p337-thermal-pivotal-preflight/REPORT.md)
   在既有L32/L64流上验证了5,242,880个pair/site callback、midpoint恒等式与
   carrier实现；有限回放只出现18个kernel-changing事件，呈稀疏、外部、
   shell-localized形态，没有观察到kernel-preserving topological事件。它没有
   finite-L rejection threshold、总体centering或full J2，因而不增加counter、
   seed或距离window。下一动作回到上述有限 minor 判伪器，不把preflight升级为
   物理显著性。

这保留了 N25 非零 J2 的有限事实，也解释了它为什么不能由 raw C 的正号
或绝对可和性决定。它没有把 remote pivotal 宣称为已证实的渐近主导项，
也没有识别 H4、Jordan 或任何连续场。

## 来源边界

- raw summability：execution-independent theory
  [`eed2190c`](https://github.com/LightChainr/Matching-One/blob/eed2190c04b67084ab5aef5827e00377853a0bca/notes/p337-critical-spatial-summability.md)；
- support、pointwise bound 与 two-component witness：PR #509
  [`baa5d33b`](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md)；
- 完整 original-U source functional：
  [`regular-pair-joint-u-functional.md`](regular-pair-joint-u-functional.md)；
- 已完成有限 J2：execution
  [`410015f5`](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md)。
- 已关闭 Issue #536 的有限双通道来源、开放 #537 的当前 gate 与已完成语义preflight：
  [`p337-thermal-pivotal-two-channel.md`](p337-thermal-pivotal-two-channel.md)；
- 四臂、六臂、pivotal measure、near-critical与cross-Fourier适用边界：
  [`p337-thermal-pivotal-literature-bridge.md`](p337-thermal-pivotal-literature-bridge.md)。

式(1)--(4)、(14)、(17)--(19)是有限配置恒等式或确定性 support 结论；
式(9)--(12)是三点良好分离块的严格上界账，式(13)明确记录尚未证明的
fused-block 充分输入；(16)与(20)是尚待证明的充分条件，未在本文中冒充
定理。文献桥进一步说明朴素的三中心四臂乘积即使代入三角格精确
`alpha_4=5/4`仍给`R^(1/4+o(1))`，所以Delivery A必须提供共同landing增益或
逐通道signed cancellation；near-critical使用还须认证root不超过相关长度。
