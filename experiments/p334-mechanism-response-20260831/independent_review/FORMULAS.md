# P334 determinant-square U-statistic：独立数学复核

设固定prefix Z下的quartet矩阵为 `X_i=(a_i,b_i;c_i,d_i)`，i=1,…,n。假定它们条件iid；同一矩阵的四个分量可以任意相关。目标为

\[
\delta(Z)=\det\mathbb E[X_i\mid Z]=\mu_a\mu_d-\mu_b\mu_c,
\quad
\delta(Z)^2=\mu_a^2\mu_d^2+\mu_b^2\mu_c^2-2\mu_a\mu_b\mu_c\mu_d.
\]

这不是 `(E[det X])²`，也不是 `E[(det X)^2]` 或样本均值行列式的平方。

## 一般Möbius公式，固定阶数时O(n)

对带标签的位置1,…,k，定义

\[
p_B=\sum_{i=1}^n\prod_{r\in B}x_{i,r},\qquad
T(x_1,\ldots,x_k)=\sum_{i_1,\ldots,i_k\ {\rm all\ distinct}}\prod_{r=1}^k x_{i_r,r}.
\]

即使两个位置是同一分量（如a,a,d,d），位置仍保留不同标签。集合分拆格的反演给出

\[
T=\sum_{\pi\in\Pi_k}\left[\prod_{B\in\pi}(-1)^{|B|-1}(|B|-1)!\right]
\prod_{B\in\pi}p_B.
\]

理由：`product_B p_B`恰将索引在每个block内相等的所有赋值求和；按索引相等模式分组后，在集合分拆格上反演即可只保留全异索引。k=4只有15个非空子集矩与15个分拆，故一趟累积矩即可，时间O(n)、额外空间O(1)。

阶4的显式式子为

\[
T=p_1p_2p_3p_4
-\sum_{\{i,j\}}p_{ij}p_kp_l
+\sum_{\text{3 pairings}}p_{ij}p_{kl}
+2\sum_{\text{4 triples}}p_{ijk}p_l
-6p_{1234}.
\]

第一负和共有6项，第二正和3项，三元组和4项。

## 可直接实现的det与det²

记`S_x=sum_i x_i`、`S_xy=sum_i x_i y_i`，并令 `(n)_k=n(n−1)…(n−k+1)`。

\[
\widehat\delta=\frac{S_aS_d-S_{ad}-S_bS_c+S_{bc}}{(n)_2}.
\]

对任意两分量x,y，记`S_rs=sum_i x_i^r y_i^s`、`A=S_10`、`D=S_01`，则

\[
T(x,x,y,y)=A^2D^2-S_{20}D^2-S_{02}A^2-4S_{11}AD
+S_{20}S_{02}+2S_{11}^2+4S_{21}D+4S_{12}A-6S_{22}.
\]

对四个不同位置，

\[
\begin{aligned}
T(a,b,c,d)={}&S_aS_bS_cS_d\\
&-(S_{ab}S_cS_d+S_{ac}S_bS_d+S_{ad}S_bS_c
+S_{bc}S_aS_d+S_{bd}S_aS_c+S_{cd}S_aS_b)\\
&+(S_{ab}S_{cd}+S_{ac}S_{bd}+S_{ad}S_{bc})\\
&+2(S_{abc}S_d+S_{abd}S_c+S_{acd}S_b+S_{bcd}S_a)
-6S_{abcd}.
\end{aligned}
\]

最终

\[
\boxed{\widehat{\delta^2}=\frac{T(a,a,d,d)+T(b,b,c,c)-2T(a,b,c,d)}{(n)_4}}.
\]

因为每个全异索引乘积的条件期望是对应均值之积，分母恰为全异有序元组数。n=72时，分母为 `(72)_2=5112`、`(72)_4=24,690,960`。不得将上述有序分子除以`choose(n,4)`或n⁴。

独立等价写法也可用于小n审核：设

\[
h_{ij}=\frac{a_i d_j+a_jd_i-b_i c_j-b_jc_i}{2}.
\]

则det估计是所有`i<j`的h均值；平方估计是

\[
\frac1{\binom n4}\sum_{i<j<k<l}
\frac{h_{ij}h_{kl}+h_{ik}h_{jl}+h_{il}h_{jk}}3.
\]

这明确要求乘积中的两对没有重叠quartet。直接平方det估计会包含共享quartet、自身平方项，通常有偏。

## 负值必须保留

有限样本的无偏平方估计可以为负。令X以各1/2概率取`(1,1;1,0)`和`(1,1;1,2)`，则总体均值是全1矩阵，δ²=0。n=4恰好各出现两次时，无偏平方估计为−1/3。全分布精确计算：

- `E[U4]=0`；
- `E[max(0,U4)]=1/8`，截零造成正偏；
- `E[U2²]=1/4`，直接平方也造成正偏。

n=72各出现36次时，U4=−1/71。不得截零、取绝对值或将负值解释为负的总体平方。若下游需要总体δ²非负的约束，应在单独声明的推断层处理，不能篡改无偏观测量。

## 验证与应用条件

`verify_ustat.py`以Fraction精确算术完成：25组随机有理数据、n=4…8的全异有序元组枚举；Möbius公式与上述显式公式；独立的无序四元组/三配对式；n=4、5、8、72常矩阵；含非零目标的三点有限分布与零目标的二点有限分布期望枚举。结果见`verification.json`。

旧8+新增64只能在给定prefix后属于同一分布、相互不重叠且条件独立时按n=72合并；旧8若因quartet观测值筛选而被条件化，不能仅靠增加64个样本恢复这个前提。若3053个prefix由其预先状态double-R0确定，这不同于按已观测quartet值筛选。跨prefix、跨方向、同一外层batch的依赖仍须在最终scorer保留，72个quartet不能代替外层prefix/批的推断单位。

本验证仅证明代数与无偏性；不验证实际raw标记的独立性、旧新批对齐或实际outer权重。生产浮点实现宜使用补偿求和/高精度累积，并检查不同尺度下的误差；公式存在正常的抵消，不能用截零掩盖数值误差。

运行命令：`python3 verify_ustat.py`（仅Python标准库）。本次实际用 `/Users/lc/python-envs/research-py311/bin/python verify_ustat.py` 执行，125次有序枚举比较、25次独立配对比较及两种有限分布的精确期望全部通过。非零目标例中 δ=5/6、δ²=25/36，精确枚举的E[U2]、E[U4]分别与之完全相等。
