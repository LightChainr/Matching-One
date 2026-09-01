# 固定 prefix 联合出生形状：预定读数（运行前）

2026-08-31。先完整阅读 dc4bb041 的协方差层级及 03603388 的无偏恒等式。它们已经完成原8 quartet 的全总体 prefix/label/suffix 分解，本轮不重复该分解，也不把旧结果当新发现。新增信息来自已经生成、尚未用于条件形状的 cell00 新64 quartet，与原8合并为72；不再生成任何 prefix、label 或 suffix。

## 问题与四个读数

局部 A 和 integral A 的两源均值 Jacobian 已有二维证据。检验它是否可与“每个固定完整 prefix 内，两出生时间仅平移而联合形状保持”相容。

在每个物理几何中令 X=K1/(N+1)、Y=K2/(N+1)、C=(X+Y)/2、W=Y-X。数值实现把 C 减去已固定的 k0/(N+1)，方差和协方差不变，不按结果选择中心。

只预定以下四种物理读数：

1. H Var(C | Z)：出生中心条件展宽。
2. H Cov(C,W | Z)：中心与间隔条件耦合。
3. H Var(W | Z)：出生间隔条件展宽；任意保留每条路径 W 的共同平移都要求其为零。
4. E_Z[1_cell00 sum_(physical orientation o, source s=plus/minus) ||H_s Cov((C,W) | Z,o)||_F²]：条件形状响应的非负总体能量，包含 1,2,1 Frobenius 权重，防止正负 prefix 响应相消。估计无偏并允许负值，绝不截零。

前三项先在物理几何算完整条件形状导数，再形成 S=(first+second)/2、D=(first-second)/delta_cos4；保留两条固定 plus/minus exact-census score，不选最佳方向。第四项用物理几何及固定 plus/minus 源的平方和，不使用带几何差分缩放的 D 平方。

若每个 Z 内有 X_t=X+t a(Z)、Y_t=Y+t b(Z)，四个目标都为零；非零的均值或能量可排除这一条件纯平移解释。均值为零不足以支持该解释；能量未分辨也不证明局部不变。这里不拟合分布、不扫 p、不搜峰、不新增 source、width 或多个形状变体。

## 无偏估计与复用

直接复用 375cd3a12b2b7a87d79148a59f62b95898f9e471 的 exact-score quartet archive，其含原 birth_k、score、完整 census 类计数，不重读/重建原8的层级分解。新64从已冻结 extension/*.csv.gz 与同一 exact census 构造相同矩。

令每 quartet 的 b_f=(f_U+f_V)/2，h_f=(s_U-s_V)(f_U-f_V)/2，其中 f_U 是两 suffix 均值。条件形状估计为 mean(h_fg)-P_Q(b_f,h_g)-P_Q(b_g,h_f)，P_Q 使用 q!=r 的有序跨 quartet 乘积。平方能量使用展开后的2/3/4阶不同 quartet U乘积；禁止直接平方有噪声的条件导数。

计算 old8、new64、combined72 三者及配对差；它们共享 prefix，combined72 还包含前两者，不作独立复现。三个版本的 cell00 贡献均保留每批1000原prefix、20原批次分母，不能叫全总体响应或cell00条件均值。结果保留每prefix值、20批向量、delete-one值和共享协方差factor。

把既有 results-extension/score.json 的 A(p_ref)/integral A 局部行列式均值逐批向量直接拼入共同factor，展示均值二维与形状响应是否并存；不重算或扩充其端点，不据显著性选读数。单N分别报告估计/SE与批次删除稳定性，不假定N325/N425共同常数，也不将源/几何/端点计作独立证据。

## 实现核验

在科学读数前，以小型固定数组的全排列 U统计核验 O(Q) 集合分拆公式，并以有限离散分布穷举证明条件形状能量估计的均值等于真值。对旧/新数据逐prefix核对 ID、quartet域0..7/8..71、原rank00、精确score类内中心、已保存来源哈希。主分析只运行固定读数，云机等待root明确分配后使用。
