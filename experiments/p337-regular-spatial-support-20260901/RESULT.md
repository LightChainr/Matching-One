# 正则 pair 的空间传递至少需要两个共享占据组件

**单一共享组件不能传递这套相互作用的首阶 Q 激活。** 对两个非相邻、edge-port 集合不重叠的 vacant 标记，固定外部占据后，若至多一个外部 equality component 同时接触两处，则双插入系数精确分解为两个单插入系数的乘积，其首阶 Q 导数为零。这个选择规则适用于保持逐配置 Q1 零值的正则、颜色置换不变补全；它没有依赖新选的 counterterm。

本次进一步把规则接到了**占据总和后的空间响应**，并给出固定 canonical 补全的有限上界：

```text
C_xy(p) = ∂logQ ∂lambda_x ∂lambda_y log Z |_(Q=1,lambda=0)
        = E_p[a_xy(A)],
|C_xy(p)| ≤ (43/16) Pr_p{x,y vacant，至少两个不同占据组件接触两处}.
```

右侧是含 vacancy 条件的无条件概率。它不是原 global U 的数值界，也没有指定距离衰减指数。

## 固定模型与本次新增内容

固定输入为执行分支 [2ba8863f 的 canonical Kreg=K2+K0](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/local-pair-two-insertion-algebra.md)，singlet 系数恒为1，端口顺序 N/E/S/W，q/E/rank 始终按原占据定义。

该分支已经完成所有有限网络 Q1 零定理、单点原 U 的非零混合响应 `W=-0.04503611397592696`，以及固定四路径占据的 `13/8`。**这些不计作本次新增结果。** 新增的是一般两点的共享组件选择规则、完整固定核、占据平均恒等式和上界，以及原 U 的正确双插入接口。

两个局部参数 lambda_x、lambda_y 独立；插入只在对应点 vacant 时出现。整个计算没有选择 N、距离、Q 或补全参数来提高得分。

## 有限核与独立核查

外部八端口连通有 Bell(8)=4140 种抽象分区。distinct exterior components 可以同色；计算不能把它们强制染成不同颜色。

| 两处共享外部分量数 | 分区数 | 激活恰为零 |
|---:|---:|---:|
| 0 | 225 | 225 |
| 1 | 1369 | 1369 |
| 2 | 1922 | 672 |
| 3 | 600 | 0 |
| 4 | 24 | 0 |

全表正1394项、零2266项、负480项，范围 `[-9/4,43/16]`。这是所有抽象连接上的包络，**没有声称每项都能由指定平面/环面实现，或43/16就是某个物理几何的最大值**。它仍给每个实际占据的合法上界。两组件条件仅必要，不保证非零，亦不保证占据总和后没有符号抵消。

[主计算](compute.py)对每个 exterior partition 求颜色 coarsenings，使用 falling-factorial 的精确 Q1 导数；合计167894个分区/coarsening对。[独立计算](review/verify_diagrams.py)直接展开56个 delta 乘积为15个 equality diagrams，sew 两核后逐项求 Q 系数、自由颜色数及归一化的导数。两种算法的4140项[逐项完全一致](review/COMPARISON.json)。单核15项、两端独立C4转动、交换标记和既有四路径13/8也全部通过。全部运算为整数/有理数。

## “两个”确实能达到：一个实际占据见证

[预先写明的16×16占据](results/witness.json)含两个 vacant 标记 `(3,7),(11,7)`、两个相互分离的占据带。它们在两处分别连接 N/E 和 S/W，八端口分区为 `00110011`；全占据图 rank0。只核查这一构型的连通性，没有枚举16×16总体。

对这两组件的颜色，若两色相同，核为零；若不同，每个核为 `(Q-2)/(4Q)`。所以其归一化双系数为

```text
B_xy(Q) = (Q-1)(Q-2)^2/(16Q^3),    a_xy = 1/16 > 0.
```

因此不能把必要阈值提高到三个或四个共享组件。这是实际原图上的条件系数，不能由它推断整个 Bernoulli 平均正号。

## 与原 U 的连接及停止动作

[证明 §§3、6](PROOF.md)给出精确接口：在有限图、非零归一化和简单共同根的局部解析域内，

```text
∂lambda_x ∂lambda_y ∂logQ U |_(Q=1,lambda=0) = W[a_xy].
```

W 是已存在的完整原 U 线性响应泛函；每几何分别中心化，保留共同根移动和分母变化四项。没有额外 `a_x*a_y` 项，`Cov(a_x,a_y)` 也不能替代 a_xy。若另行使用所有站点共同 `epsilon/N`，每个无序对对 `∂logQ∂epsilon²U` 贡献 `2W[a_xy]/N²`（裸的 Q1 epsilon 导数仍为零），不得混淆约定。

**停止把“只有单一共享占据组件”的传播解释用于这套首阶 Q 激活。** 该排除不是对全部 Matching-One 机制，也不是 rank1 sector 整体的排除。剩余问题已经具体缩成：由固定 sewing 核加权的多组件空间概率，经过同一个 W 后能否产生所需响应。下一步应约束这项概率/符号抵消或给出可失败的空间预测，不再补单点 N25 score、另造 source covariance，或扫补全系数。

本次没有求出 C_xy 的距离曲线或双点 U 数值，没有做连续场/Jordan 识别，没有新增随机样本、云作业或独立统计票。P154/P334/F4 停线不变；现有合同不自动延伸。

从本目录复算：

```bash
python3 compute.py results
python3 review/verify_diagrams.py
python3 review/compare_csv.py
```

只需 Python 标准库。输入、代码、证明和结果哈希见 [MANIFEST-SHA256.json](MANIFEST-SHA256.json)。
