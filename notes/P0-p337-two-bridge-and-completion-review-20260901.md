# [P0 建议] 两桥空间因式分解与正则补全不可消去的双插入

**建议优先级：P0（最高研究建议优先级）。** 本文提交两个有限等式张量结果及可证伪的下一步，不把 P0 当作连续场身份的证据等级，不自动授权新生产。原 canonical 补全固定为 `Kreg=K2+K0`；P154/P334/F4 与 m64 普通采样的停止决定不变。

阅读基线：`baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff` 的空间支持结果，以及 `2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb` 的局部补全、四路径几何和已完成单点响应。以下**不重报**有限网络 Q1 零定理、至多一个共享组件的零规则、canonical 单点 W、或 canonical 四路径 13/8 为新成果。

## 结论与现在应改变的判断

1. **恰好两个共享外部组件时，首阶 Q 空间激活精确等于两端条件签名的乘积**：`a_xy=kappa_x*kappa_y`。这里的端点签名取决于实际共享组件怎样接入端口，不是旧的闭合单点激活 a_x。canonical 的 62 个带两桥标签的端口模式给出七种签名，完整因式分解通过全部 1,922 个相应 Bell8 分区核验。不能再把两个普通单点源的 covariance 当成这个双插入。
2. **对同一两-sector 正则补全族，四路径反射接线有不可调到零的激活**：若 `c(1)=1`、`alpha=c'(1)` 为实数，且两端使用同一补全和固定 pair 归一化，则 `a_4(alpha)=3/2+(alpha-1/2)^2/2 >=3/2`。这反对“靠有限 singlet counterterm 消去所有双插入”的解释，不建议把 canonical alpha=0 改成最小点 alpha=1/2。
3. 下一目标应是固定核加权的**带符号多组件空间概率**及原 W 的传递；不是再扫 Q、seam、c'(1)，也不是把一个条件正系数当作占据平均的正号或连续场数。

## 1. 对象、端口和已完成结果

端口顺序为 N/E/S/W。固定外部占据，所有剩余张量为普通 equality tensors 或 constant vacant tensors。其 exterior equality components 的颜色在整数 Q>=4 下独立均匀取值；不同组件允许同色。两标记 x,y 均 vacant。非相邻且 edge-port 集合不重叠时，共享 equality component 必须包含占据连接；相邻标记的共用 edge-node 也必须计数，不能遗漏。

沿用原有

```
b_x(Q)=E_col[K_x],
c_xy(Q)=E_col[K_x K_y],
a_xy = d_logQ c_xy |Q=1.
```

对固定 canonical 和本文明确说明的 regular endpoint-invisible 补全，所有非空有限网络插入在 Q1 都为零。因此 occupation-summed logarithmic coefficient 为 `C_xy=E_{mu_1}[a_xy]`，其完整 original-U 接口是 `W[a_xy]`，不是 `Cov(a_x,a_y)`。这些接口已在基线证明，本次只向前推进它们的精确核与可识别性。

原图的 q/E/rank 始终按原占据定义，不根据 virtual diagram join 重分配 rank。

## 2. 恰好两桥的精确空间因式分解

设两端恰好共享两个 exterior components，颜色为 A,B；其余颜色私有。先对每端私有颜色求和。S_Q 不变性使所得函数只有两种值：

```
f_x^=(Q)  : A=B,
f_x^ne(Q) : A!=B,
```

y 同理。这是有限颜色和的有理延拓；Q1 的 `f^ne(1)` **不是**在一个颜色的集合上实际条件于不同颜色。

条件独立立即给出

```
b_x = [f_x^=+(Q-1)f_x^ne]/Q,
b_y = [f_y^=+(Q-1)f_y^ne]/Q,
c_xy=[f_x^= f_y^=+(Q-1)f_x^ne f_y^ne]/Q.
```

相减得到全 Q 的有理恒等式

```
c_xy-b_x*b_y
 = (Q-1)/Q^2 * (f_x^=-f_x^ne)*(f_y^=-f_y^ne).       (1)
```

regular endpoint-invisible 假设给 `f_x^=(1)=f_y^=(1)=0`、`b_x,b_y=O(Q-1)`。定义 `kappa_x=f_x^ne(1)`，则

```
a_xy = kappa_x*kappa_y.                             (2)
```

这也等于 conditional connected logarithm 的首阶 Q 系数。背景 occupation measure 的 Q 导数乘以零系数，故不会在这一阶补出另一个 `a_x*a_y`。

### 2.1 canonical 签名的完整有限集合

用 0、1 表示两个有序共享组件，2、3 表示按出现次序命名的私有组件。两个共享标签都必须到达本端；这样共有 62 个四端口模式。

| kappa | 模式数 |
|---:|---:|
| -1/2 | 2 |
| -1/4 | 4 |
| 0 | 12 |
| 1/2 | 24 |
| 3/4 | 8 |
| 1 | 8 |
| 3/2 | 4 |

几个直接可核验的例子：`0011 -> -1/4`，`0101 -> -1/2`，`0123 -> 1/2`，`0213 -> 0`，`0102 -> 1`，`0112 -> 1/2`，`0212 -> 3/2`。`0011` 在两端给 `1/16`，吻合已发布的实际两组件见证。其它模式可给负乘积；非零两组件支持并不保证正号。

因此在 exactly-two-shared stratum 上，canonical 有

```
|a_xy|<=9/4,
|C_xy^(2)| <= (9/4) P{x,y vacant, exactly two shared components}. (3)
```

这是该层的上界，不把整个空间系数的已有 `43/16` 上界改成 `9/4`。完整源仍须保留 s=3、s=4 两层。

### 2.2 不做 Q 网格也能精确求 kappa

对固定端口模式的私有组件，枚举它们与两个固定不同颜色的可能 equality patterns。若模式还使用 j 种额外颜色，其乘数为 `(Q-2)_j`；在 Q1 等于 `(-1)^j j!`。所以只需 Kreg 在 Q1 的 regular equality-pattern 值做有理求和。全部 62 个签名及其 alpha 线性系数在证书中保存。

双端验证使用另一组权：完整 exterior components 的颜色 coarsening，k>=2 个颜色块时 `(Q)_k'|1=(-1)^(k-2)(k-2)!`。它对 1,922 个恰好两共享分区逐项重算 a_xy，与两个四端口签名乘积完全相同。另在 Q=4,5,6 做 literal colour 条件和检查 (1)；不是通过有限 Q 插值推断 Q1。

## 3. 正则补全族：四路径通道不能被有限 counterterm 消去

保留基线模型族

```
K_c(Q)=K2(Q)+c(Q)K0(Q),
c(Q)=1+alpha*(Q-1)+O((Q-1)^2),
```

其中 c 解析、两端相同、alpha 实数，K2 的系数固定为 1。regularity 只要求 c(1)=1，不要求 alpha=0。canonical c identically1 是既定模型，不因以下比较而改变。

令 `H(Q)=(Q-1)K0(Q)`，则 `K_c=Kreg+alpha H+O(Q-1)`，其中最后余项更精确地为 `(Q-1)` 乘一个 regular unequal-pair tensor。闭合中 all-one-colour 项恒零；其它 colour multiplicities 在 Q1 有简单零。因此首阶 Q 双插入只取决于各 regular kernel 在 Q1 的值，即只依赖 alpha，不依赖 c''(1) 等高阶系数。

在已经具有实际 8x8 四条分离路径实现的反射接线上，两个 kernel 的颜色闭合是 Frobenius pairing。除以未插入的 `Q^4` 只做 regular normalization，不改变零值处的首阶系数。固定双线性型

```
B(L,M)=d_Q [Q^(-4)<L(Q),M(Q)>] |Q=1
```

在 `(Kreg,H)` 基底中的矩阵为

```
G = [[13/8, -1/4],
     [-1/4, 1/2]],
det G=3/4.
```

这些系数也可由基线三个 norm/cross rational formulas 直接展开：`<Kreg,K0>|1=-1/4`（可去极限），`Res_Q1 ||K0||^2=1/2`。

允许两端暂时用不同 alpha 作代数核对，得到

```
a_4(alpha_x,alpha_y)
 =13/8-(alpha_x+alpha_y)/4+alpha_x*alpha_y/2
 =3/2+(alpha_x-1/2)*(alpha_y-1/2)/2.                 (4)
```

同一微观补全使用同一个 alpha，因而

```
a_4(alpha)=3/2+(alpha-1/2)^2/2 >=3/2.                (5)
```

不等式在 alpha=1/2 达到，但**这里不改变 canonical alpha=0，不推荐调参到该点**。canonical 的 13/8 已在基线得到；新的内容是整个上述正则补全族的下界。

从代数上看，沿 H 消去交叉项的 Schur complement 是

```
13/8 - (-1/4)^2/(1/2) = 3/2.
```

这个有限 Gram 的正性不是一般解析延拓的酉性或 CFT metric 正性；不同接线或不同两端 alpha 不必正。例如 alpha_x=-3、alpha_y=3 已使 (4) 为负。把 H 商掉是一个补全敏感性诊断，**不声称 counterterm 是物理 gauge**。

### 3.1 可以排除的具体解释

可以排除：在 pair 系数保持 1、两端使用相同实解析 c、且 c(1)=1 的家族内，通过选择 c'(1) 使所有双插入首阶 Q 响应消失。四路径接线始终留下至少 3/2。

不能据此排除占据平均的符号抵消。若 E4 是这个指定反射 wiring 的事件，该受限层的贡献为 `a_4(alpha) P(E4)`，但其它 exterior wirings 可以为负。不存在从 (5) 到整个 C_xy、再到 W[a_xy] 的自动正下界。也不能把数值 13/8 或 3/2 当作尺度指数。

## 4. 对下一项研究的 P0 建议

### P0-A：预测带符号两桥统计，而不是普通单点源 covariance

在**同一 canonical 核**下，用固定共享-component incidence 定义

```
a_xy^(2) = 1_{both vacant, s=2} kappa_x*kappa_y,
a_xy = a_xy^(2)+a_xy^(3)+a_xy^(4).
```

s=3、4 的系数直接复用已有固定 Bell8 表，不新增 descriptor 或 outcome-selected kernel。三层要严格加回原源，保留同批 covariance。理论预测必须针对带符号平均或原 `W[a_xy]`，不能只预测至少两桥的概率就宣称预测了原 U。

**首份有意义的交付：**事前选定两分离尺度/几何上的一个带符号比值、符号或可控上界，并说明 s=3、4 的余项为何不改变判断；或给出一个无法压低余项的反例。这里尚未取得该空间预测，没有据本代数自动开启生产。若只是更精确地检出两桥事件，没有解决符号/余项，则不升级为机制解释。

### P0-B：以补全稳健的接线约束筛选机制，停止 counterterm 救场

候选若将 canonical 正则相互作用压缩成独立单点权重，必须解释为何已知四路径通道 (5) 可以在所声称的观测与系综层面被省略，并给相应余项；不能通过重新选择 c'(1) 把其条件双响应调为零。

**首份有意义的交付：**在同一原占据测度下，对指定反射四路径与其它 signed wiring 的竞争给出预测，再经过原 W 完整传播。无需重做旧四路径13/8、旧单点W或此前Q4/Q1 seam score；本下界也不替代 generic-Q 到尺度场的实际映射。

若局部参数统一为 `epsilon/N`，无序对进入二阶 epsilon 导数时必须带 `2/N^2`；独立 lambda_x、lambda_y 下没有该因子2。原 W 的归一化、共同 root motion 与 slope response 四项均保留。

## 5. 验证、复算与证据边界

新包仅为标准库 Python、证明与小型有理证书；没有导入未合并分支代码，没有新 occupation population、MC、服务器任务或物理结果重评分。

```
python3 scripts/p337_bridge_completion_review.py --verify results/p337-p0-bridge-review/exact.json
python3 -m unittest discover -s tests -p 'test_p337_bridge_completion_review.py' -v
```

17 项 focused tests 已实际通过；全部 1,922 个两共享分区因式分解、25 个不同 alpha_x/alpha_y 核对、literal integer-colour 求和、norm identities、C4/桥标签对称、错误输入拒绝均有检查。证书生成后做独立确定性复算；输出命令拒绝覆盖旧文件。**未运行全仓测试或宣称远端 CI 通过。** 本文的一般等式依赖解析证明，有限核验不替代证明；所有新系数属于相同 kernel algebra，不是新增独立统计票。

该结果不提供 C_xy 的实际空间衰减、全局 U 新数值、square-site arm exponent、连续场数目或 Jordan 身份；不证明整个 regular completion 空间有统一下界，只证明 §3 的明确两-sector 家族与归一化。用户授权的 P0 是建议等级，不改写科学 claim ledger。

## 固定来源

- [kernel/completion algebra @2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/local-pair-two-insertion-algebra.md)
- [actual four-path wiring @2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/local-pair-two-insertion-geometry.md)
- [completed one-mark result @2ba8863f](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)
- [at-most-one support and full W interface @baa5d33b](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/PROOF.md)
- [canonical Bell8 bounds and actual two-component witness @baa5d33b](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md)

外部背景仅作研究接口：Couvreur/Jacobsen/Vasseur, arXiv:1704.02186（S_Q 张量与多组件连接）；Vasseur/Jacobsen/Saleur, arXiv:1206.2312（渗流对数观测与算符混合）。本次没有凭这些文献给当前核赋予已验证的连续场身份；外部论文不是以上有限恒等式的证明替代。
