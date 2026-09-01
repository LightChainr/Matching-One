# P337 canonical 空间核在方格临界点绝对可求和

**结论。** 对真正的无限方格 site-percolation 临界点 `p_c`，canonical
`Kreg` 的两点首 Q 空间核不是一个可由远距离累计形成发散响应的源。存在
常数 `eta>0`、`C<infinity`，使标准平方 `L x L` 环面上

```text
E_pc |g_xy^(L)| <= C (1+d_L(x,y))^(-2-eta),
sup_(L,x) sum_(y!=x) E_pc |g_xy^(L)| < infinity,
sup_(L,x) sum_(d_L(x,y)>R) E_pc |g_xy^(L)| <= C R^(-eta).
```

因此它的每体积齐次二阶自由能系数有绝对收敛的无限体积极限。这个结果
排除的是固定 lattice normalization 下的 marginal/relevant 长程累计解释；
它没有给原 `U` 的渐近定理，也没有识别 H4、Jordan 或其它连续场。

没有为本结论新增样本、枚举、服务器任务、距离点、补全系数或拟合指数。

## 1. 已消费的最新执行进度

以下问题已完成，不能再从旧入口重复立项：

- [`a237968f`](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md)
  已用两个新依赖块拒绝 L64、距离16的有限零传递；L32/距离8与
  L64/距离16的比值只是一项冻结的两点比较，没有拟合空间指数。
- [`410015f5`](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md)
  已完成 `J2=d_logQ d_epsilon^2 U`：总量
  `-0.0055194314248394015`，其中相邻项
  `-0.0017510744544027990`、其余非相邻项
  `-0.0037683569704366022`。additive-linear closure 和
  NN-contact-only closure 均已被排除；非相邻类包含短对角，不能改称长程尾。
- Draft PR
  [#532](https://github.com/LightChainr/Matching-One/pull/532) 的提交
  [`2e1c57b`](https://github.com/LightChainr/Matching-One/blob/2e1c57b42129c507b45fd1a8111eb32e78b1f80b/notes/P0-p337-two-bridge-and-completion-review-20260901.md)
  给出恰好两共享组件的端点签名因式分解，以及指定同一两-sector 补全族中
  四路径系数的 `3/2` 下界。这些有限代数结果有用，但不能给占据平均的
  符号。该 PR 后续评论中“再计算 uniform J2”的建议已被 `410015f5`
  超过；评论中的八通道压缩和 `alpha=3/2` 两点非负补全尚无同一 PR 的
  新提交，不登记为当前已交付资产，也不据此改变 canonical `alpha=0`。

## 2. 确定性双四臂约束

令 `d=d_L(x,y)` 为环面 `l_infinity` 距离。对 `d>=16` 取
`r=floor(d/4)`。两个 `B_(r+2)` 球忠实嵌入平面且不相交；这个判断在
`d=L/2` 仍成立。近距离、相邻点和小环面只含一致有限数量的位移，稍后用
Bell8 点态界吸收。

既有精确核已经证明：非相邻 x、y 至多共享一个 exterior component 时
`g_xy=0`，且始终有 `|g_xy|<=43/16`。在非相邻端口集上，两个共享
exterior components 必须是两个不同的原 occupied NN components。若两个
这样的组件都接触 x、y，则在每个标记周围的 `A(2,r)` 中：

1. 两个不同组件各给一条互不相交的黑 NN 径向臂；
2. 它们在局部也不能由黑 NN 路连接；
3. 两条黑 crosscut 分出的两个扇区，各由 square/matching 离散对偶给一条
   白 matching 径向臂。

所以每点都出现黑、白、黑、白交替四臂。令 `pi4^square(p;2,r)` 表示这个
方格 site 四臂事件的概率。丢弃全局共享组件条件后，两个环带及 x、y
读取的顶点集合互不相交，在 Q1 Bernoulli 乘积测度下无条件独立，故

```text
E_p |g_xy|
 <= (43/16) Pr_p{x,y vacant and >=2 shared occupied components}
 <= (43/16)(1-p)^2 [pi4^square(p;2,r)]^2.          (1)
```

下文记 `C_L(x,y)=E[g_xy^(L)]`，所以同一界也控制 `|C_L(x,y)|`。
独立性只在丢弃外部 wiring 后成立；没有声称给定共享组件、K、rank 或
外部颜色后两个局部事件仍独立。

## 3. 方格 site 的严格四臂界足以推出可求和性

van den Berg--Nolin 直接研究临界方格 site percolation，黑臂用 NN、白臂
用 matching 邻接，与式(1)的模型相同。他们在第5.2节得到有限尺度不等式

```text
pi4(3n) <= C n^(-1) sqrt(pi2(n)).                 (2)
```

来源：[On the four-arm exponent for 2D percolation at criticality](https://ir.cwi.nl/pub/31186/31186.pdf)，
Theorem 1.1 与式(5.6)后的结论。RSW 在不交的几何级环带上给某个
`rho>0` 的 `pi2(n)<=C n^(-rho)`，所以

```text
pi4(n) <= C n^(-1-rho/2).                         (3)
```

他们的 `pi4(n)` 使用固定内半径1，而式(1)使用内半径2。这里需要保留一个
固定微观尺度接口：四条臂在 `partial B_2` 的着陆型只有有限多个；在
`B_2` 的有限顶点集上作有界局部修改，可把每种交替着陆延伸到
`partial B_1`。Bernoulli 有限能量和有限原像数给

```text
pi4^square(p_c;2,r) <= C0 pi4^square(p_c;1,r),
```

常数与 r 无关。这里不需要宏观臂分离，也没有把三角格的 `5/4` 指数移植
到方格。结合(1)--(3)，取 `eta=rho>0`，得到

```text
E_pc |g_xy^(L)| <= C d^(-2-eta).                  (4)
```

平方环面上与 x 的 `l_infinity` 距离为 m 的点至多 `8m` 个。于是

```text
sum_(d_L(x,y)>R) E_pc |g_xy^(L)|
 <= C sum_(m>R) m*m^(-2-eta)
 <= C' R^(-eta),                                  (5)
```

一致于 L 和 x。近距离项由 `43/16` 的点态界吸收。

## 4. 无限体积极限和齐次耦合

临界方格 site 没有无限 occupied cluster。固定一个平面位移 z，把逐渐增大
的环面与同一无限 iid 配置局部耦合；接触固定八端口的 occupied clusters
几乎必然有限，并由有限 vacant 边界与外部隔开。因此 L 足够大时，环面与
平面的端口分区及 `g_0z` 逐配置相同。点态界给支配收敛；再用(5)交换有限
位移和尾部，得到

```text
C_L(0,z) -> C_infinity(0,z),
chi_L := sum_(y!=x) C_L(x,y) ->
chi_infinity := sum_(z!=0) C_infinity(0,z),
sum_(z!=0) |C_infinity(0,z)| < infinity.           (6)
```

若每站点使用字面局部因子 `1+lambda Kreg`，则

```text
(1/N) d_logQ d_lambda^2 log Z_L |Q1,lambda0 = chi_L.  (7)
```

同点项在这一 Q 阶为零。理由还需要 canonical 单插入系数
`beta_x(1)=0`，使 `log Z` 的同点平方从 `(Q-1)^2` 才开始；不能只用局部
因子对 lambda 线性来解释。

沿 `lambda=epsilon/N` 的路径，

```text
d_logQ d_epsilon^2 log Z_L |Q1,epsilon0 = chi_L/N,
s2(A)=N^(-2) sum_(x!=y) g_xy(A)
     =(2/N^2) sum_(x<y) g_xy(A).                   (8)
```

因此 `E|s2|<=M/N`。对任意不随插入参数变化的固定有界 occupation
observable O，在真正 `p_c` 的固定热坐标下，

```text
|d_logQ d_epsilon^2 <O>| = |Cov(O,s2)|
 <= 2 ||O||_infinity M/N.                          (9)
```

## 5. 对原 U 的边界和新的停止规则

式(4)--(9)不能解释 `410015f5` 的负 `J2`，也不能把它推出为零：

- `410015f5` 是两个 N25 原几何在 pooled root 的完整 `U` 投影；它包含
  热导数、共同根移动、斜率分母及显式 `A_N=N^(13/8)/2`，`U` 不是式(9)
  的固定有界 O。
- 其非相邻项包含所有非 NN 位移，大量贡献可以来自短对角和中程；它不是
  宏观尾部。
- 两几何的 `E[s2]` 都为正，完整 `J2` 却为负，已经直接说明 raw
  susceptibility 的符号不能决定原 U 的符号。
- `a237968f` 的 `p_ref=0.592746050790` 与 N25 pooled root 都不是本证明
  已认证的精确 `p_c`。没有近临界窗口控制时，不能把本定理事后改成那两个
  有限参数的渐近证书。

据此执行以下研究取舍：

1. **停止** canonical raw pair susceptibility 的距离网格、指数拟合、更多
   completion/alpha 扫描，以及“远距离累计发散”解释。有限距离非零、完整
   J2 非零和临界绝对可和性已经分别回答三个不同问题。
2. **下一理论闸门只看 thermal transmission。** 对固定 signed `g_xy`
   写出 `partial_p E[g_xy]` 的 site-pivotal 支持，判断热导数的空间和是否仍
   一致可求和。若可求和并能同时控制原 root/slope 分母，则该局部相互作用
   不能生成所需的渐近放大，应正式降级；若唯一某个具名 pivotal/landing
   通道逃过可求和界，先给它对同一个原 U 的符号或尺度预测，再冻结一次
   新读数。
3. 不把 PR532 评论中的更多桥 irrep、三插入或 positivity completion 直接
   升为生产。只有它们对上述 thermal-pivotal 闸门给出互斥预测，才有资格
   进入下一次合同。

这把当前未决问题从“再找一个空间 descriptor”压缩为一个明确可失败的问题：
**绝对可和的 raw canonical interaction，能否只经热/pivotal 与 moving-root
操作获得原 Matching-One 的尺度增强？**

在 exact `p_c`，式(4)还排除净 `C_xy` 中非零振幅的 `r^-2` 或
`r^-2 log r` 主导尾。严格方格定理只给某个正的 `eta`，没有给它的数值；
所以 `r^-5/2`、`r^-5/2 log r`、更快衰减、变号或别的投影仍然允许。

## 6. 证据边界

- 使用方格 site 的严格 `alpha_4>1` 结论；没有使用方格 `5/4` 普适性猜想。
- 结论是绝对上界，没有 `C_xy` 的下界、符号、非零振幅或等价渐近式。
- `chi_L` 的极限不自动给 Q1 邻域自由能解析，也不授权交换更高阶导数。
- 相邻不同点保留在近距离常数中；同点 Q1 项才为零。
- PR532、a237、410 与本定理分别属于有限代数、独立有限 MC、同一 N25
  exact population 的新 contraction、以及无新数据的解析推论，不能合票。
