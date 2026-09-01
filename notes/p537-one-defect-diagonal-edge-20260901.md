# #537 full-graph one-defect 判决：contact/collision diagonal edge 存在

日期：2026-09-01

基线：`090cdc80`
状态：本地隔离分支 exact 结果；未 push、未合并、未使用服务器。

## 判决

Issue #537 的 blanket full-graph 命题已经被一个 literal physical edge 否定：同一
N25 axis 构型只把一个站点 `z` 从 0 翻到 1，就同时改变

```text
landing/rank:  q+1 = 0 -> 1
source/Bell:   9240712 -> 6848576
source value:  g16 = 4 -> 0, hence Delta a = -1/100.
```

该 edge 在完整两几何共同冻结的 pooled root、`R`、`mu_a`、`beta_diag1` 下，C4
物理轨道和 geometry pool `1/2` 都恢复以后，Eq. (10) 的 signed weight 为

```text
source-midpoint part   -1.0888815582478189e-11
root-Schur part        +2.5901918540035547e-12
full                   -8.298623728474635e-12
```

三个有理区间都严格不含 0。source-midpoint 部分自己已经非零，且 `Delta a=-1/100`；
因此结论不依赖把 `-beta B` 如何分配到 Bell cells，也不是“Bell 标签变了但 kernel
值没变”的伪 diagonal edge。

这触发 stop rule A：**完整物理图不能自动分解成两个独立空间 defects，所以不能从
双变量 interaction 直接领取 six-arm 增益。** 无需继续枚举完整 graph 来执行这个
存在性停止规则。

## 图和 signed mass 的精确定义

固定 geometry `g`、off-`z` 背景 `eta`、source pair component `lambda=(x,y)`，两个节点是

\[
v_i=(g,\eta,z,i,r_i,\mathcal B_i,C_i,B_i,W_i),\qquad i\in\{0,1\}.
\]

- `r_i=q_i+1` 是有限 rank index；
- `C_i` 是 `x4+y4+z4` 十二端口在同一全局 carrier 图上的 first-occurrence canonical
  component map；
- `\mathcal B_i` 是 `C_i` 限制到 `x4+y4` 后重新 canonicalize 的 Bell key；
- `B_i` 是四个 `z` cardinal ports 的 global NN-black partition；
- `W_i` 是删去 `z` 后四个 thermal ports 的 matching-white partition。

物理边 `e=(v_0,v_1)` 只改变 `X_z`。完整总体先冻结 `p,R,mu_H,mu_a,beta_lambda`，然后定义

\[
P_i(e)={1\over2}\nu_{g,-z}(\eta)w_i>0,
\]

\[
S_i(e)=P_i(e)\widetilde H_i
\{(a_i^\lambda-\mu_{a^\lambda,g})u_i-\beta_\lambda b_i\},
\qquad W(e)=S_0(e)+S_1(e).
\]

结果文件同时保存一个 fixed direction 的 `P0/P1,S0/S1` 和恢复四方向 C4 orbit 后的
signed masses。没有在 cell 内重估 root counterterm，也没有对零质量 cell 做密度除法。

独立 Bell key 足以描述这个**存在性** witness，因为同一个 producer 保存了唯一物理
背景、唯一 `z` flip 和两个 Bell keys；共同 `C` map 又逐端口证明 Bell 正是同一 joint
fibre 的限制。若要证明“没有 diagonal edge”或证明两个 defects 位于可分离 annuli，
Bell key 单独不够：它忘记 `z` ports 与 source ports 是否属于同一 outer component，
此时必须保留这里的 joint `C/B/W`。

## literal witness

```text
geometry       axis N25, quotient (5,0)
x,z,y          0, East(x)=1, y=(-1,-1) [vertex 24]
eta mask       12567 over the frozen 23-site order
k_minus        7
joint C        23090870354448 -> 92359816642816
Bell           9240712 -> 6848576
q              -1 -> 0
g16            4 -> 0
```

producer 扫描 12,568 个背景后命中并立即停止。Python oracle 用独立 DSU/Euler/Bell
实现重建同一构型；scorer 再从 joint key 限制重建 Bell，不信任 producer 的 Bell 字段。

## 必须保留的边界

该 witness 是 contact/collision edge：`z` 紧邻 `x`，`y` 为 diagonal source；
`arm_mask=3`，`local_contact_mask=3`，source 确实接触 global black landing component。
它不是 alternating ordinary separated four-arm，`J_B+J_W-1` identity 在这里不适用，
也不证明任何渐近下界。

因此被删除的是 **blanket full-graph two-independent-defect route**。separated sector
仍可单独研究，但先采用无自由参数的固定分解

```text
contact/collision: d_NN(z,{x,y}) <= 1
separated:         d_NN(z,{x,y}) >= 2.
```

下一 P0 对象应是 contact/collision sector 对 surviving leading four-arm signed functional
的贡献。separated sector 只有在进一步证明 row/column changes 可定位到不交 annuli 后，
才允许使用 two-defect/six-arm bound；距离条件本身还不是 annular certificate。

只读 joint-sector 聚合中的同-sector determinant 是多条 edges/cells 的 rank-two 量，不等于
单条 physical diagonal edge。本判决不继承该未提交 scorer，而使用 count-one literal
witness、共同 joint map 和 allocation-robust source part。

复现入口：

```bash
python3 experiments/p537-one-defect-gate-20260901/verify.py
```

机器可读结果在 `results/p537-one-defect-gate-20260901/result.json`。
