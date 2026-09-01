# #537 full-graph one-defect 判决：typed-carrier diagonal edge 存在

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

因此被删除的是 **blanket full-graph two-independent-defect route**。最初按 NN 距离切分
contact/separated 的方案也已被下一条固定见证否定；距离不是 carrier identity。

## 距离二仍是 joint-incidence contact

第二个冻结构型在 N5 row-major 坐标取 `x=0,y=6,z=2`，off-z occupied cells 为
`{1,3,4,5,7,9,10,12,15,16,17}`。独立 quotient 重建把三点映到 internal
`(0,4,3)`，三对 NN 距离均为 2，但一次 physical `z` flip 仍给出

```text
rank                         0 -> 1
Bell                         274568 -> 8256
joint C                      21990249529872 -> 3298535014656
joint terminal incidence     2 -> 1
g16                          8 -> 0
```

该 edge 的 `arm_mask=11`，source 接触 global black landing component，off-state 没有
non-port degree-3 branch。恢复同一 globally frozen root/counterterm 与 C4 pool 后，source
midpoint 为 `-1.0121115955209059e-10`，root counterterm 为
`+5.341390620686106e-12`，full Schur weight 为 `-9.586976893140449e-11`；三项都有
严格不含零的有理区间。

所以 `d_NN<=1` 不能定义完整 contact 类，`d_NN>=2` 也不能认证 separated sector。
这个 N5 距离二构型仍共享 source/thermal joint carrier，正确有限分类是
`joint-incidence/typed carrier`。真正的 separated sector 至少要先排除共同 carrier incidence，
再证明 row/column changes 可定位到两个不交 annuli；本见证不提供该证书或渐近下界。

## 与固定 radius-one selected sector 的关系

平行固定 commit
[`df4a64f6`](https://github.com/LightChainr/Matching-One/blob/df4a64f68232eec5aa5b8c8a5d920062aaa7808e/results/p537-one-defect-diagonal-edge/REPORT.md)
在 N25 axis/tilted、`x=West(z)` 的 alternating radius-one **selected sector** 中保留
6,846 个 kernel-changing row classes、740,950 条 physical fibres，总 Schur signed mass
为 `-4.948839916450813e-6`。12 个 rank-stage×source-orbit cells、两个 row sums 与六个
column sums均严格不含零；`0→1` 给 117.63%，`1→2` 给 -17.63% 抵消。它不是完整
physical graph 或完整 `T_N`，也来自同一个 N25 dependency block。

在同一 selected sector 内，contact mask 0 有精确的 0 classes / 0 fibres；mask 1、2
分别给约 `-4.09370e-6/-4.11267e-6`，mask 3 给 `+3.25753e-6`。这个空集不能外推到
更大 collar、其他 endpoint 或完整图。该 branch-only asset 的旧 fibres 没有单一背景或
canonical joint `x+y+z` map；本分支的两个 literal witnesses 补上构造性 joint-map 信息。
这些结果互补，但不能当独立渐近 evidence votes。

`df4a64f6` 还联合保留 contact mask、birth stage 与 source orbit。在这个 selected sector
中，one-arm masks 1/2 只支持 NN source，所有 non-NN source 都落在 mask 3。把两个并不
完全相同的 one-arm masks 仅作为 single-contact aggregate 合并后，globally frozen
pooled-root Schur signed table 为

```text
                 single-contact aggregate    double contact
0->1                 -2.8838028012e-6       -2.9372878647e-6
1->2                 -5.3225654942e-6       +6.1948162436e-6
```

determinant 的严格区间为
`[-3.3498535471290615e-11,-3.3498535471290614e-11]`。这只证明该有限 selected
contact×birth signed table 是 rank two；它不证明一个完整 operator 已经识别，也不是 exact
exchange-even projection。

下一 P0 不再枚举 contact descriptor 或距离层。当前最小理论对象是
`typed joint-incidence/contact-fusion × topological completion`，目标是证明它怎样投影、
抵消并缩放到 full original-`U`，最终仍需得到 `T_N=o(D/A_N)`。只有先冻结这样的精确
传递恒等式，才允许在一个 held-out size 上检查同一 tensor/determinant；不允许用单坐标救场。

selected-sector aggregate determinant 是多条 edges/cells 的 rank-two 量，不等于单条
physical diagonal edge。两个 literal gates 使用固定单一背景、共同 joint map 和
allocation-robust source part；它们与 aggregate 回答不同层次的问题。

复现入口：

```bash
python3 experiments/p537-one-defect-gate-20260901/verify.py
```

机器可读结果在 `results/p537-one-defect-gate-20260901/result.json` 与
`results/p537-one-defect-gate-20260901/result-nonadjacent.json`。
