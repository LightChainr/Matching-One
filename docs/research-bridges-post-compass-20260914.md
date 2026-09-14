# 研究桥接推进：从校准文件到两个可判定机制

2026-09-14。本文接在 `docs/research-compass-beyond-exactness-20260914.md` 之后，目的不是再增加一层任务管理，而是把本轮真正改变判断的推理压成两个 bridge：

1. Matching root 的 leading non-common correction 如何被拆成 microscopic regularization mismatch、spin character、thermal tangent/normal response 与最终 CFT module；
2. fixed-subcritical common-label process 如何从已有 component Poisson proof 进入 space × intensity-clock marked Poisson cloud，而不是继续求抽象 splitting kernel 的更多后果。

现有精确结果不重复证明。以下若无额外注明，新的 continuum/module 解释均为 conjecture/programme；AGG/BK 部分是对 #739 已写证明的重新组合，需独立审计后才能升级为 author-proof。

## A. Root correction：先分 regularization coupling，再谈 field 名称

固定宽度 safe transfer 的两条对象

```text
I4,w^0(p),
I8,w^0(1-p)
```

在 critical limit 流向同一个 magnetic/topological continuum sector。已有 `x_m=5/48`、level-one momentum 与 oblique magnetic-metric controls 支持这一点。因此最自然的 RG 组织不是先假定“两个 continuum state 的 matrix element 不同”，而是把它们看作同一 continuum state 的两个 microscopic regularizations。

写

```text
E_a(ell,p)
 = ell^-1 Phi(X,
              u_0^(a) ell^y0,
              u_4^(a) ell^y4, ...),
X=t(p) ell^yt,
yt=3/4,
```

其中 `a=4,8`。matching/complement map 允许把 lattice couplings 分成 common/even 与 difference/odd：

```text
u_j^+ = (u_j^(4)+u_j^(8))/2,
u_j^- = (u_j^(4)-u_j^(8))/2.
```

lower-dimensional common corrections 可以在每条 gap 里很大，却从 `Delta=I4-I8` 消掉。matching root 只对 `u^-` 敏感。

### A1. 当前最强 finite evidence 已经先确定 spin character

oblique safe transfer 在同一个 square-site NN / complementary-matching microscopic model 上只改变 period direction，得到

```text
Delta(pc) ~ cos(4 theta) ell^-17/4,
p_root-pc ~ -cos(4 theta) ell^-4,
```

且 `(1,0),(1,1),(2,1),(3,1),(3,2)` 在 physical normalization 后 collapse 到共同 amplitude；near-node `(5,2)` 同样被显著压低。

因此当前最保守的优先排序应是：

```text
leading observed irrep: spin 4 >> spin 0,
module identity: unresolved.
```

scalar `x=21/4` 仍可作为 subleading/hidden contribution，但不应再与 spin4 作为 leading observed amplitude 的等权默认解释。真正需要测的是同阶 spin-0 projection，而不是再拟合一个 exponent 4。

一个干净的 equal-circumference pair 是

```text
axis:    u=(1,0), n=10, ell=10,
oblique: u=(3,4), n=2,  ell=10,
H4(3,4)=-527/625.
```

同一 `ell` 上写 scaled critical mismatch

```text
Y(theta)=B0+B4 H4(theta)+higher.
```

两点直接解 `B0,B4`，避免把不同 circumference 的 ordinary finite-size remainder误认成 angular harmonic。

### A2. RG-improvement interpretation

若 leading difference coupling 是 thermal-family spin4，

```text
y4 = 2-(xt+4) = -13/4,
xt=5/4.
```

则 critical mismatch 的 per-length exponent自动是

```text
y4-1 = -17/4,
```

而 root shift 是 relevant thermal tuning 对 irrelevant mismatch 的补偿：

```text
p_w^*-pc ~ ell^(y4-yt)=ell^-4.
```

这把 matching root 解释成有限尺寸的 improvement condition：用一个很小的 thermal displacement 抵消两个 regularization 之间第一个 dual-odd irrelevant coupling。它也直接解释 square/kagome 的 4 vs 6：若更高旋转对称性强制 `u_4^-=0`，下一允许 coupling 才决定 pseudo-critical shift。

### A3. 更强猜想：TANGENT_SPIN4

仅有 root shift 仍不能说明 spin4 是否改变 universal charge curve 的形状。设 physical dimensionless charge function

```text
Q_{ell,theta}(X)
 = F(X)+u4(ell) H4(theta) G4(X)+...,
u4~ell^-13/4.
```

一般分解

```text
G4(X)
 = [G4(0)/F'(0)] F'(X) + G4_perp(X),
G4_perp(0)=0.
```

第一项只平移 root；把每个方向自己的 root 移回0后，只剩 `G4_perp`。

强猜想为

```text
G4(X)=c4 F'(X),
```

即 leading spin4 对这个 observable 完全沿 thermal tangent。则

```text
Q_{ell,theta}(X)
 = F(X+c4 u4 H4)+o(u4),
```

所以会同时出现清楚的 orientation-dependent root shift 与近乎为零的 fixed-b/intrinsic normal response。

最直接的检验不是更多 root，而是 root-centered oblique curves。进一步用各自 root slope 归一 thermal metric：

```text
s_theta=dQ/dp|p*,
Psi_theta(y)=Q_theta(p*+y/s_theta).
```

若 `Psi_axis-Psi_oblique` 的 cos4 component 比 raw critical/root signal 明显降一个阶，则支持 tangent mechanism；若仍保持 `ell^-13/4`，则 spin4 有真正 normal shape component，后续 module/matrix-element identification才有高价值 target。

这也是 #802 中 `(T,N)` 的一个直接 actual-site realization：

```text
T = root translation,
N = root-centered, slope-normalized shape response.
```

### A4. c=0 module identity：generic-Q 解简并比 Q=1 贴标签更干净

generic-Q Potts 表示论给一个关键事实：diagonal `(h_{r,1},h_{r,1})` Kac fields 的 null descendants 在 generic Q 上是真正 quotiented；energy `Phi_{2,1}` 因此有明确的 Kac module。到 `Q=1,c=0`，energy `Phi_{2,1}` 与 2-hull `Phi_{0,2}` 在 `(5/8,5/8)` 碰撞并形成 logarithmic pair。

于是它们的 spin4 descendants也在 `x=21/4` 碰撞：

```text
Q4 Phi_{2,1},
Q4 Phi_{0,2}.
```

所以 `Q=1 + x=21/4 + cos4theta` 仍不能决定 module。若 generic-Q 两支 amplitudes 在 `Q->1` regular，则得到 ordinary `ell^-17/4`; 若系数以 `±1/(Q-1)` 相消，极限会生成

```text
ell^-17/4 [A+B log ell].
```

因此 module/Jordan 判别应优先做 generic-Q branch tracking 或 Q-velocity，而不是自由拟合 `A+B log L`。这条路线与 TANGENT_SPIN4 正交：前者问 field/module，后者问该 field 对目标 observable 的 correction function 是 tangent 还是 normal。

## B. Common-label bridge：原 Poisson proof 可能已经包含 no-merger 与 marked process

#780 定义 early parameter `p_-` 与 final parameter `p_+`。对每个 `p_+` final essential component C，令 `n_C` 是它包含的 distinct `p_-` essential lineages 数，并定义 merger factorial density

```text
M_w(p_-,p_+)
 = E sum_{C: final anchor row 0} binom(n_C,2).
```

原任务把 `M_w/nu_- ->0` 当作新的 blocker。重新读取 #739 的 `poisson-birth-windows.md` 后，这一步似乎可以由旧 AGG/BK/localization 估计直接推出。

### B1. 短 final component：merger pair 是旧 `b2` pair 的子集

取 `H=w^2`。若 final C 的 span `L(C)<=H` 并含两个 early components，则它们的 unique anchors 行距 `<=H`；early components 本身也是 C 的子集，所以 span也 `<=H`，由 #739 的 localized anchor indicators 精确计数。

#739 §6 已证明，对任意两个重叠 anchor windows 中的不同 full winding components，各自存在 disjoint occupied winding witnesses，site-BK 给

```text
E[I_i I_j] <= U_w(p_-)^2,
```

其中

```text
U_w(p)=poly(w,H) exp[-(w-1) kappa(p)].
```

每行的 close-pair 数至多 `O(w^2 H)` 个，所以

```text
M_w^short <= poly(w) exp[-2 kappa(p_-) w].
```

这就是原 component-Poisson proof 的 `b2/m`，只是按 lineage-merger 语义重新读取。

### B2. 长 final component：旧 localization tail 是 superexponential

对 `L(C)>H`，粗界

```text
n_C <= |C| <= w L(C)
```

给

```text
M_w^long
 <= (w^2/2) E sum_{C:anchor0} L(C)^2 1{L(C)>H}.
```

#739 §4 的 strip-crossing proof 对任意 `r>=w` 给 unnormalised span-tail intensity

```text
E[# {C:anchor0,L(C)>r}] <= C_I w exp(-c_I r)
```

uniformly on compact subcritical intervals。用 tail-sum identity 得

```text
M_w^long <= poly(w,H) exp(-c_I H)
          = exp[-Theta(w^2)].
```

因此

```text
limsup (1/w) log M_w <= -2 kappa(p_-),
M_w/nu_w(p_-) ->0,
```

因为已有 `-log nu_w/w -> kappa(p_-)`。这甚至不要求 `p_+-p_-=O(1/w)`；任意固定 compact-subcritical pair 都有相同尺度分离。

若审计无缺口，no-macroscopic-merger 不再需要新的大计算。

### B3. no-merger 以后，mark process 直接继承同一个 AGG dependency graph

在 bounded intensity-clock window 上固定 top parameter `p_+`。对每个 localized final component C，用同一份 U-label filtration定义 first essential birth mark `tau(C)`。rare merger可先任意 tie-break，其总质量由 `M_w` 控制。

由于 final component 不碰 guard rows，所有较低 p 的 ancestors 也是同一 local window 内的子图；`tau(C)` 是 local label function。把 parameter window 分成有限 bins `A_r`，定义 typed indicators

```text
I_{j,x,r}
 =1{final p_+ anchor at (j,x), tau(C) in A_r}.
```

不相交 windows exact independent；重叠 windows 的 pair bound仍由两个 p_+ disjoint winding witnesses控制。因此 AGG process theorem 原样给有限 mark bins 的 independent Poisson approximation。

### B4. mark intensity 不需要 p-analyticity

对任意中间 p，每个 p-essential component 唯一属于一个 p_+ final component，故单位长度精确有

```text
nu_w(p)=E sum_C n_C(p).
```

而 birth mark `<=p` 的 final-component intensity 是

```text
mu_w((-,p])=E sum_C 1{n_C(p)>=1}.
```

逐 C 有

```text
0 <= n_C-1{n_C>=1} <= binom(n_C,2),
```

所以

```text
0 <= nu_w(p)-mu_w((-,p]) <= M_w(p,p_+).
```

若参数仅按 intensity ratio 定义

```text
nu_w(p_w(x))/nu0 -> Lambda(x),
```

则 `M_w/nu0->0` 立即给

```text
mu_w((x1,x2])/nu0
 -> Lambda(x2)-Lambda(x1).
```

这已经识别 marked Poisson 的 mean measure，不需要先证明 `nu(p)` analytic，也不必先写 `p_w=p0+x/(vw)`。

birth-anchor 与 final-anchor 的位置差至多 O(H) on localized event，而

```text
H nu0 = w^2 exp[-kappa(p0)w+o(w)] ->0,
```

因此在 rescaled longitudinal coordinate 上 anchor motion 自动消失。

若 finite-mark partition -> PRM 的标准升级及 torus/cylinder coupling无隐藏问题，目标 birth cloud应为

```text
sum_C delta_(nu0 y_C, Lambda_w(tau(C)))
 => PPP(dy dLambda).
```

此后 #785 的 exact clock kernel/Laguerre hierarchy才成为已建立 model map 的 consequence，而不再是 bridge 本身。

## C. 近期待验证顺序

按信息增量排序，而不是按文件数：

1. **#780 proof audit first**：审计 B1/B2 的四个 pathwise/counting 点；若通过，写成独立 lemma，并立即尝试 B3/B4 的 typed-AGG marked process。此路线理论收益远大于再做 common-label static tables。
2. **#802 equal-ell / centered-curve test**：axis `(1,0),n=10` 与 `(3,4),n=2` 共用 `ell=10`; 同时输出 critical mismatch、thermal derivative、root 与 5--7 个 root-centered curve points。先判 `spin0 vs spin4`，再判 `tangent vs normal`。
3. **#768 改用 RG factorization 语言**：区分 microscopic odd coupling `u_4^-` 与 universal magnetic matrix element；不要把两者都叫“sector matrix-element difference”。
4. **#586 generic-Q reserve**：只有 actual-site centered-curve 说明 spin4 normal component确实需要 field identity时，再投入 generic-Q energy/2-hull branch tracking。若 tangent-only 已解释当前 observable，则 module研究仍有数学价值，但不是 root 机制的前置。

## D. 当前大胆但可失败的统一猜想

一个较强的工作图景是：

```text
finite matching topology
    -> fixes which global channels can differ;
microscopic regularization mismatch
    -> supplies a small dual-odd spin4 coupling;
thermal RG response
    -> converts it into the L^-4 pseudo-critical root shift;
rare complete-component AGG/BK structure
    -> supplies a Poisson space×intensity birth cloud;
record/splitting kernels
    -> are consequences of that cloud, not independent mechanisms.
```

在这个图景中，项目近期真正未知的低维对象很少：root 侧是 `spin irrep + tangent/normal + module` 三层；common-label 侧是 `local marked birth map` 一层。若后续验证继续支持这种压缩，应主动减少派生模型任务，而把资源集中到能击穿这两个 bridge 的反例上。

## Claim boundary

- oblique root/free-energy values、existing AGG/BK/localization inequalities 是现有 branch 的 deterministic/author-proof inputs，本文不重新认证。
- `regularization mismatch` 与 `TANGENT_SPIN4` 是新的 RG interpretation/conjecture。
- energy–2-hull generic-Q collision 是文献支持的 representation-theory input；其对 Matching-One safe-sector correction 的具体投影仍未识别。
- `M_w/nu_->0` 与 marked-PPP route 是从现有 #739 proof 组合出的 candidate author-proof；在逐行审计写入 branch 前不得升级成已证 theorem。
