# 反向审计：从单一机制叙事转向 angular irrep、source quotient 与 topological action spectrum

Date: 2026-09-14

Status: correction / synthesis note.  本文专门审计近期 #771 分支（包括本代理自己写入的若干统一化猜想）在新出现的 #47/#61/#802/#807/#808 与独立 oblique 复核之后哪些仍成立、哪些必须降级、哪些应撤回。除明确标为 exact 的代数外，不把新的组织方式升级成 continuum field identification。

## 1. 先给裁决：四个需要纠正的地方

### 1.1 `spin4 × x≈4 dressing` 不能解释 H4 projector 后的 residual

此前提出过：square root 的 `ell^-4, ell^-6, ...` 可以由一个 H4/spin-four sector 被共同的 `x≈4` scalar correction 以 `ell^-2` dressing 产生。这个说法只可以保留在 **H4 coefficient 自己的 radial expansion** 中：

```text
P4(ell)=a4 ell^-4 [1+c2 ell^-2+...].
```

如果一个 angular projector 精确满足 `L[H4]=0`，那么任何仍然正比于 H4 的 scalar dressing 也被同一个 projector 精确杀掉。故 projector 后仍存活的量不能再由“同一 H4 被 scalar dressing”单独解释。

因此旧的强解释：

```text
post-H4 residual ~= dressed H4
```

撤回。保留的弱版本只是：

```text
H4 coefficient 本身可以有 ell^-2, ell^-4,... 的 radial dressing。
```

#808 正在测的正是 **angularly orthogonal** residual；若它是真信号，应归入 H0/H8/H12/... 中至少一个，而不是重新塞回 H4。

### 1.2 “matching-even/odd continuum field” 语言必须降级

有限 matching/complement identity 精确给出的是原图与 matching 图之间的 observable/source pairing。它并不自动构造一个作用在同一 local CFT/OPE space 上的 involution。

因此本文以后区分：

```text
finite pair-exchange S/D projection : exact / empirical lattice sector;
RG tangent parity                    : hypothesis until a local map J is constructed;
OPE/interchiral parity               : stronger hypothesis.
```

`alpha<->gamma, beta fixed` 等 pathwise statement 仍可作为有限 pair-exchange 事实使用；不得仅据此称某 continuum field “matching-even/odd”。尤其 `V_<1,3>` 等 lower competitor 不能只靠未经证明的 OPE parity 排掉。

### 1.3 #802 已证明 source normalization / source redundancy 不能省略

新勘误给出逐配置恒等式

```text
H_W = N - 2K + H_B,
```

所以两个原先被当作独立的 black-pair / white-pair source 实际只差 constant + thermal score。它们属于同一个 normalized source class，不是两个独立 nonthermal directions。

同时，若 transfer 用未正规化行权 `exp(gH)`, physical free energy 必须包含 row normalizer `f(g)`：

```text
I(g)=f(g)-log lambda_tilde(g).
```

漏掉 `f_g` 会把 pure normalization/thermal motion 错当 shape response。

### 1.4 单一 `epsilon=e^{-w kappa}` witness fugacity 过强

#780 的 BK double-witness theorem 的确给 merger 一个相对 `<= exp[-w kappa+o(w)]` 的上界，但它没有证明等号；#800 的 slow doublet 是 spectral object，不是一个已经识别的概率事件；#760 的 wrap-sensitive correction 也尚未证明其最小 defect action 就是 `kappa`。

因此统一猜想

```text
all one-extra-witness effects have exponent kappa
```

降级。更稳健的对象是 defect-specific **incremental topological action spectrum**，见 §5。

## 2. 正确的分析顺序：angular irrep first, radial field second

对同一 microscopic square model、固定 physical circumference `ell`，charge root 作为 orientation 的 D4-invariant scalar 可写成离散角向展开

```text
p_ch(ell,theta)
 = P0(ell)
 + P4(ell) cos(4theta)
 + P8(ell) cos(8theta)
 + P12(ell) cos(12theta)
 + ... .
```

这里 `P0,P4,...` 只是 lattice angular irreps / Fourier coefficients；没有先赋予 CFT field、matching parity 或 radial exponent。

研究应分成两个阶段：

1. **angular spectroscopy**：尽量在 same-circle / same-modulus / same-Smith geometry 上提取 `P4,P8,...`；
2. **radial spectroscopy**：再研究每一个已经分开的 coefficient 随 `ell` 的幂、log、mixing；最后才匹配 continuum map/field。

这避免把 “一个 power” 同时解释成 angular spin、CFT dimension 与 pair-exchange parity。

### 2.1 两角 H4-null projector 的精确 leakage 公式

令

```text
h_i = H4(theta_i)=cos(4theta_i),
L[f] = (h1 f(theta2)-h2 f(theta1))/(h1-h2).
```

则 exact：

```text
L[1]  = 1,
L[H4] = 0.
```

用 Chebyshev 恒等式

```text
H8  = 2 H4^2 - 1,
H12 = 4 H4^3 - 3 H4,
```

直接得到

```text
C8  := L[H8]  = -(1+2 h1 h2),
C12 := L[H12] = -4 h1 h2 (h1+h2).
```

因此一个两角 projector 要同时近似 notch H8/H12，理想条件是

```text
h1 h2 = -1/2,
h1+h2 = 0,
```

即

```text
h1=+1/sqrt(2), h2=-1/sqrt(2).
```

N377 的 `(4,19)` / `(11,16)` 恰好接近这个双条件：它不是任意挑的两角，而是一个近似的 **H4/H8/H12 triple-notch geometry**。其当前 `C8≈0.0036146`、`C12≈-0.1378` 与上述公式一致；H16 并未被同时 notch，因此仍须保留更高 harmonic adversary。

### 2.2 命名纠正：`p_perp4` 不是自动的 `p_H0`

只要 `C8,C12,...` 非零，两角 H4-null combination 严格只能叫

```text
p_perp4 = L[p_ch],
```

而不是已经纯化的 `H0/scalar root`。

#808 的 double-notch 使 scalar 解释更可信，但若最终把结果命名为 “H0” 或 `V_<1,4>`，仍需要：

- 数值上控制 H8/H12/higher 的允许贡献；或
- 理论上给这些 harmonic 的 radial 下界/selection rule。

否则正确结论只能是 `post-H4 angular-orthogonal residual`。

## 3. 新的 working spectrum：不要再把 square 的 `4,6,...` 当成一个序列

当前最合理的分层是：

```text
P4(ell) : leading root harmonic ~ ell^-4, 已有强 deterministic evidence;
P0(ell) : post-H4 scalar candidate，#808 检验 ~ ell^-7;
P8(ell) : 独立 harmonic，radial exponent 尚未由数据确定;
P12...  : 更高 adversaries，优先作为 leakage bounds 而非默认新任务。
```

`P4` 内部仍可能有 `ell^-6` 等 dressing；但这些项在 exact H4-null projector 中全部消失。

因此 #47 的历史 `L^-7` 与 #808 的 scalar `ell^-7` 候选，不应再用“next spin4 descendant”解释。普通 thermal spin4 quasiprimary tower的下一项也并不自然给 relative `q=3`。若 scalar `ell^-7` 被确认，`x=33/4` (`x-x_t=7`) 是一个结构上匹配的 candidate，但 field identity、log collision 与 pair-exchange tangent action仍未解决。

更低的 interchiral/scalar competitor（例如 #61 保留的 `V_<1,3>` adversary）必须继续保留，直到真正构造 RG tangent map 或 radial data 排除。

## 4. Source space 的 exact 修正：Bernoulli chaos grading

这里有一个可以逐配置定义、但不冒充 RG parity 的微观 source grading。

在 primal model 参数 `p` 上定义标准化 Bernoulli variable

```text
xi_v = (n_v-p)/sqrt[p(1-p)].
```

matching/complement model 使用 `q=1-p` 与 `n_hat=1-n`，则

```text
xi_hat_v
 = (n_hat_v-q)/sqrt[q(1-q)]
 = -xi_v.
```

因此对任何有限 site set `A`，Hoeffding/Walsh chaos basis

```text
Psi_A = product_(v in A) xi_v
```

在 exact complement pairing 下满足

```text
Psi_A -> (-1)^|A| Psi_A.
```

这是 microscopic source-function space 的精确 grading，不是 continuum OPE parity。

### 4.1 #802 source redundancy 被这个 basis 一眼解释

对一条相邻 pair `(i,j)`：

```text
n_i n_j
 = p^2 + p sqrt(pq)(xi_i+xi_j) + pq xi_i xi_j,

(1-n_i)(1-n_j)
 = q^2 - q sqrt(pq)(xi_i+xi_j) + pq xi_i xi_j.
```

两者的 degree-2 chaos **完全相同**；差异只有 degree 0 normalization 与 degree 1 thermal score。故 modulo `{1, K}` 后 black-pair 与 white-pair source 必然相同。这正是独立核验发现的 `H_W=N-2K+H_B`。

### 4.2 新的 source design rule

在花算力比较两个“物理源”之前，先做 exact decomposition：

```text
source = constant + thermal/one-site part + genuine higher chaos.
```

- constant：只改 normalizer；
- uniform degree-1：thermal reparameterization；
- higher chaos：才可能提供新的 normalized source direction。

因此 #802 下一份真正独立 source，至少应在 quotient `source / span{1,K}` 中线性独立。

最简单的 exact pair-exchange-odd nonthermal source来自 **odd chaos degree >=3**，而不是另一个 pair-count source。它可以作为 #61 RG-tangent programme 的 microscopic basis element；但它是否流向某个 continuum “odd field”仍须另外证明。

## 5. 从 witness number 改成 incremental topological action spectrum

固定严格亚临界 `p`。对一个 baseline topological event `B_w` 与额外 defect `E_w`，定义增量 action（若极限存在）

```text
sigma(E|B)
 = lim[-1/w log P(E_w and B_w)/P(B_w)].
```

baseline one-winding 的 action 是 `kappa(p)`。

### 5.1 已有 rigor 只给部分 inequality

#780 对 merger 使用两个 vertex-disjoint essential witnesses与 BK，得到绝对双见证概率至多 `exp[-2 kappa w+o(w)]`；除以 baseline `nu~exp[-kappa w+o(w)]` 后得到

```text
sigma_merger >= kappa
```

（作为 liminf / upper-probability statement）。**没有证明 equality。**

### 5.2 #800/#760 应各自有自己的 action

定义概念上：

```text
sigma_split : slow-doublet tunnelling action;
sigma_wrap  : cylinder/plane mass-locality wrap-defect action;
sigma_merge : lineage-merger incremental action.
```

目前只有经验上 `Delta gamma/(w nu)=O(1)` 暗示 `sigma_split≈kappa`；这不是证明。`gamma_w-kappa` 的 leading wrap defect 也可能被 symmetry cancellation 推到更高 action，或由不同 constrained network 控制。

故不再默认

```text
sigma_split=sigma_wrap=sigma_merge=kappa.
```

### 5.3 更大胆但更稳健的统一猜想：Wulff network action ratios

对每一种 topological defect type `T`，猜想存在一个 constrained Wulff/network action

```text
sigma_T(p),
```

并且在 near-critical isotropic limit

```text
sigma_T(p)/kappa(p) -> c_T,
```

其中 `c_T` 只依赖 defect topology/network geometry，而不必是整数。

于是 massive tail 的自然 expansion 应是

```text
sum_T exp[-c_T s] * P_T(s),
s=w kappa(p),
```

而不是简单的整数 powers `exp[-j s]`。

#758 的 loop/branch variational rate本身已经提示网络代价不必等于整数 witness count。

`w kappa -> infinity` 仍然是 #780 pure-fragmentation 的一个 **充分** 远尾条件，因为已有 merger upper bound会消失；但不能据此宣布 `s=w kappa` 单独控制所有 spectral/locality corrections。

## 6. Clapeyron/source-response 公式的 normalization-safe 版本

若 safe transfer 使用未正规化微观行权 `W_tilde(g)`，行配分 normalizer为

```text
f(g)=log Z_row(g),
```

则 physical free energy是

```text
I(g)=f(g)-log lambda_tilde(g).
```

令 raw row score

```text
H_raw=partial_g log W_tilde,
mu=E_Q H_raw.
```

exact：

```text
I_g=f_g-mu.
```

所以两个 topological safe phases 的 coexistence差满足

```text
Theta_g
 = (f4_g-mu4) - (f8_g-mu8).
```

只有当 `W` 从一开始就是 normalized physical row transition weight，才可以把它简写为纯 `-mu4+mu8`。

因此 normalization-safe Clapeyron equation是

```text
dp_ch/dg = -Theta_g/Theta_p.
```

这保留了之前的 finite Perron/Feynman--Hellmann结构，但纠正“raw score difference就是 physical numerator”的过强写法。

### 6.1 source gauge quotient

若两个 source 相差

```text
constant + a * thermal score,
```

则它们描述的是同一 normalized probability family的 gauge/reparameterization。fixed-b normal response应不变，而 root tangent按 thermal坐标改变。

这与 #773 的 `source modulo span{1,K}` exact quotient一致，也是 #802 勘误的正确抽象形式。

## 7. 计算路线也需要改：N1105 当前不是默认下一步

独立 `dpfloor` 审计说明，统一 `p_c` 与 tight solver 后 numerical floor 足够低；但随后真正独立 oblique implementation 的 N325 rehearsal 表明状态空间成本指数爆炸，N1105 在现有 automaton上远超当前资源。

因此：

```text
precision gate : pass after fixes;
algorithmic cost gate : fail for generic N1105 tomography.
```

不要把“精度上 GO”误写成“计算上 GO”。

当前更高信息/成本比的是 #808 N377 specialized double-notch Phase 0；若它也撞墙，下一步应改 state representation / contraction algorithm，而不是硬扩 N。

## 8. 下一轮最值得做的四件事

### A. N377 只叫 `post-H4 residual gate`

Phase 0 先给 state/resource cost。若可算，输出 `p_perp4` 与精确 H8/H12 leakage budget。只有当 higher harmonics在允许 radial规模下不足以解释 residual，才升级 “angular scalar”。

### B. 构造 microscopic exchange tangent basis

用 Bernoulli chaos degree 1/2/3 的少数 translation/D4-symmetrized local sources：

- degree1 thermal control；
- degree2 pair source control（应验证 black/white quotient redundancy）；
- degree3 genuine nonthermal exchange-odd source。

计算 normalized safe-phase source-response矩阵，而不是先给 continuum field命名。这是 #61 的一个可执行 lattice-side tangent map。

### C. 对 topological defects 测 action，不测统一 prefactor

对 #800 slow splitting，优先拿 direct tagged-operator eigenvectors并估计

```text
-w^-1 log Delta gamma_w
```

与 independent mass interval比较。不要再先拟合 `Delta gamma/(w nu)` 为常数。

#760 类似：先识别 minimal wrap defect，再决定其 `sigma_wrap`，不要预设是 `kappa`。

### D. angular coefficient 分开后再做 continuum field map

若最终确认：

```text
P4 ~ ell^-4,
P0-pc ~ ell^-7,
```

再去比较 `x=21/4` spin4 与 `x=33/4` scalar/log-block候选；在此之前不需要继续扩大 operator label list。

## 9. 对本代理此前几份 note 的明确状态

- `equal-circumference-spin4-projector...`：**保留设计思想**；full tomography 的计算可行性已被 cost-wall改写。
- `sector-even-dressing-of-spin4-tower...`：**降级**为 H4 coefficient 内部 radial dressing；不再解释 post-H4 orthogonal residual。
- `dimension-21over4-resonance...`：总维数简并与 angular decomposition **保留**；“matching-odd field”措辞降级为 pair-exchange response，且 scalar `x=21/4` 是否存在需数据决定。
- `rare-topology-fugacity-and-mass-clock...`：mass-clock/Palm-score桥 **保留为证明目标**；统一 witness powers降级为 defect-specific action spectrum。
- `topological-clapeyron-response...`：Perron coexistence结构 **保留**，但 raw-source公式必须使用本节 normalization-safe版本。
- Krushkal/Euler bond-FK notes：finite embedded-bond algebra **保留**；不迁移为 square-site local identity。
- `topological-source-master-scaling...`：作为组织图保留；其 rare tail应从 integer witness series改成 `sum_T exp[-c_T s]` action spectrum，operator parity语言按 #61 降级。

这次审计的目标不是减少猜想，而是把可失败的位置重新放对：先测 lattice irrep / normalized source / network action，再谈 continuum operator identity。