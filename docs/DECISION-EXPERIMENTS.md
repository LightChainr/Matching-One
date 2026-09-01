# Mechanism decisions / 独立实验与精确传输判决

2026-09-01 · #334两项不同的独立实验、#154的165M固定新块和次级时钟线均已完成。
优先级只分配注意力，不锁Issue、不合并PR、不把一次参数化失败扩展为整条
研究路线失败。[最终科学交接](../notes/independent-decisions-final-20260831.md)
给出定义、数据独立性和精确来源。

## 最新完成：联合 Q 激活已传入原 U，且不限于最近邻接触

[实际结果](../notes/regular-pair-joint-u-result.md)，`open_pr #267`
`f8e30859f05e86ef35d257fc900f97e74f41e21c`，给出固定canonical Kreg的
`J2=∂logQ∂epsilon²U|1,0 = −.0055194314248394015`。
预先固定的NN部分为`−.001751074454402799`，非NN部分为
`−.0037683569704366022`；三个精确有理区间均排除零。

因此，**first-Q线性可加global closure**和**只有四个NN位移传入原U**
两个具体零模型均被排除。非NN在N25仍可相距很近，这不排除所有有限范围
contact/OPE机制，也不证明宏观分离场的尺度律。J2已移出待算。

定义/约简`7557da5271f85a69ea5426b61ce7e67b94ee8ff2`，预数据NN分解、
producer/scorer `99b58fc18666cfa6d35b96b52bb84c78dec43a55`。每几何只枚举
原点为空的`2^24`构型，三个joint源矩除`16N`，分母仍是旧全`2^25`总体。
相邻空点共享真实edge ID；完整q/E中心化、移根和斜率项全部保留。
总耗时4.63秒、一次编译及一次向量评分；没有旧源重算、找根、MC或科学测试。
这仍是同一N25精确总体，不能增加独立统计票。

执行分支随后在`branch_only`[410015f5](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md)
用另一套完整`2^25`遍历实现得到完全相同的total/NN/nonNN，并补齐相邻edge
kernel与distinct-site坐标不变性。它是同一精确总体上的实现交叉核查，不是
第二份统计证据；本总览自己的平移约简结果不因此重复计票。

另一个问题也已由新数据完成：execution `branch_only`
[a237968f](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md)
给出L64/r16空间均值`C64=6.85546875e−6`，99% MC区间
`[5.2033972758e−6,8.5075402242e−6]`。L32/L64各200k配置、200batch，
两尺寸独立；每配置32个pair是相关平均。它排除有限非接触空间零假设，
不是J2，也不把两个尺寸的ratio变成已识别指数。

PR #509（OPEN Draft live head `ac5f5fe7`；本线科学
commit `ec3941b03b2694e827db1cba34766a82e6146a5a`）的
[baa5d33b选择定理](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md)
进一步证明：非相邻两点至多共享一个外部占据组件时，首阶Q激活精确为零；
实际两组件外部给`a_xy=1/16`，且
`|C_xy| <= (43/16) Pr{两端vacant且至少两个共享组件}`。这删除了单组件
传播模型，却不是global U的数值界或距离指数。

PR #509的后续[eed2190c解析结果](https://github.com/LightChainr/Matching-One/blob/eed2190c04b67084ab5aef5827e00377853a0bca/notes/p337-critical-spatial-summability.md)
又证明exact `p_c`下raw canonical `g_xy`绝对可和，排除“raw红外尾自身发散”
这一具体解释。它没有控制`∂p E[g_xy]`、带`q/E`的三点pivotal influence、
near-critical一致性或原U的baseline thermal jets；因此不能把raw定理误写成
thermal/original-U尺度机制已经关闭。

## 宏观 window pilot 已完成，并按合同停止

旧“下一问题转向宏观尺度比”已经实际执行，状态统一改为
`completed_stopped_unresolved`。[冻结合同](../analysis/regular_pair_macro_joint_u_contract.json)、
[合同说明](../notes/regular-pair-macro-joint-u-contract.md)、
[producer说明](../notes/regular-pair-macro-joint-u-producer.md)、
[producer](../scripts/regular_pair_macro_joint_u_sampler.cpp)、
[scorer](../scripts/analyze_regular_pair_macro_joint_u.py)、
[结果](../results/regular-pair-macro-joint-u/REPORT.md)、
[机器结果](../results/regular-pair-macro-joint-u/latest.json)和
[run receipt](../results/regular-pair-macro-joint-u/run.json)已经齐备；raw `N100/N400.csv`
及四张固定window表保存在`results/regular-pair-macro-joint-u/raw/`。

同一canonical `Kreg`、Euclidean `1/4–2/5`窗口和完整moving-root U下，
N100的`T_N=+11.877846`，同时95%区间`[-2.2556162,+26.011308]`；
N400的`T_N=−542.50382`，区间`[-1307.6425,+222.63486]`。两区间均含零，
点估计异号，且N400投影需求`2,119,100`超过冻结上限`2,000,000`。
所以`D17/D21`没有评估；不top-up、不换window、不追加local/GPU/cloud生产。
这不是对`x=17/4`或`21/4`的否定，而是该固定窗口、该预算下无法进入场比判决。

N100与N400分别属于两个独立依赖组（每个50k paired configurations）；同一N内
axis/tilted共用occupation与anchor流，`total/s2/s>=3`也是同一源的相关加和，
不能拆成独立票。100个paired batches才是推断单位。pilot到此停线；其它探索
仍可并行，优先级不是许可或任务锁。

## 最新完成：contact carrier已传输但只占小份额；P0转向原始观测量

已关闭的Issue #536保留一项**issue-only、非main事实**的N25 midpoint精确分解。
仓库团队已把开放的[#275](https://github.com/LightChainr/Matching-One/issues/275)
调为P0，把[#537](https://github.com/LightChainr/Matching-One/issues/537)
保留为开放P1理论线。下面的有限分解仍是P1机制来源：
`observer-pivotal=−0.005935841948386444`、
`kernel=+0.0004141042799437205`、`root=+0.000002306243603322`，合计
`−0.005519431424839401`。它说明下一步不是首次发现通道，而是把已有有限分解
变成空间pivotal支持判决。

冻结的prospective readout已经
[完成](../results/p337-thermal-pivotal-preflight/REPORT.md)：复放既有L32、
L64数据块各自前32个counter，总计64个；不换seed、不新采样。5,242,880个
pair/site callback全部通过结构恒等式。有限子集内只出现18个kernel-changing
callback（L32 11、L64 7），没有kernel-preserving topological callback；
L64非零项属于external shell 2/3 shared-transition。这个结果完成语义与成本
preflight，不估计中心化总体或full J2，也不授权追加counter。

PR #509 `ec3941b0`给出精确压缩：取`u=M`为热坐标，则
`partial_u partial_epsilon Yhat=J_N/A_N`。其显式N25 six-block
clean-two-bridge临时合同的P4/root minors全部非零，故该端口级rank-one模型被
排除。修正后的`branch_only`
[`a3bc80c8`](https://github.com/LightChainr/Matching-One/blob/a3bc80c86585220f96dfeff022bd575b6c21d29f/notes/p537-finite-landing-transfer-definition.md)
同时明确：N9的`1001/536870912`与N16的`-chi/2`只是raw supporting controls。
它们把同一covariance derivative拆成kernel/readout两列，缺少独立归一化的
axis/tilted P4和完整`-beta H B`双线性Schur项，不能冒充canonical ordinary
landing falsifier。

因此Delivery A先执行
[complete transfer contract](../notes/p537-signed-landing-transport-contract.md)：
在global root/means/`R`/`beta_lambda`下形成每个`(pair,z)` fibre。最新
已关闭但未合并的PR #544（recorded head `e8e9c7cf`）所含精确定理
[`139c7e58`](https://github.com/LightChainr/Matching-One/blob/139c7e5850c5be1e8e17c3d58be806080b1e2b73/notes/p537-global-four-arm-emptiness.md)
给出`Delta q=1-b_z-w_z`及`b_z+w_z<=1`：两条全局occupied分支和两条全局
vacant separator不能同时独立。因此必须先把四臂运输到固定finite collar，再把
outer reconnection/rank attachment作为独立状态；global-component表不能作为
canonical矩阵。`near_block`虽只是放宽控制，但
[`ddbbad64`](https://github.com/LightChainr/Matching-One/blob/ddbbad6462abb9fcebdaae1f7162e287ed3bebed/results/p537-siteflip-aggregate/preferred-minor.json)
已经完成它的full-site-flip评分：`15 x 12` P4/Schur矩阵中，rank
`0->1`/`1->2`行与axial2 source absent/present列的det为
`-4.039642418513639e-18`，exact interval排除零，且degree-branch与extra-contact
都为0。已关闭但未合并的[#544](https://github.com/LightChainr/Matching-One/pull/544)
在科学commit [`e21f6220`](https://github.com/LightChainr/Matching-One/blob/e21f6220775d51450d73d4258c3edf2dbdf0a910/results/p537-finite-collar/REPORT.md)
完成coarse-marginal collar value test：
的`B_inf(z,1)`、outer `J_B=J_W=1`、all-corner、full root/Schur与axis/tilted
P4四格矩阵，在唯一声明的minor给出
`det=+2.6904188461441777e-14`，exact interval严格为正。这否定该coarse
signed-mass矩阵的rank one；producer未保留source+z joint 12-port identity，
而选中landing rows的local raw `a/qa/Ea` sums全部为零，所以不能把它命名为
直接source transmission或canonical physical closure。

四格的row margins为`(+4.3555592e-7,-1.7359553e-7)`，column margins为
`(+1.0122295e-7,+1.6073744e-7)`，均非零；这已经停止自动two-defect/extra-arm
升级。`60.7191%`、`2.24641%`、比值`27.0294`、cross-ratio
`chi=.0369967`和interaction residual `+9.9275800e-8`仅是endpoint-
availability下的global-centering/componentwise-beta Schur cell-mass统计。
branch-only [`c958b30b`](https://github.com/LightChainr/Matching-One/blob/c958b30b/results/p537-finite-collar-joint/REPORT.md)
已完成exact joint-incidence refinement。`1,390`个12-port sectors中`118`个
exact nonzero；`106/106`个pooled rectangles和`30/30`个完整two-geometry
rectangles均exact非零。第一个完整几何sector的
`det=-4.005552609094306e-18`；所有sectors精确重构coarse matrix及
`+2.6904188461441777e-14` determinant。same-sector determinants之和只占
`0.2101274%`，cross-sector bilinear terms占`99.7898726%`。rank two因此在
共同物理identity内存活，而coarse幅度不可归给单一landing sector。

branch-only
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-one-defect-diagonal-edge/REPORT.md)
已经完成edge gate。broad Bell/rank图有`135,253`个exact row classes和
`8,250,462`个physical fibres；kernel-changing+`g16`-changing子图保留
`6,846`个classes、`740,950`个fibres，总signed mass为
`-4.9488399165e-6`。冻结判决为`TWO_INDEPENDENT_DEFECT_GAIN_REJECTED`。
first birth `0->1`给`-5.8210906659e-6`，second birth `1->2`给
`+8.7225074943e-7`，形成`17.63%`反向抵消。12个rank-stage×source-orbit
cell全部精确排除零。

在上述fixed-anchor、NN thermal-`z`、source-present kernel-changing选择内，
contact mask给出更强的结构选择律：mask 0没有任何row class或fibre；mask
1/2/3承载全部kernel-changing边。radius-one剩余对象因此是contact-local、
first-birth-dominated的signed four-arm/OPE functional，不是两个分离defect。
`c958b30b`与`22f01e33`复用同一N25 exact population，前者是value decomposition，
后者是edge decomposition，不能算两票。

`df4a64f6`进一步把birth stage与single/double contact同时保留，得到exact N25
矩阵`[-2.88380e-6,-2.93729e-6; -5.32257e-6,+6.19482e-6]`和
`det=-3.3498535471e-11`。因此幸存对象不是scalar contact loading。
commit
[`f46c38c3`](https://github.com/LightChainr/Matching-One/blob/f46c38c3088c1a9f4df8ab0f256b88639f0b34a3/results/p537-contact-stage-n65/REPORT.md)
随后完成冻结的20M N65 held-out block：四个marginal同时保持`-- -+`，
`Delta_cs=-8.6882161e-14`，95% CI
`[-1.1459612e-13,-5.9168198e-14]`，判决
`CONTACT_FUSION_COMPLETION_TRANSMITS`。这些格子是pooled-root original-U
allocation；最终`J_N`还包含全局`A_N/M_t`归一，不能把raw格子随意改名为完整
`J/U`。

branch-only `95a695c7bd05eeb6ce1b24cb5299a13823b40fcb`保存P45与新MC的完整
协方差，并明确标注post-hoc shape检查：single-amplitude N25→N65模型
`chi2=10.3371/3`、`p=.01591`。它与`f46c38c3`使用相同seed和byte-identical
shards，是同一次production的依赖重分析，不是第二票。传输问题到此已完成；
all-z N25不能再作为它的前置门槛。

branch-only `bab37f21`把finite determinant写成noncommuting completion
commutator，并给出post-reveal的`N^-3/N^-29/8` fingerprint；这是机制候选，不是
prospective field identification。current branch head
[`f9ba1ff6`](https://github.com/LightChainr/Matching-One/blob/f9ba1ff6/results/p537-full-t-transport/REPORT.md)
用同一N65 sufficient statistics重建complete response：
`J65=-.00162250989+/-.00018553008`，
`J65/J25=.29396323+/-.03361398`。typed carrier占full `T`的份额由N25
`5.892%`降到N65 `2.551%`；pooled-root相对reference-p的位移只改变
`|J65|`约`.0155%`。carrier真实传输但不饱和原始响应，且这些N65结果仍只有
一个dependency block。

[#537最新gauge audit](https://github.com/LightChainr/Matching-One/issues/537#issuecomment-5490511026)
仍是issue-only：它在exact N9控制中证明complete `T_t`对共同thermal-coordinate
shift不变，而筛选后的contact-stage allocation一般不变；附加N65 sidecar尚未进
Git。这个结果加强#275的normalizer要求，但不覆盖f46冻结分账或冒充新的N65 score。

已关闭但未合并的[#545](https://github.com/LightChainr/Matching-One/pull/545)
记录head `e1f19e4c`；其L4 exact det、matching-root det
`-2.501463041122436e-6`和all-scale族是独立local-collar回归，不是第二个
canonical vote。PR #544关闭后的同名分支当前branch-only head为
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-aggregate-wedge-l8-mc/REPORT.md)：
L8 `G4=.000680204061±.000012697002`，冻结判决为
`UNRESOLVED_MODEL_BOUNDARY`；`L^-35/8`是post-hoc发现，不是prospective证明。
不再增加reduced-`G4`尺寸；existing-data `Xi_L=N^2 G4 sqrt(P_L)`也已计算：
L4/L5为`8.00214/7.99820`，L6为`7.72546+/-0.08469`，L8为
`8.03977+/-0.15006`。它是reduced normalization线索，不是指数。

下一P0只交付一份**#275 finite forward-prediction contract**，默认不新产数据。

1. 原始观测字典必须区分local charged insertion、torus propagating sector、
   homology-marked trace、scalar composite和normalized expectation；保留original
   `q/E`、同一physical normalizer和pooled moving-root `U`。
2. 至少两个候选机制分别给出带单位的预测向量，覆盖至少两个已有相关raw
   coordinates。保留每个模型允许的amplitude，报告covariance-weighted design
   rank与separation。
3. 若不可识别，证明rank degeneracy并只命名一个缺失的physical rotation或
   modulus relation；不以新descriptor掩盖退化。

一个unnormalized numerator character zero不能穿过twist-dependent normalization
后继续写成零。声称common vacuum/Ward differential operator的候选，必须证明
该operator在rank0/rank2 restricted traces上以同一规则作用；否则必须明确给出
defect/contact/normalization项。Q4/Jordan同样需给restricted-sector预测，而非只
借用同spin场名。

#537保持P1：非局域remainder和uniform near-critical transport未解。all-`z`
N25、N85 normalized tensor和N130 `F=C+R`都是条件式后续；只在#275冻结原始
observable/normalizer与相反预测后再启动。

同一N25 provisional端口合同的same/reversed bridge-order `(F,dF/dM)`
Wronskian已精确排除零，所以finite-collar列必须保留cyclic order。该signed
contrast尚未因而成为物理reflection-parity sector，也不是独立数据票。

| 预设结果 | 事前判决条件 | 机制更新 |
|---|---|---|
| distinct vectors | 两候选在保留allowed amplitudes后的weighted design rank增加，且分离超过冻结阈值 | 允许在旧数据上消除其中一类；明确标作retrospective model selection，不写成prospective验证。 |
| rank-degenerate | 合并design不增rank或加权分离为零 | 停止在同一观察上继续拟合；只请求事前命名的第二旋转/modulus relation。 |
| normalizer-dependent | numerator projection为零但各twist分别归一后非零 | 排除“character zero自动排除完整normalized observable”；在physical normalizer上重写预测。 |

[PR #530 open-pr support @0e1917f9](https://github.com/LightChainr/Matching-One/blob/0e1917f983982f411c9c72147cce3f90ce7330fa/experiments/p337-original-u-influence-20260831/README.md)
提供含共同root、均值和分母响应的中心化`phi`，作为候选共用的原U权重接口。
这些PR只提供有限代数/估计接口，尚未是main结论。若现有坐标不可识别，则交付
rank certificate与唯一缺失输入，不新增seed、window、额外臂或第四个descriptor
救场。文献
[适用边界](../notes/p337-thermal-pivotal-literature-bridge.md)已经说明，三处四臂
事件的朴素乘积不足以关闭bulk，六臂只接管已证明的fused collision annulus；
还必须给common-landing增益或signed cross-observable cancellation，并认证
pooled root处于near-critical窗口。#539只保留P2 exact-N25复现支持，并应优先
导入#536原附件而非重复生产。该tail是#537 P1理论对象；#275的observable/
normalizer/identifiability是当前attention target。优先级不锁其它任务。

另一个独立前沿也已前移：PR #533当前head
[`5aa929a6`](https://github.com/LightChainr/Matching-One/blob/5aa929a6ef4ff92d71f6f36097b763c05ad71ab8/notes/closed-source-axis-sector-closure-gate.md)
已经在`beta=L/m^2->0`下闭合声明的完整axis/tilted强源符号律：all-reversal
kernel、two-cloud压力、额外carrier、transverse/companion winding和sector odds
都只给`o(1)`余项，`Ustar/A_N<0`。此前`L^2/m^3` crossover与nonlocal
west-step候选是证明伪影；仅fixed-m热力学仍在该joint limit之外。
#534/#535现已由仓库团队关闭。该强源解析结果与临界canonical macro pilot分开，
且保持`open_pr`，不写成main已集成事实。

另一个不同源的有限任务也已结束：PR #509
[ef3b2c68](https://github.com/LightChainr/Matching-One/blob/ef3b2c68f824e29421747c805ea7a505aca41908/experiments/p337-homogeneous-n50-20260831/RESULT.md)
以完整状态合并覆盖每几何`2^50`配置，得到齐次N50 `Sstar`的
`U=1.0615603876876551`、`V_S=+0.0543457826695583`并严格排除有限零传递。
约49.85 CPU秒；合同结束，不自动追加N100、t/epsilon网格、尺度拟合或场认定。

## 已完成的前一步：正则单点直接传输被整类排除，双点 Q 激活保留

执行的`branch_only`
[2ba8863f结果](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)
已经构造并评分固定`Kreg=K2+K0`：直接`∂epsilon U|Q1=0`，完整
混合`∂logQ∂epsilon U|Q1=−.04503611397592696`，精确区间排除零。
它使用同一N25总体的新规定源矩，不是新的独立随机证据，也不是旧纯K2响应。

本总览`open_pr #267`的[单点定理](../notes/regular-one-site-q1-thermal-quotient.md)
和[统一counterterm判决](../notes/regular-pair-counterterm-gram.md)，科学提交
`21563da4b0cf721a2aa512901f6ffc966ffa8384`，进一步压缩机制：

- 所有entry-regular、homogeneous、one-original-binary-site equality张量
  在Q1仅改变共同Bernoulli参数；允许同时扰动occupied/vacant张量。
  完整移根/斜率U消去该参数，全部直接coupling导数为零。因此这个类别内
  “消pole且保留旧非零直接V”的搜索已结束。多格点、奇异confluent、
  占据connectivity重新加权及Q激活不在该排除中。
- 对`K2+c(Q)K0`、`c=1+alpha(Q−1)+...`，相同counterterm的first-Q
  双点Gram为`3/2+(alpha−1/2)^2/2 ≥ 3/2`。高阶Taylor项不影响它。
  两孔求和后除以`(1+v_x)(1+v_y)`仍正，排除此物理外部条件中的独立可加源。
  不同alpha的交叉配对没有正下界；条件正值也不保证global U为正或非零。

**当时命名、现已完成的比较：** 固定canonical `Kreg`、原N25方向对、site-average
归一化，取得真正联合闭合的`J2=∂logQ∂epsilon²U|1,0`。
保持Q1基线不变而first-Q有效log weight仅对epsilon线性可加的模型预报`J2=0`。
非零则排除此global closure；为零则放下“条件双点正值必然进入global U”的
主张，不增加可拟合counterterm挽救它。完整tensor的无条件符号尚未推导，
不能预填为正。`Cov(a_x,a_y)`不能代替联合tensor，旧V的尺寸比也不能套用。
这只规定一个可改变判断的输出，不展开新descriptor、source或certificate目录。

该前一步仅一次0.168秒有理化简及解析证明，无新occupation枚举、U评分、MC、
根搜索或科学测试套件。旧bounded occupation tangent的固定尺寸比作为
不同源并行保留，完整口径见[Next Targets](NEXT-TARGETS.md)。

## 已完成的前一步：纯 K2 局部传输非零，其双插入有 Q1 极点

执行[923f66b9](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md)
已交付固定局部四端口的site-average响应`V_av(25)=+.0018155512845251097`。
这是完整原U响应，local与full seam的构型级等同性也已由双向反例排除；
“首次局部插入接口/评分”移出待办，已有Q1/Q4 trace结果不重算。

此前[5864de49的精确结果](../notes/local-pair-two-insertion-obstruction.md)
给出`G(Q)=Tr(Kbar²)=Q(Q−3)(3Q²−9Q+8)/[8(Q−2)(Q−1)]`。
在17×17方格torus的两孔、四条互不连接占据路径中，它就是实际双插入闭合。
Q1留数为`1/2`；连通条件响应的留数是`1/[2(1+v_x)(1+v_y)]`，
公共归一化和单点二次counterterm都不能消掉此混合项。

**被排除的具体机制：** 原Kbar无需额外处理便可在所有物理外部条件中形成
正则有限强度Q1局部张量族。**仍成立：** 单次插入的有限非零U响应。
不据此断言对所有外部构型求和后的均匀partition发散。
新[固定cut分解](../notes/local-pair-crossing-sector-resolution.md)同时显示
四个colour块非零，裸thermal零overlap来自singlet/standard抵消，不能当RG零耦合。

当时提出的“正则单点补全保留旧linear响应”现已由上节整类定理解答，
不再列为下一任务。另一条对固定bounded占据tangent比较`W_N=N V_av(N)`，其[预先推导比例](../notes/local-pair-size-response-predictions.md)
`W_(4N)/W_N`在明确single-field/nonzero-loading假设下分别趋向2（x17/4）
或1（x21/4）。后者不唯一识别thermal-Q4；`Cov(t_x,t_y)`也不能替代joint tensor。
此前0.114秒符号代数与一个显式图，没有新MC、occupation枚举或测试套件。

## 已完成：稳定 Q1 通道的原 U 传输与完整颜色导数

[5c1f9d3b 的结果](../results/n25-stable-colour-q1/REPORT.md)已从既有seam计数
与已有Q1根算出`B1=∂epsilon U=−.001904836180602413`及
`B1_logQ=∂logQ∂epsilon U=+.005036496028411871`。两有理区间均排除零。
第二项包含beta的显式Q导数、测度、移根与热斜率；它不是总U的Q切向。
一次局部有理计算6.759秒，没有新枚举、随机生产、Q4重算或测试。

[N25 packing证明](../notes/n25-stable-colour-completion.md)给出`c|u|≤2`，
因此已有S4计数足以定义全stable通道：
`beta(Q)=I12(Q−3)/(2Q)+I21(Q−3)/2`。B1中的两部分分别为
`−.001945570733316785`和`+.00004073455271437206`；单独两essential簇的
限制会漏掉第二项。完整颜色导数中显式beta项为`+.0028979888236179917`，
测度/移根/斜率合计`+.002138507204793879`，均非拟合份额。

[更大tori的精确反例](../results/colour-specialization-gap/REPORT.md)仍排除
“任意closure只乘维数因子便可从Q4延拓”的一般说法；N25的几何证明解除
本例的这个障碍。局部pair的实际构造、正响应及其与full stable trace的
不同支撑现已完成，下一步见上节；不再把Q1或首次local计算列为待办。

## 已完成：Q4 有限归一化传输接口及其非零评分

[041257b8 的解析推导](../notes/four-leg-trace-denominator-interface.md)
在固定m=2、同一空间恒等颜色缝给出`Z22=(T+3D2−4C3)/6`，精确投影
完整S4的`[22]`通道。其逐构型支撑仅在rank1，所以原q/E未归一化分子
严格为零；但原分母中的份额为`z22=2Z22/(T+R)`，通过完整四项root/slope
公式传入`J22=d_x U`。x是固定分量归因，不是拟合的新颜色数或局部正源。

轴L≥5时`Z22(a)=L(L−3)4^(2L)a^(2L)+...`，由两个不相邻竖直列精确
计数得到；它证明trace分量非零，**不证明matching根上的J22非零**。
执行[54352b2e](https://github.com/LightChainr/Matching-One/blob/54352b2eefa651ca482ca84837053c792e82c71e/results/p337-s4-trace-transmission/score/score.json)
现已给出`J22=+.000005440121494634842`，严格正界排除该finite null。
新的rank1 seam信息已取得；Gaussian90°旋转、rank分母因子及热变量的
对应证明使它与我们的J22相同，不需要重复score。执行来源为branch_only。
定义及完成指针见[机器接口](../analysis/four_leg_trace_interface.json)。
Q4颜色表示不等于唯一连续四腿场；本例Q1延拓与评分已另见上节。

## 已完成：N25 的统一有限 m 符号区间

[85d5e44b 的统一余项](https://github.com/LightChainr/Matching-One/blob/85d5e44ba8aed471470373f972c670dc7c82bdcf/notes/closed-source-uniform-projection-tail.md)
和[有理评分](https://github.com/LightChainr/Matching-One/blob/85d5e44ba8aed471470373f972c670dc7c82bdcf/results/p337-uniform-projection-tail/score.json)
已对所有实m≥64证明各律自身原pooled根上`Ustar<0<Udrop`。
归一化后有向外取整界`−.618102m^−11<Ustar/A<−.454124m^−11`及
`1.376734m^(−42/5)<Udrop/A<1.844309m^(−42/5)`。
计算复用同一N25整数总体，耗时0.444821秒，零枚举、零MC、零耦合点求值；
不是新独立证据。最小阈值和可采样幅度仍未知，但“有限m窗口待建立”已过时。

执行的轴fixed-m winding和轴/斜Poisson同步极限亦已完成，后者使实际pooled U
超多项式小；它们适用不同极限。fixed-m斜几何仍需几何order25的twist penalty
及rank2/rank0 odds mismatch两个界，见[统一注意力表](NEXT-TARGETS.md)。
这些执行来源均为branch_only；没有因此合入代码或启动新的生产。

## 已完成：弱 Q 路径异号，正则端点的 Q 激活被排除

[e87d5de2 的精确结果](../results/weak-q-path-comparison/REPORT.md)在同一N25
基线上得到闭合源Q路径响应+0.063082681707085、保留rank投影的普通site-RC
Q路径响应−0.269828026713487。其差是已指定的局部B控制+0.332910708420572。
旧有理包围按`V_Sstar/2=V_(CB-r/2)+V_Bvac`直接变换，耗时0.0202秒；
没有重算根、旧scorer、枚举或随机数据。这排除两条命名路径的有限U切向等同，
不构成跨尺寸普适性判决。

[整族选择推导](../notes/weak-q-paths-and-regular-selection.md)给出有理恒等式
`ell P_[2](Q)=0`，所以此正则未标记一插入端点的所有正则Q导数都为零。
排除的是这个机制定义下的四腿sqrt(N)激活。上面的Q4 trace/归一化接口
与J22评分现已完成，指定stable Q1延拓及响应也已完成。canonical正则局部
补全及其混合Q激活也已完成，不能把它写成尚无任何Q1附近传输结果。
不能把这些对象等同于正则端点或自动算成一个已经识别的替代模型。

## 已完成：两套固定微观律给出相反的 U 强耦合尾

原闭合源Sstar与已定义的去投影对照Sdrop=Sstar+r，在固定t、N增大时
pressure-density差不超过2t/N；这不保证有限拓扑观察者相同。
[本次精确首项提取 `fbbaa2aa`](../results/projection-drop-tail/REPORT.md)给出N25
`U_drop/A=+(625/384)lambda^(42/5)+...`，原模型为
`U_star/A=−(625/1152)lambda^11+...`。去掉m^(-r)会打破最低条带的
归一化占据对称抵消，改变符号和幂次。这里只比较两个已指定模型，没有拟合源系数。

对轴L×L与ell1>=L+2的同面积伴随几何，已推导
`U_drop/A~+(L−2)lambda^(2L−2+2/L)/Delta`；执行的原模型为
`U_star/A~−(L²−6L+6)lambda^(2L+1)/Delta`。
N100/N225为理论预测，N25新首项只使用旧整数直方图，耗时0.0366秒。
渐近符号已经分开；后续85d5e44b已给出N25的m≥64统一区间。
统一尺寸余项和可采样幅度不由这个固定N窗口确定。

## 已完成：循环/rank分离与一孔传输

[本总览的精确分离 `5483aa82`](../results/decimation-cycle-rank/REPORT.md)得到
`V_q=.053247535115`、`V_Sstar=.126165363414`、
`2V_beta_null=.072917828300>0`。显式单位系数q不能解释完整响应；
不扩大为任意重拟合rank源或固定K/rank内中心化源的排除。

[执行的一孔结果 `f5c4a74a`](https://github.com/LightChainr/Matching-One/blob/f5c4a74a20bad8589c39e1034cfb209462110dbe/results/p337-endpoint-defect/score/REPORT.md)
已完成：`R=U U_st−U_s U_t=27.766563581230>0`，
`Xi=−U_st=−10.755718407564`。源无关整体gain预报的U_st为
`.531368026778`，实际为`10.755718407564`。主gain模型与次级thermal-only
mixed null都失败，硬端点闭合仍成立。两判决来自同一exact packet，不能合票。

[固定算子分离 `e1b96895`](../results/defect-reweight/REPORT.md)也已完成：
baseline重新加权贡献 **+4.550327123237**，加权跳变 **−15.306045530801**，
合计 **−10.755718407564**。省去 `Cov(w,O_intact)` 的jump-only模型被精确排除。
重新加权包含rank-changing和rank-preserving构型，不能误称仅后者份额；
两项异号、固定系数，不拟合份额或补第四个自由源。只补了交替邻点1/8构型，
旧完整统计及根包围复用，共2.036秒；不重复“首次Xi”或“首次分离”。

## 同一正源的转折、负尾与尺寸律已完成

执行的 `branch_only` [a70eeff0四耦合结果](https://github.com/LightChainr/Matching-One/blob/a70eeff09f51ce2fa0fea5ae637e9191efbf2e1f/results/p337-closed-source-finite-coupling/score/REPORT.md)
在事先固定m=2,4,8,16全部给出U_t<0；结合U_t(0)>0，排除“始终单调增强”。
至少一局部峰位于(0,log2)，不声称峰唯一或这些点确定衰减指数。
固定N强耦合的双态与rank1耗尽提供U→0的理论解释；完整证明提交晚于score，
结构预测已在b70dc4bd冻结时提出。随后e3c8d3a的精确负尾和762dbaf4
尺寸律进一步保证至少一次后续过零、负谷和从负侧回零；不定位或声称唯一。
跨尺寸首项及N25的m≥64窗口已给出，下一fixed-m斜几何问题已缩成
twist penalty与sector odds，详见[Next Targets](NEXT-TARGETS.md)。
不自动增加峰附近采样点。

## 新独立F4块：80M已结束，未分辨

[固定结果 `25ca3635`](https://github.com/LightChainr/Matching-One/blob/25ca3635ea64655923c32adee4b62d683579cdcd/results/p337-f4-transmission-20260831/scored/REPORT.md)
在N65/85/130/170各20M、100批。四个共同95%区间均含零，也均未全部
落入±.5实用带：零模型`NOT_EXCLUDED`，决定为
`INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP`。这不是零耦合或成功传输；
N25精确正结果不覆盖此判决。各N独立seed，同N ordinary/forced-face与
两个方向共享排列并配对，不能按四条流独立计票。

## 新增已完成判决：抽稀强制项确实进入原global U

[固定N25完整枚举](../results/decimation-plaquette-u/score/REPORT.md)给出
bulk F4源的 `V_F4=+.1944146864609`；移根和热斜率修正全部包含。
`V_F4/(25^(13/8)/2)` 的有理数包围严格为正。因此“F4只是热重参数化、
原U端点传输可只保留簇数C”的预测在这对有限商上失败；N50端点漏项
为 `2^(13/8)*V_F4=+.5996568681566`。这是定理传输的端点值，
不是另一份N50随机生产。合同先于枚举固定，完整整数系数保留。

[随后完成的源字典](../notes/decimation-closed-source-and-global-u.md)进一步给出
`S_hat=C+F4+T_NN−4K+2N` 的精确端点不变性。这里没有拟合第五个描述量，
也没有挽救#154/#334的失败模板：端点变换本身强制了全部系数。
新结果范围是degree4 square-site、合法checkerboard端点和固定N25方向对；
不同Smith组不能当作同一生产lineage的连续场证据。

源note导出的单A-vacancy原U计算已由上述f5c4a74a完成。
其非零Xi排除了共同thermal-only延伸，非零R进一步排除了源无关gain。
不重复接口、首次评分或新采样；固定算子分离已见本页顶部。

## 1. 已完成：#154原U弱传输触发停止优先投入规则

[正式结果 `f4999e29`](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md)
已交付全部九片、主六坐标和协方差。B进入主导、C完成主导均被排除；
四通道区间完全位于W的±.30带，W保持not_excluded而非物理机制确认。
两个net同时区间为N85 [−.071640,.158581]、N340 [−.157394,.278745]，
均位于冻结±.50带。**这个lag1源退出当前主要H4解释的默认优先投入。**
不追加样本、不换lag/source、不拟合第四模板；下面保留原定义及判决规则。

正式源为对每个最终K在K−1注入、早期rank内中心化的bulk CB+CW源，lag=1。
固定几何N85/N340，新样本5M/160M，各200批，共九分片。
[合同及生产代码冻结于 `0820b8d2`](https://github.com/LightChainr/Matching-One/blob/0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f/experiments/p154-prospective-transmission-20260831/CONTRACT.json)，
后续 `14b2c98e` 提交的authorization已给实际freeze SHA；合同中保留的“待授权”描述
不是新的团队许可要求。生产与一次正式主评分已完成，不新开重复run。

新块内重估条件均值、共同matching根和所有配对delete-one。
p微分只作用于Binomial权重，源定义与K内条件均值函数保持固定。
根括号[.55,.65]，原归一化指数13/8、两方向精确角差不变。
旧发现数据仅用于事前预测和预算，不进入新评分。

### 六个主坐标与不可混用的通道分解

记F1/F2为两激活CDF，q=F1+F2−1，E=1−F1+F2，
D=bar(q')，A=N^(13/8)/2。合同的**读出**分解为

```text
U_entry = −A P4(F1')/D
U_completion = +A P4(F2')/D
U = U_entry + U_completion
v = v_entry + v_completion
```

两通道均使用完整源的同一个pdot和Ddot。这与分别开启早rank0/rank1
**源开关**后归因的V0/V1不同；不替换旧数值。每N保留entry、completion、
net共六坐标及完整共同协方差，不把可加的net当作独立证据。

### 固定数值模板与失败规则

- W：entry与completion均在[-.30,.30]。
- B：entry≥.60，completion在[-.30,.30]。
- C：completion≥.60，entry在[-.30,.30]。

每项要求两个N同时成立。六个共同区间用
`z=Phi^-1(1−.05/(2*6))=2.6382572735`及配对delete-one SE。
某项必需区间与允许带无交集即判模板被反驳；未排除不等于唯一胜者。
三个模板互斥、非穷尽，是有限数值预报，不是三套完整物理理论。

净U的独立主目标：两个N的共同区间均完全进入±.50，就停止将该lag1源
在这个分辨率下作为当前主要H4解释；这不证明精确零或完整clock闭合。
两N的净U下界都>.50则实质正传输目标存活。负、混合、跨N不一致或未分辨
照实报告。固定预算结束，不加样、不换源、不拟合第四个救场矩形。
这些是声明的渐近共同区间，不是有限样本exact certificate。

## 2. 同一新块的次级解释：M10/M11均未排除，但没有识别时钟

[次级结果 `612df8ec`](https://github.com/LightChainr/Matching-One/blob/612df8ec1cbe3be3938ee2e1f6183a1aefc6510b/notes/p154-clock-line-secondary-result.md)
四个同时区间均含零，两条线均not_excluded。两条线都允许零幅度，因此不能用
弱响应选择一个时钟，不能覆盖主停止规则，也不把二幅度混合拿来填满读出平面。
这是读取结果前、生产期间/之后的次级注册，不是第二份生产前主合同。

[完整传输图 `c2828e34`](https://github.com/LightChainr/Matching-One/blob/c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847/notes/p154-birth-clock-transmission-map.md)
定义 `alpha_jg=N J_jg/F'_jg=a00+eta*a10+tau*a01+eta*tau*a11`。
整个a00(p)共同函数在U中抵消；其余contrast值及p导数共同决定U。
N340的a10名义值gain246.4、a01约.2694；这是未标记曲线灵敏度，
不是已测源幅度，也不自动给出大效应。

[次级冻结 `83f3eba8`](https://github.com/LightChainr/Matching-One/blob/83f3eba88d7f1290704f82610c28669dc5e12f3c/notes/p154-clock-line-secondary-freeze.md)
已在读取新源结果之前固定两项纯、局部平坦first-jet限制：
M10只保留a10值，M11只保留a11值，允许任意共同a00及其导数。
两者不等同于我们已在旧数据排除的纯cos4相对位移。

现已只用完整官方六坐标及同批未标记q/E计算

```text
d_(N,m) = (C_entry*v_completion−C_completion*v_entry)
          / hypot(C_entry,C_completion)
```

基线gain在相同delete-one内重估；不除以随机entry，不拟合源幅度。
四残差顺序N85.M10、N85.M11、N340.M10、N340.M11。
各边际采用 `Phi^-1(1−.05/(2*4))`，一N排除零即反驳该两尺度限制；
否则只记not excluded。保存与原六坐标的共同协方差，不新增独立票数。
两条都失败就停止此scorer，不补拟合混合线；原主模板失败不能由次级结果改写。

## 3. 已完成：#334 source-normal完整闭合比较

[结果 `1164ba91/d0a9daf1`](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md)
使用N325/N425各500k新prefix，20批/N，与旧20批分离。
原冻结M0为完整条件label均值两-score闭合；M1预报四own-C均值
T_forecast=3.6565e−8。二阶源normal phi保持class质量和两个原score均值，
以合法非负有限政策采样，max|phi|权重留在期望中。

新T=(3.0852005663±.3918738407)e−8，3SE区间
[1.9095790443,4.2608220883]e−8，下端>delta1e−8。
**M0按原规则淘汰；正响应迁移预报未被排除。**
第一Jacobian近似J=BG并未被该结果整体否定，也没有识别global异常的起因。
固定块已结束、原始结果已提交，不再登记为producer待封存或待采样。

## 4. 已完成：#334四-contact残差固定幅度比较

[另一独立实验 `14b2c98e`](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md)
使用原一阶H源 `q±=(1±H/8)/d`，每N300k新prefix、60批。
主量为完整人口加权的固定四-contact预测残差clock-loading投影R，
不是normal C响应T。旧系数/均值/R_old点预报在4b3c21b7冻结。

R_new/R_old在N325为.4988857 [.4360616,.5617098]，
N425为.5169035 [.4506760,.5831311]。两者均避开C0[-.25,.25]和
C1[.75,1.25]，**两预报均淘汰**。每N为t59的97.5%区间，双N共同95%。
推断条件于旧固定训练值，不证明未知总体参数发生减半；
“约一半”不转成新的确认模型。

## 收口与下一问题

#334两实验共1.6M新prefix、四个分离N-domain；内部所有配对读数保持依赖。
二者不能合为一个残余效应或统一场计数。它们把下一科学问题推向具名
global observer的传输，不自动批准更多contact描述量。

#154旧scalar预报.01/.018对零难分辨的预算只约束那个候选，不能取消
新的165M较大效应区间实验。
现在两套评分均已完成；原规则确实给出了有限尺度排除。新的checkerboard
闭合源给出一个具名微观候选及已完成的原U端点计算；下一主实验仍需它的
内部传输定量预报，不能从本次弱结果或contact约一半反向拟合。
未形成此类预报时不自动占用生产队列；探索仍开放。

执行与结果交接只留仓库。十机均可按更新Huawei Skill直接使用，仍不覆盖
他人输出或操作不明进程。本次新增本地有限全枚举和精确源推导，未开关机、
追加随机采样或运行科学测试套件。
[此前决策稿 `f670d26e`](https://github.com/LightChainr/Matching-One/blob/f670d26e8056116ec1787e0a2f29009b5db24a7d/docs/DECISION-EXPERIMENTS.md)
保留未冻结候选、旧预算与历史映射，不再充当当前合同。PR #267保持Draft。
