# P334 / P429 / P487：成果与生产档案上下文恢复

本轮参照 Draft 工作树 `7cdf1d575f077cebec5cdfb85a25695cc03d1728`。
任务为只读恢复已有研究、识别重叠，再整理交接；**下文候选方向不启动任何研究任务**。
没有运行重放、可靠性计算、测试或 Monte Carlo，没有服务器/GitHub 写操作。
只读 CSV 的行数、字段与 counter 覆盖属于档案清点，不是新增总体分析。

## 一句话现状

两 N425 选定 checkpoint 已从 pair-only 时钟推进到**完整物理 birth clock**，
又得到 canonical K2 profile 与可消去的条件 suffix covariance。
后续**另外 12 个预先冻结选择的真实前缀也已全部求解**，其 47 个并行通道的 winning-time
和原始 direct/collective 分解也已完成。固定 PR484 的一般 two-port theorem 已把表示
推广到所有声明范围内的 embedded rank-one prefixes，并证明多通道多项式因子化。
PR #491 的 cut 网络和 PR #492 的白 matching blocker 是相关但不同的表示。
现在未完成的是**更广总体上的几何机制分布、配对方向响应、求解成本与噪声收益分布**，
不是首次证明一般 two-port、首次计算 12 个前缀或首次写出 channel-winning 公式。

## 固定 PR484 与已捕获关联分支：新增完成结果

这次补充只读固定 PR484 head `705819e95d1146fdedb06e9c7628344f108b80af`，
不追更高 head。相关 12-prefix/race 来源以下分别固定为已存在的本地 source refs。
它们不是固定 705819e tree 中的文件，不得把跨分支的证据链接误写成同一次集成。
7401 thermal 文件也不在该平行 tree；两点 git diff 中的缺席不表示研究成果被删除或推翻。

| 来源 | 完整 commit | 路径 | 状态与新增范围 |
|---|---|---|---|
| PR484 固定 head | `705819e95d1146fdedb06e9c7628344f108b80af` | `notes/p334-general-two-port-birth-theorem.md` | 一般 embedded rank-one → 多个并行 two-port 网络；factorization 与 winning-channel 导数公式已给出证明。 |
| 12-prefix selection freeze | `b9cbe13e05da778f385d633c2d6a716e080b06c9` | `analysis/p334_twelve_prefix_selection_manifest.json` | 在新映射/求解前固定选择，不是结果文件。 |
| `analysis/p334-twelve-prefix-clocks-20260831` | `bd95f2a048d5780568b689bd42e0a684daf74315` | `notes/p334-twelve-real-prefix-clocks.md`；`results/p334-twelve-prefix-clocks/full_clocks.json`；`scientific_summary.json`；`selected_prefixes.json`；`maps/*.json` | 12/12 `solved_full_physical`；完整条件时钟、成本与 crossing 已保存。 |
| `analysis/p334-component-race-20260831` | `e3d978216220b5c55ff8da3062473e713c1246c7` | `notes/p334-exact-component-birth-races.md`；`results/p334-component-birth-race/summary.json` 与每 counter JSON | 已完成 12 前缀、47 factors 的实际 winning-time 与 direct/collective 分配。科学实现 commit 为 `d519b71696616e7d19581580d440273c3553cef5`；e3d9782 只把文档 factor 总数纠正为 47。 |

### 一般表示和并行因子化已不是缺失理论

固定 PR484 的证明选任一 occupied essential connected component K，以
`q(v)=det(P ell,v)` 表示横向增益。对真正 embedded 图：完整 G\K 的 cycle q-gain 为零；
每个 complement component 最多两个 attachment addresses，差值大小为 `abs(det P)`。
rank-two birth 是任一 two-address component 的两 terminal 连通；不同 components 的
随机顶点集互不相交，共享的只有固定 terminals。

故完整安全多项式为 `F=(1+z)^r product_i F_i`。固定总插点数 k 的 component survival
不能直接相乘，必须取乘积多项式系数进行卷积。对第 i 个含 n_i 随机点的通道，
`B_i=n_i F_i-(1+z)F_i'`；乘其余安全因子得到 winning-channel/birth-time 联合法则。
n_i 是随机点数，不是安全多项式次数。该公式不需要 sitewise forcing 扫描。

仍未得到的是所有 prefixes 都有小 treewidth/低求解成本；证明也不直接适用于未连接的
交叉 matching-diagonal 抽象边。它不自动证明 PR492 的 white cycle-packing 充分条件普遍存在。

### 12-prefix 初步适用结果：完成，不是全 population 的估计

freeze 在同一 N425 second、k0=252、age10、ell=(12,-19) 的归档 rank-one rows 中，
排除旧两例，再取最早 12 个 counter。共有 147 eligible rows；H2、b2、clock、treewidth、
runtime 和 port shape 都不是选择变量。没有替补、额外 MC 或改变筛选。

结果已经是：

- 12/12 完整求解；mean T 为 5.772687743062006…24.89885088052643；
  P(T>40) 为约 0.00006809…0.11284544。
- 四个 prefix 单通道，八个含 2–14 并行 two-port components，共 47 factors；
  原 single-factor initial pass 的八个未完成项已被多因子方法完成，不能继续列成未解。
- core 随机点 87–151，实际 treewidth 上界 4–7；记录完成阶段总 wall time 约 0.56 秒，
  单 prefix 最大约 0.28 秒，均未达到冻结的 20 秒/200,000-state 限额。
  这些是所选条件小集合的成本观察，不是全 20k/100k population 吞吐保证。
- 66 对固定 prefix 中 11 对完整 survival crossing；不是总体 crossing 比例。
  例 counter43042500083 相对43042501006，在 k1…10 survival 更低，k11…153 更高，
  mean T 反而长 1.1725270404465133 步。这个展示对是在全冻结结果出来后描述性选择。
- 因此“原两例全程同序”不是同 geometry/age/line 下的一般排序规则；
  H2 或通道数量都不能单独排序完整时钟。

### Winning-channel 数值分配也已经完成

component-race source 已保存各通道 F_i、B_i、全系统 winning polynomial、实际 site lists、
精确有理 win probability、条件时间及尾部份额。第一次导出此前未保存的 factor coefficients
确实额外调用过小 DP；不能把那一步描述成 zero-DP 归档读取。后续 direct/collective 代数
则直接复用已保存系数，新的使用无需再求解同一网络。

counter43042500083 的 71/5/11-site 通道分别赢得 45.9074%/25.7058%/28.3869%；
counter43042501006 的 107/2-site 通道分别赢得 79.0869%/20.9131%。
83 的整体 collective share 为 32.2122%，1006 为 47.7172%；1006 的 collective 完成全部
来自 107-site 通道。它们在 T>40 births 中的 collective share 分别为 60.5392%/71.6711%。
这些是已求出的**条件竞争通道**份额，不是跨前缀均值差的因果比例，也不是总体 H4 分解。

仍未完成的范围因此须更精确：12 个全为同一 second orientation、同 age/line 条件集合，
没有给两方向 population covariance；也没有 sitewise pivotal support/pi_v 普查，
没有全体 20k/100k 的 cost/noise-weight map。下一步不能重新命名已完成的 12-prefix
或 winning-channel 任务，应先完成总仓库事实整理再决定是否扩展这些范围。

## 已完成链：pair → triple → quartic → full → thermal

以下均给出实际读取的 source ref/完整 commit；引用分支结论不等于本轮合入 main。
本轮未做任何合并。完整/thermal 两个 head 不在当时本地 `origin/main` 或 Draft branch 的祖先中。

| 读取分支 | 完整 commit | 主报告 / 结果路径 | 已完成范围 |
|---|---|---|---|
| `analysis/p334-pair-only-clock-20260831` | `ad6c595a70c66ea4421c816b4c65b1cfe3d9c803` | `notes/p334-pair-only-survival-clock.md`；`results/p334-pair-only-clock/pair_only_survival.json` | 两实际 pair 图的完整 independent-set 时钟；不是完整 rank-two event。 |
| `analysis/p334-pair-triple-clock-20260831` | `d5d2cc89e77ebb2ec6252df75dc858e9c240e6ce` | `notes/p334-pair-triple-survival-clock.md`；`results/p334-pair-triple-clock/pair_triple_survival.json` | 全部 583/509 triples，完整 ≤3 截断时钟。 |
| `analysis/p334-quartic-clock-20260831` | `1614a17e10997656fdf2d5520846fff2a228a5cd` | `notes/p334-quartic-birth-clock.md`；`results/p334-quartic-clock/quartic_survival.json` | 全部 1,178/2,866 minimal quartics，完整 ≤4 截断时钟。 |
| `analysis/p334-contracted-full-clock-20260831` | `6358ba49ef390c10a3f501b589ba7ba1d4e05b09` | `notes/p334-full-physical-birth-clock.md`；`results/p334-contracted-full-clock/full_physical_birth_clock.json`；`whole_event_networks.json` | 两 checkpoint 全部 173 个未来变量的完整物理 event；无有限 trigger 截断。 |
| `analysis/p334-conditional-thermal-averaging-20260831` | `7401c931117b693250139d7523406ba181decb24` | `notes/p334-conditional-second-birth-thermal-response.md`；`results/p334-conditional-thermal-averaging/score.json` | 消费完整 T 分布得到 canonical K2 曲线、18 个 p 加积分量、完整 19×19 条件 covariance。 |

共同对象：N425 second，period `[[425,268],[0,1]]`，seed `20260831430425`，
A counter `43042514269`、B counter `43042505280`，k0=252，age=10，ell=(12,-19)，
H2=0、d=b1=173、b2=14770。它们始终是**同两例选择性条件状态**，不是五批独立证据。

| 时钟 | A 平均等待步数 | B 平均等待步数 | B−A |
|---|---:|---:|---:|
| pair-only | 21.60624182 | 31.75099733 | 10.14475552 |
| ≤3 | 18.00285614 | 23.11369450 | 5.11083837 |
| ≤4 | 17.75453001 | 21.17202754 | 3.41749753 |
| 完整物理 | 17.73237780 | 20.77877866 | 3.046400854498077 |

full source 已给出：B 在非平凡共同区间 k=3…154 survival 更高、hazard 更低；
两例最大物理安全集合均为 154，k=155 hazard 均为 1。
实际 terminal cores 含 122/146 个随机顶点，treewidth 上界 4/6；
既有 DP 报告运行时间约 0.01/0.21 秒。这些只是两个网络的已记录成本，不能外推总体吞吐。
因此旧 pair/triple/quartic 笔记的“下一步求更高阶/完整时钟”已经被后续 source 消费。

thermal source 在 p_ref=0.59274605079 给出 f2(A/B)=0.1016281483/0.08617025757，
实际 binomial-tail suffix variance=0.01185405546/0.01171741881。
不能使用 Bernoulli `f2(1-f2)` 代替这个实值生产读出的方差。
它没有估计总体平均 suffix variance、总方差下降比例或总体可解前缀比例。
其条件-law 输入未显式提供 K1，故未产完整 A/E；**这不是整个生产档案缺少 K1**：
下述两类 geometry CSV 都保留 `k1`，未来若使用应明确按 replica 关联，不能从 age 猜测。

## PR #491 与 #492：互补，但不是同一项重复计数

### #491：黑 NN occupied cut 的一般网络表示

- open PR #491，head `research/p487-cut-network-mechanism-20260831`，
  commit `ab90201e88409310632812727e0138c56b455644`。
- 全文读取 `notes/p487-cut-network-theorem.md`、
  `scripts/p487_rank_one_cut_network.py` 与 `results/p487-cut-network/inputs.json`。
- 定理范围：honestly embedded torus graph，rank-one prefix，沿 occupied simple essential
  cycle 切开并收缩旧 occupied components；所有未来 vacant subsets 的 rank-two event
  等价于两个边界 terminal 的 **vertex** connectivity。
- code 在 HNF square-NN 上排除 coincident darts/loops；保持最初 cut，插点/收缩/吸收形成
  update-closed network。随机变量是原 vacant sites，保留固定基数无放回采集，不是独立 edge。
- pair 图从几何先给出的 L/R ports、直接 vacant 邻接和 neutral-component bicliques 得到；
  非 Ferrers 并未被撤回。两例 W2 的差 540 已解释为 within-biclique 472 加 shared-site 68。
- `essential_cycle()` 按 occupied labels/BFS 顺序选 cut，独立于 future pair outcome；
  但不能据此称 raw cut/component 编号已经是 D4/translation-covariant 物理角向 mark。
  定理的 event 图 cut-independent，不代表每个更细 component presentation 都 cut-independent。

### #492：白 matching essential cycles 的 blocker 描述

- open PR #492，head `research/p429-dual-cycle-blocker-20260831`，
  commit `0e52dbaeed53dfffa94592e53e38129c179c5078`。
- 全文读取 `notes/p429-dual-cycle-blocker.md`；证书为
  `results/p429-dual-cycle-blocker/certificate.json`。
- 在 digital-Alexander rank identity 的定义域内，minimal black completion hypergraph
  是所有 white-matching essential cycle 顶点集合族的 blocker。
- 若两条 essential cycles 仅可在 singleton triggers 上交叠，则 minimal pair 图 bipartite。
  这是充分条件；普遍存在这种 cycle pair、必要性、annular packing/minmax 尚未证明。
- 已有两个 N425 cycle 证书：长度 A=20/25、B=19/43，仍是相同两例、相同 108 边/W2。
- #491 的 embedded 黑图 theorem 不自动解决 white matching 对角连接的全部拓扑语义；
  不能把普通平面 edge-Menger theorem 无说明移植成这里的 vertex/cycle-packing 结论。

两 PR 均已给出新的几何解释。再跑两例的 pair/c3、安全更新或 cycle certificate 不会形成
新的 population 结果；也不应把两份报告重复算作方向响应证据。

## 档案实际可用范围

### 22 个图：选择性结构集合，不是总体

source `experiment/p334-cooperative-closure-pilot-20260830`，
commit `1b5a9dea07e1c62f69798fddbf4899ff986c0b72`。

- `results/local-20260831/P334-cooperative-closure/trigger_graph_raw/*.json` 共 22 图。
- 选择规则见 `scripts/p334_trigger_graph_structure.py::selected_rows`：两个 N425
  same-scalar-state witnesses，加每个 size/orientation 最早五个含 minimal trigger pair 的 counter。
- witness 原选择又是 age/line 匹配组中 degree-square range 最大的碰撞；
  见 `scripts/p334_checkpoint_scalar_collision.py`。它服务于精确反例，不服务于 population 比例。
- 图 JSON 保存 N/h12/k0/seed/counter/ell、safe site labels 与完整 pair edges，
  **没有全体 occupied masks**。两例的 mask/prefix 另在 `scalar_state_collisions.json`、
  后续 `full_triples.json` / `full_quartics.json` 中；其余可由原 counter 恢复。

### 20k cooperative 全档案：已有 exact cooperative rows

同一 `1b5a9dea07e1c62f69798fddbf4899ff986c0b72` source，路径：
`results/local-20260831/P334-cooperative-closure/raw/N325.geometry_pilot.csv` 与
`N425.geometry_pilot.csv`，配套 `.metadata.json`；
freeze 为 `analysis/p334_cooperative_closure_freeze.json`。
实际 runner commit `6712ec5b00bdf3dc0c6f8733ef85eda58b86cb2f`。

| size | base permutations | first/second rank-one rows | unique risk counters | both-risk counters | seed | counter 区间（右端不含） |
|---|---:|---:|---:|---:|---:|---|
| N325 | 20,000 | 8,997 / 8,966 | 13,954 | 4,009 | 20260831430325 | [43032500000,43032520000) |
| N425 | 20,000 | 8,910 / 9,081 | 13,967 | 4,024 | 20260831430425 | [43042500000,43042520000) |

k0 为 193/252。HNF shears：N325 first/second=57/18；N425=132/268。
共 35,954 个 rank-one checkpoint rows；没有理由只筛双方向同时 risk 的交集。
保存字段包括 `batch,replica,k0,k1,k2,age_steps,ell_u,ell_v`、当前 morphology/H2、
`next_site`、common/clone sites/outcomes，及
`checkpoint_b1_safe_count,checkpoint_b2_safe_pairs,checkpoint_sum_child_b1_sq,`
`branch_q_after_safe_count,branch_q_after_denominator,q_after,q_after2`。
没有全 population 的 mask/port/component-incidence 列。

既有 scorer `scripts/score_p334_cooperative_closure.py` 在每 size 的 20,000 base clusters 上
形成两方向各 11 个坐标的完整 22D covariance；Draft 的
`results/p334-fork-directional-allocation/score.json` 消费此联合矩阵，方向结论仍有限。
不能用 22 个图代替这批全档案。

### 100k branching production：独立档案，但未存 exact b2/Rao–Blackwell rows

archive source `analysis/p429-branching-continuation-pilot-20260830`，
commit `751f8b384883b3ce92e5efa77c35f45a86afa84d`。
路径：`results/server-20260830/P429-branching-continuation/production/N325/N325_100k.geometry_pilot.csv`
与 `N425/N425_100k.geometry_pilot.csv`；freeze `analysis/p429_branching_continuation_100k_freeze.json`。
实际 runner commit `c7134775cf6712b6a659022962c9a3ad9efab6e4`。

| size | base permutations | first/second rank-one rows | unique risk counters | both-risk counters | seed | counter 区间（右端不含） |
|---|---:|---:|---:|---:|---:|---|
| N325 | 100,000 | 45,490 / 45,692 | 70,366 | 20,816 | 20260830429326 | [42932600000,42932700000) |
| N425 | 100,000 | 44,761 / 44,613 | 69,437 | 19,937 | 20260830429426 | [42942600000,42942700000) |

共 180,556 rank-one rows。保留相同种类 k1/k2/geometry/H2 与 common/clone 信息，
但 CSV 止于 `branch_both_survive`：**没有 b2、sum_child_b1_sq 或 exact q_after**。
它不是已有 20k 22D covariance 的更多同型 rows；如果未来增加新观测，必须从同一
base-cluster 构造新的联合 influence，不能复用不对应的误差矩阵。

一个重要角向语义差别：100k CSV/metadata 的 `a,b`/`first_rep` 记录的是
`(325,57)/(325,18)` 或 `(425,132)/(425,268)` 这类 HNF 参数；
不能把它们当作 Gaussian 几何代表直接代入 cos(4θ)。对应物理 reps 已在 20k freeze
明确为 N325 `(17,6)/(18,1)`、N425 `(16,13)/(19,8)`，应沿用声明的几何和方向映射。

## 可重用程度、成本与仍未完成的边界

两 runner 的 `splitmix64`、`SplitMixStream::below` 和 `counter_permutation` 源码已逐段读取，
其 uint64 counter-derived stream 与 unbiased Fisher–Yates 相同。
已有 `scripts/p334_checkpoint_scalar_collision.py::archived_permutation` 是对应恢复实现。
每 unique counter 只须恢复一次长度 N 的 permutation，两方向共享 label prefix；
前 k0 个 labels 给 occupied mask，并可使用已有 `next_site` 标识来源。
**档案不存在丢失随机种子造成的不可恢复障碍。** 本轮没有执行这种恢复。

20k 两 size 的工作量是最多 27,921 次 prefix 解码、35,954 个切网络构造；
100k 是 139,803 次 prefix 解码、180,556 个切网络构造。
这给出数据工作量，不是实测耗时。35,954 个网络的 Python/C++ throughput 尚未测量，
不能用两例 DP 的 0.01/0.21 秒承诺全总体时长。

完整 reliability 的 **6358 旧版本**不是已完成的 population adapter：
`scripts/p334_contracted_birth_network.py::build` 固定 N425、shear268 和 transverse covector
(19,8)，并对两例检查删除 essential root 后 full available graph balanced；
`p334_full_birth_reliability.py` 要求观察到的唯一 two-port component。
通用 `safety_polynomial(network,all_sites)` 的 DP 可重用。后续 bd95 已在固定 12-prefix
实现中处理多通道，并全部求解；固定 PR484 也已证明一般几何表示。因此这些旧代码限制
是版本/实现范围说明，**不是当前尚待证明“多通道是否存在/因子化”的科学障碍**。
尚不能据此承诺所有 rows 的低成本求解；PR491 cut、PR484 一般 gain 表示和已有
12-prefix 实现覆盖仍须按各自范围记录。

### 保留为候选，不在本轮执行

- population incidence/overlap 必须先从 occupied geometry 定义 component/port 来源，
  再读取其 response；不能从已揭晓 trigger graph 上挑出恰好解释结果的分组。
  原 BFS cut 的辅助编号不直接作为物理 H4 mark。一个待明确的候选是使用原 occupied
  rank-zero components 的接触集合及其重叠，或另给 covariant/cut-invariant 统计定义。
- 配对统计沿既有规则：每 size 保留全 n 个 base clusters；方向条件均值的 influence
  为 `I_risk*(z-mu)/risk_rate`，非 risk 的该方向为零；两方向在同一 replica 联合。
  不能只取 both-risk intersection，不能把方向/clones 当独立样本。
  两 size、20k 与原 100k counter domains 的依赖关系分别标明。
- 若目标是更广总体完整时钟/conditional thermal readout，12-prefix 的首次 cost map、
  并行 factorization、11 对 crossing 和 47-factor winner allocation 均已有结果。
  尚缺的是更广前缀集合中的成本分布、配对方向差异及条件 suffix covariance 总体权重；
  不是给两例再求 higher triggers/full clock，也不是重做这 12 个条件时钟或 channel derivative。
- white essential-cycle packing 的普遍存在性仍为理论候选；22 例证书扩充也不自动回答
  population 机制。不得把对角 matching 语义省略后宣布 Menger/minmax 已解决它。

## 团队交接

1. 总览/ledger：把 pair、triple、quartic 的旧“下一步”链接到 full/thermal 完成节点，
   再吸收固定 PR484 一般 theorem、bd95 的 12-prefix 和 e3d9782 的实际通道分配。
   保留各次结果原本的适用范围、条件选择与独立 source refs，不把它们计作新 MC replication。
   将 PR491 与492 标成不同描述、同一旧两例依赖块。
2. 生产档案维护：保留上述 22/20k/100k 三种范围与 SHA，特别是 100k 缺 exact cooperative
   列但不缺 counter/k1/k2 的区别；目前不启动重放。
3. 后续科学分工只在完整仓库审读后确定。若选 population 路线，应由一个实现拥有
   geometry-defined marks、显式来源关联和配对 covariance，避免分别重建相同前缀。
4. 本轮不改 Issue 状态、不锁任务、不设置新的许可顺序；这里的候选不是已派发任务。

审读覆盖：#491 完整当前 body 与上述主实现/证明/输入；#492 完整证明笔记；
pair/triple/quartic/full/thermal 五份完整科学笔记；full network/DP 两份完整实现；
20k/100k freeze/metadata/CSV 字段和选择/协方差实现；两个 runner 的 RNG 核心。
补充全文读取固定 705819e 的 general-two-port theorem、b9cbe13e selection manifest、
bd95 的 twelve-prefix 科学笔记和结果状态/精确均值、e3d9782 的 component-race 科学笔记
及结果文件范围；这些均只读保存 artifact，不重算系数。
未运行其附带 tests/verify/replay/census，未声称本笔记替代全仓库上下文审读。
