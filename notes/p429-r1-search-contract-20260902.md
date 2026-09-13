# 本地执行提示词：Matching-One P1 有界摘要非压缩搜索

> Provenance, 2026-09-02, additive, does not rewrite the frozen protocol.
> The two search-protocol tokens in §0 and §7 remain locked:
> `NO_COMPRESSION_WITNESS_FOUND` and `BOUNDED_SEARCH_CLOSES_UNDER_SUMMARY`.
> Manuscript and certificate wording for the same mathematical outcome uses
> `BOUNDED_SUMMARY_INSUFFICIENT`. Size claims in the certificate are “smallest
> witness found in the declared enumerated families”, not global minimality.
> See `notes/p429-r1-claim-wording-erratum-20260902.md`.

把本文件从「角色」到「停止规则」整段复制给本地 agent。不要改冻结合同。不要新开 Monte Carlo，不要写论文，不要合并 PR，不要开新 Issue。不需要 GitHub 权限：全部在本地仓库目录完成。

---

你是 Matching-One（GitHub: LightChainr/Matching-One）P1 论文的定理加强执行 agent。工作目录是本地仓库根目录。用确定性精确计算，不要猜。不申请、不使用 GitHub token。

## 0. 任务是什么

P1 已有闭合定理链：生存律非闭、cut-network 表示、更新闭包、平行 two-port 分解、无界预测类族。本题**不是**投稿前提，只回答一个结构问题：

> 是否存在一对有限平面两点顶点网络，它们在**预先冻结**的摘要下完全一致，但在**预先冻结**的深度-2 复合语言上至少一个实验概率不同？

只许返回下面两个裁决之一：

- `NO_COMPRESSION_WITNESS_FOUND`
- `BOUNDED_SEARCH_CLOSES_UNDER_SUMMARY`

禁止第三种措辞（“几乎闭合”“倾向非压缩”“需要更大 n”）。

## 1. 冻结合同（搜索前锁死，禁止自适应）

### 1.1 摘要（禁止事后加描述子）

对每个有根两点网络记录且只记录：

1. 完全安全子集多项式 S(z)：按占据基数 k=0..n 计安全子集个数，S_k = #{mask : |mask|=k, occupied(mask) 不把 L 连到 R}
2. 可开关顶点数 n（L、R 是永久端点，不进 n）
3. H2 = n − S_1（单点终端连接数；由 S,n 决定，冗余）
4. b2 = C(n,2) − S_2（最小对触发数；由 S,n 决定，冗余）
5. 端点局部有根邻域：先 r=1，若 r=1 在声明搜索类上闭合，再单独做一次 r=2。邻域 = {L,R} 的半径 r 诱导球，顶点类型 (dist_L, dist_R, touches_L, touches_R)，在类型保持的置换下取典范。

**后继 hazard 矩（含 delayed-fork 用的 ∑x²）不进摘要。** 冻结在外。不要看见反例后再加进去。

H2、b2 与 (S,n) 等价，分组键实质是 (n, S, neigh_r)。

### 1.2 冻结实验语言（深度 ≤2）

观察只许是 L–R 不连通 / rank-one 生存。顶点占据均匀：长度为 p 的共享前缀是 n 个可开关点上的均匀有序 p-元组；吸收的前缀贡献 0。之后每个克隆独立走长度 c≤2 的延续。

必须精确计算的实验，报告分裂时按此顺序取**第一个**不同者：

```
E1_c2, E1_mix, E2_c1, E2_c2, E2_mix, E0_c2, E0_mix, E0_c1, E1_c1
```

含义：Ep_cq = 共享前缀长 p，两克隆延续长都是 q；mix = 一克隆 c=1、另一克隆 c=2。
E1_c1 是 delayed-fork，已知一旦供给后继二阶矩就被决定。它排在最后：若更深实验也分裂，不要把它当作第一区分器。

普通延续 P(survive k) 由 S(z) 完全决定，只作健全性检查。

概率用 `fractions.Fraction`，禁止浮点、禁止 Monte Carlo。

### 1.3 图范畴

- 简单无向图；L、R 永久；可开关点 0..n−1 为随机变量
- 连通载体（L、R 与每个可开关点在同一连通分支）
- 存在 L–R 路；无 L–R 直边
- 平面：G ∪ {L–R} 平面（L、R 可放在同一面）
- 连通性是**顶点占据**两点连通，不是独立边可靠度

### 1.4 声明的搜索类 G（显式有界生成族）

协议允许“n≤12 穷举平面图”或“等价的显式有界生成族”。采用后者，必须写进证书：

| 代号 | 定义 | 硬上限 |
|---|---|---|
| EXH | 穷举连通两点简单图，n≤5，G∪{LR} 平面 | n≤5 |
| HID | 对 EXH 的每个核做走廊隐藏：两侧 hop∈{1,2,3} | n≤12 |
| SP | 两点 series-parallel，从单位边与 path(1) 用 series/parallel 生成 | **n≤12，本机必须跑到 12** |
| W | Wheatstone（n=2 加对角）与 SP 的 series/parallel 复合 | n≤12，SP 伙伴可用 n≤6 |
| MP | 多路（平行若干 path） | n≤12 |
| GR | 2/3/4 行梯子，L 接左列、R 接右列 | n≤12 |

走廊约定（必须统一，否则 n 对不上）：

```
series(A,B)：把 A.R 与 B.L 粘成新可开关点 z，n = nA+nB+1
hops_hide(core, Lh, Rh) = series(path(Lh-1), series(core, path(Rh-1)))
path(0) = 单位边（n=0, lr_edge=True）
因此 Lh=Rh=1 给出 n = core.n+2，L 只邻 zL、R 只邻 zR
```

不要只复现 PR #549 的平行放大。成功结果必须是：在**本摘要合同**下丢失了冻结语言的信息。

## 2. 算法（确定性）

对每张图：

1. `safe_table[mask]`：占据 mask 是否仍 L–R 不连通。2^n DP/BFS，n≤12 可承受。
2. S_k = 按 popcount 累加 safe。
3. 冻结实验：前缀用整数计数，最后才化成 Fraction。p=2 只枚举无序对再 ×2。n≥4 时剩余点数固定，分母是常数，内环禁止 Fraction。
4. r=1（必要时 r=2）邻域典范键。
5. 按 (n,S,H2,b2,neigh_r) 分组。只对**多成员摘要类**算实验。
6. 类内行为向量（上列 9 个实验）若多于 1 个，记一次分裂。
7. 找到分裂后：用删边/删点约化，保持“同摘要、不同行为”；输出近最小对 + 独立证书。

同构：WL 着色 + 胞腔置换/回溯，L、R 固定。穷举 EXH 与 SP 生成需要典范键去重。HID 狩猎阶段可以不去重（不同核 + 相同走廊几乎必不同构）；报告见证时再比 canonical_key。

平面性：n≤5 穷举必须用禁示子/收缩（K5、K33 + 边收缩，v≤8）。HID/SP/W/MP/GR 由构造保证平面。v>8 不要对一般图做收缩搜索。

## 3. 本机性能（必须按此组织，否则会跑死）

- 语言：Python 3.10+，**只用标准库**（无 networkx）
- EXH n=5：约 2^20 边掩码，连通 + 典范 + 平面，约 1 CPU 分钟，约 6775 张平面图；EXH n≤5 合计约 7398。立刻 pickle。
- SP 增长约 ×3.2/点：n=8≈1401，n=9≈4798，n=10≈16750，n=11≈4.5 万，n=12≈14 万。按 n 分层生成并每层 pickle。n=11–12 的 series 笛卡尔积大，但本机应当跑完；不要在 n=10 停。
- 狩猎两遍：先 S+邻域，只对多成员摘要类跑实验。
- HID L1R1：7398 张 n=core+2，实验很便宜。L2R2 同样。
- 不要对 HID 全量 canonical_key。
- 实验整数核：n=10 全量 SP 多成员类约 1 万张，整数核大约数十秒；Fraction 内环会慢一个数量级。

独立健全性（实现后先过再搜索）：

- path(1)：S=(1,0)
- path(2)：S=(1,2,0)；E1_c1=0；E0_c1=1
- 两个平行 path(1)：S=(1,0,0)
- Wheatstone：S=(1,0,0)，但 r=1 邻域 ≠ 两个平行 path(1)
- series(path(1),path(1))：n=3，S=(1,3,3,0)

## 4. 已知可复现见证（当作回归，不是答案的替代）

本地实现正确时，**必须**收回下面这对。允许找到更小的；不允许“找不到就改摘要”。

n=5 穷举中唯一在忽略邻域时实验就分裂的 S-类：

```
S_core = (1, 5, 8, 4, 0, 0)
40 张：E1_c1 = 13/20，E1_mix = 1/3，succ_m2 = 4/5
16 张：E1_c1 = 27/40，E1_mix = 41/120，succ_m2 = 6/5
E1_c2、E2_* 在该类上仍相同
```

两端 1-hop 隐藏后 n=7，S 变为

```
S = (1, 7, 21, 35, 33, 15, 2, 0)
H2 = 0, b2 = 0
```

主见证（无悬挂、各 10 边，r=1 同构，r=2 不同，只分裂 E2_c2）：

**G_A** 核 `exh_n5_4313` 隐藏。顶点 0..6，ℓ=5，r=6。

```
L — 5 — 2 — 1 — 6 — R
         \ /
          3
L — 5 — 4 — 0 — 6 — R
```

边：`[0,4],[0,6],[1,2],[1,3],[1,6],[2,3],[2,5],[4,5], L-5, R-6`
adj_L=1<<5，adj_R=1<<6
adj = `(1<<4|1<<6,  1<<2|1<<3|1<<6,  1<<1|1<<3|1<<5,  1<<1|1<<2,  1<<0|1<<5,  1<<2|1<<4,  1<<0|1<<1)`
4-点割：`{5,2,1,6}` 与 `{5,4,0,6}`，只交于 `{5,6}`
P(E2_c2)=**937/1050**

**G_B** 核 `exh_n5_2451` 隐藏。

```
L — 5 — 2 — 0 — 6 — R
         \     /
          1   4
          |
          3 — 6 — R
```

边：`[0,2],[0,4],[0,6],[1,2],[1,3],[2,5],[3,6],[4,5], L-5, R-6`
adj = `(1<<2|1<<4|1<<6,  1<<2|1<<3,  1<<0|1<<1|1<<5,  1<<1|1<<6,  1<<0|1<<5,  1<<2|1<<4,  1<<0|1<<3)`
4-点割：`{5,2,0,6}` 与 `{5,4,0,6}`，交于 `{5,0,6}`
P(E2_c2)=**313/350 = 939/1050**
间隙 **1/525**

其它冻结实验两边相同：E0_*=1，E1_c1=E1_c2=E1_mix=E2_c1=1，E2_mix=33/35。

机制：两边都恰好 2 个连通 4-集（S_4=33）。S(z) 看不见割的相交；E2_c2 是 2-前缀上 p_surv(剩余,2)² 的均值，能看见。Delayed-fork 两边都是 1，故**不是**后继二阶矩可观测。

更小边数变体：悬挂核 `exh_n5_1056` 隐藏 vs G_B，9+10 边，同一 S/r=1/E2_c2 间隙。删点得不到更小见证。

SP 旁证（n≤10 已见）：至少 29 个 r=1 摘要类分裂，最小 n=8、12 边、E2_c2 间隙 2/1575。典型形态是 2-hop 走廊串两个同 S 的 SP 核。本机把 SP 推到 n=12 后，把新增分裂计入证书，不要改主见证，除非找到 n<7 或同 n 更少边且无悬挂。

r=2 窥视：同一 56 元 S-类做 L=R=2 隐藏后 r=1 与 r=2 都对齐，但 E2_c2 间隙消失。不要把主对写成 r=2 反例（它们的 r=2 邻域本就不同）。

## 5. 搜索流程（按序，写日志）

1. self-check
2. EXH n=1..5，pickle，r=1 狩猎；闭合则再 r=2。预期：两半径都闭合。
3. HID Lh=Rh=1 对全部 EXH 核，r=1 狩猎。预期：收回 n=7、E2_c2、间隙 1/525。取近最小对，删边/删点约化。
4. HID Lh=Rh=2 对全部核，r=1；再 r=2 窥视。预期：该构造上无分裂。
5. SP 生成到 **n=12**，每层 pickle，r=1 狩猎。n≤10 预期已有分裂；n=11–12 是本机新增。
6. W、MP、GR，r=1。
7. 主见证取最小 n，其次少边、无悬挂、每个可开关点在某条简单 L–R 路上。
8. 写证书。r=1 已分裂则裁决 A，**不要**为了“再冲 r=2”改合同或加描述子。

## 6. 输出（必须落地到本地文件）

写三个文件：

- `artifacts/bounded_summary_search.json`  schema: `matching-one/bounded-summary-search/v1`
- `artifacts/bounded_summary_search.md`
- `research/summary_search/verify_witness.py`  硬编码主对关联表，不依赖搜索，跑完打印 `VERIFY_OK`

JSON 至少含：

```
verdict
frozen_summary（写明后继矩未包含）
frozen_experiments（上列顺序）
enumerated_class：各族规模、n 直方图、max_n、SP 实际到达的 n
r1 / r2 统计
witness 或 closing_radius
```

见证字段：A/B 关联表、n、S、H2、b2、split_experiment、p_A、p_B、gap（有理数 a/b）、全部冻结实验有理数、planar(G∪{LR})、same_S、same_neigh_r1、same_neigh_r2、4-点割、约化记录、环面嵌入一句（见下）。

环面嵌入只许写：对在已闭合的平面两点 cut-network 范畴；cut-network 表示把 rank-one 环面状态映到该范畴；two-port 演算可把 gadget 当两点块。**不要**声称已构造具体环面占据，除非你真的构造了。

## 7. 裁决模板

### A. `NO_COMPRESSION_WITNESS_FOUND`

给出：两张图的关联表；平面两点验证；相同 S(z)；相同冻结摘要；第一个冻结分裂实验；精确有理概率与间隙；短解析机制；机器证书路径；环面嵌入句。然后只对**测过的表示类**给一段手稿定理。禁止升格为“任何有限标量表示都不行”。

手稿定理骨架（数字换成你的）：

> 存在两个连通平面两点顶点网络 G_A、G_B，各有 N 个可开关顶点且两端点共面，使得完全安全子集多项式相同、H2 与 b2 相同、半径-1 端点局部有根邻域作为带类型图同构，但 P(E; G_A)=p_A ≠ p_B=P(E; G_B)，其中 E 是冻结深度-2 语言中的指定实验。因此 (S(z),n,H2,b2,r=1 邻域) 不是该冻结语言的充分统计。不断言 r=2 失败，不断言欧氏隐维下界，不含连续极限/CFT。

若收回第 4 节那对：N=7，E=E2_c2，p_A=937/1050，p_B=313/350，间隙 1/525。

### B. `BOUNDED_SEARCH_CLOSES_UNDER_SUMMARY`

仅当声明类在 r=1 **以及**（若 r=1 闭合）r=2 上都无分裂。写清类定义与规模、到达的最大 n、闭合半径、关闭信息是否有更小状态含义、一段“cut-network 仍是充分但非极小表示”的手稿段。负结果禁止升格为全图定理。

若你收回了第 4 节那对还报 B，视为实现错误，先修算术再搜索。

## 8. 禁止声称

- 欧氏隐空间维数下界
- 一个实标量不可能编码状态
- 微观占据过程非 Markov
- 连续记忆、CFT、LCFT、场计数
- 对每个最近邻 HNF 环面普遍成立
- r=2 非压缩（除非你真的有同 r=2 邻域、不同冻结实验的一对）
- 全图极小充分统计

## 9. 停止规则

遇到下列之一立刻停，写证书，回到投稿模式：

1. 有效见证已约化且独立脚本复算通过；或
2. 声明类在 r=2 上穷尽到合同上限（EXH n≤5 全；HID hops≤3 n≤12；SP **n≤12**；W/MP/GR n≤12）

不要加顶点上限，不要加描述子，不要开 Monte Carlo，不要开新 Issue，不要改手稿正文（本任务只出证书与定理段落）。不要把 PR #549 平行 gadget 当作本题成功。

## 10. 仓库纪律

- 只新增：`research/summary_search/*`、`artifacts/bounded_summary_search.json`、`artifacts/bounded_summary_search.md`
- 不改 `notes/p429-*`、不改 `docs/PUBLICATION-PORTFOLIO.md`、不 merge、不 force-push、不连 GitHub API
- 已有相关 PR（#549 生存律 no-go，#491 cut-network，#435/#434 N16 反例，#484 two-port）当作只读背景。本题不重做那些定理。若本地没有 clone，在当前工作目录新建 `research/summary_search/` 即可，不强制访问远程仓库。

开始执行。先 self-check 与 EXH pickle，再 HID L1R1 回归第 4 节有理数，再 SP 到 n=12。
