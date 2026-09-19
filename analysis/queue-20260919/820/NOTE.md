# #820 — R2 先例核查：收窄否定面，降级“新状态族”主张

日期2026-09-19。消费者：#813 的 R2 定位及 #144 分类合同。

## 结论

**发现实质先例，R2 的“独立新状态族”一级定位应降级。** 顶点诱导子图状态和、四端点 all/none site 生成元、按整条超边取舍的嵌入多项式均已有文献。保留下来的研究问题是：仓库的有限商格点 rank-image/双侧 genus 读出，究竟是哪个已有不变量的特化，或需要什么最小增广。**本轮没有证明完整 matching 多项式等于某一个现成标量多项式，也没有证明其新颖性。** 不改写历史文件或直接修改总控分级。

## 1. 六个指名家族逐项判断

这里把“字面保留局部端点状态”的合同和“最终标量多项式相等”分开。单看缺少 J4 能证明前者的障碍，不能自动证明任意自然编码下的后者。若要把局部论证升级成标量否定定理，还须把端点标记、胶合兼容及局部状态保持写成合同的一部分。

| 家族 | 正文核对 | 本轮判定 |
|---|---|---|
| multivariate Tutte / random-cluster | [S] §2 式(2.1)、FK展开式(2.10)，A⊆E | 一个 ordinary edge/site 的局部状态保留实现确被阻挡；不是对超边/顶点版本的否定。 |
| Bollobás–Riordan | [H] §4.6，边子集和及 vertex-partitioned 推广 | 即使预先划分顶点，随机 ground set 仍是 ordinary edges。固定预接线不随 site 开关；单边只能增加一个连通合并，不能把四个原本分离端点一次全接通。 |
| Las Vergnas | [LV] Definition 3.2、Proposition 3.3 | matroid perspective 的自然 ground set 仍是图边，局部障碍同上。它的表面秩读出不自动增加 J4 状态。 |
| Krushkal | [H] 式(4.9)，A⊆E 与两侧 regular-neighbourhood genus | 双侧 genus 是非常近的读出先例，但 ordinary edge alphabet 仍受限。不能据此宣布 typed site 版已是它的直接特化。 |
| topological transition | [C] §4.2 对四价子类的解释 | **四价**版本的三种配对确不等于四个 cluster ports 的 J4；但广义 Eulerian/hypermap transition 允许八价顶点，二态边界平滑的先例直接绕过这个四价论证。必须标明子类，不能把 cluster ports 与 interface half-edges 混为一谈。 |
| 2025 embedded-graph vertex polynomial | [V] Definitions 2.4、3.1、Remark 3.2、Lemma 3.3 | 名字中的 vertex 是 partial/twisted dual 的顶点数；状态索引仍是 A⊆E 或边集合的3/6划分，并非 site 子集。仓库对此项的识别正确，但它不代表全部 vertex-indexed 文献。 |

六项的窄局部障碍不因下节先例而消失；不能再用这六项的名单证明“state family is new”。

## 2. 被遗漏的顶点/超边方向

**顶点可靠性 [N]。** Q(G;x,y)=Σ_{X⊆V}x^{|X|}y^{k(G[X])} 已按 site 子集求和；§8 式(11)通过变量代换和 y 导数给出独立顶点存活后的分量分布。这直接否定“边指标是已有状态和的共同限制”这种泛化。Q 本身没有 ambient H1 标记：同一个抽象圆环图在环面上可作可缩或本质嵌入，Q 不变而满占据 rank 不同。因此不能把 Q 单独说成完整 R_site 的已知答案。

**四端点 site 生成元 [F]。** §2 的端点划分以及 Fig.3(c)、§4 的 generator III 明确使用单个概率 t 的 site，使四个内端点全部连接或全部不连接。J4 的 Bernoulli 局部元不是本仓库首次出现；它也不是所有内部边独立的普通 bond gadget。没有把该文的特定自对偶临界条件移植到当前方格。

**整条超边多项式 [C]。** §3.2 Definitions 4–5 以超边子集 A 求和，局部因子为 u_e^{d(e)-1}；§4.2、Theorem 2 将它映到 medial 的广义 transition polynomial。degree-4 超边对应 degree-8 medial 顶点，非零的两个局部平滑是

    (01)(23)(45)(67),  (07)(12)(34)(56).

它们编码 junction 的边界，不是把四个 cluster 端点两两连接。取平方格 edge-midpoint 为 hypervertices，每个原 site 的四臂星为一条超边；保留/删除整条星与 site 占据逐一对应。孤立 midpoint 仅多出可缩分量，不改变黑 ambient H1。这个**黑连通/同调状态字典是本轮的推导**，不是声称 [C] 已写出 Matching-One 的双侧公式。白侧细胞化及所有退化 quotient 仍需另核对。

重要限制：[C] 的面/边界数权重不是自动的 g_B−g_W 源，未经证明不能通过重新命名变量取得仓库的 charge。八价实现也不是原四价 topological-transition 子类的直接特化。故先例足以撤去“全新二态 junction 状态族”的卖点，但不足以断言完整 rank-image 商已被同一篇文献包含。

**更细超边划分 [W]。** Definition 2.2 对 β≤α 的 permutation refinements 求和，允许比二态更细的局部划分。这也是应比较的先例；本轮不声称它与 [C] 或 R_site 相同。

## 3. 独立边端点支撑引理：成立，专门出处无法确认

合同必须包括：有限无向图、至少三个**互异**端点、每条边独立且0<q_e<1、没有把端点预先识别。若全部端点可连通，取一个连接全部端点的极小树 T。极小性保证每个叶子是端点。取叶端点 t，令唯一开边集为 T 去掉其叶边，其余全闭；则端点划分恰为 {t} 与“其余至少两个端点”两个块。每个指定边状态的概率都严格正，故该真非平凡划分有正概率。此证明不依赖平面性。

[S I] 式(15)提供标准独立边配置概率，Definitions 25–26、32 明确端点与划分概率的语言。**没有在本次已读段落找到上述叶树推论被单独命名的出处；写“无法确认专门出处”，而不是“确认文献无此结果”。** 本轮自证仅消费于 #820 的 gadget 合同，不作为新颖引理。

另一个合同边界：允许多边 block 后，“只能源求导”不等于“不能投影”。三边树令 D=t∂_t，则 (D−1)(D−2)/2 对占据数k=0,1,2,3的权重为1,0,0,1，可以删除混合状态。它违反 one ordinary edge/site，也不是保持正独立 bond 概率的实现；但说明禁止 block projector 必须明确，而不能仅靠“用了导数”作否定。

## 4. 检索记录与读取分级

实际查询（2026-09-19）：`subgraph component polynomial`; `vertex topological reliability polynomial`; `hypermaps dichromatic polynomial`; `2506.07522`; `hypermaps dichromatic polynomial topological transition polynomial`; `hypermap Krushkal`; `The enumeration of vertex induced subgraphs`; `homological reliability polynomial`; `independent terminal partitions percolation`; `all-or-nothing hyperedge independent`; `network reliability partition support`; `percolation three terminal generator connectivity partitions independent bonds all connected`; `hypergraph percolation all none connectivity partitions site bond representation`; `network reliability splitting formula terminal partition lattice independent edges`。

覆盖：普通/嵌入边多项式、vertex-induced/subgraph-component、node reliability、hypermap coarse/refined deletion、generalized transition、percolation terminal generators、K-terminal splitting。专门同调可靠性与所有语言/数据库未穷尽；部分查询返回大量无关结果，未拿它们证明不存在先例。完整 rank-image 先例最可能位于 hypermap 双侧 genus 增广、homological node reliability、partition-algebra transfer 和关联 bond/site 状态表示中。

正文已读如下（均仅指所列段落，不是整篇）：

- [S] A. D. Sokal, *The multivariate Tutte polynomial (alias Potts model) for graphs and matroids*, arXiv:math/0503607v1，§2 式(2.1)–(2.10)。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/math/0503607
- [H] S. Huggett、I. Moffatt, *Types of embedded graphs and their Tutte polynomials*, arXiv:2212.14233，§4.6 的 BR/Krushkal 状态和及式(4.9)。**PRIMARY_TEXT_READ**。没有假装已读两家的原始论文。https://arxiv.org/abs/2212.14233
- [LV] J. A. Ellis-Monaghan、I. Moffatt, *The Las Vergnas polynomial for embedded graphs*, arXiv:1311.3762，Definition3.2、Proposition3.3。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/1311.3762
- [V] Qi Yan、Qingying Deng、Metrose Metsidik, *Introducing a vertex polynomial invariant for embedded graphs*, arXiv:2506.07522，§2–3 上述定义及 p6 截图。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/2506.07522
- [N] P. Tittmann、I. Averbouch、J. A. Makowsky, *The Enumeration of Vertex Induced Subgraphs with respect to the Number of Components*, arXiv:0812.4147，§2.1–2.2、§8式(11)，定义页截图。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/0812.4147
- [C] J. A. Ellis-Monaghan、I. Moffatt、S. Noble, *A coarse Tutte polynomial for hypermaps*, arXiv:2404.00194v2，§2、§3.1–3.2、§4.2 Definitions8–9/Theorem2；印刷p11、p21截图成功，p20截图请求失败，那里只读解析正文、不声称看清其图9。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/2404.00194
- [W] R. Cori、G. Hetyei, *A Whitney polynomial for hypermaps*, arXiv:2311.06662，§2 Definition2.2及其紧接解释。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/2311.06662
- [F] O. Khatib Damavandi、R. M. Ziff, *Percolation on hypergraphs with four-edges*, arXiv:1506.06125v2，§2、Fig.3(c)、§4 generatorIII；印刷p4、p7截图。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/1506.06125
- [S I] F. Simon, *Splitting the K-Terminal Reliability*, arXiv:1104.3444，式(15)、Definitions25–26、32。**PRIMARY_TEXT_READ**。https://arxiv.org/abs/1104.3444
- R. K. Wood 1985 *A factoring algorithm using polygon-to-chain reductions for computing K-terminal network reliability*，doi:10.1002/net.3230150204。**ABSTRACT_ONLY**；不作支撑引理的专门出处。

## 本轮执行范围

实际读取两份 p144 分支笔记及上述正文，执行 `python check.py`：3/4/5顶点的全部8/64/1024个简单图中，12,888个满足全端点连通的端点集合均构造出严格中间划分，零违反；检查八半边105种配对只剩两种 checkerboard 非零状态；Fraction 验证三边 block 投影权重。没有执行格点/环面枚举、旧 oracle、临界点计算、完整标量特化证明、父PR测试或全仓CI。没有修改旧结果、状态文件或 frozen scorer，没有合并PR。
