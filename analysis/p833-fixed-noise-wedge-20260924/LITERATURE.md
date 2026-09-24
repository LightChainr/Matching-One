# #833 固定噪声查询：定向扩展检索与先例边界

检索日：2026-09-24。消费者：NOTE.md 的楔形充分构造、接口擦除下界和下一步路线。13项主记录：8项读到下列指明正文，5项仅摘要；另记一项1992会议论文的元数据入口。不是13篇全文逐页审完，也不是系统原创性认证。前轮PR #836已有一般正实现/PSR矩阵，本次不重复把它算作新增阅读。

检索覆盖：site/matching分隔、wedge/sponge窄条渗流、可靠连接与路由、信息收缩和擦除比较、带噪声公式/电路、可靠细胞自动机、随机访问码。定向词包括 percolation wedge logarithmic width, corner connection high density, fault tolerant percolation routing, expanding wedge site crossing, noisy formula amplification, information erasure cut, random access code lower bound；沿最近的作者/期刊来源回溯原文。不根据搜索无命中认证不存在先例。

## 核心结论

可靠长条连接、Peierls路径计数、对数厚度及随机访问码下界都已有直接前例。这里需要独立补足的是：在一个含旧微观比特的NN site环面内，用不读取该比特的非均匀独立新站点源，既让1接通，又阻止0被错误旁路。正文的充分构造没有把任何外部bond定理冒充这个site定理。

## L1 — Duminil-Copin：初等Peierls计数

**来源：** Hugo Duminil-Copin, *Introduction to Bernoulli percolation*，所取讲义封面日期October 7, 2018。作者IHÉS文件：
https://www.ihes.fr/~duminil/publi/2017percolation.pdf

**阅读层级：PRIMARY_TEXT_READ（指定部分）。** Theorem 1.1及其二维对偶/Peierls证明、式(1.1)、site练习入口；相关PDF页截图已核。文件名含2017，不把它当正文版本日期。

**实际输入与限制：** 用路径数乘独立坏点/坏边概率给高密度连接界，是经典方法。主讲义的bond对偶不直接给本文带退化端点的4/8有限楔边界引理。本文另证有限分隔，并用8邻接简单路径至多8·7^(n−1)的粗上界；不主张这套方法新颖。

## L2 — Grimmett–Li：site matching与面中心

**来源：** Geoffrey Grimmett, Zhongyang Li, *Percolation critical probabilities of matching lattice-pairs*, arXiv:2205.02734v3，2024-02-20：
https://arxiv.org/abs/2205.02734v3
https://arxiv.org/html/2205.02734

**阅读层级：PRIMARY_TEXT_READ。** §1 matching图定义、§2.1嵌入、§2.2增补面中心的平面构造。

**实际输入与限制：** 方格NN的matching伙伴含面内对角连接；白面中心可避免把交叉对角线解释成独立无交路径。本文只用分隔的几何部分，不调用该文临界点相等结论来证明查询。没有声称该文已经处理我们的两个底角、尖点及有限噪声常数。

## L3 — Damron–Lam：对数子图与窄条穿越

**来源：** Michael Damron, Wai-Kit Lam, *Asymptotics for first-passage percolation on logarithmic subgraphs of Z²*, arXiv:2502.18235v2，2025-03-05：
https://arxiv.org/abs/2502.18235v2
https://arxiv.org/html/2502.18235v2

**阅读层级：PRIMARY_TEXT_READ。** §1.1的楔形/临界海绵尺寸背景、Proposition 2.5和Corollary 2.7的陈述与证明。

**实际输入与限制：** 长窄矩形中位置数量与指数连接代价竞争，提供h与log w相配的直接尺度先例。对象是bond/FPP，不是非均匀site概率下连接到指定旧站点的查询。下一步若要提高允许噪声，需要把实际site相关长度与端点接入分开；不能简单替换文献中的p符号。

2026年的机构新刊记录中也检索到此题名，但本次推导锁定可核v2，不把平台/新刊日期当成2026新增定理，也不依赖未获取的最终刊本改动。

## L4 — Naor–Wool：渗流式可靠quorum

**来源：** Moni Naor, Avishai Wool, *The Load, Capacity, and Availability of Quorum Systems*, SIAM J. Comput. 27(2), 423–447 (1998)，DOI 10.1137/S0097539795281232。作者页面及其PDF入口：
https://www.avishaiwool.sites.tau.ac.il/publication/journals
https://epubs.siam.org/doi/10.1137/S0097539795281232

**阅读层级：PRIMARY_TEXT_READ。** §5.1 Definitions 5.1–5.3、Lemma 5.5、Propositions 5.7–5.8、Appendix A.1。读取作者PDF抽取正文；本工具对相关图页截图失败，未宣称图像已核。

**实际输入与限制：** 平面路径结构可以获得很强的可用性保证，是“用渗流实现可靠连接”的直接旧先例。其quorum元素与成对的primal/dual边失效机制不同于本题独立site噪声；连接一个已有微观比特、并在0时保持隔离，不由其quorum可用性自动推出。

## L5 — Angel–Benjamini–Ofek–Wieder：巨型分量中的路由

**来源：** Omer Angel, Itai Benjamini, Eran Ofek, Udi Wieder, *Routing Complexity of Faulty Networks*, Random Structures & Algorithms 32(1), 71–87 (2008)，DOI 10.1002/rsa.20163；预印本2004：
https://arxiv.org/abs/math/0407185
https://arxiv.org/pdf/math/0407185

**阅读层级：PRIMARY_TEXT_READ。** §1.1的访问/路由合同、§1.3 Theorem 4、§2 Lemma 5及证明。PDF截图工具失败，结论基于可读正文，不基于未核图像。

**实际输入与限制：** 超临界格点中的路由复杂度结论含“端点已经连通”的条件，路由器还能探查边。我们的查询源不能看到故障图；错误率必须无条件计入微观端点失效。把它的条件性路由成功当成我们的无条件位查询成功，会直接丢掉NOTE.md §6的擦除底。

## L6 — Polyanskiy–Wu：信息渗流与擦除比较

**来源：** Yury Polyanskiy, Yihong Wu, *Application of information-percolation method to reconstruction problems on graphs*, Math. Stat. Learn. 2(1), 1–24 (2019)；本次读arXiv:1806.04195v2，2020-05-21：
https://arxiv.org/abs/1806.04195v2
https://arxiv.org/html/1806.04195v2

**阅读层级：PRIMARY_TEXT_READ。** §4 Theorem 3的less-noisy/erasure比较、Remark 3–4及证明。

**实际输入与限制：** 通道收缩可以与擦除图比较，但less-noisy、more-capable和显式Blackwell降级不是同一关系。本文接口下界只用一个直接公共随机数耦合：目标第一行端口被强制关闭后，两旧准备从此相同。没有把一般通道比较定理加强成其未声称的后处理关系。

## L7 — Nayak：随机访问码

**来源：** Ashwin Nayak, *Optimal lower bounds for quantum automata and random access codes*, arXiv:quant-ph/9904093 (1999)：
https://arxiv.org/abs/quant-ph/9904093
https://arxiv.org/pdf/quant-ph/9904093

**阅读层级：PRIMARY_TEXT_READ。** Theorem 2.3和§4.4证明；相应PDF页截图已核。

**实际输入与限制：** 任意比特可可靠读取时的熵下界是已有方法；本文只使用其经典Fano版本，并给出自包含链式推导。新输入是此NN格点中固定噪声、对数高度的可靠位查询，不是信息不等式。潜在标签数、存储比特数、普通非负矩阵秩保持区分。

## L8 — Duminil-Copin–Tassion：site亚临界指数衰减的下一步入口

**来源：** Hugo Duminil-Copin, Vincent Tassion, *A new proof of the sharpness of the phase transition for Bernoulli percolation and the Ising model*, arXiv:1502.03050v3，2018-01-21；相关期刊DOI 10.1007/s00220-015-2480-z：
https://arxiv.org/abs/1502.03050v3
https://arxiv.org/html/1502.03050

**阅读层级：PRIMARY_TEXT_READ。** §1.1 Theorem 1.1、§1.2的site适配说明。

**实际输入与限制：** 有限程亚临界连接概率的指数衰减可以替换粗路径数来改善大噪声区间，但仍需对非均匀楔域的端口失败和旁路另作界。当前证明在极小固定ε下完全使用初等路径计数，不把未算的相关长度常数称为显式保证。

**日期核查：** 本次HTML页顶部出现2026渲染日期，但arXiv提交史明确最后修订为2018-01-21。不得把它写成2026年新sharpness论文。

## L9 — McDiarmid：格点子域渗流

**来源：** Colin McDiarmid, *Percolation on subsets of the square lattice*, Journal of Applied Probability 17(1), 278–283 (1980)，DOI 10.2307/3212947：
https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/percolation-on-subsets-of-the-square-lattice/3196989E4DD7FBF7ABA7D3C4AC551A85

**阅读层级：ABSTRACT_ONLY。** 期刊元数据与摘要。平台2016上线日期不是发表年份。

**相关性与不能用的部分：** 子区域中的正概率无限连接是楔域问题的重要历史先例；本次没有读到足以直接引用的端点误差或有限域定量常数，不能把“正概率存在”写成“以99%概率读取指定旧比特”。

## L10 — Grimmett：临界海绵尺寸

**来源：** Geoffrey Grimmett, *Critical sponge dimensions in percolation theory*, Advances in Applied Probability 13(2), 314–324 (1981)，DOI 10.2307/1426686：
https://www.cambridge.org/core/journals/advances-in-applied-probability/article/abs/critical-sponge-dimensions-in-percolation-theory/03E257AD0A038B9BF97B07EB3B38C604

**阅读层级：ABSTRACT_ONLY。** 原期刊摘要/书目信息；全文相关结论的后继比较另见已读L3。平台2016上线不作为原年代。

**相关性与限制：** 宽度log长度的几何尺度不是本项目新发明；不能仅凭该摘要确定原定理所有site/bond、边界和概率常数约定。

## L11 — Evans–Schulman：带噪声电路中的信息传递

**来源：** William S. Evans, Leonard J. Schulman, *Signal propagation and noisy circuits*, IEEE Trans. Information Theory 45(7), 2367–2373 (1999)，DOI 10.1109/18.796377；Caltech作者存档：
https://authors.library.caltech.edu/records/mcye7-j6k13

**阅读层级：ABSTRACT_ONLY。** 作者机构记录摘要可读，链接PDF访问失败，不把它算成正文审阅。

**相关性与限制：** 含噪门的多路径放大/信息衰减是应检索的先例，但额外计算层通常能够再生信号；旧site位只有一次向第一新行的接口，不能在未编码冗余前就等同可反复读取的电路输入。

## L12 — Dubiner–Zwick：read-once放大

**来源：** Moshe Dubiner, Uri Zwick, *Amplification by Read-Once Formulas*, SIAM J. Comput. 26(1), 15–38 (1997)，DOI 10.1137/S009753979223633X：
https://epubs.siam.org/doi/10.1137/S009753979223633X

**阅读层级：ABSTRACT_ONLY。** 出版者摘要、作者发表列表。另沿作者页定位1992会议论文 *Amplification and Percolation*；其旧PDF链接不可访问，仅METADATA_ONLY，不额外算作一篇读过的研究：
https://www.cs.tau.ac.il/~zwick/online-papers.html

**相关性与限制：** 放大功能与相关性控制应作为准备编码的潜在工具，而不是对当前单出口位复制输入后仍声称原合同没有改变。

## L13 — Gács：可靠细胞自动机的真实冗余模型

**来源：** Peter Gács, *Reliable Cellular Automata with Self-Organization*, J. Stat. Phys. 103, 45–267 (2001)；2024-01-25的修订加强版arXiv:math/0003117v2：
https://arxiv.org/abs/math/0003117v2

**阅读层级：ABSTRACT_ONLY。** 231页正文未作逐定理审阅；摘要明确该版是对2001刊本的修正加强。初稿2000、刊本2001、修订2024，不能只按最近网页时间标注年代。

**相关性与限制：** 主动局部修复和分层编码能研究长期保存信息；我们的旧图不被持续修复，信息只经第一行接口传出。下一步若引入冗余准备，必须计算新增出口与独立旧位数量的代价，而不是借用该结果跳过端口条件。

## 路线取舍

1. 对当前微观准备，NOTE.md已给足够小固定ε的O(logw)查询作者证明；与PR #836下界匹配高度阶数，不再把这个目标保持为纯猜想。
2. 后续最有区分力的问题是单出口最优误差：当前证明给O(ε)上界与补件中的[1−(1−ε)^5]/2=(5/2)ε+O(ε²)下界，差别是常数和可否达到擦除底。可先在不更改旧准备的前提下分析，不购买更大宽度机器。
3. 若追求固定ε下误差趋零，端口擦除下界迫使改变准备或允许独立重新制备。那是新的编码合同，应把出口数d、总宽度和保留的独立信息位数一起计量。
4. 原自然iid三行中该准备族的质量仍指数小；本次的最坏情况结论没有自动成为自然平均复杂度定理。

## 检索透明性

所有声称读到正文的项目都只指上述明确段落。本工具对L4、L5的PDF图页渲染失败；文本可读，但图形未独立核查。对L11及1992旧链接的PDF读取失败已明确保留。外部PDF不随研究包分发。没有发现一份可直接覆盖“旧微观位+非均匀site源+同一环面rank读出”的现成定理；这句话只描述本轮命中范围，不认证未发现的先例不存在。
