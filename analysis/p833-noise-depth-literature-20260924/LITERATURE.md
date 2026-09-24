# #833 文献定位矩阵与后续路径

检索日期：2026-09-24。下面列17项有实际内容访问的一手文献/作者机构记录：11项读到指定正文，6项仅摘要。`PRIMARY_SECTION_READ` 是读了写明的章节或定理，不等于从头到尾审完论文。PDF中的相关公式页已经视觉核对；未取得的原典不以二手引用替代其正文。此清单不是系统综述，也不能认证本仓库定理新颖。

## 1. 检索范围与实际阅读

|编号|文献与一手入口|阅读状态与定位|对本项目的实质输入|不能直接搬来的部分|
|---|---|---|---|---|
|L01|Benvenuti–Farina, *A tutorial on the positive realization problem* (2004), DOI 10.1109/TAC.2004.826715；[作者机构记录](https://iris.uniroma1.it/handle/11573/232104)|ABSTRACT_ONLY，机构摘要|正实现阶数可大于普通阶数是既有问题；不能把“间隙存在”本身作为本仓库创新点|正文定理未审；不从摘要认领任意控制族、任意准备族的最小阶数判据|
|L02|Benvenuti, *Minimal positive realizations: A survey* (2022), DOI 10.1016/j.automatica.2022.110422；[机构原记录](https://dis.uniroma1.it/publication/25841)|ABSTRACT_ONLY。另一机构仓储列有预印本，但本轮正文访问返回403|整理一般正系统、最小性、谱/锥与经典困难的参考入口|这是LTI传递函数综述，不是已核的受控概率字串模型通用最小化定理|
|L03|Anderson, *The realization problem for hidden Markov models* (1999)；[ANU作者仓储](https://openresearch-repository.anu.edu.au/entities/publication/9c5fdbdc-c7d6-4845-8450-f15ca78ac466)|ABSTRACT_ONLY|HMM实现、Hankel与多面体结构的历史背景|不把搜索到摘要写成已核全部证明，也不等同观测空间维数与正状态数|
|L04|Jaeger, *Observable Operator Models for Discrete Stochastic Time Series* (2000), DOI 10.1162/089976600300015411；[MIT原刊](https://doi.org/10.1162/089976600300015411)|ABSTRACT_ONLY，官方摘要；正文受访问限制|线性依赖过程、OOM与HMM的区别早已明确|未读其正文构造，不能把本次具体图模型与其全部假设自动等同|
|L05|Littman–Sutton–Singh, *Predictive Representations of State* (NIPS 2001)；[原文PDF](https://papers.nips.cc/paper_files/paper/2001/file/1e4d36177d71bbb3558e43af9577d70e-Paper.pdf)|PRIMARY_SECTION_READ，§1–2、Theorem1、Lemma1；PDF对应页面6、8、10（零起始索引）|受控action-observation tests、概率向量、算子闭合和递归更新的直接原典。基搜索和“线性PSR不大于POMDP”已有|这不是正最小性证书；其条件更新含归一化，不能把有符号坐标当潜在概率。论文还讨论组合结构可能导致的指数紧凑性，因此一般“指数分离”也不能仅凭名称认领|
|L06|Monràs–Winter, *Quantum learning of classical stochastic processes: The Completely-Positive Realization Problem*, arXiv:1412.3634；[正文](https://arxiv.org/html/1412.3634)|PRIMARY_SECTION_READ，§1–2，尤其Theorem6及其后的几何说明|等价正实现与共同不变的尖凸多面锥关联；本文定理6明确归属Dharmadhikari1963|主要形式是单一过程的quasi-realization；本项目全部准备、多个可控输入需分别携带约束。量子后半篇未审，不据此主张量子/上下文性|
|L07|Taghavian–Sjölund, *Minimal positive Markov realizations*, arXiv:2502.21102v3 (2025)；[正文](https://arxiv.org/html/2502.21102)|PRIMARY_SECTION_READ，§I–III以及§V的界的解释|在Markov规范形内用线性规划最小化；某些谱类中上下界匹配|作者明确说所得最小Markov型阶数一般只是无限制最小正阶数的上界。这里“Markov”是实现规范形，不等于本题所有受控HMM；不能把LP未找到较小规范形当不可能性|
|L08|Hsu–Kakade–Zhang, *A Spectral Algorithm for Learning Hidden Markov Models*, arXiv:0811.4413v6 (2012)；[正文](https://arxiv.org/html/0811.4413)|PRIMARY_SECTION_READ，模型、可观表示及误差定理的rank/奇异值假设（§2–4）|精确rank与可学习性之间有条件数、奇异值和样本误差的独立接口|原发射矩阵满列秩假设不被“只有2或3个输出、数百隐状态”自动满足；可考虑输出块，但必须重新计成本和条件数|
|L09|Faenza–Fiorini–Grappe–Tiwary, *Extended formulations, nonnegative factorizations, and randomized communication protocols*, arXiv:1105.4127；[正文](https://arxiv.org/html/1105.4127)|PRIMARY_SECTION_READ，§1–3的定义与Theorem2|非负秩对应在期望意义计算矩阵的随机通信模型；明确精确的是哪种操作合同|不是本题统一归一化的潜在编码与二值响应；不能把RAC记忆下界直接重命名为普通非负秩或扩展复杂度下界|
|L10|Nayak, *Optimal lower bounds for quantum automata and random access codes*, arXiv:quant-ph/9904093v3 (1999)；[PDF](https://arxiv.org/pdf/quant-ph/9904093)|PRIMARY_SECTION_READ，Theorem2.3及§4.4；相应PDF页1、5截图|任意位可读出时的信息下界；本项目用其经典Fano版本|必须先构造在真实格点噪声下可靠的位查询。不能从精确线性秩直接跳到有噪声RAC；“状态数”与“比特数”也不同|
|L11|Barrett, *Information processing in generalized probabilistic theories*, arXiv:quant-ph/0508211，发表2007；[正文](https://arxiv.org/html/quant-ph/0508211)|PRIMARY_SECTION_READ，§IV.3–IV.4的方形状态及归一化结构|四种确定两测量结果、三维线性坐标的方形操作空间，是既有结构背景|本模型仍为经典site图。不能由方形/张量语言认领量子性、非局域性或物理上下文性|
|L12|Marzen–Crutchfield, *Circumventing the Curse of Dimensionality in Prediction: Causal Rate-Distortion for Infinite-Order Markov Processes*, arXiv:1412.2859；[正文](https://arxiv.org/html/1412.2859)|PRIMARY_SECTION_READ，§II–IV，Theorem1(Causal Information Bottleneck)|用未来相关信息和损失定义近似预测复杂度；有限窗口和无限记忆的区别|其平均信息失真/因果态设置不等于本题sup控制、sup准备的TV损失。先定准备分布，不能把典型性和最坏性混用|
|L13|Polyanskiy–Wu, *Application of information-percolation method to reconstruction problems on graphs*, arXiv:1806.04195v2 (2020)；[正文](https://arxiv.org/html/1806.04195)|PRIMARY_SECTION_READ，§2 Theorem1、Remark1以及§4的擦除比较|用通道收缩和图连接控制可提取信息，是本轮从代数rank转到实验可读性的直接方法启发|Theorem1是顶点iid标签、边XOR噪声等指定图观测；不等于环面同调输出。本轮擦除界由显式空列耦合自证|
|L14|Evans–Schulman, *Signal propagation and noisy circuits* (1999), DOI 10.1109/18.796377；[Caltech作者仓储](https://authors.library.caltech.edu/records/mcye7-j6k13)|ABSTRACT_ONLY，PDF下载失败。与1993会议版区分|可靠计算须有噪声、深度、大小预算；反馈/编码不能绕开信号衰减的模型假设|普通噪声门电路有明确门、扇出和布线；NN site连通图还没有实现这些元件。不能直接援引容错电路定理解决本查询|
|L15|Duminil-Copin–Tassion, *A new proof of the sharpness of the phase transition for Bernoulli percolation and the Ising model*, arXiv:1502.03050v3 (2018)；[正文](https://arxiv.org/html/1502.03050)|PRIMARY_SECTION_READ，§1.1 Theorem1.1、§1.2 site适配说明|有限程亚临界连接指数衰减是处理白matching故障路径的基础工具|正文主定义先是bond；site适配要明确。不能把square-bond pc=1/2移植到square-site；本文不求任何新pc。按arXiv版本日期，不把HTML动态排版日期当发表年代|
|L16|Damron–Lam, *Asymptotics for first passage percolation on logarithmic subgraphs of Z²*, arXiv:2502.18235v2 (2025)；[PDF](https://arxiv.org/pdf/2502.18235)|PRIMARY_SECTION_READ，导言与§2.1，Proposition2.5/Corollary2.7及其证明；PDF页0、2、12、14已截图|矩形中“机会数×指数连接代价”控制穿越，直接指向log尺度和罕见横向切断|对象为bond percolation与FPP；不是带指定旧端点的inhomogeneous site查询。可参考的具体命题不是直接可用的完整查询构造|

|L17|Johnson–Marijuán–Pisonero, *The INIEP: Irreducible and Positive Realizations*, arXiv:2605.23552v1 (2026-05-22)；[作者预印本摘要](https://arxiv.org/abs/2605.23552)|ABSTRACT_ONLY，作为2026更新筛查|不可约非负矩阵的谱实现是持续发展的邻接问题|其对象是给定谱/特征多项式的矩阵实现，不是同时保持准备、所有受控字串概率和归一化的随机实现；同名“positive realization”不能混用|

Heller1965、Vidyasagar2011的原文先例链在本轮未完成直接正文核查，不能把它们加入“已读定理”计数。其他检索命中的二手聚合页、无关临界指数表和量子扩展，没有用于科学结论。


## 2. 阅读真正改变了什么判断

### 正实现线：整理已有成果，而非重复证明一般现象

L01–L07显示，“普通阶数与最小正阶数不同”“用未来实验表示状态”“组合结构产生较大差距”都已有成熟背景。本仓库可主张的范围应限定为具体格点的可实现几何、全宽度条件和可复算证书，而不是发现这些一般原理。

L07尤其提醒：一个针对规范形的LP只能证该类的最小值。没有找到小模型，不能替代本仓库的正下界证明。新一轮无需重新跑大规模负权重优化。

### 信息线：先证明查询可靠，再使用记忆下界

L10给出下界方法，但只有物理查询错误率已被控制，它才有内容。L13把注意力转到通道可分辨性：本轮的共同擦除界同时给出TV、互信息、样本轮数和单状态近似。L09阻止把“统一归一化记忆”误写成普通非负秩。

### 几何线：固定噪声的困难是锚定门与远端输运同时可靠

L15–L16表明高密度长条的可靠性应考虑横向故障的指数代价与线性数量的位置。由此本轮先证明一般的空分隔列阻塞，再用两行矩形的最小故障集合完成平方根噪声律。这两个结果是明确的局部增量，而非新的渗流普适指数。

### 近似线：必须先选最坏或平均合同

L08与L12都不支持“精确rank大，所以有限数据下难学”的直接推断。本轮给出同一个准备族上精确rank指数大、固定高度固定噪声下一个状态却能统一近似的实例。若要研究自然样本而不是人工准备，另需先规定历史先验、允许源、h、ε、TV/信息失真和资源预算。

## 3. 建议的下一项实际数学工作

首选：不扩宽度自动机，构造并证明一个**带旧端点的固定噪声查询**，目标是存在固定ε₀>0，在h=C log w内将所有旧比特的最坏读错率压到1%。本轮已经证明h必须至少是对数级。

可试结构是扩宽的低密度隔离楔与外侧高密度导线。需要同时得到旧比特0的防旁路上界、旧比特1的端点接入下界、以及O(w e^(−ch))远端故障界。只证明普通矩形穿越，或只找到一次数值成功，都不完成该目标。

若上述构造成立，则标准RAC给固定噪声下指数初始记忆下界，而局部tensor关系仍把准备响应维数限定在3^t；这会消除当前ε=O(w^(-1/2))的限制。它依然是受控、最坏准备的定理，不自动成为自然临界采样或自主模型的结论。

备选不是“再数一个更大宽度”，而是证明任何这类有限密度锚定查询都受更强的噪声瓶颈，或把平均/典型合同单独立成问题。两条都需要具体的新增观察信息，不能通过更换术语获得。
