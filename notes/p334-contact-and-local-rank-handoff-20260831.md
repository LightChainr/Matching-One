# P334：接触结构调节出生时钟，局部二维响应已经出现

## 科学交接摘要

P334 已经从“即时 Euler/rank 看不见的控制是否影响未来”推进到两个具体问题：**哪些 prefix 结构组织出生中心的响应，以及这些结构能否预测新的条件续接读数。** 目前有三项互补交付：

- 联合出生协方差响应的约96–98%由 prefix 条件均值之间的关系承载，其中约85–90%位于固定 rank-cell 均值之间；固定同一 rank-cell 后仍留下可分辨的 prefix 响应差异。
- 四个预先声明的接触描述量捕获本几何、本源的层内有符号协方差 loading 的80–99%；共同安全标签质量在去掉两个基线时钟均值后仍携带中心响应信息。这提供了“接触结构调节时钟易感性”的具体线索。
- PR #509 的一次定向64-quartet增量已经支持**局部 A 条件均值响应存在二维成分**。它使用 `E_Z det J(Z)`，不再以 `det E_Z J(Z)`代替局部问题。E/出生间隔的局部行列式和四阶平方量仍未分辨。

**Held-out 预测尚未在本次固定提交中交付。** 执行分支的80–99%是描述性 loading 份额；新64局部秩结果也不等于接触模型已通过预测。下一份有价值的结果应消费现有数据，直接比较已认领的 `J=B G_score` 与固定接触特征预测，而非重做层级分解、首次局部秩或再次追加普通全域样本。

本交接固定执行分支 `43a30e49fdfcd38cdc9e085346c9e7eb49fa7650`，状态 `branch_only`；PR #509 固定 `8ad30617b0a3076a5c01a208eb213096d8879b32`，状态 `open_pr`。局部秩数值交付提交为 `7beb99ce86df903ec41ee2ec2a1de59184cf5b6b`。这些状态均不表示已经集成进 `main`。

## 1. 联合协方差变化主要是均值输运，rank-cell 身份仍不充分

这里使用不含额外 uniform-order 时钟的归一化出生秩
`X=K1/(N+1), Y=K2/(N+1)`，以及 `C=(X+Y)/2, W=Y−X`。
`Z`为完整原始 prefix，`U`为下一标签，`R`为其后的 suffix，
`G_rank=(rank_first,rank_second)`。先在各实际几何中计算，再做
`S=(first+second)/2`、`D=(first−second)/Δcos4`。

共同 next-label 政策不改变原 prefix 分布和即时 rank-cell 质量，
也不改变给定 `(Z,U)`后的 suffix law。其物理 score 为
`s_o(Z,u)=π_a[L_o(u)−mean_a L_o]`，在联合安全 degree 类 `a`内严格中心化，
其余标签为零。因此 `E_U[s_o|Z]=0`，但选择到不同标签后未来时钟可以变化。

完整分解已经交付：

```text
Cov(X,Y)
 = E_Z E_U Cov_R(X,Y|Z,U)
 + E_Z Cov_U(E[X|Z,U], E[Y|Z,U])
 + Cov_Z(E[X|Z], E[Y|Z]).
```

最后一项的源导数又拆为固定 rank-cell 内的 prefix 均值项，
以及不同 rank-cell 均值之间的输运。**rank-cell 质量不变，不意味着这些均值关系不变。**

主要数值如下，所有误差为原20批 delete-one SE：

- N325，plus→S：总协方差响应 `4.45778e−7±3.56657e−8`；固定 rank-cell 内的 prefix 项 `3.69856e−8±9.35027e−9`；固定完整 prefix 内的项 `9.32057e−9±6.40755e−9`。
- N425，plus→S：总响应 `4.21911e−7±3.39756e−8`；同 cell prefix 项 `3.19177e−8±8.33132e−9`；固定 prefix 内项 `1.23916e−8±8.22064e−9`。
- N325，minus→D：总响应 `−1.30906e−6±9.20886e−8`；同 cell prefix 项 `−1.33437e−7±2.75886e−8`；固定 prefix 内项 `−4.97139e−8±2.31287e−8`。
- N425，minus→D：总响应 `−9.73307e−7±5.82399e−8`；同 cell prefix 项 `−1.06430e−7±2.12545e−8`；固定 prefix 内项 `−3.89734e−8±1.88738e−8`。

同 cell prefix 项约3.8–5.0 SE，给出“条件均值响应只由 rank pair 决定”之外的证据。
完整 prefix 内的协方差导数较弱，minus→D 仅约2.1 SE；不能把强总体响应描述为同样强的固定 prefix 内 shape 变化。
分解中的 suffix 项是**选择具有不同剩余噪声的标签**，没有改变某个固定标签后的动力学。

此前弱的 lifetime 方差同样隐藏抵消。plus→S 的层内项为
`5.50527e−8±2.81101e−8` / `7.62338e−8±2.63454e−8`，
层间均值项为 `−6.30451e−8±2.01832e−8` / `−6.50410e−8±1.79515e−8`，
N325/N425 的总和仍弱。这不支持“每条路径的 lifetime 都保持不变”；也没有唯一确定每条路径如何移动。

精确来源：

- `44dc9e3396e39105cae85a29d04b39d0afc82d84:results/p334-birth-covariance-hierarchy/score.json`。
- `2bc3529468fbcba589182acaf98fa4855eb0a85e:results/p334-rankcell-covariance-transport/score.json`。
- [固定快照中的解释与全部公式](https://github.com/LightChainr/Matching-One/blob/43a30e49fdfcd38cdc9e085346c9e7eb49fa7650/notes/p334-birth-covariance-hierarchy.md)。

## 2. 接触结构已解释大部分有符号 loading，但尚未完成预测

四个固定特征是 joint-safe 标签质量、本物理源的精确 score energy、
safe contact degree 和 safe R0-loop count。安全标签的均值仍使用**原全部 vacant 标签分母**，
不是在 safe 子集中重新归一。两基线时钟均值为 `mu_C,mu_W`。
各 receiving geometry 只取自身 rank0 的原 prefix，在三个对应 rank-cell 内分别中心化，
保留原总体权重：有各 cell 的截距，但同一组斜率。

分析的主量为
`2 Cov(mu_C,H_C | within G_rank) − Cov(mu_W,H_W | within G_rank)/2`，
即本几何、本源的层内出生协方差响应 loading。按 N325 first、N325 second、
N425 first、N425 second 排列：

- source-energy 单特征捕获份额：`51.68±8.05`、`49.09±9.25`、`51.86±7.63`、`53.10±7.60`%，误差单位为百分点。
- 四接触特征捕获份额：`90.01±9.29`、`80.28±8.08`、`97.02±6.97`、`99.00±10.05`%。
- 剩余 loading：`8.9578e−9±9.2003e−9`、`1.6999e−8±7.7924e−9`、`2.4642e−9±5.9149e−9`、`7.6293e−10±7.7432e−9`。不能把四行都写成已闭合。

更有机制价值的独立问题是：接触结构是否只重述两个基线时钟均值？
移除 `mu_C,mu_W`的线性投影后，joint-safe mass 与本源中心响应的部分协方差仍为
`9.54781e−8±3.92033e−8`、`8.29040e−8±3.38796e−8`、
`9.01373e−8±3.72281e−8`、`1.07475e−7±2.36335e−8`，四行同号，约2.4、2.4、2.4、4.5 SE。
相应 lifetime 方向更弱。

这定位了**接触结构对中心易感性的额外线性信息**。它排除不了非线性的时钟解释，
没有识别唯一因果接触特征，也不是80–99%的响应方差被解释。
仅用时钟特征精确重现上述 clock-response loading 是投影代数身份，不能充当闭合发现；
去卷积后的响应方差估计仍有负值和较大误差，本交付没有可解释的响应 R²。

数值来源为
`9022659843ff0e9c2919c37e9468b0e7b5307268:results/p334-prefix-response-projection/score.json`；
读数实现为
`011f50e3efff5cffecd6171e497c9b879d5eb465:scripts/p334_prefix_response_projection.py`；
解释见
[`413da0e05c5c6a7a1f5feb843a24e8d4929f0d51:notes/p334-prefix-response-projection.md`](https://github.com/LightChainr/Matching-One/blob/413da0e05c5c6a7a1f5feb843a24e8d4929f0d51/notes/p334-prefix-response-projection.md)。
精确描述量来源是
`1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd:notes/p334-exact-prefix-structure.md`。

## 3. targeted64 支持局部二维出生中心响应；间隔与平方目标仍弱

PR #509 已完成对原 cell00 的1502/1551个 prefix 各新增64个 quartet，
共3053个既有 prefix、781568条新 suffix；未新增独立 prefix。
主读数合并旧8+新64为72组。`J(Z)`的行是两个实际几何的条件均值，
列是两物理 loop-score 源。所有局部目标仍用原每 N 20000个 prefix 分母，
非00 cell 的局部行列式贡献严格为零。

N325 / N425 的最终 paired20-batch 结果为：

- `E_Z det J_A(p_ref)(Z)`：`1.15937e−8±2.14282e−9` / `1.49767e−8±1.97301e−9`。
- `E_Z det J_integral A(Z)`：`5.32779e−10±7.33706e−11` / `5.26807e−10±6.77299e−11`。
- `E_Z det J_E(p_ref)(Z)`：`−3.79564e−10±5.97081e−10` / `3.65034e−10±6.62209e−10`。
- `E_Z det J_integral E(Z)`：`2.10345e−11±2.11663e−11` / `4.23907e−12±2.44714e−11`。
- `E_Z[(det J_A(p_ref)(Z))²]`：`7.32068e−15±5.05121e−15` / `4.05057e−15±5.61509e−15`。

新64单独的 A 局部行列式为 `1.24176e−8±2.43619e−9` /
`1.51007e−8±2.00470e−9`。这是新的条件随机流对照，不是独立总体复现。
其正的局部均值给出“每个 prefix 的这个2×2均值 Jacobian 都 rank≤1”之外的数值证据。
它不说明每个 prefix 都 rank2，不排除任意改变潜变量分布形状的标量模型，也不是连续场数目。

积分提供直接物理解释：在上述归一化约定下，
`integral A=1−2C`，`integral E=1−W`。
因此已分辨的是两个几何的出生中心响应，间隔响应尚未给出相同判断。
平方目标使用四个互异 quartet 的 U 统计，有限估计允许为负，不能截零；
其当前精度不足不撤销二阶局部行列式均值的正证据。

精确来源：

- [`7beb99ce86df903ec41ee2ec2a1de59184cf5b6b:experiments/p334-mechanism-response-20260831/REPORT.md`](https://github.com/LightChainr/Matching-One/blob/7beb99ce86df903ec41ee2ec2a1de59184cf5b6b/experiments/p334-mechanism-response-20260831/REPORT.md)。
- 同提交：`experiments/p334-mechanism-response-20260831/results-extension/score.json`、`results-extension/prefix_statistics_N325.npz`、`results-extension/prefix_statistics_N425.npz`、`analyze_extension.py`。
- `8ad30617b0a3076a5c01a208eb213096d8879b32:docs/NEXT-TARGETS.md`已把层级分解及首次局部检验移出下一步。

## 4. 协方差沿革：新坐标与新尾部均保留原 prefix 依赖

原 fork/contact block 为 `e32a8593`/`959a7fa2`，N325/N425各20个原批、
每批1000个 prefix，每 prefix 旧8个 quartet；同一 N 的两个几何保持配对。
精确 score 档案
`375cd3a12b2b7a87d79148a59f62b95898f9e471:results/p334-exact-score-quartet-moments/`
仅重读已有 suffix，并没有增加样本。

共同协方差沿革为：`b582015e`的16948坐标，扩展到
`e2ef9983f426890a299f5a6e1a2eba8b6d072855:results/p334-euler-dipole-connected-clock/`
的17866坐标，再由
`ce20158a5928e55b67324cba7ed3a18a5c163b39:results/p334-birth-covariance-hierarchy-joint/score.json`
与同目录 `N325.complete_common_factor.json.gz`、`N425.complete_common_factor.json.gz`
扩展至每 N **20754坐标**。后一包保留旧/新 estimator 的配对差、层级、rank-cell、响应份额；
这些因子不能作为独立证据相加。

旧 matched-mask 与完整 census score 针对同一物理导数，后者使用跨 quartet、
跨 prefix 的无偏乘积；最大旧新总响应差约1.95个配对 SE。
这属于估计器变化，不能称新物理效应。
接触投影另存每 N 726个原始坐标、358个派生坐标及20行 LOO/factor；
在固定执行快照中，追加到20754坐标共同因子的工作仍是交接项，不能虚称已集成。

targeted64 增加的是条件尾部信息。新随机键域与旧域分离，仍使用同一 prefix、counter、
原20批及两几何。`results-extension/score.json`把 old8/new64/combined72 保存为
每 N **20×1200**的共同 factor，协方差为 `factor.T @ factor`；
新64与合并72不是两个独立总体实验。局部秩包的生成4.164秒、聚合2.220秒均已实际完成，
不应再次登记为待运行。本笔记没有启动计算或服务器。

## 5. 下一份机制结果：从接触 loading 走向真正的预测

固定 `43a30e49…`的
[`notes/p334-prefix-response-projection-scientific-card.md`](https://github.com/LightChainr/Matching-One/blob/43a30e49fdfcd38cdc9e085346c9e7eb49fa7650/notes/p334-prefix-response-projection-scientific-card.md)
与 `notes/p334-projection-work-allocation-20260831.md`明确把 held-out `J=B G_score`、
固定接触特征预测和 cell00 conditional-shape 分析交给 PR #509 团队。
`8ad30617…`中它们也仍作为下一项，未出现完成报告或分数。
本交接的准确状态是**已认领、固定提交未交付**；不能据此断言现实任务尚未启动。

可立即承接的科学问题是：**两个源方向产生的中心易感性，能否由精确源 Gram
`G_score(Z)=E[s sᵀ|Z]`和固定接触结构预测，还是仍需要额外的条件时钟形状？**
这比继续报告一次正的总体 det 或80–99%描述份额更能区分机制。

- 先消费已认领的 held-out 预测，报告两源→两几何 J 的有符号分量及残差，
  不只报告平均 loading。对 `J=B G_score`，B必须在声明的训练部分固定；
  否则对逐 prefix 可逆 Gram 自由选择 B 会成为代数重写。
- 原8基线时钟与新64响应的跨流乘积可改善对共同 suffix 噪声的处理，
  现有 NPZ 已包含 `counter,batch,labels,old8,new64,combined72`及 J 分量。
  它仍共享同一批 prefix，单独做该乘积是描述性跨流读数；
  held-out prefix/batch 的预测才回答跨 prefix 泛化，二者应分别命名。
- 若固定特征能预测新读数，接触调节的易感性就从相关性进入可迁移的有限机制模型；
  若剩余有符号响应或条件 shape 系统偏离，下一靶点是该残差的具体 prefix 结构，
  而非重新宣布“没有机制”或重新做已完成的局部秩实验。

即时 Euler/rank 不变、局部均值 rank2、固定 prefix 协方差变化弱、
总体响应由均值输运主导可以同时成立。它们描述不同层次；
目前尚未给出唯一微观生成元、完整条件 law、普适性或原 norm-4 能量算符身份。
