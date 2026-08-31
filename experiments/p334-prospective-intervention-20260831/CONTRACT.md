# P334：新 prefix 上的 contact 干预预测合同

**状态：实现与旧输入 QA 完成，等待唯一 F0 commit；尚无本实验的新 prefix、census 或 tail。生产只在主任务提交本合同和代码并给出 full SHA 后开始。** 本轮从同一 archive 的描述性扩展，转到独立 prefix 群体的预定干预预测；不增加第五、第六个 descriptor。

## 1. 真正剩下的问题与两个主候选

contact 总交付固定于 `f225d2e25cbdf37d655a3a8b2c9515569db4dedb`。original00/new64 的四变量拟合和实际数值来自 `323de7d5ee4a980b3c77e1a972cb6c812a9f88e5`；原特征与 source 定义见 `1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd`。#334 最新 contact 讨论是 [5477178424](https://github.com/LightChainr/Matching-One/issues/334#issuecomment-5477178424)。这些是同原 prefix 的 descriptive signed loading，不是新群体预测。

只保留两个主候选，针对同一个 residual functional：

- **C0：这一个 clock-loading 残差投影消失**。新群体的 R 接近 0；可能发生抵消，不能称为四 contact response 闭合或充分状态。
- **C1：archive 剩余 loading 可传递**。新群体的剩余 loading 接近固定的旧值 R_old。

两者是可以被新干预杀掉的操作性限制，并非声称物理机制只有两种。原 80% loading 不自动提供两套接触理论，也不是 response-variance R²。固定旧点值（完整人口分母）为：

| N / receiver | R_old | 旧 20-batch SE |
| --- | ---: | ---: |
| 325 first | 4.976e−9 | 1.727e−9 |
| 325 second | 5.190e−9 | 1.851e−9 |
| 425 first | 4.651e−9 | 1.676e−9 |
| 425 second | 3.885e−9 | 1.351e−9 |

实际 scorer 从固定 score JSON 取未四舍五入数值。每 N 唯一主读数为 paired 两 receiver 的平均：N325 R_old=5.083315314795407e−9；N425 R_old=4.2680835963154025e−9。固定科学容差：C0 为 R_new/R_old ∈ [−0.25,0.25]；C1 为 [0.75,1.25]。不随新结果调整。两个尺寸的主区间各用97.5%双侧Student-t（60批、df59），Bonferroni共同覆盖率95%。置信区间完整落在一侧才称与该容差相容；横跨两侧或太宽即信息不足，落在两侧之外即两者均失败。**不把“未拒绝精确零”当作残差投影消失。** 四行及共同协方差完整报告但不按行挑赢家；两个尺寸不合成一个普适参数。

## 2. 预测器原样冻结，不看新 tail 调模型

对 receiver o，四变量严格保持：

`T_o(Z)=(joint_safe_mass, own_score_energy, own_safe_degree, own_safe_loop)`。

固定 323de7d5 的 **new64/original00** β、训练 μ_T、训练 μ_response：

`rhat_{o,s,F}(Z)=mu_response_old + beta_old^T (T_o(Z)-mu_T_old)`。

β 已在旧输出保存；均值由该输出的 `raw_batch_means` 和 `cell00.mass` 直接恢复，不再拟合。每 N 单独用本 N 模型，禁止 cross-size 重标定。训练时特征已见、原尾部已见；新预测阶段只允许见 prefix 本身、rankcell、exact census 与 T。新 baseline clock μ_C、μ_W 不作输入；它们仅在预测文件封存后用独立 baseline tail 组估计，供最终评分。

现成 β 数值预测 **C、W**，因此也预测 `X=K1/(N+1)=C−W/2`、`Y=K2/(N+1)=C+W/2` 与积分 `A=1−2C`、`E=1−W`。它没有训练 `A(p_ref),E(p_ref)` 的 β；只知道出生一阶矩不能推出固定 p 的 CDF。不能靠纯平移假设补这个缺口（旧 shape 分析本就限制该假设）。固定 p 的 A/E 仍保存为同实验的预定粗状态零响应控制，**不宣称已有其 contact 数值预测**。

## 3. 保持什么粗状态，改变什么 contact

新 prefix 按原完整采样协议产生，不拒绝/替换任何 prefix：N325 几何 `(325,57;0,1)` / `(325,18;0,1)`、k0=193；N425 几何 `(425,132;0,1)` / `(425,268;0,1)`、k0=252。前缀是均匀排列的前 k0 个位置，两取向共用 label 顺序；prefix RNG 起始 key 域与旧域分离并记录，但不把不同 PRNG key 说成数学独立性证明。

只在双 R0 上施加干预，其他八 cell 响应置零，**分母仍是全部新 prefix 数**。记录九 cell 计数以核验新群体，而不是抽够若干00就停止。取 d=N−k0。joint-safe label 按 `a=(e_first,e_second)` 分类，n_a 为完整 census count，π_a=n_a/d。物理源为

`H_s(u|Z)=pi_a*(L_s(u)-mean_a L_s)`，`L_s=e_s−c_s`，类外 H=0。

沿单一物理 source s 定义真正的正概率干预：

`q_±(u|Z)=(1±epsilon H_s(u|Z))/d`。

本几何每个 vertex primal degree=4，0≤L≤3，因此 |H|≤3。固定 epsilon=1/8，密度比至少5/8；census 逐 label 核验界限，不在新结果后改 epsilon。每类 ΣH=0，故 q± 正、归一化且每类总概率 **精确** 等于 π_a，类外概率仍1/d。每个 source 的正 arm 都提高相应 loop mark 的均值，除非其类内方差恰为零。

两 arm 使用同一完整 prefix，按同 a 配对 label，joint-safe 使添加后两 rank 逐对相同，e_first/e_second 逐对相同，k=k0+1 逐对相同，两图 V−E 的 Euler 值逐对相同。因此明确保持完整向量 `(k,rank_first,rank_second,Euler_first,Euler_second,e_first,e_second)` 以及所有原 prefix 函数。**不声称保持组件数、添加后的 contact census、未来 clock law 或任意尚未声明的 coarse summary。** 若“coarse state”另含这些量，必须先补可行性证明，不能改名蒙混。

## 4. 有限效应严格等于已训练的一阶 source，免去 Taylor 外推

对任何完整未来观察量 F，有精确恒等式

`E_q+[F|Z]−E_q−[F|Z]=2 epsilon E_uniform[H_s F|Z]`。

这是仿射密度恒等式，对有限 epsilon 成立；不使用 exp 源的二阶或高阶截断。定义干预单位响应 `tau_F=Delta_F/(2epsilon)` 后，训练 response 与新实验的 estimand 相同。

推荐执行 **同 class maximal-coupling 的 Rao–Blackwell 化 contrast**：令整数 `hnum(u)=n_a L_s(u)−sum_a L_s`，则 H=hnum/d；正负整数质量在每类相等。令 `W_s=sum_u max(hnum,0)/d²=E_uniform[H_+]`，抽 class 权重为其正整数质量，类内正 arm 按 hnum_+，负 arm 按 (−hnum)_+。两 arm 用共同的均匀剩余-label排列（各删除已激活 label），再分别运行两取向 rank engine。此时

`tau_F=W_s E_contrast[F_+−F_-]`。

共同概率质量的两 arm label 与 suffix 可完全相同，其差严格为0，故可解析消去；不是稀有非零样本的结果后筛选。保存 W 和每个 contrast。若 W=0，该 source 响应精确0，不重抽 prefix。要如实称为正政策耦合的 contrast 抽样；不能把所有生成尾部说成直接 iid 来自完整 q±，也不复用旧 uniform 尾部做新权重。估计 Δ 时乘2epsilon。

## 5. 直接评分剩余20%，不另加 descriptor

每个新00 prefix、每 receiver/source，独立 uniform baseline 组估计 μ_C、μ_W；另一组干预 contrast 估计 τ_C、τ_W。两组起始 key 域不重叠。已封存预测仅由 T 得到。主 own-source 量为

`R_new = pi00 * [2 Cov_00(mu_C, tau_C−rhat_C) − 0.5 Cov_00(mu_W, tau_W−rhat_W)]`。

baseline 与干预组条件独立，故它们的交叉乘积无同-tail偏差；同一个 baseline 组可共同评分所有 physical source/receiver，并保留其相关性。有限 prefix 协方差用 m00−1，外乘实际 m00/M；不能将00重新归一化成人口。预测器 T 精确、没有 tail 测量误差；旧 β 固定，不在新批拟合或改截距。

新的原批分组删除每次重算00质量、均值和 covariance。**主检验条件于冻结 β、训练均值和旧点预测，old 数据不进入验证得分。** 另保留旧20个 β/μ/R 的 LOO，但不据此增加主任务；若后续呈现无条件新旧差，必须联合计算 `D_loo=R_new(beta_old,loo,mu_old,loo)−R_old,loo`，不得只加旧 R 方差。旧候选±25%带针对“固定 archive 数值可否预测”；不称为未知旧总体 R 的精确带。

预定粗状态控制：若未来条件 label 均值只依赖上述被保持的 coarse vector（可附原 prefix 的四 contact 量），则所有 ΔA、ΔE、ΔK1、ΔK2 精确为0。这个限制比“first Jacobian 可由 contact 预测”强；后者本身不蕴含零。零响应只作控制，不升为第三个主实验。

## 6. 有界信息预算与一次冻结后的自动分片

固定 **300,000 个新 prefix/尺寸**，每尺寸60批、每批5,000；60个 shard 可均分到2–4台机器。每新00 prefix 固定32条独立 uniform baseline 尾部，以及 **每物理 source 64个正负 coupled contrasts**（两source合计256条 arm 尾部，两取向在同一尾部共同求值）。其他 cell 只保存 rank/census身份并按零响应纳入分母。预算不依 p 值、模型哪边胜或00计数而追加。

旧8quartet baseline / New64 tangent 精度下，SE按 sqrt(20000/300000) 缩放约为旧值的0.258，四行约0.35–0.48e−9。这只是设计参考：新 coupling 的条件方差不同，不能称为已证明power；旧训练误差也不会随新样本消失。拟定最小可分辨差为0.5 R_old，容差宽0.25 R_old；若最终区间分不开两候选，就报告预算内未分辨并停止。不得根据新 tail 方差追加预算。

封存顺序：

1. **唯一 F0**：代码、训练参数JSON/哈希、两个候选及容差、几何/clock定义、N/M/Q、所有prefix/label/tail key域、shard集合、统计规则一起提交。F0前连新prefix+census也不生成。
2. 每shard先生成完整5,000个prefix与00的exact census，保存所有预测、输入hash和时间；不调用未来rank paths，也不允许按新features调整已冻结规则。
3. 该shard预测文件落盘并记录SHA后，自动进入固定baseline/干预阶段，无需第二次人工审批或全局commit。独立baseline用于评分，不回流预测。
4. 回收60×2完整shards、校验预算/哈希，执行唯一scorer，报告固定预算下的相容/排除/未分辨，再结束对应机器任务。

具体counter方案：N325原prefix seed `20260831430325`，新counter `[53032500000,53032800000)`；N425 seed `20260831430425`，新counter `[53042500000,53042800000)`。旧counter分别从43032500000/43042500000开始的20,000个，故本任务两新块与这两旧块不相交。每prefix index=shard×5000+offset；不拒绝偶然重复的配置。原fork seed保持`202608311638334`，新条件流地址置bit63，结构为`bit63 | N<<48 | index<<20 | kind<<16 | draw<<2`。baseline kind=1、两source label kind=4/5、两source共同剩余排列kind=8/9，各槽固定。旧8/旧prefix-new64的地址在bit41以下，因此在同一bijective key映射下起始key命名域不同；不宣称随机独立性的数学证明。

静态分片：HZsCM6承担每尺寸shards0–29，TV2N0X承担30–59；各14workers、BLAS线程1。两尺寸共120个shards，调度顺序不改变样本/预算。编译命令：

```sh
g++ -std=c++17 -O3 src/prospective.cpp -lz -o prospective
```

F0给定后，HZ示例（TV替换为30:60）：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run_shards.py --freeze-commit F0_FULL_SHA --shards 0:30 --workers 14
```

`run_shards.py`校验FREEZE.json里的文件hash，逐shard运行prefix/census→保存预测→tails→sufficient统计。`score_prospective.py --input production --output results`只在全120shards取回后运行一次。拒绝覆盖已有产物，失败必须查明阶段/随机域，不能静默重抽。机器二进制不提交，提交源码、训练快照、参数、合同和QA。当前没有独立验证结果。

## 7. 最新 source-normal 结果的边界

`20a95bf253752a02a5bd5b4ca1153494488f17dd` 已区分“class间不同斜率”和“class内高阶mark结构”，尚不能在两者间裁决。其 source-normal positive center response 来自 `2c3a5ca2`，不是这里四contact loading的20%残余。成功预测 J、未零 normal response 和未零 R 可以同时成立。本实验不给该分解新增坐标、不重跑Hessian，也不把 R 的可传递性解释为其中某个唯一微观机制；下一座理论桥仍是把被实测拒绝/保留的 residual constraint 对应到明确的 contact 动力学限制。
