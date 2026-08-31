# P334：保持即时Euler/rank分布的有限强度空间扰动

**有限t=±1的两种共同标签策略均已计算。** 既有一阶响应延续为可测的有限未来birth响应，同时两取向即时rank及Euler增量的联合分布按定义完全保持。此次采用旧标签与旧续接的配对importance估计，没有按新策略直接重采样，也没有增加独立prefix。

主读数单位为原A=F1+F2−1；S为两取向平均，D为原Δcos4归一化方向差。误差是原20个配对batch估计的1SE，多个视图共享数据。

| N | 源→读出 | t=−1响应 | t=+1响应 | 同估计器精确t=0导数 |
|---:|---|---:|---:|---:|
|325|g_plus→S(A)|5.47057592e-05 ± 4.32e-06|-5.42530839e-05 ± 4.23e-06|-5.44956191e-05 ± 4.27e-06|
|325|g_minus→D(A)|-0.000142572163 ± 8.28e-06|0.00014259443 ± 8.4e-06|0.000142632686 ± 8.34e-06|
|425|g_plus→S(A)|5.9366526e-05 ± 4.93e-06|-5.86576298e-05 ± 4.84e-06|-5.902831e-05 ± 4.89e-06|
|425|g_minus→D(A)|-0.000132256641 ± 9.45e-06|0.000131577085 ± 9.7e-06|0.000131948945 ± 9.57e-06|

两尺寸g_plus→S(A)为负、g_minus→D(A)为正，反向t改变主要响应符号。g_plus→D仍弱，不能由这次有限扰动宣称新的交叉方向已分辨。

## 这次实际新增了什么

从原N325/N425各20000个prefix精确枚举所有空位，补齐同时保持两取向rank的joint-degree类 a=(e_first,e_second) 内，(L_first,L_second)的完整整数计数。N325枚举2640000空位，输出470087个非零类/loop计数行；N425枚举3440000空位，输出507735行。没有回放后续序列或调用DP。

令d为空位数，π_a=|A_a|/d，L_o=1{old_rank_o=0}(e_o−c_o)，g±=(L_first±L_second)/2。沿用此前已定义的策略：

```text
q_t(u|Z) = π_a exp(t π_a g(u)) / Σ(v∈A_a) exp(t π_a g(v))
outside the joint-safe classes: q_t(u|Z)=1/d
w_t(u)=d q_t(u|Z)
H_g(u)=w′_0(u)=π_a(g(u)−mean_a g), outside=0
```

每类的概率仍为π_a，故任何有限t均保持两即时rank与两个Euler增量1−e的联合分布。完整census使归一化无需从16个已采标签估计。读取原每prefix8个独立quartet，各标签两条续接取平均Ȳ；主估计为

```text
ΔY(t|Z) = E[(w_t(U)−w_t(V))(Ȳ_U−Ȳ_V)/2 | Z]
Ydot(0|Z) = E[(H_g(U)−H_g(V))(Ȳ_U−Ȳ_V)/2 | Z].
```

因为U,V条件独立且E_uniform[w_t|Z]=1，半差期望正是E_qt[Y]−E_uniform[Y]。允许U,V来自不同类；与旧同类mask一阶估计器期望相同，单次样本值无需完全相同。未来CDF用Binomial生存函数Pr(Binomial(N,p_ref)≥K)，p_ref=0.59274605079；p积分用1−K/(N+1)。F1/F2、A/E、两取向的续接与标签配对保持不变。

每个batch的分母始终为1000prefix×8quartet×2半差，cell00等分层保留原总体分母。两尺寸各20batch，存下全部768维batch向量及covariance factor，N之间不伪造配对。两条续接没有增加next-label样本数。

## 有限偏离可测但很小

| N | 通道 | [Δ(+1)−Δ(−1)]/2 − Ydot(0) | 有限偶项[Δ(+1)+Δ(−1)]/2 |
|---:|---|---:|---:|
|325|plus→S(A)|1.61976e-08 ± 2.31e-09|2.26338e-07 ± 9.42e-08|
|325|minus→D(A)|-4.93899e-08 ± 4.91e-09|1.11339e-08 ± 1.97e-07|
|425|plus→S(A)|1.62321e-08 ± 3.14e-09|3.54448e-07 ± 1.01e-07|
|425|minus→D(A)|-3.20824e-08 ± 5.82e-09|-3.39778e-07 ± 2.04e-07|

若只保留t的一阶项，主通道的有限响应已有良好数值近似；上表依赖的微小差值可有较小SE，并不意味着大幅非线性。t乘以各类π_a，t=1在当前分布上是温和扰动，不能外推到任意大t。

所有实际标签census键匹配成功；每类归一化余量≤1.43×10⁻¹⁴个标签计数单位。全部有限权重范围约0.7424–1.3307。理论uniform-proposal的有效样本比例1/E[w²|Z]最小0.99596；这只是importance权重稳定性，不是prefix独立性或实际SE的替代。

## 运行、复现与下一步

Huawei TV2N0X实际运行：C++ census分别2.06161/2.76299秒，Python分析6.38408秒；包安装与传输另计。ARM64、GCC10.3.1、Python3.9.9、NumPy1.26.4、SciPy1.13.1。主程序、census与结果的远端/本地SHA256均一致，见RUN_RECEIPT.json。

```bash
g++ -std=c++17 -O3 src/finite_census.cpp -lz -o finite_census
./finite_census 325 census/N325
./finite_census 425 census/N425
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analyze_finite_source.py --output output
```

现有census/output目录受防覆盖检查保护；重现时在新副本中给新的输出目录，或先保留旧结果。全部冻结源码/旧输入的commit和SHA在SOURCES.json。EXPERIMENT.md说明本轮是看过一阶结果后的分析设计，不是预注册验证。

这批完整census同时交给独立子代理，用于同prefix局部rank检验的exact-score估计；它不再需要同类label对的稀疏筛选。下一步应先消费这个检验的结果，再决定是否在支持两独立score的prefix定向增加续接。无需再扩大全部普通尾部。

统一边界：有限两尺寸、原随机prefix群体、同20batch依赖；已有数据的精确importance估计与直接干预采样分开；有限Euler不可见空间响应不等于连续场身份，也不证明每个prefix有两个响应方向。
