# 全孔连通核：哪些变量能积分，哪些信息仍须保留

本页补足从已完成的全 epsilon 面核到总体计算的两个具体步骤；不另开描述量搜索。**非循环 A 可以严格积分，原 U 的符号仍依赖外部同调连通。** 齐次 N50 的实际计算与资源记录归入[唯一计算包](../experiments/p337-homogeneous-n50-20260831/CONTRACT.md)。

## 精确积分全部非循环开关

固定 B 占据配置，以全部 A 潜在存在建立原 NN 图。令 A_c 为处于至少一个图循环上的 A，L=25-|A_c|。a 不在 A_c 当且仅当其全部关联边都是桥；无论其它 A 如何关闭，添加 a 都不会增加 beta 或完整环境同调像。于是 q/E 不变，同一源 S*=2 beta-r-3K+2N+1 每添加此 A 只减3。

注意：判据是**图循环支撑**，不是2-core，也不是只保留非收缩循环。两个循环之间的桥链会留在2-core；收缩循环虽不改变 rank，却改变 beta，不能删掉。

加占据数记账场 u，alpha=1-epsilon*(1-p)，w=exp(u-3t)，F=1-alpha+alpha*w。所有非循环开关的精确 message 是 F^L。给定保留状态，设 rho=alpha*w/F，则

```text
K = K_c + J,              S = S_c - 3J,
J | retained ~ Bin(L,rho),
E K = K_c+Lrho,            E S = S_c-3Lrho,
Cov(K,S | retained) = -3Lrho(1-rho).
```

L 依赖 B，不能在 B 平均中当常数消去。源响应需保留 F 的导数；例如

```text
partial_p^j F^L = (L)_j epsilon^j (w-1)^j F^(L-j),
partial_p partial_t log(F^L) = -3L epsilon w/F^2.
```

一般 epsilon 下，u 只是 K 记账场，不能偷换原 p 方向。epsilon=1 才能统一以总 K 的 Bernoulli score 表示热导数。

## 初始潜在图允许精确的块卷积

在本题 honest 方格图上，一个循环 A 不会同时属于多个循环双连通块。相邻两个占据 B 端口可经另一个潜在 A 角点绕过中心 A；m>=3 的端口连通，m=2 相邻也连通。m=2 对向若无外部绕路则 A 本来就在桥上，否则不是 articulation。因此初始潜在循环块中的 A 集合互不重叠，共享 articulation 只能是固定 B。

beta 按块相加，B 占据数只在全局计一次；每块随机 A 可分别积分。但各块同调像必须合并后只取一次 rank，不能直接把各块 rank 相加。令

```text
Z_b(H;u,t) = sum_(z_b:H_b=H) Pr_alpha(z_b)
                * exp[(u-3t) K_Ab + 2t beta_b],
Z_x(u,t) = exp[u K_B+t(2N+1-3K_B)] F^L
             * sum_(H_b) exp[-t rank(sum H_b)] product_b Z_b(H_b;u,t).
```

收缩潜在块的 H=0，可先标量积分。只求当前 q/E/S 时保留 H 的有理 span 足够；只保留单块 rank 不够。若先删孔、再按动态图分块，同一个 A 可能变成两个循环块的共享开关，不能独立复制抽取。

指定真实 N50 `(5,5)` 图，固定 B 在 `(0,-1),(1,0),(1,2),(2,3)`，25 个 A 中只有4个循环 A；桥 `(1,1)` 留在2-core但可积分。两收缩四环各有两个 A，故整个2^25条件 A population 直接得到闭式

```text
B(u,t)=(1-alpha)^2+2alpha(1-alpha)exp(u-3t)+alpha^2 exp(2u-4t),
Z_x(u,t)=exp(4u+89t) F^21 B(u,t)^2.
```

这是一个固定 B 的准确例子，不宣称典型 B 也能同样压缩。两固定 B、七指定真实构型、十六对抽象桥状态及一个有理 message 点已核验，见计算包中的 `cycle-check.py/json`；没有 B 扫描或随机采样。

## 既有 Russo 结构能给出的最强局部比界

此处复用仓库既有 Russo/birth identity，新增的是它对全孔原 U 的限制。齐次 t=0，令 a,b,c 分别为全部站点环境的 rank `0→1,1→2,0→2` pivotal 概率之和（不含 p(1-p)因子）。则

```text
q_p=a+b+2c,          E_p=b-a,
|E_p| <= q_p-2c,
|U|/A_N <= (2/|Delta|) [1-2(c1+c2)/(q1_p+q2_p)].
```

该界不会单独证明 H4 小量或确定取向符号。实际4×4周期方格给出尖锐反例：在 `(0,2)` 开关的上下端口占据、左右端口为空时，只占据其余竖列产生 rank0→1；另加完整横行产生 rank1→2。两者局部端口、外部端口连通、Delta beta=1、Delta S=-2 完全相同，但 E_p/q_p 分别为−1与+1。两整行交叉去掉交点还实现0→2，故不能删去 c。

源响应 J_O=Cov(S,O) 不是非负 Russo 和。固定揭示顺序，令 d_i O 为给定过去、对未来平均的开关差分，则

```text
J_O=p(1-p) sum_i E[d_i S * d_i O].
```

不能替换成同一个完整环境中未经条件平均的差分乘积。度4与原源恒等式给 |Delta S|<=3，因而

```text
|Jq| <= 3p(1-p)q_p,
|JE| <= 3p(1-p)(q_p-2c).
```

这些界仍不提供原 U 源响应中的 mixed 热导数、共同根和分母变化。缺失的是带外部同调的 pivotal 分配及其与条件 source score 的联合权重。现有 birth archive 已可恢复的量不重新立项；没有新连通比较定理时，本通用界路线到此停止。

## 用充分状态直接计算总体

另一条有限图途径已经实现：[商空间状态证明](../experiments/p337-homogeneous-n50-20260831/STATE-PROOF.md)只保留活跃黑 NN 连通性、全局已知同调 H 与路径位移在 R²/H 中的像。rank1时保留规范方向 h 和 det(h,d)，rank2后丢弃位移。初 S=2N+1；占点−3，冗余边+2，新增长全局 rank 减相应量。count/sumS 按 K 精确合并，终点由 q=r−1 恢复原六矩。

这是原 q/E/S 所需的充分状态，不支持事后恢复完整整数 twist、每组件 winding 或有限 t 分布。两种几何的 N25 完整表已与既有独立2^25枚举逐整数相同。动态规划合并部分配置的连通等价类有既有方法先例，参见 [Akhunzhanov 等的环面 site-percolation 精确计算](https://arxiv.org/html/2204.01517)；其单方向 wrapping polynomial 不等于本项目的 rank/S 输出，也不替代本次资源实测。
