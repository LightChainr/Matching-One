# Mechanism decisions / 独立实验与精确传输判决

2026-08-31 · #334两项不同的独立实验、#154的165M固定新块和次级时钟线均已完成。
优先级只分配注意力，不锁Issue、不合并PR、不把一次参数化失败扩展为整条
研究路线失败。[最终科学交接](../notes/independent-decisions-final-20260831.md)
给出定义、数据独立性和精确来源。

## 最新已完成：Q4 有限归一化传输接口，下一判别只剩指定 J22

[041257b8 的解析推导](../notes/four-leg-trace-denominator-interface.md)
在固定m=2、同一空间恒等颜色缝给出`Z22=(T+3D2−4C3)/6`，精确投影
完整S4的`[22]`通道。其逐构型支撑仅在rank1，所以原q/E未归一化分子
严格为零；但原分母中的份额为`z22=2Z22/(T+R)`，通过完整四项root/slope
公式传入`J22=d_x U`。x是固定分量归因，不是拟合的新颜色数或局部正源。

轴L≥5时`Z22(a)=L(L−3)4^(2L)a^(2L)+...`，由两个不相邻竖直列精确
计数得到；它证明trace分量非零，**不证明matching根上的J22非零**。
下一计算固定为N25原方向对、m=2的J22：严格非零界将排除这个有限通道的
complete normalizer-neutrality；零则保留该有限null。需要三缝p-jets或
rank1的`(K,g,n mod6,c)`信息，旧`(K,g,q)`档案不能替代。
定义见[机器接口](../analysis/four_leg_trace_interface.json)。本轮没有新枚举、
随机生产或测试。Q4颜色表示不等于唯一连续四腿场，也不自动给出Q→1延拓。

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
现已建立；J22仍待计算，Q→1及有限confluent延拓仍需分别定义。
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
