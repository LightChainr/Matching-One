# Independent decisions / 独立实验判决与模型退出记录

2026-08-31 · #334两项不同的独立实验、#154的165M固定新块和次级时钟线均已完成。
优先级只分配注意力，不锁Issue、不合并PR、不把一次参数化失败扩展为整条
研究路线失败。[最终科学交接](../notes/independent-decisions-final-20260831.md)
给出定义、数据独立性和精确来源。

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
现在两套评分均已完成；原规则确实给出了有限尺度排除。下一主实验需由具名
微观模型给出实际source激发的contrast值/导数和原U预报，不能从本次弱结果
或contact约一半反向拟合。未形成此类预报时不自动占用生产队列；探索仍开放。

执行与结果交接只留仓库。十机均可按更新Huawei Skill直接使用，仍不覆盖
他人输出或操作不明进程。本总览仅做指定目录只读检查，未开关机或追加采样。
[此前决策稿 `f670d26e`](https://github.com/LightChainr/Matching-One/blob/f670d26e8056116ec1787e0a2f29009b5db24a7d/docs/DECISION-EXPERIMENTS.md)
保留未冻结候选、旧预算与历史映射，不再充当当前合同。PR #267保持Draft。
