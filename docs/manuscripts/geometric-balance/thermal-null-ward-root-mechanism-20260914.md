# 热零向量、E4 形状律与 Matching-root：一个可失败的机制

2026-09-14。服务 #768/#802。主要读取快照为 #771 `8e1282f6d4397f03c80a2883b917c2d31a4b93a3`；写入前分支已继续推进，不覆盖其他文件。本篇不再增加一套宽度阶梯。

**核心新目标**：把“thermal spin4”从尺寸幂和旋转标签，推进到可计算的一点响应与整个环面形状函数。以下 Virasoro/Ward 计算是明确条件下的推导；实际 square-site 根的模函数规律仍是猜想。标准 Ward/Zhu 技术和热通道的 level4 选择已有文献，不重新计作通用方法发现。

## 1. 热候选的具体代表，而非 spin/momentum 排除

取 c=0，热主场 epsilon 的 h=bar h=5/8。局部平坦坐标中的二级奇异向量为

    chi=(L_-2-(2/3)L_-1^2)epsilon.

当前 `verify-768-momentum-exclusion-20260914.md` §2.1 仍把局部场当成满足 `[P,phi]=(2pi/L)(h-hbar)phi` 的平移动量本征算符，遗漏了插入位置导数。[GL](2.3)明确给出 `[L0,V(v,z)]=(z partial_z+h)V(v,z)`。空间零模不因局部 spin 非零就自动消失。这里直接计算热模块，不只用 identity-family I3 作反例。

定义

    Q epsilon=[L_-4+(20/11)L_-3 L_-1-(160/561)L_-1^4]epsilon,
    U4=-(11/29)[Q epsilon-(80/33)L_-2 chi].

有理 PBW 计算给出 L1 chi=L2 chi=0，L1 Q epsilon=(80/11)L_-1 chi，L1 U4=0。因此 U4 是全 Verma 空间里的准初级代表。在商掉 null 后裔与总导数后，它的类恰为 L_-4 epsilon；level4 的五维空间中，被商子空间秩4，L_-4 不在其中。

更有用的是不商 null 时的精确关系：

    U4=L_-4 epsilon+(80/87)L_-2 chi+L_-1 V3.       (1)

V3 是一个三级态。故“null 是否影响本观察”可用一个具体振幅检验，而不是先替整个 c=0 理论贴模块标签。

## 2. 圆柱零模确实允许非零

若 `<m|epsilon(z)|m>=C_m z^-h`，平面 Ward 给

    <L_-2 epsilon>/C_m=m+h,
    <partial^2 epsilon>/C_m=h(h+1).

在 C_m 非零、该三点通道 null 解耦时，m+h-(2/3)h(h+1)=0，得到 m=5/96。这与热—磁性融合相容；该融合选择已有 [JPS] §4.2。

周长2pi圆柱、z=exp(t)中，正规化应力张量插入是

    <T_cyl(t)epsilon_cyl(0)>/C_m
       =m-c/24+h exp(t)/(exp(t)-1)^2.

Laurent 展开 `exp(t)/(exp(t)-1)^2=t^-2-1/12+t^2/240+...` 因而给

    <L_-2 epsilon>=(m-c/24-h/12)<epsilon>=0,
    <L_-4 epsilon>=(h/240)<epsilon>=<epsilon>/384.    (2)

null 后裔解耦时 `<U4>=<epsilon>/384`。恢复周长 L，实的两手征和系数为 pi^4/(12 L^4)。这并不保证实际格点 C_m 非零；它已经排除了“平移普遍杀掉该 spin4 零模”的错误推理。真空、磁性外态与 rank 投影 trace 仍须区分。

## 3. Torus Ward 给出完整 E4 因子

取复周期 u,v，tau=v/u，Im(tau)>0。设一个线性插入泛函 F 对主场 epsilon 满足普通周期 Ward 恒等式、插入位置平移不变、总导数期望为零，且无额外 seam/contact 项。普通 CFT trace 满足相应标准关系；**rank/同调投影是否满足，是实际模型的待证接口**。

定义

    G4(u,v)=sum_(lambda in Zu+Zv,lambda!=0) lambda^-4
           =(pi^4/45)u^-4 E4(tau),
    E4(tau)=1+240 sum_(n>=1) sigma_3(n) exp(2pi i n tau).

Weierstrass 函数展开为 `wp(z)=z^-2+3G4 z^2+...`。单主场 Ward 恒等式使 F(T(z)epsilon) 的位置依赖为 h wp(z)F(epsilon) 加一个与 z 无关的项。取 z^2 系数：

    F(L_-4 epsilon)=3hG4 F(epsilon).                (3)

无需除以可能为0的 F(epsilon)。同一式可从 [GL](2.10) 的 Zhu 递推得到。如果 F(L_-2 chi)=0，由(1)

    F(U4)=(pi^4/24)u^-4 E4(tau)F(epsilon),
    F(U4+bar U4)=(pi^4/12)Re[u^-4 E4(tau)]F(epsilon).  (4)

**全模块 null 解耦只是充分条件，不是必要条件：目标一点振幅 F(L_-2 chi)=0 就够。** 若其不为0，精确差额是

    F(U4)-3hG4 F(epsilon)=(80/87)F(L_-2 chi).        (5)

对于普通未正规化手征 trace f、u=1，右侧可写为 `(80/87)(2pi i)^4 D_(h+2)D_h f`，其中 D_h=q_mod partial_q_mod-h E2/12。广义对数 trace 或投影接缝异常应另列，不能机械复制普通主场式。

当 null 在 trace 中解耦，D_h f=0，故 f=C eta_D(tau)^(5/4)，其首项为

    q_mod^(5/96)[1-(5/4)q_mod-(35/32)q_mod^2+(45/128)q_mod^3+...].

这里 eta_D 是 Dedekind eta，q_mod 不是占据概率。本轮精确核对了乘积展开与微分递推。[JPS] 已讨论热通道二级后裔消失、四级后裔先贡献；它同时指出 Q=1 的无条件能量一点函数为0。**不能把该零函数当成 rank 条件热响应，更不能凭它取消本任务的源映射。**

## 4. 形状、方向和周期基是同一个协变量

周期换基 `u'=u(c tau+d), tau'=(a tau+b)/(c tau+d)` 下，E4 的权重4变换给

    u'^-4 E4(tau')=u^-4 E4(tau).

若 u=ell exp(i theta)，几何因子为

    ell^-4 Re[exp(-4i theta)E4(tau)].              (6)

一般斜模参数中不能随意更换相位符号。固定面积 N=|u|^2 Im(tau) 时，N^2 倍因子为 `(Im tau)^2 Re[exp(-4i theta)E4(tau)]`。

它包含：长圆柱 E4(i infinity)=1；方形 E4(i)=1.45576289226870932246242200359886929...；长宽比2矩形 E4(2i)=1.00083698843473765919291512747422264...；六角 rho=exp(i pi/3) 上 E4(rho)=0，且 E4'(rho)=-(2pi i/3)E6(rho)不为0。

既有 double-null 方案已经有角向/六角零点。本篇增加的是零点之间的完整形状律，而非再增加一个零检验。

## 5. Actual-site 猜想：一个系数控制整个根形状

所需条件明确为：

H1：实际 rank 扇区具有上述临界 Ward 接口，使用同一个热场正规化；
H2：leading matching-odd correction 是同一微观系数乘 U4+bar U4，且目标 null 异常消失，其他同阶 scalar/module/seam 项已排除或分开；
H3：根局部非退化，余项对固定模参数紧集一致，并另外控制 finite-aspect 到 cylinder 的对应。

这时每个 Z_r 的 first-order insertion 都满足 `partial_g4 Z_r=C4(u,v)partial_t Z_r`，总 Z 同样满足，故正规化 P_r 也同样满足。将热度量、符号和微观系数收进一个未知常数 A，预言

    p*_Lambda-pc=A Re[u^-4 E4(tau)]+o(|u|^-4).     (7)

**(7)是 square-site 机制猜想，不是(3)证明的格点定理。** A未计算，每个形状不能另拟合一个A。一个后裔可在临界一点响应上等价于热移动，却不等于全理论的冗余算符；不能对临界恒等式直接再微分，忽略二次热插入与接触项。

因此特别预言，同一周长L的方形环面根与无穷长轴向圆柱 charge 根，若共同A非零，则

    (p*_square(L)-pc)/(p*_cylinder(L)-pc) -> E4(i). (8)

这把根指数、角向spin4、六角零点和此前的切向盲点，连接成同一个可失败机制。

## 6. 复用已返回数据，不生产新宽度

团队L5/L6整数rank-by-K表，配仓库已有同宽圆柱根，给：

|L|方形环面根|已有圆柱根|用pc_ref形成的位移比|相对E4(i)偏差|
|---|---|---|---|---|
|5|0.591988256518333844610968680211928879|0.5922358232050258|1.48520835204458070455|+2.0226824%|
|6|0.592395070817704237693855807642505434|0.592507356205638|1.47041447230994514445|+1.0064537%|

pc_ref=0.59274605079210 是已有参考值，不是严格pc区间。圆柱根也是已有打印值，不能把位移比的全部显示位数当认证精度。每个k的三rank和等于binom(L^2,k)已复核，但这不认证实际rank分类。未重新枚举2^25或2^36，未重跑Perron。

两个小尺寸趋近目标不是渐近证明，也不唯一识别thermal模块；这是同一旧数据的新比较，没有拟合振幅/指数来选目标。参考误差灵敏度为 `R'(c)=(a-b)/(b-c)^2`，当前约-951、-1971。今后若已有三个同面积不同形状的根，可用根差之比同时消去pc与A；不为此自动开启形状阶梯。

## 7. 现在只追真实接口

#768：计算 rank 投影下应力张量运输是否产生seam/contact项，以及 F(L_-2 chi) 是否消失。可只证明一个实际扇区/源的一点关系，不要求先解决整个对数模块。零范、奇异与在指定关联函数中解耦不是同一句话。

#802：先修正上一轮已指出的正规化与 H_W=N-2K+H_B 冗余，保留两压力与各自导数。复用现有 finite-aspect 根检验(7)，一个不符合E4的法向/形状响应比更多同轴根更能改变判断。非E4残差可能是null或投影异常，也可能是同阶scalar项，不能直接命名。

#780：新 `intrinsic-birth-clock-construction-20260914.md` 已用最终簇出生CDF beta 构造单调时钟；本轮吸收该进展，不重新购买时钟存在性或更多Poisson核矩。压力时钟是另一实现，不是第二套独立生产。

本轮只加本文件、独立控制脚本、测试和结果。没有修改原始U、STATUS或团队文件，未启动机器。12项局部测试核对PBW、null异常、模变换、eta系数和原始表算术；高精度模函数/求根不是区间认证，未运行全仓CI。

## 来源

[GL] Gaberdiel–Lang, arXiv:0810.0106，PDF印刷页3–4、6，(2.3),(2.10),(2.14),(3.4)–(3.7)：标准torus Ward/Zhu和null微分方程。只取实际读取的递推，不将一般有理性假设移植到渗流。https://arxiv.org/abs/0810.0106

[JPS] Javerzat–Picco–Santachiara, arXiv:1907.11041v2，§2.2、4.2，(4.11)–(4.18)：热/磁性融合、二级null及四级后裔，Q1无条件能量零与连通函数的奇异正规化。其Potts/FK连通不等于本项目rank投影。https://arxiv.org/html/1907.11041v2

[HS] He–Sun, arXiv:2004.07486v2，§2.1(4),(5)：未变形CFT的周期Ward恒等式。这里不引入T Tbar变形模型。https://arxiv.org/html/2004.07486v2

内部输入：主读SHA上的L5/L6 rank-sector-C JSON，fixed-width-charge-spectrum-derivatives-w4-w8 JSON，double-null、two-observable-response、topological-clapeyron及intrinsic-birth-clock四份说明。完整输入路径和blob见结果JSON。

复现：`python scripts/thermal_ward_root_control.py --output results/research-dispatch/thermal-null-ward-root-20260914.json`；`python -m unittest discover -s tests -p 'test_thermal_ward_root_control.py' -v`。依赖mpmath，精确代数用标准库Fraction。
