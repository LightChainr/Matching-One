# Matching-One 长期探针:12 对话执行计划(Plan)

**状态**:研究路线图 · 2026-09-06 制定 · 归属 #599 框架 · **v2 已按 C1 快照(2026-09-06 下午)修正**

## 执行状态(2026-09-06 晚)

| 对话 | 状态 | 主交付 |
|---|---|---|
| C1 快照校正 | ✅ | `probe-baseline-snapshot-v2.md`(#600/#601/#603 吸收,Round1/2 重新定位) |
| C2 分叉根计数代数 | ✅ | `probe-fork-root-count-algebra-20260906.md`(单根秩≤2;深度=代数次数;d*=k) |
| C3 w9 记忆度 | ✅(acquisition) | `probe-w9-memory-acquisition-20260906.md` + 成本探针 |
| C4 语言形式化 | ✅ | `probe-experiment-languages-formal-20260906.md`(κ_L/r_L、最小区分深度、签名类表) |
| C5 选择规则+L2π/纤维 | ✅ | `probe-markov-selection-rule-and-l2pi-fibre-20260906.md` |
| C6 #582 判别器 | ⏸ 数据依赖挂起(需 #582 档案 Q_N 数据/授权);协议与决策判据已备 | 判别器协议见 threshold note Part 2;C10 排名 #7 |
| C7 阈值 no-go 定理化 | ✅ | `probe-threshold-no-go-theorem-20260906.md`(偶任务对 g 不变 1e-16;gap=g·m→0) |
| C8 可辨识性目录 | ✅ | `probe-invariant-catalogue-20260906.md`(不变:选择零/商恒等式/扇区惰性) |
| C9 阶梯定理 v1 | ✅ | `probe-complexity-ladder-v1-20260906.md`(全部箭头 P/F/C/U + 反例档案 X1–X9) |
| C10 猜想排名+廉价测试 | ✅ | `probe-bridge-conjectures-ranked-20260906.md` + T1/T2 实测(T1:动态签名饱和于轨道商;T2:标记通道见 c) |
| C11 提交包 v2 | ✅ | `submit-pack-v2/` + `matching-one-probe-full-20260906.bundle`(单分支,不合并) |
| C12 总评 | ✅ | `probe-mission-update-20260906.md`(北星命题获支持并加精确边界) |
**已完成基线**(不再占用对话):Round1 = #598 Phase C 商因子化数值(P398,谱一致至 1e-16,内积陷阱 37–61% 定量);Round2 = 一般有限群响应选择定理(#244+#598 统一,完整证明,C3 数值验证,isotypic vs 轨道商纠正)。

**交付约定**:每个对话产出 1 个主 Markdown 文档(+可选脚本/JSON);无 GitHub 写权限时以 md 存档于 `/workspace/matching-one-round1/notes/`;授权恢复后统一打包提交。

---

## 战略地图

研究分三条并进主脉,12 个对话按"信息增益优先 + 依赖约束"排序:

```
脉 A  预测/线性实现层级理论   (C2→C4→C9→C10)
脉 B  表示论/平衡/内积验证     (C5, 依赖 C1 快照)
脉 C  RG/阈值可辨识与判别器   (C6→C7→C8)
包装层 快照校正(C1)、综合(C9/C10)、提交(C11)、复评(C12)
```

对话依赖(C1 → C5; C2 → C4 → C9; C6 → C7; C1/C2/C5/C6/C7 → C9; C9 → C10 → C11 → C12)。C3 独立,可并行或后移(计算重型,受预算约束)。

---

## 对话卡片

### C1 · 仓库前线快照校正(基线重映射)
- **目标**:研究开始于准确快照。重读 #596/#595/#599 及其引用分支的最新评论/结果,核对 round1/2 结论是否与仓库最新状态冲突;校正 atlas 的"已知未知"。
- **输入**:issue #596/#595/#599 + 相关 open PR 头、最新结果 JSON。
- **成功**:atlas v2 + 一份"与仓库现状冲突清单"(若为空也写明);无未声明冲突进入下一步。
- **失败退出**:仓库冲突需要owner裁定 → 停在该对话并上报。
- **输出**:`notes/probe-baseline-snapshot-v2.md`。

### C2 · Program B 深化:分叉概率的代数次数与秩增长
- **目标**:回答 #599"最小 future 语言深度使响应秩必增"。获取并精读 #549 协议脚本(`p429_parallel_gadget_lower_bound.py`)与其 notes;本地重建单轮/双轮共享更新-双克隆协议。
- **关键问题**:(1) 双轮概率在类坐标 a 中是否恰为二次(→秩 3 可达)?(2) 组合(同时触及多 gadget)测试的代数次数是否等于语言深度?(3) 若 #549 族内秩无法增长,给出精确的秩上界定理。
- **成功**:给出显式测试集与响应矩阵(精确/有理),证明秩 3(或证明秩 2 上界)。
- **输出**:`notes/probe-fork-polynomial-degree-<date>.md` + 脚本 + JSON。

### C3 · P398 w9 记忆度饱和检验与重计算预算评估
- **目标**:把记忆度曲线(数值阶 4,9,12,13,14)推进到 w9(4862 态),区分"饱和"与"亚线性增长";复核 #588 修正结论。
- **方法**:稀疏均匀化+快 Krylov;先做 5 分钟成本探针,预算不足则产出 heavy-test acquisition 文档(按任务 §4 纪律)。
- **成功**:w9 的数值阶/能量阶 + 增长法则判定;或规范 acquisition 文档(预测结果、算力、每分支结论)。
- **风险**:纯 Python w9 均匀化矩阵指数可行但耗内存;上限保护。
- **输出**:`notes/probe-w9-memory-saturation-<date>.md`(或 acquisition)。

### C4 · 最小严格分离族:预测类 vs Hankel 秩 vs 非负秩 vs 正实现
- **目标**:构造/定位一个原生对象族,使 预测类计数、普通 Hankel 秩、非负秩、正实现维、强 lumping 维 出现**可任意放大的严格分离**,并用已有工具(C2/C3、#401 三层、P398 字典过滤)生成"最小分离表"。
- **方法**:优先仓库原生(#549 体系 / #401 q2/q3 矩);否则声明性合成对象并注明。
- **成功**:一张表(全部箭头标 PROVED/FALSE/UNKNOWN)+ 至少一个最小反例文件。
- **输出**:`notes/probe-strict-separation-family-<date>.md`(Program A 半成品)。

### C5 ·(重定位)马尔可夫生成元选择规则正式陈述 + isotypic 纤维数值(填补 #601 Q3 空缺,避开 #603)
- **目标**:(i) 将 Round2 Theorem 1 改写为对准 #601 Q3 空缺的正式陈述(明确引用 Schur/Diaconis/Wonham,只保留马尔可夫生成元逐点形式 + ℓ 阶张量幂 + C2-accident 说明);(ii) 补 #603 未做的部分:非 C2 对象上 Theorem 2(c) isotypic 纤维实现数值(C3)、L2(π) 内积下的全空间=商谱一致性。
- **成功**:陈述被 #601 Q3 语言精确承接;纤维实现数值闭合;两种内积一致 ≤1e-10。
- **输出**:`notes/probe-markov-selection-rule-formal-<date>.md` + 脚本 + JSON。

### C6 · Program E 判别器实验(#582 档案重分析)
- **目标**:在 #582/#584 已有统计量(results + notes)上跑"光滑 Taylor null"判别器清单(cocycle/步长一致、二阶曲率 vs 一阶导数、跨 scale-word/lineage 留出);判定主导一方向流是否需要超光滑结构。
- **方法**:不新增采样,只用档案充分统计量 + 合成 smooth family 作为 null 参照系。
- **成功**:对第 5 个跃迁的"超曲率余项"给出可否被光滑族复现的结论;输出判别器清单的实证栏。
- **输出**:`notes/probe-smooth-null-discriminator-run-<date>.md`。

### C7 · Program F:阈值不可辨识小定理族(正式化)
- **目标**:把 Round1 范围笔记升级为**可证明的小定理**:在声明的任务语言(有限源/读出/时域/分支深度)下构造两族系统,任务数据全同而"无限体积"可判据不同;列正定理所需条件(#594 no-go 目标)。
- **方法**:以 #435/#549 的 p-参数化 gadget 为基底构造显式族;若在库内做不出严格版本,给抽象定理 + 一个合成精确实例。
- **成功**:至少一个自足精确的 no-go 构造(定理+证明+数值)。
- **输出**:`notes/probe-threshold-no-go-theorem-<date>.md`。

### C8 · Program G:可辨识性目录(observer/quotient/内积下不变 vs 可变)
- **目标**:系统目录"什么随 observer 字典/内积/商约定/归一化/lift 变化、什么不变"。分类库中候选(选择零、投影响应子空间、谱、minimal distinguishing language、表示重数、扩展类)。
- **成功**:目录 md(每项给:定义、变/不变、证据或反例、依赖)。
- **输出**:`notes/probe-invariant-catalogue-<date>.md`。

### C9 · Program A 完成态:复杂度阶梯定理文档(成熟 v1)
- **目标**:整合 C2/C4/C5/C6/C7 全部结果,产出 `complexity-ladder-theorems.md` 成熟版:复杂度概念精确定义、全部箭头(PROVED/FALSE/TRUE-UNDER/UNKNOWN)、每个反例文件化;修正/升级 round1 atlas。
- **成功**:每个"未知"箭头有最便宜验证的描述;无未标注推测混入事实。
- **输出**:`notes/probe-complexity-ladder-v1-<date>.md`。

### C10 · 桥接猜想排名与 top-3 廉价测试
- **目标**:更新"未解决桥接猜想"排名(#599 遗留),每个猜想给最便宜 falsification;对 top-3 实际执行廉价测试(如可离线完成)。
- **成功**:排名表 + 3 个廉价测试的实际结果/结论。
- **输出**:`notes/probe-bridge-conjectures-ranked-<date>.md`。

### C11 · 提交包 v2(双速度沟通)
- **目标**:把 C1–C10 的 exact findings 整理为紧凑 issue 注记(#598/#593/#582/#599 相关)与统一 md 提交包;尝试经授权通道推送 PR/评论;若不可行则归档可一键提交的包。
- **成功**:所有材料按仓库现状正确归类;提交动作完成或被明确阻塞(原因与解锁方式写明)。
- **输出**:`/workspace/matching-one-round1/submit-pack-v2/`。

### C12 · 总评与方向修订(长期复评)
- **目标**:以 10 个对话的证据重估核心命题("复杂度相对实验语言与对称定义")——支持/反驳/需修正;产出新的长期方向集(若继续)。
- **成功**:一份可辩护的 mission-update 文档,含遗留最大开放问题的排序。
- **输出**:`notes/probe-mission-update-<date>.md`。

---

## 执行约束与纪律(每对话适用)

1. 证据分级强制:FACT/DERIVATION/CONJECTURE/COUNTEREXAMPLE/IDENTIFIABILITY-LIMIT/CHEAP-TEST/HEAVY-TEST 标注;禁止 conjecture 静默升格。
2. 计算纪律:确定性、无新采样、无宽度无谓升级;重计算先探针成本,超预算改 acquisition 文档。
3. 仓库对象优先;合成对象必须显式声明。
4. 语言:研究文档英文(与仓库一致);用户沟通中文。
5. 每对话结束更新 `/workspace/matching-one-round1/notes/PROGRESS.md`(一页:结论、产出、开放)。

## 风险与预案

| 风险 | 预案 |
|---|---|
| #549 协议脚本缺失/定义不清(C2/C4) | 降级为"由发布闭式可证的秩上界定理"+ 标注需 owner 提供脚本 |
| w9 计算超预算(C3) | acquisition 文档,纳入 HEAVY-TEST 台账 |
| #582 档案统计不足(C6) | 退到"判别器清单 + 合成参照"并可复现的流程文档 |
| GitHub 写权限持续缺失(C11) | md 归档提交包,给出授权后的精确步骤 |
| 发现与仓库新结果冲突(任意) | C1 式冲突协议:停下、写对照、请求裁定,不硬推 |

## 判定"计划成功"的总指标

12 个对话后应产出:复杂度阶梯定理 v1(箭头图全部标注)、≥3 个自足精确反例/定理(其中 ≥1 个 no-go)、#582 判别结论、可辨识性目录、排名桥接猜想;并明确回答北星问题:experiment-relative 复杂度框架在仓库对象上被证实、证伪或被如何修正。
