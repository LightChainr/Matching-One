# 团队机制收敛建议：把 #154 / #334 统一为 birth-current transmission 判决问题

日期：2026-08-31

状态：团队建议 / 决策框架，不改写既有科学结果，不自动改变 Issue 优先级或停止已冻结任务。

## 结论先行

仓库已经完成了一次有效的任务数量收敛：当前 P0 只保留 #154 与 #334。下一步不应继续主要通过“增加描述性 observable / contact feature / residual coordinate”推进，而应把这两个 P0 合并成一个共同的机制判决链：

```text
microscopic source
    -> ingress / egress birth currents
    -> J_grad = j_in - j_out
    -> d_p E_top
    -> original pooled-root / slope-normalized U
```

核心问题不再是“还有什么结构与 residual 相关”，而是：

> **哪一个已冻结的微观模型，能够在未见数据上预测 `j_in-j_out -> E_top' -> U` 的变化？**

这是当前最直接的 hypothesis-killing 接口。

## 1. 为什么 #154 与 #334 已经是同一个问题

#334 已经给出 projective continuity 结构。对 rank-one plateau line `ell`：

```text
J_grad,ell = j_in,ell - j_out,ell,
J_act,ell  = j_in,ell + j_out,ell.
```

并且

```text
sum_ell J_grad,ell = d_p P(r=1) = f1-f2.
```

另一方面

```text
E_top = P0 + P2 = 1-P1,
```

因此

```text
d_p E_top = f2-f1 = -sum_ell J_grad,ell.
```

而 #154 的原始 norm-4 `U` 正是对 `d_p<E_top>` 做 pooled-root、H4 projection、有限尺寸归一化和 matching-slope normalization 后得到的量。

所以真正的机制闭环不是“local source -> E_top 相关性”，而是：

```text
source
 -> birth-current imbalance
 -> E_top thermal response
 -> U.
```

这也解释了为什么以下两类结果可以同时为真：

1. fixed `(K,rank1)` 内部 O4 / winding / contact source response 很强；
2. 同一中心化 source 对 global U 的直接响应可以严格为零。

强 `J_act` 或 rank-one 内部重分布，不必产生 `J_grad`。因此不能把“局部 response 很强”自动升级为“解释 global norm-4 residual”。

## 2. 当前最需要避免的机制误判

团队目前最容易发生的误判是：

```text
测到强 local / common-clock / contact response
        ↓
把它当作 global U 的候选来源
```

但 global `E_top'` 看的是 ingress/egress 的差，而不是总 activity。

因此后续所有 #334 microscopic mechanism 结果都应至少明确回答：

```text
它主要预测 J_act？
还是主要预测 J_grad？
还是二者都有？
```

如果一个结构只能稳定预测 `j_in+j_out`、两个 birth center 的共同移动或 rank-one 内部重新分配，那么它仍是有效的局部/过程机制，但不能被称为 Matching-One global anomaly 的主要生成机制。

## 3. 暂停扩充 contact feature；先检验 transmission

#334 当前 contact-regulated clock susceptibility 是一个具体而有价值的机制 lead：已冻结 contact structure 能解释大部分 signed loading，在高精度 original00/new64 population 上仍保留稳定同号的约 20% remainder。

建议：

- **不要因为这 20% remainder 继续自动添加第五、第六、第七个 contact descriptor。**
- 固定当前 descriptor span。
- 下一项核心问题改成：这些 contact coordinates 能否预测真正进入 `J_grad`、`E_top'` 和 `U` 的 response？
- 若当前 span 只能解释 birth-center / common-clock loading，而不能预测 held-out `J_grad/U`，则应把它降级为局部机制，而不是继续扩 feature 直到 residual 变小。

同理，固定 prefix 内 `E[det J_A(Z)]>0` 说明局部二维 A-response 确实存在，但这不是“第二个 global norm-4 state 已找到”。必须进一步证明二维 response 的第二方向能传播到 `J_grad -> U`。

## 4. 将机制空间压成三个嵌套模型

### M0 — lifecycle/current-only

假设：global U 的 source response 只需要 entry/exit、rank population、K1/K2、lifetime/current imbalance；不需要额外 prefix spatial state。

已有支持：fixed `(K,rank1)` 内部中心化 source 对 U 严格零，与“内部强响应不直接进入 global U”一致。

淘汰条件：冻结 contact/geometry 后，在 held-out 数据中它们对 `J_grad` 或 U 仍有稳定额外预测力。

### M1 — contact-regulated current

假设：entry/exit current 还需要当前冻结的一组 contact/loop state；这些 microscopic coordinates 调节 `j_in-j_out`，从而进入 global U。

已有支持：contact structure 能组织大部分 signed clock-response loading。

必须通过的检验：

- descriptor 事前冻结；
- 预测 held-out `J_grad` / E / U，而不只预测 local clock；
- 至少跨一个独立 N、geometry 或新 counter block 验证；
- validation 失败后不能在同一 block 上增加新 feature rescue。

### M2 — extra transfer coordinate

假设：完整 lifecycle + 冻结 contact state 仍不足，需要额外 transfer coordinate（例如 Jordan/even-mode/一般 matrix state）。

当前状态：原 q2 已被拒；Jordan + one-even-mode family 可存活，但额外 eigenvalue 尚未识别。这只是 survivor，不是物理识别。

升级条件：只有 M1 在独立 validation 上失败，且 #275 / theory side 能给出对同一 target vector 的冻结跨尺度/跨几何预测时，M2 才获得新的 P0 资源。

## 5. 当前正在做的任务如何处理

### #154 temporal-source：继续，而且应提高为主 transmission test

当前固定唯一 lag `ceil(sqrt N)`、不扫 lag；早 source 在 `(K',r_early,g)` 条件下中心化，使即时 rank response 为零，再测未来 rank / entry-exit / U。

这是当前最直接测试“早期隐藏微观状态能否通过 lifecycle/current 传到 global U”的设计，应完成。

但揭晓前必须冻结至少：

- M0 / M1 对主要承载 channel 的预测（entry、exit、lifetime 或其组合）；
- 预期符号或至少相对方向；
- 什么结果算 M0 失败；
- 什么结果只算“lagged correlation 存在”，而不能升级成 field/mechanism identification。

否则即使得到高显著度，也只会新增一个相关量。

### #334 held-out contact prediction / conditional shape：继续

这些任务至少含 held-out / cross-quartet 评价，比继续做描述性 projection 更有价值。

完成后设置硬暂停点：

> 不根据同一 validation residual 再新增 contact descriptor。

### mixed-source Hessian：可完成廉价 exact/reuse 分析，但降为 exploratory

即使 mixed response 非零，也最多说明两个已声明、可交换 source coordinates 的 response 非加性；不能自动解释成 memory、noncommutation、新场或 global U mechanism。

它不应自动开启新的 P0 分支。

## 6. 统计纪律：当前主要风险已经不是 covariance，而是 adaptivity

仓库当前对 covariance / dependency group 的处理是正确方向：同一 random block 的多个 observer 保持一个 covariance block，不把多个 derived view 当独立证据。

但仍有更高层的 adaptive-search 风险：

```text
看 residual
 -> 设计 contact
 -> 看 contact residual
 -> 设计 shape
 -> 看 shape
 -> 设计 Hessian
 -> ...
```

每个单独 p-value 都可以计算正确，但 hypothesis generation 本身一直在消费同一个 archive。

建议把数据明确分为：

```text
discovery/training archive
vs
decision/validation block.
```

下一块独立 counter-domain 的中等规模 validation，信息价值高于继续把当前 block 扩大十倍。

在 validation 之前冻结：

- mechanism class；
- descriptor / lag；
- primary receiver；
- sign / ratio / vector prediction；
- covariance score；
- rejection rule。

如果 validation 失败，不允许在同一 block 上新增 feature 来把 M1 修回来。新机制可回到 discovery 阶段，但要等下一块独立数据再验证。

## 7. 下一次新 production 应是联合机制验证，而不是默认更大 N

先做一个零新样本的统一步骤：

1. 把 #154 与 #334 的相关量统一写到 `j_in / j_out / J_grad / J_act` current coordinates；
2. 写清原 U 的 moving-root + slope normalization 在 current language 下的完整公式；
3. 判断现有 contact result 主要落在 activity common mode 还是 gradient mode；
4. 冻结 M0/M1/M2 对下一块数据的不同预测。

之后才启动新 validation production。

新 block 的成功标准不是“又显著看见一个结构”，而是：

> 冻结的 contact/lifecycle 模型能否在未见数据上同时预测至少一个 birth-current quantity 和一个 global E/U quantity？

做不到，就不能说 #334 已经解释了 norm-4。

## 8. 三队职责建议

不需要再进行组织重构，保留现有三队名称，但硬化职责边界。

### 总览队

负责 model compression + prospective prediction。

- 维护科学 atlas；
- 压缩 M0/M1/M2；
- 在 validation 前写冻结预测；
- 不在 reveal 后立即发明 rescue feature。

#275 主要服务这一层：给同一 target vector 的理论预测，不再平行增加 source 候选。

### 执行队

只负责 frozen calculation / production。

输入必须是明确的 target vector、dependency block 和 decision rule。执行队不负责根据 reveal 发明第四个 mechanism class。

### 俯瞰队

改成真正的 red team。

负责：

- chronology / post-reveal 标记；
- dependency / covariance；
- 检查模型是否 reveal 后扩容；
- 判断 claim boundary；
- 需要时调用 #370 对明确模型做 certificate/no-go。

俯瞰队不再同时充当第二 exploratory-analysis team。

## 9. 仓库治理建议：不建第三套战略文档

#267 保持历史 scientific atlas；#509 保持本周期 delivery。不要再创建一个新的长期导航体系。

只建议把 `NEXT-TARGETS.md` 的 P0 表从“下一分析问题”升级成“机制判决表”，每个 P0 增加五个字段：

```text
model being tested
primary receiver
frozen prediction
rejection condition
data status: discovery vs decision
```

当前文件可先作为建议 PR 审阅，不在本提交中自动改写既有 priority board。

## 10. 资源建议

下一个决策周期建议：

- 35%：#154 current / temporal transmission；
- 30%：#334 frozen-contact held-out prediction；
- 20%：#154/#334 共同 current-language、model freeze 与独立 validation 设计；
- 10%：#14 严格界独立路线；
- 5%：维护与必要的 #370 concrete verification。

以下路线不分配新的主线 production 预算，除非出现明确的新判别目标：

- #13 相邻 serial-algebra census；
- #370 generic framework expansion；
- #419 generic high-pass expansion；
- #398 自动扩 width / 扫参数。

它们可以保留廉价独立 exact/support 工作，但不应竞争 P0 的认知预算。

## 最终 stop rule

以后一个新的 P0 分析在开始前必须能回答：

> **什么结果会让我们停止或降级一个明确机制？**

如果答案只是“这个量是否非零”“这个 residual 还能解释多少”“是否还能找到另一个相关 descriptor”，默认降为 exploratory/support。

项目现在缺的不是新的 observable 数量，而是更高的 mechanism-elimination rate。

## 建议的下一次正式判决

冻结 M0 / M1 / M2，然后把下一次 validation 的核心问题写成：

```text
Can a frozen microscopic lifecycle/contact model
predict the unseen birth-current imbalance J_grad
and the resulting original norm-4 U simultaneously?
```

这是目前最可能真正减少 Matching-One 机制不确定性的团队级实验。