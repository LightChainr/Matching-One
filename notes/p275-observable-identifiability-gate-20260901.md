# #275：先固定原始观测量，再决定 H4/H8 是否可识别

## 当前判决

#275 是唯一 P0。任务不是恢复旧 E_top、C3 或 Q-lift 扫描，而是把

```text
raw observable -> normalizer -> pooled moving-root U -> candidate prediction
```

写成同一个可检验对象。#537 保持 P1；#154、#334、#337 保持 P1，均无自动续跑。

## 已完成的前向设计

[精确相位合同](../experiments/p275-c3-phase-contract-20260901/REPORT.md)已经把
“再加一个物理旋转”拆成三个互斥条件：

1. 若两几何共享同一未知复振幅，`sin(6 delta) != 0` 即可区分 H4/H8。
2. 若允许未知非零实增益及符号，`delta=7.5°` 可区分，`15°` 有精确符号别名反例。
3. 若允许任意复增益，任何两个旋转都不可识别。

所以现阶段不能仅凭形式 character 正交、spin 标签或一个“最优角度”申请生产。
真正缺少的是 amplitude/phase transport：原始 trace、normalized expectation 与
moving-root quotient 中，哪一项在两个物理几何间共享，哪一项允许独立变化。

## 唯一下一交付与停止规则

下一份 P0 交付只允许以下两种结果之一：

- **可识别：** 给出两个候选在至少两个相关原始坐标上的单位、normalizer、允许振幅类和数值/符号预测；先用已有资产计算协方差加权设计秩。只有已有资产缺少唯一必要坐标时，才冻结一次新采集。
- **不可识别：** 证明在物理允许的振幅/normalizer 自由度下两个候选列空间相同，并明确指出唯一能解除退化的输入；该 C3 识别线随即降为 P2，不再增加角度。

结果出来后不得增加第三种振幅合同挽救失败。有限 character 非零、另一个 source、
另一个 minor 或两点幂律都不算 P0 完成。

## #537 的交接边界

N65 的 frozen contact-stage 四符号与 `Delta<0` 是一次真实前瞻结果，但只对
canonical selected carrier 分账成立。该分账在共同 thermal gauge
`a -> a+cK+d` 下一般会移动，尚未成为坐标无关物理算符；完整 full-T quotient、
两点指数和 N145 外推不能替代 #275 的原始观测量字典。#537 后续只接受解析
等变性/归一化闭合或明确反例，不自动排新尺寸。
