# #537 N65 审计、full-T 边界与停止线

## 1. 前瞻顺序与可复现性

GitHub PushEvent 与固定提交顺序为：

```text
45ca59ba  frozen contract
76e2d82e  producer/scorer
f46c38c3  first result
```

合同先于 20M 新 counter block，四个 shard、100 个 production batches、seed
`20260901537`，没有 top-up。独立完整重放的四个 TSV 与原始 shards 逐字节一致；
冻结 scorer 的中心结果也被复现。原 `result.json` 保持不改，新增
[`AUDIT.json`](../results/p537-contact-stage-n65/AUDIT.json)只完成合同已经要求的
positive exposure 与 covariance 保存。

## 2. 安全的 primary 结论

固定 canonical `g16/(16N)` source、同一 selected carrier 与同一 Schur 分账下，

```text
                  single contact        double contact
birth 0->1       -1.57857291e-7        -9.22765780e-8
birth 1->2       -3.09129846e-7        +3.69680209e-7
```

预注册符号 `[-,-;-,+]` 与 `Delta<0` 均通过：

```text
Delta = -8.6882160551e-14
95% CI [-1.1459612320e-13, -5.9168197900e-14]
```

完整 6×6 covariance 为对称 PSD；collapse covariance 等于 `T C T^T`。六格和为

```text
-1.8958350633e-7 +/- 6.1808292092e-8
95% CI [-3.1072775883e-7, -6.8439253831e-8], z=-3.0673.
```

忽略 off-diagonal covariance 会把总和 SE 低估约 28.36%。六个边际区间虽都排除
零，仍不能算六张独立 evidence votes。`theta=-1` 在这个开符号象限内是代数恒等，
全部 delete-one replicate 均为 -1；零 SE 不表示无限精度或“最大耦合”。

## 3. thermal-gauge 边界

selected contact-stage cells 不是独立物理算符。在共同热坐标变化

```text
a -> a + cK + d
```

下，完整 full-T quotient 不变，但固定 class `C` 的分账一般变为

```text
T_C -> T_C + c C_C,
C_C = p(1-p) sum_z <I_C Htilde>.
```

真实 N65 raw 的 horizontal representative 保留 canonical 中心附近的同一符号和
determinant；然而 unit gauge `c=+/-1` 足以翻转多个 cells 与 determinant。故 primary
只拒绝**冻结 canonical 坐标内**的 scalar/separable law，不能命名为 gauge-invariant
CFT operator 或 physical commutator。

## 4. 完整 full-T secondary

[`f9ba1ff6`](https://github.com/LightChainr/Matching-One/commit/f9ba1ff690b07beefcc71e669f1f29581d4e264e)
只重用同一 N65 sufficient statistics，没有新样本。遗漏的 `+e1` NN source column
可由另外三个 NN 方向作 C4-unbiased 重建：canonical kernel 对全部 4,140 个 Bell8
分割一致，16,560 次共同 90-degree rotation 无违例，两套 N65 quotient tori 均有
对应自同构。该重建在 expectation level 成立，不是样本内逐条恒等；三方向共同
jackknife covariance 已传播。

安全的有限结论是

```text
J25 = -0.00551943142484
J65 = -0.00162250989 +/- 0.00018553008
J65/J25 = 0.29396323 +/- 0.03361398.
```

因此完整 canonical response 在这一个 exact-to-MC 尺寸方向明显收缩；selected
carrier 不饱和 full response。由两个尺寸计算的 power 1.281 只是描述，不证明
`J_N -> 0`、little-o、5/4 指数或 CFT 身份。

## 5. 不吸收与下一停止线

- `bab37f21` 的 `N^-3`、`N^-29/8`、commutator/triangular/CFT 解释来自同一
  N25/N65 dependency block，统一列 post-hoc P2 hypothesis。
- N145 合同与 producer 在 #537 降为 P1、#275 升为 P0 后出现；合同直接输出 full-J，
  并不直接测 contact-subtracted remainder，scorer 也未强制 frozen samples/seed/
  shards/decision。任何已运行 raw/result 只作 parked P2 audit asset，不 top-up、不
  改 priority、不用于再拟合指数。
- #537 只保留 contractible-collar quotient、bounded normalized pivotal 与
  near-critical uniform transport 三项证明或反例。任一 physical counterexample
  退休 `J_rem~N^-5/4` 路线；三项全闭合并推出 `T_N=o(M_t/A_N)` 才能关闭 #537。
