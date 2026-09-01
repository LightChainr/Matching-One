# 下一步：先固定原始观测量，再做一次可识别性判决

**更新至 2026-09-01 16:00（Asia/Shanghai）。** 当前 GitHub 生命周期为
#275=`priority:P0`、#537/#154/#334/#337=`priority:P1`。本文件只列现在允许继续的
工作；旧路线、完整结果与固定提交见 [STATUS](STATUS.md)。

## P0：#275 original-observable identifiability

当前核心不确定性是：现有 H4/H8、trace/source 与 normalized original-`U`
结果是否在说同一个物理观测量。下一交付必须把

```text
raw q/E coordinates
  -> physical normalizer
  -> pooled moving-root U
  -> two candidate prediction vectors
  -> covariance-weighted profile rank
```

闭合为一个包。候选至少覆盖两个相关原始坐标，逐项写清单位、共同/独立振幅、
相位、root 与 normalizer；不能用相同 spin 标签代替这张 map。

[C3 精确相位设计](../experiments/p275-c3-phase-contract-20260901/REPORT.md)
已经实际完成 66 项 exact checks，并给出第一个决策：

- 共享未知复振幅时，第二旋转满足 `sin(6 delta) != 0` 即可区分 H4/H8；
- 允许未知非零实增益及符号时，`delta=7.5°` 可区分，`15°` 有精确别名反例；
- 允许任意复增益时，任何两个旋转都不可识别。

所以“再测一个角度”本身不是 P0。理论必须先证明原 observable 属于哪种
amplitude/phase transport 合同，再把现有 C3、rho-child、P43/P57 或 modulus
资产接入同一 covariance。详见[当前判决与停止树](../notes/p275-observable-identifiability-gate-20260901.md)。

**Stop rule：** 若保留全部物理允许振幅后两个候选张成同一子空间，登记
`UNIDENTIFIABLE_WITH_CURRENT_ASSETS`，把该候选对降为 P2；只允许理论指出的一个
缺失输入。若设计满秩，只做一次冻结评分并降级失败模型；不得新增第三模型、
descriptor 或事后 amplitude class。已有资产确实缺少唯一必要坐标时，才可冻结
一次采集。

## P1：#537 保留 finite facts，停止尺寸外推

N25 one-defect/contact 判决已否定 automatic two-independent-defect/six-arm gain。
N65 20M 新块又在冻结的 canonical selected-carrier 分账中复现
`[-,-;-,+]` 与严格负的 `Delta`；这是一次真实前瞻判决。

合同补充审计保留完整 6×6 covariance、positive exposure 与 selected total。
它同时限制解释：六格相关，不能算六票；`theta=-1` 由符号象限恒等强制；
selected contact-stage cells 对共同 thermal gauge 一般会改变，不能直接命名为
坐标无关物理算符。

`bab37f21` 的两点指数、commutator/CFT 解释只列 post-hoc hypothesis。
`f9ba1ff6` 的 full-`T` N65 scalar 已通过 C4/kernel/normalization 审计，登记为
同一依赖块上的 finite original-`U` secondary；N25/N65 只给描述性收缩，不给
little-o、指数或独立 evidence vote。

N145 200M 合同/producer 是 P0/P1 重置后出现的分支资产。它依赖 reveal 后的两点
幂，scorer 还未强制样本数、seed、shards 与 frozen decision，且输出测的是 full-J
而非直接测 contact-subtracted remainder。该线不得 top-up、不得据此改 priority；
任何已落盘结果先隔离为 P2 audit asset。

**Stop rule：** #537 只接受 contractible-collar quotient、bounded normalized
pivotal domination、near-critical uniform transport 三项的证明或物理反例。任一
反例立即退休 `J_rem~N^-5/4` 路线；三项全成立并推出
`T_N=o(M_t/A_N)` 才达到关闭条件。第三尺寸拟合不能替代证明。

## 其他队列

- #154/#334 的独立生产已触发各自降级，不换 lag、不补 prefix、不在验证块救场。
- #337 保留 exact/paradigm P1，但只服务 #275 的具名 observable map。
- #419/#370/#398/#539/#542 与代数全族搜索保持 P2；不新增 generic certificate、
  低次多项式高度/次数、距离/动量网格或描述性坐标。
- 十台 Huawei 容器按需可用，但机器数量不决定任务数量。当前没有经过上述
  stop rule 的远程生产任务。

结果继续交付 Draft #509；维护 PR #528 只同步导航与生命周期。**不合并，
不删除历史数据、冻结合同或分支。**
