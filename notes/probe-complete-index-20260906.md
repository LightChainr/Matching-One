# Matching-One 长期研究探针 · 完整交付

本目录整合两个自主研究程序的全部成果,对应 GitHub 仓库 `LightChainr/Matching-One` 的长期研究探针任务(#599 框架)。

- **项目一 · 预测复杂度探针**(`matching-one-round1/`)——围绕 P398 对称商、平衡实现、有限群响应选择规则、实验语言形式化、阈值 no-go、复杂度阶梯,共 12 个对话单元(C1–C12)全部推进完成。
- **项目二 · 切割网络 / 分支预测程序**(`matching-one-cutnetwork-round1/`)——围绕"分支预测状态 / Hankel 秩 / 语言深度",共 10 个程序(A–J),可精确闭合方向全部闭合。

---

## 验证结果(2026-09-07 重新运行)

全部 **12 个可执行脚本 + 1 个共享库**(`p398core.py`)逐一遍历运行,**全部通过**,核心数值复现到机器精度:

| 脚本 | 关键输出(验证值) |
|---|---|
| `cutnetwork/scripts/rank_lower_bound.py` | 秩升级:K=2..8 恒 `rank(M)=k+1`、`degree(E_j)=j` |
| `cutnetwork/scripts/radius_insufficiency.py` | 无有界半径商:k=1 时 fork gap 精确 `1/98`(#435 见证) |
| `cutnetwork/scripts/slack_family.py` | 正/有符号分离:正方形 rank=3、非负秩=4(精确) |
| `round1/scripts/probe/verify_gate1_independent.py` | `GATE1_VERDICT: PASS`(轨道商=lumping、固定点、奇读出) |
| `round1/scripts/probe/rank_vs_classes_demo.py` | F_{k,a} 仿射、Vandermonde 秩=k+1 |
| `round1/scripts/probe/group_selection_demo.py` | C3 选择规则:零峰值 2.6e-17、二阶斜率 2.006 |
| `round1/scripts/probe/c2_root_rank.py` | 根计数=深度=代数次数;rank=k+1 |
| `round1/scripts/probe/c7_no_go.py` | 偶任务对隐藏耦合 c 不变 ~1e-16,谱隙随 c 变化 |
| `round1/scripts/probe/c3_cost_probe.py` | w9/w10 成本探针(e^{tG}v = 0.06 s / 0.245 s) |
| `round1/scripts/probe/phase_c_quotient_factor.py` | full≡quotient 至 2.3e-16;PLAIN 内积偏差 44%(w8) |
| `round1/scripts/probe/c5_l2pi_fibre.py` | L2(π) 下 full≡quotient 至 2.5e-16;C3 纤维 1.4e-15 |
| `round1/scripts/probe/c10_cheap_tests.py` | T1 动态签名饱和于轨道商(750);T2 奇通道可见 c |

验证期间修复了两个脚本的输出路径 off-by-one(`rank_lower_bound.py`、`phase_c_quotient_factor.py`),并清理了误生成的缓存/误放文件。

---

## 项目一 · 预测复杂度探针(C1–C12)

### 已闭合的定理与数值事实

1. **商因子化定理**(#598 Phase C,独立复现 + 精化):P398 全空间与 C2 轨道商的有限时域任务奇异谱在宽度 4–8 上一致至 ~1e-16;商上误用普通欧氏内积会改变谱 37–61%(内积诊断)。
2. **有限群响应选择规则**(#601 Q3 空缺):任意有限群 K 下,源/读出/微扰的表示论张量积判据给出一阶(及 ℓ 阶)响应的精确选择零;`偶函数=轨道常值函数` 是 C2 特例,一般群需 isotypic 分解。C3 合成对象数值验证。
3. **实验语言形式化**:静态/测试语言、类数 κ 与响应秩 r 分离;P398 字典静态签名类 w4–8(w8:D0=32→D1=156→D2=209 vs 轨道 750 vs 微观 1430);完整声明字典在 w≥5 无法分辨 C2 轨道。
4. **阈值 no-go 定理**:两副本 C2 族上,偶任务数据对隐藏跨副本耦合 c 不变至 1e-16,而整链谱隙随 c 变化、c=0 处隐藏扇区丧失遍历性——有限任务无法约束隐藏扇区临界点位置。
5. **#549 根计数代数**:单根分叉语言对类坐标仿射(秩≤2);d 个独立根给 d 次多项式;深度 d 封顶秩 d+1。
6. **复杂度阶梯 v1**:全部箭头标 P(证)/F(反例)/C(条件)/U(未知),附反例档案 X1–X9。
7. **可辨识性目录**:不变(选择零/商恒等式/扇区惰性)vs 随字典/内积/语言变化的量。
8. **C10 廉价测试实测**:T1(P398 偶语言动态分辨饱和于轨道商)、T2(标记通道可见隐藏耦合)、T3(#549 协议代数)。

### 交付文件

- `notes/` —— 17 篇研究笔记(含总评 `probe-mission-update`、阶梯 `probe-complexity-ladder-v1`、猜想排名、语言形式化、no-go、选择规则正式表述等)
- `scripts/probe/` —— 10 个确定性脚本(含共享库 `p398core.py`)
- `results/` —— 8 组结果 JSON
- `PLAN-12-turns.md` / `PROGRESS.md` —— 执行计划(全状态)与进度日志
- `PR-BODY.md` / `COMMENT-599-BODY.md` —— 提交用 PR/issue 正文草稿
- `matching-one-probe-full-20260906.bundle` —— 完整 git 分支包(拿到 GitHub 写权限后可直接 push/开 PR)

> 说明:C6(#582 判别器实证)因需要仓库 #582 档案原始数据、无授权而挂起,判别器协议与决策判据已在笔记中备好。

---

## 项目二 · 切割网络 / 分支预测程序(Program A–J)

### 已闭合的定理(✅)

| Program | 成果 |
|---|---|
| **B** 类计数→秩下界 | 构造 k+1 个实验,响应矩阵**精确秩 k+1** |
| **C** 深度 vs 秩 | `r_d(k)=min(d+1,k+1)`——固定深度秩有界而类数无界 |
| **D** 无有界半径商 | 任意固定半径 r 的终端邻域摘要都不足(把 #550 的 r=1 证书升级到 ∀r) |
| **E** 正 vs 有符号 | slack 矩阵族:普通秩≤3、非负秩无界(正方形 n=4 自证分离) |
| **F** 区分维度 | 分离/仿射/线性维度 = `{1, 2, k+1}` 可任意分离 |
| **H** 记忆 vs 分支 | 有界记忆 + 无界分支类(单个 L0 类含 k+1 个分支类) |
| **I** 可靠性矩代数 | 深度 d = 二阶后继矩的 **d 次张量幂**(代数次数,非矩阶) |

### 综合回答(北星问题)

在 #549 族上,"预测状态"是一个按深度索引的**矩状态**:语言 L_d 的最小线性预测状态维数 = `d+1`,且 `1 = r_0 < 2 = r_1 < … < k+1 = r_k` 严格细化、稳定于深度 k。**"状态"不是网络身上的标量,而是一条由实验深度索引的阶梯。**

### 开放项(有明确下一步目标)

- A(一般树组合语义)、D(最粗精确商存在性,正向)、G(一般树自动机 Myhill–Nerode 对应)、J(串联/粘合同余)、E 的 native 嵌入(把 slack 结构实例化进切割网络分支代数)。

### 交付文件

- `notes/` —— 7 篇研究笔记(秩下界、区分维度、矩代数、语言/商/组合、无有界半径商、slack 分离、综合层级)
- `scripts/` —— 3 个确定性脚本
- `results/` —— 3 组结果 JSON

---

## 文件结构

```
matching-one-probe-complete-20260906/
├── README.md                         ← 本文件
├── matching-one-round1/              ← 项目一(预测复杂度探针 C1–C12)
│   ├── notes/  scripts/probe/  results/  PLAN-12-turns.md  PROGRESS.md
│   ├── PR-BODY.md  COMMENT-599-BODY.md
│   └── matching-one-probe-full-20260906.bundle
└── matching-one-cutnetwork-round1/   ← 项目二(切割网络/分支预测 A–J)
    ├── notes/  scripts/  results/  README.md
```

## 如何复现

所有脚本为纯 Python(标准库 + numpy),无采样、确定性、可复现:

```bash
# 项目二(切割网络)
cd matching-one-cutnetwork-round1/scripts
python3 rank_lower_bound.py       # 秩升级 + r_d(k) 定理
python3 radius_insufficiency.py   # 无有界半径商
python3 slack_family.py           # 正/有符号秩分离

# 项目一(预测复杂度)
cd matching-one-round1/scripts/probe
python3 verify_gate1_independent.py     # Gate-1 复现
python3 phase_c_quotient_factor.py      # 商因子化 + 内积诊断
python3 group_selection_demo.py         # C3 选择规则
python3 c5_l2pi_fibre.py                # L2(π) 一致 + 纤维实现
python3 c7_no_go.py                     # 阈值 no-go
python3 c10_cheap_tests.py              # T1/T2/T3 廉价测试
```

> 研究笔记以英文书写(匹配仓库语言);本索引与用户沟通以中文。
