# Matching-One cut-network / branching-predictive probe — full round

Programs A–J of the cut-network brief, at their current state of completion.

## Completion matrix

| Program | State | Result |
|---|---|---|
| A 组合实验语言 | 部分 + 开放 | 并联族的精确树语义(L0..Lk 严格细化、稳定于 Lk);顺序共享更新的一般语义待闭 |
| B 类计数 → 秩下界 | ✅ 定理 | 构造 k+1 个实验,响应矩阵精确秩 k+1 |
| C 深度 vs 秩 | ✅ 定理 | `r_d(k) = min(d+1, k+1)`,固定深度秩有界而类数无界 |
| D 切割网络最小商 | ✅ no-go 定理 | **任意固定半径 r 的终端邻域摘要都不足**(#550 的 r=1 升级到 ∀r);coarsest 商本身仍开放 |
| E 正 vs 有符号 | ✅ 分离族 + 必要特征 | 显式 slack 族:普通秩 ≤3、非负秩无界(4@n=4 自证);刻画 #549 为何不分离;native 嵌入 open |
| F 最小区分实验集 | ✅ 精确 | 分离/仿射/线性维度 = `{1, 2, k+1}`;一举证伪"一个见证=秩" |
| G Myhill–Nerode 对应 | 猜想 + 族内成立 | 分支 Hankel 秩 ⇔ 最小线性预测表示(并联族已证) |
| H 记忆 vs 分支 | ✅ 精确 | 单 L0 类含 k+1 个分支类;共享前缀测的是反事实未来相关 |
| I 可靠性代数 | ✅ 精确 | 深度 d = 二阶后继矩的 d 次张量幂(代数次数,非矩阶) |
| J 组合代数 | 部分 + 开放 | 不相交并是同余;串联/粘合待测 |

## The central result (programs B + C)

For the #549 parallel-gadget family under the declared composition of the
published fork:

* `P_j(a)` = success probability of the `j`-group fork = polynomial of exact
  degree `j` in the class coordinate `a` (nonzero leading coefficient);
* response matrix with columns `E_0..E_k` has **exact rank k+1** (major
  outcome A);
* `r_d(k) = min(d+1, k+1)` — fixed depth has bounded rank `d+1` while classes
  grow as `k+1` (major outcome B): **depth is the complexity resource.**

## Files

```
notes/synthesis-predictive-state-hierarchy-…           north-star answer (state ladder)
notes/branching-hankel-rank-lower-bound-20260906.md    rank upgrade + r_d(k) theorem
notes/no-bounded-radius-quotient-20260906.md           Program D (∀r no-go)
notes/distinguishing-dimensions-and-positive-rank-…    Program F + E
notes/memory-branching-and-moment-algebra-…            Program H + I
notes/language-quotient-nerode-composition-…           Program A + D + G + J
scripts/rank_lower_bound.py                             exact Fraction evaluator
scripts/radius_insufficiency.py                         bounded-radius no-go evaluator
results/rank-vs-depth/latest.json                       k=2..8 degrees and ranks
results/radius-insufficiency/latest.json                S(z)^k vs fork gap
```

Reproduce everything with `python3 scripts/rank_lower_bound.py` (stdlib only).

## Boundary

Exact statements on the declared composition of the published #549 fork
protocol (per-group independence from future-vertex disjointness).  The lift
to the full N16 network needs the N16 site data and is not claimed.  Ordinary
real rank; nonnegative rank equals it on this family and is open elsewhere.
