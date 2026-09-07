# 集中理论分析 · bounded task rank vs thermodynamic threshold

对 `LightChainr/Matching-One` issue #594 Q3 的正式理论回答。

**一句话结论(verdict)**:`BOUNDED_TASK_RANK_HAS_NO_THRESHOLD_CONTENT_WITHOUT_COMPLETENESS`
——有限时域任务阶有界,在不加"临界扇区完备性(CSC)"假设时,对任何热力学阈值 `p_c` **零约束**。

## 交付文件

| 文件 | 内容 |
|---|---|
| `notes/bounded-task-rank-threshold-no-go-20260907.md` | **主文档**:setup/定义、Main verdict、no-go 定理(完整证明)、direct-sum 是否平凡的三个层次、逃逸条件 A–E 必要性判定、P398 后果、square-site 后果、#609 companion |
| `scripts/no_go_theorem.py` | Theorem 1 的精确构造 + witness(两个 family 相同响应、不同 threshold,全部可复现) |
| `scripts/irreducible_masking.py` | 第二目标:irreducible 下 slowest-mode 精确屏蔽(fine-tuned)的机制验证 |
| `results/no-go/latest.json` | 谱隙闭式、witness、task-order 有界性的数值 |
| `results/irreducible-masking/latest.json` | irreducible 精确屏蔽的 partial-fraction 系数 |

## 核心成果(可引用层次)

1. **Proposition 1(Kalman)**:隐藏扇区 ⟺ 落在不可控不可观子空间。
2. **Theorem 1(no-go, direct sum)**:存在两个解析 Markov 族,finite-horizon 响应对所有 `L,p,t≤T` 完全相同,但 closing-gap threshold 分别在 `1/2` 和 `1/3`。`sup_L r_task = 2` 有界而 `p_c` 不同。
   - 隐藏扇区 = 有偏近邻游走,谱隙精确闭式 `γ_L = a+b−2√(ab)cos(π/L)`,参数化 `a=e^{p−p_c}, b=e^{−(p−p_c)}` 得 `γ_∞ = 4sinh²((p−p_c)/2)`,singularity 只在 `L→∞`。
3. **§3 分层**:direct sum(干净定理)→ 对称保护暗扇区(P398 reflection-odd,非 block-diagonal 但 reducible)→ irreducible(fine-tuned 单模正交屏蔽,generic 下不成立)。弱耦合在 irreducible/local/positive 下给 balanced-rank 版的 no-go。
4. **§4 逃逸条件**:唯一**必要**条件是 A(critical-sector completeness);B/E 是 A 的放大器,C 是更强的充分条件,D 是独立于 rank 的机制(observable 符号变号定理)。
5. **§5 P398 后果**:`r_balanced≈3–4` 只说明 declared task 低有效 I/O 复杂度,不说明低维 critical state;给出可直接入论文的 scope statement,并区分 P398 对 threshold estimation(无直接贡献)vs observability/state-semantics/experiment-design/selection-rules(有用)。

## 复现

```bash
cd scripts
python3 no_go_theorem.py        # Theorem 1 + witness
python3 irreducible_masking.py  # 第二目标 irreducible 屏蔽
```

全部纯 numpy、确定性、无采样、无 Monte Carlo、无 width 9/10。
