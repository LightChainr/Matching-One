# note-reproduce —— 第一层：逐字复现 #771 包的自述

机器 `DevEnvC_TVVfoB`（账号2，16 vCPU，Python 3.9.9，`PYTHONPATH=/workspace/mo/compat`）。
工作目录 `/workspace/genealogy/`。**本机 Mac 未做任何计算。**
被核验件：`genealogy-noise-20260914/matching_one_genealogy_noise_20260914/repo/`（原样上传）。

## 1. 18 项测试

命令（按包内 README 口径，脚本 `sys.path` 自己插入 `scripts/`）：

```
cd /workspace/genealogy/in/repo
PYTHONPATH=/workspace/mo/compat PY39COMPAT_MARKER=1 \
  python3 -m unittest discover -s tests -p test_coalescent_bulk_clock.py -v
```

结果：**`Ran 18 tests in 2.845s` / `OK`，18/18 通过。**
测试名与 `TEST_OUTPUT.txt` 逐条一致（18 行，`test_01…test_18`，顺序相同）；
唯一差别是计时（本地 1.080 s vs 本机 2.845 s），非判据。

（垫片生效证据：每步都打印
`py39compat: int.bit_count backport active (route=installed, stdlib-only, selftest=OK)`，
自检 `(255).bit_count()==8`。）

结论：**通过**。证据等级：有限枚举事实 + 有限测试（不构成渐近证明）。

## 2. `--output` 再生 vs 远端结果 JSON

接口核对（先做）：脚本用 `argparse`，只有 `--output`（`type=Path`），**没有硬编码绝对路径**；
`__main__` → `main()`；写盘用 `json.dumps(report(), ensure_ascii=False, sort_keys=True, indent=2, default=encode)`，
即**排版（indent=2）**；无 `--help` 障碍（`-h` 正常）。

比对（脚本 `scripts/compare_regen.py`，独立于包）：

| 比对口径 | 结果 |
|---|---|
| (A) 字面字节：脚本 stdout/`--output` 产物 vs 提交件 | **不相同**（11747 B vs 7482 B）|
| (B) 规范化字节 `dumps(load(...), sort_keys=True, separators=(',',':'))` | **逐字节相同** ✓ |
| (C) 解析后结构相等 | **相同** ✓ |

规范化后两边 canonical sha256 都是
`d93c99e5728636b51b93846e9296b7b06b140835125721dd2472c96765ad49aa`（各 7481 B）。

⚠️ 口径澄清：`EXECUTION.json` 的两个字段名是
`json_structure_regenerated_identically` 与 **`json_compact_canonical_bytes_regenerated_identically`**
——措辞**明确限定在 compact/canonical**，与 (B)(C) 相符，**没有过度宣称**。
但要知道：脚本本身写的是**排版 JSON**，因此**它无法在字面字节层面再生成提交件**；
"逐字节相同"只在规范化口径下成立。规格 §2.2 担心的"排版 vs 紧凑"在这里是**已声明**的，
不是隐瞒。（对照：#780 包的 `EXECUTION.json` 主动写明 "Exact identity is after JSON compact
normalization, not literal pretty-vs-compact bytes." —— 同一种做法。）

## 3. `EXECUTION.json` 四个计数能否从产物复算

四个数**都能**从再生的结果 JSON 直接复算（`report()` 的输出里就有）：

| 字段 | 声称值 | 复算来源 | 复算值 |
|---|---|---|---|
| `circle_combinatorial_realizations` | 89438 | `Σ_n genealogy[n].realizations` = Σ_{n=2..6}(n−1)!·n! = 2+12+144+2880+86400 | **89438** ✓ |
| `full_history_transition_checks` | 4047 | `Σ_n genealogy[n].full_history_transition_checks` = 1+4+25+251+3766 | **4047** ✓ |
| `partition_probability_checks` | 277 | `Σ_n genealogy[n].eppf_equalities` = 2+5+15+52+203 | **277** ✓ |
| `four_coordinate_generator_checks` | 768 | `report()['four_coordinate_generator_equalities']` = 4(η)×4(a)×4(b)×3(m)×4(n) | **768** ✓ |

⇒ **不是"只出现在 EXECUTION.json、无法复算"的孤立数**；四者都是产物字段的直接聚合。

⚠️ 一处措辞注意：`four_coordinate_generator_checks=768` 是**网格检查条数**（768 条恒等式），
不是"四坐标"的维度数；`768` 与 `4×4×4×3×4` 的分解见于包内脚本 `closure_stationarity_checks()`，
`EXECUTION.json` 未给分解（可接受，但读者易误读为"四维生成的 768 个分量"）。

## 4. 未发现的问题（第一层）

- 脚本未硬编码路径；未在 Python 3.9 触发语法错误（`from __future__ import annotations` + 无 3.10 语法）。
- 18 测试与 `TEST_OUTPUT.txt` 逐条匹配；`Ran 18 tests ... OK` 属实。
- 四个计数均可复算。

## 5. 复现命令与产物

云上：`/workspace/genealogy/{scripts,out}/`
- `scripts/compare_regen.py` → `out/layer1-compare.json`，`out/regen.json`
本机：`genealogy-out/{regen.json,layer1-compare.json,scripts/compare_regen.py}`
耗时：测试 2.845 s；再生+比对 < 1 s；pip 装 numpy 2.0.2 / sympy 1.14.0 约 11 min（一次）。
