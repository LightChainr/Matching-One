# PROVENANCE — #814

本轮新增文件；日期 2026-09-19。

基线：#813 head `11104385b239be4894254772e83e61ceeabae634`，tree `9d3c1b40b001039fe227b92346d880153ab2a838`。PR 以 `research/control-20260919` 为 base，避免携入既有研究分支 diff。

读取来源（Git blob SHA，不是 SHA-256）：
- control 文档 `58adaeb010edf95bdee660f899c305c2aa980f90`。
- RESEARCH-ATLAS `66fce42649987ddcf358ee6accbe54cb0fda26b2`。
- PR771 分支 `analysis/geometric-balance-round2-integrated-20260914`：bridge `f5b3513a7047c8c9a92c435d7b3211c38d373b6e`；thermal Ward `73868e4d7b247d49909c175a41eb12e9a1c4b961`；two-stage `5e1e7ce1ebde3be2e02ea8af816f20734e3d3494`。
- 外部文献的版本/读取段落见 NOTE.md；未取得原始PDF字节，**不伪造 PDF SHA-256**。

命令：`python check.py`。结果：6 项有理数检查通过。

## SHA-256

```text
6055428ebe9fc857277a16ddf35344be0c2e381b4e62347411ce13c8acda9d16  NOTE.md
a6f0ba4b36b0d97fb2c27d28a6eca25daaa1abc18c6a08fbf53dfb6862d242bc  check.py
f1035c8a87fa9b6fb354509e204f1e736a5448f0303aa5be35bb0e8f850a4932  result.json
```

清单不包含自身，避免自指哈希。

## 本轮执行范围

只验证本票有理数记账并记录来源；没有重跑引用数据、父PR验证、转移矩阵或全仓CI。
