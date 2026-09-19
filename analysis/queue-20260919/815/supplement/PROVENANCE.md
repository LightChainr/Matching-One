# PROVENANCE — #815 supplementary execution

日期2026-09-19；同一分支追加，不覆盖初次交付。
追加基线 #823 commit `476cc029f17c3611db280aff4823a45ff40654e5`。
原参考文件：PR773 head `9ef22d1cc531c09dedf9b1195ab6aad8ef79668f` / `scripts/rank_jump_two_spine_control.py`；Git blob `9fc042535aa0707651a5e38a2e2121669d99d59a`。本目录 reference773.py 字节的Git SHA1与之相同。
正式命令 `python run_exhaustive.py`；结果与执行边界见NOTE。#775内核未获得，不能借本补件声称已复跑#775。

## SHA-256

```text
19ad7eb947cde083a40404548fdb10168f4556208a288d228cb7e0d02e359a20  NOTE.md
b6a4c719bb8d2cc8476a1f8a1c5f0d647351413b549676694fe211096573ef8a  reference773.py
a0421beb3a3490705b41504a4f777318683d931f0a49d7415ce46edce6db8bc4  exhaustive_check.cpp
278e498357a5bf1539c03cf0711f95a5482a162bd00e0d8a1b3ce8d1009cc3a2  run_exhaustive.py
5f490af563633c3e89c0db86756393a5170a02d7faa2c7c173fdad50cb15551e  exhaustive-result.json
```

清单不包括自身；编译临时文件和Python缓存不是交付数据。

## 本轮执行范围

只复跑本票有界整数内核与臂检查、原Python参考小尺寸交叉对照；没有父PR验证继承、全仓CI或长期生产。
