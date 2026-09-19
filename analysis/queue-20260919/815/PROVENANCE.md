# PROVENANCE — #815

日期2026-09-19。基线 #813 `11104385b239be4894254772e83e61ceeabae634`；文献/仓库读取范围与缺失路径见 NOTE.md。

执行 `python check.py`：三个矩形共66,064种着色，4/8二择一零违反；错误4/4负控制通过。本次不是完整jump2验收。计时属于本次单次运行，不是内核速度。

## SHA-256

```text
a11eb1795e1a32b5dac846ce3590279991b94dab57e6e53c8d325935e9f99ae3  NOTE.md
35cb77fe238a43edca1817e1346099ab77c6f2c8677e507ccae5bf2f8fc3a60c  check.py
8333df1e46e17b4fc39633ddc8a15b9169b39e04c0eb4a6488e88118060dbfaf  result.json
```

不对来源404文件编造哈希，不自指哈希。

## 本轮执行范围

只有局部矩形邻接回归；全jump2检查未执行，未跑全仓CI或父PR测试。
