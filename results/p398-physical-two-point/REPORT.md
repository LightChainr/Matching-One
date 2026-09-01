# P398：width-5 真实正权两点传播

## 科学结果

同一局部 A/L 读出在 width-5 **不再闭合为二维传播**。
精确 Krylov 秩（最大传播次数 0…7）为 `[2, 4, 6, 8, 8, 8, 8, 8]`；只用正间距
`C(1)…C(7)` 组成的 4×4 个 2×2 block Hankel 矩阵秩为 **8**。
完整 charge-one 块的维数为 8，特征多项式是否 square-free：是（与导数 gcd 为常数）。

这是使用已有物理接口得到的有限宽度科学结果，不是再造 rooted closure 工具。
正间距 Hankel 满秩且特征多项式 square-free，因此两个已指定读出的整个矩阵序列
共同探测到 8 个不同的有限宽度传播本征成分；这不等于 8 个连续极限场。

## 固定的测度与自然延拓

- Q=1 独立 square-bond cylinder，width=5，h=v=1/2。
- 正概率的 42 个 circular-noncrossing 过去连通状态，T=H V，先 vertical 再 horizontal，
  在 H 后读取；每一步是 10 个独立 Bernoulli bond bits 的精确求和，没有 Monte Carlo。
- `A=Σ_j ζ_5^j I(j~j+1)`；
  `L=Σ_j ζ_5^j I(j 为 singleton)`，j=0…4。
  两个指示函数原样保留，仅以当前圆周群的第一个非平凡 Fourier character 加权。
  width-4 的 ζ₄=i 正好还原已有 A/L；C4 和 C5 的 charge-one 并非同一群表示。
- 8 个长度 5 的轨道贡献 charge-one；两个旋转不动状态不贡献该 charge。
  平稳分布满支撑、两个读出的均值精确为零。
- 全过程在 Q(ζ₅) 中精确计算。JSON 的每个场元素依次保存
  `(1,ζ₅,ζ₅²,ζ₅³)` 的四个有理系数；小数仅供阅读。
- 这是普通的正权 42-state 模型，不把旧 signed 23-state retained-mark module
  当作随机 transfer，也没有完成旧模块的 full-Q 物理提升。

`C_ab(d)=E[O_a(X_0) conjugate(O_b(X_d))]` 是 connected 两点量。
正间距 Hankel 秩直接来自两点序列，因而不会把 8 维环境空间自动当成可观测秩。
自然延拓本来不保证二维闭合；本计算检验的是该特例能否随宽度保持。

## 实际两点读数

`C(0)`（行列次序 A,L）：

```text
  1.01777621707-1.38777878078e-17i | -0.656657811235+0.477089826209i
  -0.656657811235-0.477089826209i | 0.779103140598-1.38777878078e-17i
```

`C(1)`：

```text
  0.10368930974-8.67361737988e-19i | -0.0649014407291+0.0471536568185i
  -0.0681370535473-0.0495044671351i | 0.0678619201577-1.73472347598e-18i
```

`C(8)`：

```text
  1.86711816298e-08-4.13590306277e-25i | -1.14178174686e-08+8.29552996793e-09i
  -1.21941983252e-08-8.85960367822e-09i | 1.13933838478e-08-2.06795153138e-25i
```

`det C(0) ≈ 0.134138463813`。
完整 d=0…8 矩阵、正间距 Hankel、精确 transfer、平稳分布和特征多项式都在
[latest.json](latest.json)，无需凭幅度拟合来认定传播秩。

## 消费的已有结果与解释边界

1. `e38fe7634354b0cb2201fa55fd9b4d37ccedeef2`（`branch_only`），
   `notes/p398-positive-cylinder-propagation.md` 已完成 width-4 的正权两点 A/L 矩阵，
   两个传播本征值为 `(3±√5)/64`。本次不重算它。
2. `b35e100a3903c706dceba57c4667386eb4510ac3`（`branch_only`），
   `notes/p398-anisotropic-cylinder-spectrum.md` 已处理 width-4 完整 h/v 正权族、
   边界和 signed Jordan 特例。本次不再做 anisotropy 扫描。
3. 本结果只回答 width-5、h=v=1/2 的上述 first-character 两点传播。
   不把非二维结果写成 width-4 失效，不把 simple spectrum 泛化为所有宽度无 Jordan，
   不据此识别 Matching One 的 norm-4 场或 E_top 的微观能量投影。
4. 真正新增的信息是：固定同样的微观读出后，宽度变化是否迫使有限传播模型
   增加成分；无需增加任意新 mark 来制造该结果。

## 复现

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p398_physical_two_point.py
```

本次实际运行：Python 3.11.15 / arm64，
0.169 秒，单一 width、单一参数点、零 MC 样本；没有全仓库测试。
checkout HEAD：`bcfacf571393ac0975bfd2aa3300ad37903ec2fa`（输出为本 Draft 工作树新增结果）。
脚本 SHA256：`2885571dae7c88ef393203cf3ec3e52d7988e5911e5bf44e78c5e8275478ba78`。
JSON SHA256：`f5e33ec1574607c7acab9fbb19c6941bb473cabffbb82396d606c482f824aa81`。
