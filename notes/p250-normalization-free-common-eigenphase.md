# P250 normalization-free common eigenphase

## Deduplicated question

The compact 1M charged-cubic archive already established that raw `C113/C122`
closure survives, but the complete cubic vector was compatible with zero and
the archive had no same-operator pair denominators.  The later projective-leg
program established a nontrivial charged propagator and then explored its
71/21/8-type state dimension, projective transfer, multi-state/Hankel and
rank-eight completions.  None of those results is another normalization-free
three-point OPE test.

The fresh projective-leg 10k archive is the one existing block that closes the
required interface.  At `d=1` it has:

- the primitive charged cubics `C113` and `C122` in both hands;
- exact conjugate controls `C244` and `C334`;
- same-operator neutral pairs `G1/G2`;
- 50 shared batches retaining their complete covariance;
- a previously frozen pair/support gate that passes.

`d=2` was excluded before the new contrast because its cubic-support gate had
already failed; `d=3` was a declared tail diagnostic.

## Frozen normalization and contrast

For each hand,

\[
\Omega_{113}=\frac{C_{113}}{G_1\sqrt{G_2}},\qquad
\Omega_{122}=\frac{C_{122}}{\sqrt{G_1}G_2}.
\]

These ratios cancel positive field-amplitude normalizations while retaining
the complex phase in the exact transported Z5 deck basis.  Writing
`A=Omega113` and `B=Omega122`, a common unit-modulus hand eigenphase exists iff
the supported two-vectors obey

\[
A_+B_- - A_-B_+=0,
\]

and

\[
|A_+|^2+|B_+|^2-|A_-|^2-|B_-|^2=0.
\]

The frozen primary is therefore the three-real vector consisting of the real
and imaginary cross-product plus the norm difference.  Every nonlinear ratio
and contrast is rebuilt inside each shared-batch delete-one.

## Result

The eligibility gates reproduce:

- weakest pair denominator: 15.397 sigma;
- eight-real cubic support: chi-square 24.918/8, p=0.001605;
- maximum DFT conjugacy residual: 1.42e-15.

The primary common-unit-eigenphase score is

```text
contrast = (-0.48762, 0.21966, 1.34338)
chi-square = 3.13206 / 3
p = 0.371709
```

Thus the completion **survives** at alpha 0.01.  It is not identified.  The
descriptive common projection is `q=0.2228+0.3686i`, but its amplitude is only
`0.431 +/- 2.536` and its phase `1.027 +/- 5.891` radians.  The present block
therefore cannot resolve a common charged eigenphase even though the joint
cubic support gate is nonzero.

This is a charged projective-leg/OPE-like ratio, not ordinary global `A_top`.
No H4/H8/H12 model was rescored.  The surviving relation does not make the
topological insertion a local primary or its ratio a universal structure
constant.

