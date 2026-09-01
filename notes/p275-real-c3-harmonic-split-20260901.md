# Primitive real-C3 harmonic split after the paired Gaussian gate

Status: mechanism interpretation of the frozen N65 result; local field identity
remains open.

The complex primitive character has one nontrivial `C3` charge.  In the spin
normalization used by the project it only resolves

\[
s\equiv4\pmod 6,
\]

so `s=+4` and `s=-8` lie in the same character tower.  Real probabilities make
the two nontrivial characters conjugate, which is why a single geometry cannot
separate the aliases.

The paired physical rotation does separate them.  A global conjugation sends
`+4` to `-4`, a cyclic relabeling multiplies both geometries by one common
cube root of unity, and a signed-real gain absorbs only a phase of `0` or
`pi`.  None changes a continuous `+4 delta` slope into `-8 delta`.  At the
frozen N65 angle the two model lines are almost orthogonal, and the data select
the latter.  The result is therefore not a notation or conjugation artifact.

The minimal joint description is

\[
z_{\rm prim}(\theta)
 =a_4(\tau,N)e^{4i\theta}+b_{-8}(\tau,N)e^{-8i\theta}+\cdots,
\qquad
A_{\rm top}(\theta)=c_4(\tau,N)e^{4i\theta}+\cdots .
\]

The paired `tau=i`, N65 result says that `b_-8` dominates this primitive
projection, while the quotient-prism/global result says that `c_4` dominates
`A_top`.  Earlier norm-2 cover chains selected the negative rank-4 transfer at
other parent moduli.  The two statements can coexist if the primitive
coefficients are modulus-dependent and the square point suppresses or
overwhelms `a_4`; they should not be pooled as votes on one amplitude.

This does **not** require a local spin-8 primary.  A local H8 contribution with
the older empirical `N^-1` radial law conflicts with `x>=|s|`, which would make
a first-order local spin-8 term decay no slower than `N^-3`.  The economical
candidate is a topological observer form factor: the same underlying rank-4
sector dressed by a `C3`-scalar, spin-12 homology/readout tensor.  That changes
the measured harmonic from `+4` to `-8` without creating a new local primary.

## Two prospective decisions

1. **Persistence at fixed modulus and angle.**  The N130 pair
   `(11+3i, 9+7i)` has the same `delta=atan2(5,12)` and is a common Gaussian
   dilation/rotation of the N65 pair.  Persistence of the H8 phase there rules
   out a one-size phase accident; an H4 return marks a finite-size crossover.
2. **Convention-invariant angle test.**  For the N145 pair `(12+i,9+8i)`,
   `exp(i delta)=(4+3i)/5`.  With the same complex character define

   \[
   X={\Re[(z_2\bar z_1)^2]\over |z_2z_1|^2}.
   \]

   Squaring removes the signed-gain sign, and `X` is invariant under common
   character phase and conjugation.  Passive H4/conjugate-H4 predicts
   `X=+0.42197248`; a true `-8` observer harmonic predicts
   `X=-0.643878452245`.  Only the existing complex C3 coordinate and its paired
   covariance are needed.
