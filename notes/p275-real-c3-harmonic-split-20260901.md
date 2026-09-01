# Primitive real-C3 harmonic split after the paired Gaussian gate

Status: updated after the held-out N145 three-way gate.

The complex primitive character has one nontrivial `C3` charge.  In the spin
normalization used by the project it only resolves

\[
s\equiv4\pmod 6,
\]

so `s=+4` and `s=-8` lie in the same character tower.  Real probabilities make
the two nontrivial characters conjugate, which is why a single geometry cannot
separate the aliases.

The first paired physical rotation separated H4 from its H8 alias and selected
the H8 line within that frozen two-model set.  Post reveal, however, the same
N65 data also pass a signed-scalar H0 line (`p=0.2505`), because the chosen
angle puts `-8 delta` within one degree of `pi`.  N65 therefore rejects H4 but
does not identify H8.

The minimal joint description is

\[
z_{\rm prim}(\theta)
 =a_0(\tau,N)+a_4(\tau,N)e^{4i\theta}
  +b_{-8}(\tau,N)e^{-8i\theta}+\cdots,
\qquad
A_{\rm top}(\theta)=c_4(\tau,N)e^{4i\theta}+\cdots .
\]

The held-out `tau=i`, N145 result below shows that `a_0` dominates this
primitive projection, while the quotient-prism/global result says that `c_4`
dominates `A_top`.  Earlier norm-2 cover chains selected a negative transfer
at other parent moduli, but their simultaneous size change made the radial
gain sign an assumption.  They should not be pooled as H4 votes without an
independent positivity argument.

This stops the local spin-8 interpretation.  It also rules out passive H4
conjugation as the dominant transport at the square point.  The surviving
object is instead a nontrivial homology character whose response to physical
embedding rotation is scalar: topology is carried by the observer label, not
by local geometric spin.

## Held-out resolution

The N130 repetition was not run because the original angle leaves H8 almost
aliased with H0.  The higher-information N145 pair `(12+i,9+8i)` instead used
`exp(i delta)=(4+3i)/5`.  With the same complex character define

\[
X={\Re[(z_2\bar z_1)^2]\over |z_2z_1|^2}.
\]

Squaring removes the signed-gain sign, and `X` is invariant under common
character phase and conjugation.  Passive H4/conjugate-H4 predicts
`X=+0.42197248`; a true `-8` observer harmonic predicts
`X=-0.643878452245`; an embedding-even/signed-scalar character predicts
`X=+1`.

The held-out result is

```text
p_H0 = 0.968628
p_H4 = 0.000135506
p_H8 = 2.67966e-13
X    = 0.999926 +/- 0.007664.
```

Thus the primitive signal has a nontrivial `C3` homology character but
spin-zero embedding transport.  The N65 H8 selection was the accidental
near-alias of H8 with the real line at that angle.  The local-H8 and passive-H4
interpretations are both stopped.

The surviving form factor is exact: multiplication by the Gaussian unit `i`
acts on period homology as `(m,n)->(-n,m)` and exchanges the unoriented lines
`l0,l1`.  Hence `P_l0=P_l1` and `z=omega(P_l2-P_l0)` for every Gaussian ideal
quotient.  The real amplitude may change sign, but the complex character line
cannot rotate.  The complete proof and stop rule are in
`notes/p275-gaussian-unit-rotation-character-theorem.md`.
