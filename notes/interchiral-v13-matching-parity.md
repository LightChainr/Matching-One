# Interchiral V_<1,2> / V_<1,3> matching parity on the critical Potts branch

Status: corrected selection-rule note.  The previous version used the dual Kac branch for the label `V_<1,3>`; repository history preserves that error.  Use this version for all future interpretation.

## 1. Correct critical-Potts Kac convention

For the critical Potts branch at percolation,

\[
Q=1,\qquad \beta^2=2/3,
\]

and

\[
h_{r,s}=\frac{(2r-3s)^2-1}{24}.
\]

Therefore

\[
\boxed{V_{\langle1,2\rangle}: h=\bar h=5/8,\quad x_t=5/4,}
\]

which is the thermal primary, and

\[
\boxed{V_{\langle1,3\rangle}: h=\bar h=2,\quad x=4.}
\]

The Potts space of states contains all `R_<1,s>⊗[]` in the trivial `S_Q` representation.  Moreover `V_<1,2>` is the singlet degenerate generator of the Potts interchiral algebra.

This is a major simplification: the thermal field itself is the interchiral generator.

## 2. Matching parity of the thermal generator

The normalized matching/complement coordinate reverses the thermal coupling,

\[
t\mapsto-t.
\]

Under the working scaling-limit involution `T`, this gives

\[
\boxed{T V_{\langle1,2\rangle}=-V_{\langle1,2\rangle}.}
\]

If `T` commutes with Virasoro generators, all ordinary descendants of the thermal module inherit this odd parity.  In particular the surviving level-4 spin-4 `h=5/8` quasiprimary has

\[
x=5/4+4=21/4,\qquad \eta=-1,
\]

which is exactly the parity/dimension required by the observed central matching residual.

## 3. Self-fusion forces V_<1,3> even

The degenerate fusion includes

\[
V_{\langle1,2\rangle}\times V_{\langle1,2\rangle}
\simeq V_{\langle1,1\rangle}+V_{\langle1,3\rangle}.
\]

If matching is an OPE/fusion automorphism, the product of two odd thermal generators is even.  The identity is even, hence—absent an equal-quantum-number channel rotation—the `V_<1,3>` scalar is also matching-even:

\[
\boxed{T V_{\langle1,3\rangle}=+V_{\langle1,3\rangle}.}
\]

Thus the second-thermal scalar at

\[
\boxed{x=4}
\]

is a concrete matching-even spin-0 correction field.

This is precisely the kind of scalar `S0` whose possible mixing with the leading odd spin-4 `T4` is discussed in issue #58.  A product `T4*S0` remains matching-odd/H4 and carries a relative `q=2`, hence can produce a post-leading accelerated-root contribution near `w_ann=6` if its nonlinear response coefficient is nonzero.

## 4. Distinguish the two dimension-4 sectors

At radial exponent `L^-2`, the current theory can contain at least two conceptually different singlet corrections:

### Scalar second-thermal sector

\[
V_{\langle1,3\rangle},\qquad x=4,\ s=0,\ \eta=+1.
\]

This can contribute to orientation-independent finite-size drift and logarithmic/mixed corrections.

### Spin-4 anisotropy sector

A dimension-4 spin-4 descendant/quasiprimary compatible with square `C4` symmetry, empirically probed by

\[
P_4[S]\sim N^{-1}.
\]

Its precise percolation operator identity remains open.  Ising `T^2+\bar T^2` is an analogy, not a proof.

Feng--Deng--Blote's simultaneous power and logarithmic `L^-2` corrections should therefore not be collapsed into one operator label: scalar/logarithmic and angular spin-4 components can coexist at the same radial exponent.

## 5. Recursive interchiral parity

The fusion rule

\[
V_{\langle1,2\rangle}\times V_{\langle1,s\rangle}
\simeq V_{\langle1,s-1\rangle}+V_{\langle1,s+1\rangle}
\]

combined with `eta_1=+1` and `eta_2=-1` gives, under the same OPE-automorphism/no-channel-rotation assumption,

\[
\boxed{\eta_s=(-1)^{s-1}.}
\]

So odd `s` diagonal singlets are matching-even and even `s` are matching-odd.

This is a testable interchiral selection rule, not merely a labeling convention.

## 6. Correction to the invalid N^-4/3 competitor

The previous exploratory note treated `V_<1,3>` as `h=1/3,x=2/3` and constructed an `x=14/3` level-4 spin-4 competitor.  That used the wrong Kac branch and is invalid for critical Potts percolation.

Consequently:

- `predictions/v13_spin4_parity_failure_competitor_20260828.yaml` is now an explicit invalidation tombstone and must not be scored;
- the erroneous `h=1/3` checker was deleted;
- issue #43 should not include an `N^-4/3` V13 target.

The correct role of `V_<1,3>` is the matching-even scalar `x=4` sector described above.

## 7. Current matching-parity hierarchy

A more coherent low-dimension picture is now:

```text
V_<1,1> identity                 : even
V_<1,2> thermal x=5/4            : odd
V_<1,3> second thermal scalar x=4: even
thermal level-4 spin4 x=21/4    : odd (inherits V12 parity)
```

Together with a separate matching-even dimension-4 spin-4 anisotropy field, this naturally produces the two leading center sectors already seen numerically and the q=2/q=4 composite mechanisms of #58.

## 8. Claim boundary

The Kac labels/dimensions and Potts singlet/interchiral content are standard critical-branch facts.  The **matching parity** assignments follow from the working assumption that matching induces an OPE/interchiral automorphism that reverses the thermal generator and commutes with Virasoro descendants.

Prospective derivative parity (#48), self-dual/self-matching controls (#42/#44), and a future explicit FK/Temperley-Lieb realization of the matching map are the appropriate falsifiers.

References: Jacobsen--Ribault--Saleur, arXiv:2208.14298, Eqs. (1.3)--(1.4), (3.10)--(3.15); Feng--Deng--Blote, arXiv:0901.1370.
