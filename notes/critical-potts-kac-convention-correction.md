# Critical Potts Kac convention correction

Status: mandatory theory correction before any further operator ranking.

## 1. Source of the error

Some exploratory notes temporarily used the dual Kac polynomial

`((3r-2s)^2-1)/24`

while interpreting labels as if they belonged to the **critical Potts branch**.  For the Potts CFT convention used by Jacobsen--Ribault--Saleur, the critical branch has

\[
Q=4\cos^2(\pi\beta^2),\qquad 1/2\le\beta^2\le1.
\]

At percolation `Q=1`, therefore

\[
\beta^2=2/3.
\]

With their momentum convention this gives

\[
\boxed{h_{r,s}=\frac{(2r-3s)^2-1}{24}}.
\]

The dual polynomial with `3r-2s` corresponds to exchanging the Coulomb-gas/Kac branch convention.  Virasoro statements that depended **only on the numerical value h=5/8** remain valid; label-based operator rankings must use the corrected formula.

## 2. Thermal field label

On the critical Potts branch,

\[
h_{1,2}=\frac{(2-6)^2-1}{24}=\frac58,
\]

so the thermal primary is

\[
\boxed{V_{\langle1,2\rangle},\qquad x_t=2h=5/4.}
\]

This is especially important because `V_<1,2>` is precisely the singlet degenerate field that generates the Potts interchiral algebra in the JRS description.

The exact Virasoro level-2 null relation and the level-4 `h=5/8` quasiprimary calculations already committed in the repository are unchanged, since they used `c=0,h=5/8` directly rather than the wrong Kac label.

## 3. Second thermal scalar

The next diagonal singlet is

\[
h_{1,3}=\frac{(2-9)^2-1}{24}=2,
\]

hence

\[
\boxed{x_{1,3}=4.}
\]

This is the familiar second-thermal radial exponent `X_t2=4`, not a relevant `x=2/3` primary.

Moreover the degenerate fusion

\[
V_{\langle1,2\rangle}\times V_{\langle1,2\rangle}
\supset V_{\langle1,1\rangle}+V_{\langle1,3\rangle}
\]

means that if matching acts as an OPE automorphism and the thermal generator `V_<1,2>` is matching-odd, then

\[
\boxed{V_{\langle1,3\rangle}\text{ is matching-even}.}
\]

This provides a concrete candidate for the matching-even **scalar** `x=4` sector invoked in issue #58.  It is distinct from any matching-even **spin-4** dimension-4 descendant responsible for angular anisotropy; both can occur at the same radial power `L^-2`.

## 4. Invalidated exploratory V13 model

The temporary claim

```text
V_<1,3>: h=1/3, x=2/3
level-4 spin4: x=14/3
DeltaM ~ N^-4/3
```

was a branch-label error and is invalid for the critical Potts CFT studied here.

Therefore:

- `predictions/v13_spin4_parity_failure_competitor_20260828.yaml` must be marked invalid and must **not** be scored on issue #43;
- the associated `h=1/3` checker is not a critical-Potts V13 checker;
- issue #43 reverts to its original `13/8` primary, then the separately frozen formal `x=17/4` / `9/8` adversarial model, then zero before flexible fits.

Repository history should preserve the erroneous exploratory commit as provenance; execution artifacts must explicitly mark it superseded rather than silently pretending it never existed.

## 5. Non-diagonal spin-4 primary sequence also changes numerically

For a critical-branch non-diagonal module `(h_{r,s},h_{r,-s})`, conformal spin remains

\[
h_{r,s}-h_{r,-s}=-rs.
\]

Setting `rs=4`, the scaling dimension is now

\[
\boxed{x=\frac{4r^2+9s^2-1}{12}}.
\]

Thus the first rows are

```text
r=2, s=2     : x=17/4, Xi=[2]       (non-singlet)
r=3, s=4/3   : x=17/4, Xi=[21]      (non-singlet)
r=4, s=1     : x=6,    Xi(4,0) contains []
```

The formal `x=17/4` adversarial law itself is unchanged for `W(2,2)`, but the statement that the first non-diagonal primary singlet occurs at `x=49/4` was wrong: under the correct critical branch it occurs already at `x=6`.  This is still above the thermal level-4 descendant `x=21/4=5.25`.

## 6. Corrected low-dimensional singlet spin-4 hierarchy

The current ordinary generic-Q candidates should be organized as follows:

1. dimension-4 sector:
   - matching-even scalar `V_<1,3>` (`X_t2=4` candidate);
   - matching-even spin-4 identity/interchiral descendant candidate;
2. thermal-family level-4 spin-4 descendant:
   - `h=5/8`, total bulk `x=21/4`, matching-odd under the commuting-involution hypothesis;
3. first simple non-diagonal primary spin-4 sector whose `Xi` contains the singlet:
   - `W(4,1)`, `x=6`.

This corrected ordering strengthens rather than weakens the empirical `x=21/4` mechanism: the dangerous diagonal `V_<1,3>` field is naturally matching-even, while the first lower formal non-diagonal primary `x=17/4` remains non-singlet.

## 7. References

- J. L. Jacobsen, S. Ribault, H. Saleur, *Spaces of states of the two-dimensional O(n) and Potts models*, arXiv:2208.14298.  See the critical Potts relation `Q=4 cos^2(pi beta^2)`, the Kac momentum convention, spectrum Eq. (3.10), and interchiral generator `V_<1,2>`.
- X. Feng, Y. Deng, H. W. J. Blote, arXiv:0901.1370, for the percolation `X_t2=4` correction sector.
