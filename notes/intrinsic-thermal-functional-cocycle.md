# Intrinsic thermal-coordinate and full-curve cocycle test

Status: **C0 analysis/design**. This is a zero-new-compute analysis route for existing and already-planned threshold-rank full curves.

Related: #48, #101, #119, #125, #138, #145, #57.

## 1. Why `P4[S']` should not be interpreted in bare `p` coordinates alone

The current prospective state is unusually clean: the new N=185/265 block supports the pure laws for `P4[S]`, `P4[D]`, and `P4[D']`, while `P4[S']` is the unique clear pure-law failure. Both the frozen analytic `1/N` correction and the rank-2/Jordan-log correction survive.

But `S'` means a derivative with respect to the microscopic probability `p`. The conversion from `p` to the finite-size thermal coordinate carries the same thermal metric whose small but highly significant finite-size drift is already resolved by the center-slope data.

Therefore a first separation should be

```text
operator/scaling-function correction
versus
bare-p -> intrinsic-thermal-coordinate correction.
```

Do this before adding another operator to explain `S'`.

## 2. Intrinsic derivative normalization

Let the orientation-averaged matching curve be

\[
\bar M_N(p)=\frac12[M_{N,1}(p)+M_{N,2}(p)],
\]

and let `p0,N` solve

\[
\bar M_N(p_{0,N})=0.
\]

Use the matching curve itself as a local intrinsic thermal coordinate:

\[
u=\bar M_N(p).
\]

At the center,

\[
\frac{\partial}{\partial u}
=\frac{1}{\bar M'_N(p_0)}\frac{\partial}{\partial p}.
\]

For the matching-odd spin-4 block define

\[
\boxed{
J_N
=\frac{P_4[S'_p](p_0)}{\bar M'_N(p_0)}
=\left.\frac{\partial P_4[S]}{\partial u}\right|_{u=0}.
}
\]

If the central `D` and derivative `S'` are Taylor coefficients of one matching-odd scaling function of dimension `x=21/4`, then

\[
P_4[D](0,N)\sim N^{-13/8},
\]

and **after intrinsic derivative normalization**

\[
\boxed{J_N\sim N^{-13/8}}
\]

rather than the bare-p law `P4[S']~N^-5/4`.

The ratio

\[
\boxed{
\Xi_N
=\frac{J_N}{P_4[D](0,N)}
=\frac{P_4[S'_p](0,N)}{\bar M'_N(0)P_4[D](0,N)}
}
\]

cancels the overall lattice coupling of that block. For a single dominant scaling function it should approach the dimensionless Taylor ratio `f'(0)/f(0)`.

This is not asserted to be universal across arbitrary observables without a normalization argument; it is first a stringent **internal closure invariant** of the proposed single matching-odd block.

### Interpretation

- if most of the scaled `S'` drift disappears in `J_N` or `Xi_N`, the anomaly is substantially a finite-size thermal-metric/coordinate effect;
- if the drift survives essentially unchanged, the correction is intrinsic to the matching-odd scaling block and the q=2/Jordan competition remains genuinely operator-level;
- if `Xi_N` cannot be described by one smooth correction while `D` itself remains pure, multiple matching-odd fields or matrix mixing are implicated.

All quantities should be recomputed inside the synchronized delete-one replicate because `p0` and `Mbar'` are estimated from the same full curve.

## 3. Reconstruct the whole matching-odd scaling function

The stronger test uses the full curve instead of a center derivative.

For each frozen intrinsic level `u>0`, solve inside every jackknife replicate

\[
\bar M_N(p_+^u)=+u,
\qquad
\bar M_N(p_-^u)=-u.
\]

For the same-N spin-4 projectors define the thermal-parity pieces

\[
D_e(u,N)
=\frac12\left[P_4[D](p_+^u)+P_4[D](p_-^u)\right],
\]

\[
S_o(u,N)
=\frac12\left[P_4[S](p_+^u)-P_4[S](p_-^u)\right].
\]

For a matching-odd field with underlying scaling function `F`, these are the even/odd thermal pieces. Therefore

\[
\boxed{
T_N(u)=D_e(u,N)+S_o(u,N)
}
\]

reconstructs one branch of the same projected scaling function, while

\[
D_e-S_o
\]

reconstructs the reflected branch.

This is useful because it makes `D(0)` and `S'(0)` two local shadows of one object rather than two independently fitted channels.

## 4. Scale out the candidate x=21/4 radial character

Define

\[
\boxed{
Z_N(u)=N^{13/8}T_N(u).
}
\]

For a pure matching-odd x=21/4 eigenfield,

\[
Z_N(u)=f(u)+o(1).
\]

The two live correction mechanisms predict different multiplicative-scale cocycles but do **not** require us to know the correction scaling function `g(u)`.

### Analytic relative-q=2 correction

If

\[
Z_N(u)=f(u)+N^{-1}g(u)+\cdots,
\]

then for the same parent N

\[
Z_{5N}(u)-Z_N(u)
=\frac85\,[Z_{2N}(u)-Z_N(u)]
\]

at first correction order.

### Rank-2 Jordan correction

If

\[
Z_N(u)=f(u)+\log N\,g(u)+\cdots,
\]

then

\[
Z_{5N}(u)-Z_N(u)
=\frac{\log5}{\log2}\,[Z_{2N}(u)-Z_N(u)].
\]

Thus for every frozen `u` define the functional residual

\[
\boxed{
R_c(u;N)=Z_{5N}(u)-c Z_{2N}(u)+(c-1)Z_N(u).
}
\]

with

```text
q=2 analytic: c = 8/5 = 1.6
rank-2 Jordan: c = log(5)/log(2) = 2.321928094887...
```

The unknown leading function `f(u)` and correction function `g(u)` both cancel. This is a parameter-free **function-level** discriminator.

## 5. Use the already frozen u-grid as one correlated vector

Issue #119 already identifies the low-u levels used by the full-curve pipeline. Preserve the exact pre-existing level grid; do not add points after the norm-5 target reveal.

For each parent lineage form the residual vector

\[
\mathbf R_c(N)
=[R_c(u_0;N),R_c(u_1;N),\ldots].
\]

Score it with the full synchronized covariance across

- all u values;
- `D_e` and `S_o` ingredients;
- parent/norm-2/norm-5 sizes when shared random fields are used.

The important question is not whether one u point has a good chi-square. It is whether the **entire scale-increment function is collinear** with the fixed multiplier predicted by q=2 or Jordan.

The two prospective triples are already determined by #57:

```text
65 -> 130 (norm 2), 325 (norm 5)
85 -> 170 (norm 2), 425 (norm 5)
```

No additional large simulation is required.

## 6. A rank test if both fixed mechanisms fail

Define scale-increment vectors

\[
\Delta_2(u)=Z_{2N}(u)-Z_N(u),
\qquad
\Delta_5(u)=Z_{5N}(u)-Z_N(u).
\]

A **single** correction scaling function, whether ordinary power or rank-2 log, predicts that these two vectors are collinear across `u`.

If no scalar multiplier fits them within covariance, then one correction function is insufficient. This is direct evidence for at least one of

- two ordinary correction fields;
- an ordinary field plus a logarithmic partner;
- finite-size operator mixing;
- intrinsic-coordinate corrections not removed by the `u=Mbar` definition.

Only after rejecting the rank-1 scale-increment model should the analysis move to the larger mixing basis of #125.

A covariance-aware generalized singular-value / likelihood-ratio formulation is preferable to an unweighted SVD.

## 7. Connection to the norm-10 Jordan-rank square

If a future norm-10 child is generated, apply the same construction to the entire function:

\[
H_{2,5}Z(u)
=Z_{10N}(u)-Z_{5N}(u)-Z_{2N}(u)+Z_N(u).
\]

Then

- rank-2 `f+log(N)g` gives zero for every u;
- a quadratic-log rank-3 term gives a nonzero function proportional to `2 log2 log5`;
- an ordinary power gives its own nonzero multiplicative curvature.

So the commuting square can eventually measure Jordan rank as a **functional** property rather than from one center coefficient.

## 8. Recommended order

Before new theory fields or new compute:

1. add `Mbar'(p0)` to the P48 derived output and score `J_N`, `Xi_N` on all existing full curves;
2. reconstruct `T_N(u)=D_e+S_o` on the already frozen u grid for P49 and N185/N265 as retrospective/development evidence;
3. freeze only the implementation/sign/channel semantics before #57 target scoring;
4. after N325/N425 reveal, score the q=2 and Jordan functional residual vectors without fitting `f(u)` or `g(u)`;
5. if both fail, quantify the minimum correction-function rank before adding named operator candidates.

## 9. Why this is a stronger mechanism test

The present q=2 and Jordan fits can both mimic one scalar sequence over a short size range. The full-curve test asks them to explain the same correction at several thermal coordinates with one fixed multiplicative law.

A Jordan interpretation becomes persuasive if a whole nontrivial function of `u` transforms with the logarithmic cocycle `log Q`; an analytic correction becomes persuasive if the same function transforms with `Q^-1-1`. Either outcome is much harder to obtain accidentally than a good two-parameter center fit.
