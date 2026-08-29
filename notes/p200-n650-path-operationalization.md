# N650 path memory: exact operationalization and obstruction

Issue #200 proposed one N650 endpoint with two retained intermediate path flags:

```text
65 --(2-i)--> 325 --(1+i)--> 650
65 --(1+i)--> 130 --(2-i)--> 650.
```

The endpoint geometry is exact.  With the repository period convention

\[
P(a+ib)=\begin{pmatrix}a&-b\\b&a\end{pmatrix},
\]

both orders multiply by `3+i`.  They produce exactly the same integer period
matrices, not merely isomorphic graphs:

```text
(8+i)(3+i) = 23+11i,
(7+4i)(3+i) = 17+19i.
```

There is therefore one unmarked N650 configuration stream and one endpoint
histogram per orientation.

## What the intermediate character flag can mean

The final-to-parent deck group is

\[
K_{3+i}=\mathbb Z[i]/(3+i)\cong\mathbb Z/10.
\]

Reduction modulo `2-i` and `1+i` gives an exact CRT coordinate

\[
K_{3+i}\longrightarrow K_{2-i}\times K_{1+i}
\cong \mathbb Z/5\times\mathbb Z/2.
\]

The machine oracle enumerates all ten fiber sites and all `2^10=1024`
binary fiber configurations.  It performs the exact character transform in
both orders, representing every coefficient as an integer combination of
tenth roots of unity.  Every configuration gives

\[
\Pi_5\Pi_2=\Pi_2\Pi_5.
\]

Thus a path difference made only from transported deck-character projectors
is an exact zero control.  It cannot be a random rank-one memory signal.

## Why an intermediate homology flag is not automatically defined

A configuration on the final cover descends to an intermediate quotient only
when it is constant on every kernel fiber.  In the exhaustive ten-site oracle:

```text
descends to the degree-5 intermediate: 32 / 1024
descends to the degree-2 intermediate:  4 / 1024
descends to both:                        2 / 1024.
```

For a generic Bernoulli configuration there is no binary configuration that
*descends* to the intermediate torus.  A pushdown convention is needed before
speaking of an intermediate Bernoulli mask.  The set-theoretic direct image is
fiber OR (the universal image is AND); a one-occupied-site example makes those
two conventions disagree.

This ambiguity does not rescue an order mark.  Honest functorial images still
compose: nested OR and nested AND agree in both orders on all 1024 masks.
Likewise, if `Pi` is a connectivity partition and `R2,R5` are the two fiber
equivalence relations, then

\[
(\Pi\vee R_2)\vee R_5=(\Pi\vee R_5)\vee R_2
\]

by associativity and commutativity of equivalence-relation join.  Hence the
induced wrapping/homology direct image also has no antisymmetric path memory.

Therefore the old notation

```text
C_mark = H_[5 then 2] - H_[2 then 5]
```

is identically zero for linear character transport and for honest functorial
connectivity images.  A different value requires an explicitly nonfunctorial
truncation, representative, order-statistic, or orientation rule.  Treating
the two names as independent samples would be pseudoreplication.

## Revised acquisition semantics

Do not add the old `C_mark` field to the N650 production runner.  If N650 is
run for its unmarked radial endpoint, retain one histogram and the usual full
batch covariance.

If a charged cover experiment is desired, freeze it separately:

1. label every final site by its full `Z/10` deck coordinate, equivalently the
   exact `(Z/5,Z/2)` CRT pair;
2. keep `Pi5 Pi2 - Pi2 Pi5 = 0` as a configurationwise implementation gate;
3. for a nonzero linear response, pair `S_chi` with a declared
   opposite-character marked/covariant row `O_chibar`;
4. for an invariant susceptibility, use the full Bernoulli Hessian in
   `chi tensor chibar`, including its diagonal likelihood correction.

The smallest natural nonlinear alternative is symmetric rather than
antisymmetric.  For a topology functional `h`, acquire

\[
\Delta_{25}h=h(\Pi\vee R_2\vee R_5)-h(\Pi\vee R_2)
              -h(\Pi\vee R_5)+h(\Pi).
\]

This mixed inclusion-exclusion defect can be nonzero without inventing path
order, and is the appropriate one-shot target if #200 wants an interaction
between the degree-2 and degree-5 coarse relations.

There is an exact ten-point witness.  Start from the discrete partition on
`C2 x C5` and set `h(Pi)=10-#blocks(Pi)`.  The four ranks are `0,5,8,9`, so
the two joins commute while `Delta25 h=9-5-8+0=-4`.  In production the same
Möbius four-term construction should use transported typed ambient-H1 rank or
matching charge `q`, not this toy partition rank, and all four correlated
rows must be retained before forming the difference.

This is the same finite-volume selection rule already proved in #244.  It
does not rule out nonlinear coarse-graining memory; it shows that such a
claim needs an explicit nonlinear map and cannot be obtained by renaming the
two commuting Gaussian factors.

Reproduce the exact gate with:

```bash
python3 scripts/p200_n650_path_oracle.py \
  --output results/exact-p200-n650-path/latest.json
python3 -m unittest tests.test_p200_n650_path_oracle
```
