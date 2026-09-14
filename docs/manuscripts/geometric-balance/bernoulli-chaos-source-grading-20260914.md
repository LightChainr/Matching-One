# Bernoulli-chaos source grading: an exact microscopic pair-exchange tangent basis

Date: 2026-09-14

Status: exact finite product-measure/source algebra.  It is proposed as a microscopic basis for #61/#802; it does **not** by itself establish RG/OPE parity of continuum fields.

## 1. Exact paired standardized variables

For the primal Bernoulli site law with parameter `p`, define

```text
xi_v = (n_v-p)/sqrt(pq),  q=1-p.
```

Under the matching/complement pairing

```text
n_hat_v=1-n_v,
q=1-p,
```

the standardized variable of the paired model is

```text
xi_hat_v
 = (n_hat_v-q)/sqrt(qp)
 = -xi_v.
```

Therefore for every finite site set `A`,

```text
Psi_A = product_(v in A) xi_v
```

obeys the configurationwise exact relation

```text
boxed: Psi_A -> (-1)^|A| Psi_A.
```

For a finite Bernoulli product space the `Psi_A` form the usual Hoeffding/Walsh orthonormal basis.  Thus any square-integrable microscopic source has a unique decomposition by chaos degree.

This is a grading of **source functions under the exact complement pairing**.  It is not a theorem that the RG scaling fields reached by these sources carry the same scalar parity.

## 2. Degree 0 and degree 1 are nuisance directions for the present question

- degree 0: constant normalization;
- uniform degree 1: thermal/logit score;
- nonuniform degree 1: spatial one-site probability field, whose first-order response may vanish by translation/symmetry but is still a one-site measure deformation.

For a source-identification problem where the thermal coordinate is profiled out, the relevant source space should first be quotiented by

```text
span{1, K}
```

(or the exactly corresponding normalized thermal score).  This is the same quotient already present in #773 source calculus.

## 3. #802 black/white pair-source redundancy is a degree-2 identity

For one adjacent pair `(i,j)`, write `x_i=n_i-p=sqrt(pq) xi_i`.  Then

```text
n_i n_j
 = p^2
 + p sqrt(pq)(xi_i+xi_j)
 + pq xi_i xi_j,
```

while

```text
(1-n_i)(1-n_j)
 = q^2
 - q sqrt(pq)(xi_i+xi_j)
 + pq xi_i xi_j.
```

Their genuine degree-2 component is **identical**.  Summing over a regular row/cycle gives exactly the observed relation

```text
H_W = const - 2K + H_B.
```

Hence, after quotienting constant + thermal directions, black-pair and white-pair sources are the same vector.  Their apparent first-order distinguishability before normalization was necessarily spurious.

This turns the #802 erratum into a general source-design rule rather than a one-off bug.

## 4. The first exact exchange-odd nonthermal source occurs at odd degree >=3

After removing degree-1 thermal motion, an odd chaos source of degree three is the smallest local source carrying an exact minus sign under complement without being merely thermal.

A generic translationally symmetrized example is

```text
H_3 = sum_x sum_a c_a
      xi_x xi_(x+r_a) xi_(x+s_a),
```

with coefficients chosen to form a declared spatial/D4 representation.

Under the exact primal/matching complement pairing:

```text
H_3 -> -H_3.
```

This gives a genuine lattice-side `D` source direction.  It should be compared with an even degree-2 control and the degree-1 thermal control.

Again, the conclusion is only about microscopic source exchange.  An RG-compatible continuum tangent map must still be demonstrated before assigning a local scaling field a parity label.

## 5. A concrete finite-width tangent experiment

Choose a small, preregistered source basis, for example:

```text
H1 : uniform degree-1 thermal score;
H2 : one D4-symmetrized degree-2 pair/motif source;
H3 : one D4-symmetrized degree-3 motif source.
```

For each safe phase and width compute physical normalized Perron score responses

```text
mu_G,a = partial_(g_a) log lambda_G^physical,
```

or equivalently the normalization-safe free-energy derivatives.

Then record

```text
D_a(w) = mu_4,a - paired(mu_8,a),
S_a(w) = mu_4,a + paired(mu_8,a),
```

with the exact microscopic exchange grade carried separately.

The useful question is not “does H3 prove an odd CFT field?” but

> does the large-width response matrix become block-compatible with the exact microscopic complement grading, and which continuum scaling powers live in each block?

This is an operational version of the RG-tangent programme in #61.

## 6. Why degree alone will not identify spin four

Square-lattice spin zero and spin four both lie in the trivial `C4` representation.  A D4-invariant degree-3 source therefore does not by itself separate a scalar continuum field from spin four.

Angular geometry / same-circle projectors remain the correct analyzer of the H0/H4/... decomposition.  The source-chaos grading solves a different ambiguity: normalized pair-exchange direction versus thermal/gauge redundancy.

The two analyses should be crossed only after each is separately typed:

```text
angular irrep    : H0/H4/H8/...
source exchange  : microscopic chaos grade / S-D response
radial exponent  : measured only after the first two are separated.
```

## 7. Relation to original-U

Any candidate physical source for #275 should first be decomposed in the same way.  If two proposed source implementations differ only by degree 0/thermal directions, no amount of extra sampling can distinguish them after nuisance profiling.

Conversely, a candidate forward map that predicts different higher-chaos components supplies a real new observable direction.

This is a finite probability statement and does not replace the original-U continuum normalization contract.

## 8. Claim boundary

Exact:

```text
xi_hat=-xi,
Psi_A -> (-1)^|A| Psi_A,
black/white pair sources share the same degree-2 component,
source quotient modulo normalization/thermal directions.
```

Programme/hypothesis:

- large-width RG tangent map preserves or asymptotically respects the microscopic grading;
- particular continuum fields can be assigned to the observed response blocks;
- a chosen degree-3 source has useful overlap with the sector correction of interest.

The main correction is conceptual: construct parity first on the **source space where it is exact**, then test whether RG transports it.  Do not infer the reverse direction from a finite matching identity.