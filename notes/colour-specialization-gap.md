# A finite colour-sector label does not fix its generic-Q continuation

## The completed Q4 decision, and the next actual distinction

The [execution result](https://github.com/LightChainr/Matching-One/blob/54352b2eefa651ca482ca84837053c792e82c71e/results/p337-s4-trace-transmission/score/score.json)
has completed the [specified J22 calculation](four-leg-trace-denominator-interface.md):
`J22=+5.440121494634842e-6`. It rejects complete normalizer-neutrality
of that finite N25/m2 component. Its raw q/E numerator is exactly zero;
the geometric difference of the fractional trace's thermal derivative
supplies the dominant positive term. No repeat of that score is needed.

The next distinction is precise: **taking a full finite-colour central
projection and specializing a stable generic-Q sector need not commute**.
Below are two actual torus connectivity counterexamples to a fixed
dimension-only continuation. They do not dispute the Q4 result, or rule
out a specified continuation with its full multiplicity structure.

## One class sum on actual homology-decorated closures

Use the [explicit physical closure](https://github.com/LightChainr/Matching-One/blob/4a4390f2aeff3e79855fb11ef1185ba52c80a43a/notes/closed-source-finite-torus-pair-closure.md).
For c essential rank1 components with common first deck winding u,
factor out activity, contractible-component colour factors and the
common `Q^(-1/2)` rank projection. The remaining one-seam class function is

```text
F_(u,c)(pi)=Fix(pi^u)^c.
```

For integer Q≥4, the irreducible character and dimension are

```text
chi_[Q-2,2] = binom(X1,2)+X2-X1,
d(Q) = Q(Q-3)/2,
m_Q(u,c) = (1/Q!) sum_pi chi_[Q-2,2](pi) F_(u,c)(pi).
```

Here Xj counts j-cycles of pi. The full central trace coefficient is
`d(Q)m_Q`, not just d(Q). The sum is exact over conjugacy classes, whose
probabilities are `1/prod_j(j^Xj Xj!)`. A winding-power class function
can be virtual: a negative coefficient is not a negative probability.

The candidate under comparison is now fixed:

> After factoring the common physical weights, the coefficient m_Q of
> every fixed closure is already constant for all integer Q≥4; hence its
> Q4 central trace can be continued merely by multiplying by d(Q)/2.

This is a definite assertion about the physical closure family, not a
claim that every researcher has adopted it. The following counterexamples
exclude it exactly.

## Three fixed geometries and their exact coefficients

| Actual occupied pattern, seam along first period | (u,c) | m4 | m5 | stable m | trace at Q4 | trace at Q5 | dimension-only Q5 forecast |
|---|---|---:|---:|---:|---:|---:|---:|
| Axis5×5, only rows0,2 occupied | (1,2) | 1 | 1 | 1 | 2 | 5 | 5 |
| Axis6×6, only rows0,2,4 occupied | (1,3) | 5 | 6 | 6 | 10 | 30 | 25 |
| Axis7×7, occupied x=3y+j mod7, j=0,1,2,3 | (3,1) | −1 | 0 | 0 | −2 | 0 | −5 |

The first two patterns are separated straight essential cycles. In the
third, each row is a four-site path, and adjacent rows meet only at the
successive endpoints. The induced NN graph is one cycle: traversing it
changes the lifted position by `(3L,L)`, hence primitive winding `(3,1)`.
The same construction works for L≥7. These are individual physical
occupation configurations, not arbitrary representation matrices.

For the first row stability begins at Q=4. For the other two it begins
at Q=5. The counterexamples concern the full family, **not a claim that
these topologies occur inside the scored N25 packet**.

### Derivation without a lattice enumeration

For a uniform permutation of Q colours, factorial moments obey

```text
E prod_j (Xj)_(aj) = prod_j j^(-aj), if sum_j j*aj≤Q;
                    0, otherwise.
```

For u=1, expand `X1^c=sum_k S(c,k)(X1)_k`. Direct multiplication gives

```text
E[chi (X1)_k]
 = 1_(Q≥k+2)+(k-1)1_(Q≥k+1)+k(k-3)/2 * 1_(Q≥k).
```

In the stable range Q≥k+2 this is `binom(k,2)`. Thus c=2 gives m=1
already at Q4. For c=3 the k=3 term is2 at Q4 and3 at Q≥5; adding
`3 E[chi (X1)_2]` gives5 versus6. The omitted class moment at Q4 is
exactly one unit of multiplicity, not a change in d(Q).

For u=3,c=1, `Fix(pi³)=X1+3X3`. At Q≥4, `E[chi X1]=0`, while

```text
3 E[chi X3] = 1_(Q≥5)-1_(Q≥4).
```

It is−1 at Q4 and0 for every integer Q≥5. This proves the last row,
including the disappearance of a nonzero finite S4 component in the
stable sector. No extrapolation from noisy data is involved.

## What is and is not fixed by the generic component

For these diagrams the stable central coefficients are respectively
`d(Q)`, `6d(Q)` and0. Their rational/polynomial continuations are fixed
by those stable values, and at Q4 they give2,12,0. Projecting the actual
S4 closures instead gives2,10,−2. Thus the operations do not commute in
general. At Q1 the stable coefficients formally give−1,−6,0; these
numbers are not literal one-colour sector dimensions or probabilities.

The original **total occupation law** already has its declared real-Q
completion. This result concerns decomposition of that law into colour
components, not an ambiguity in the total measure. A particular local
pair insertion and a full isotypic seam trace can have different
multiplicity contractions. The correct next object is therefore an
explicit marked-port/intertwiner closure in the same physical family,
followed by its *complete* finite Q1 jet. One must not freeze its Q4
multiplicity and attach only a dimension derivative.

The [removable twist-jet formula](https://github.com/LightChainr/Matching-One/blob/f43b3674ce29e12629dd790bcbb7370abc5cefbc/notes/closed-source-removable-twist-jet-interface.md)
then supplies the already-derived U functional. If working with
`J=R/(sqrt(Q)-1)`, retain `J0=2R_Q|1` and
`J_logQ|1=R_QQ|1+R_Q|1/2`. A finite colour value alone supplies neither
jet. The regular endpoint's all-Q zero remains intact.

## Reproduction and scope

The [fixed calculation contract](../analysis/colour_specialization_gap_contract.json)
contains only the three physical closures above and Q4/Q5. Run

```bash
python scripts/analyze_colour_specialization_gap.py \
  --contract analysis/colour_specialization_gap_contract.json \
  --output-dir results/colour-specialization-gap
```

The script sums exact S4/S5 conjugacy classes once and stores all rational
terms. It is a reproducible algebraic consequence of the stated formulas,
not a blind numerical forecast or independent evidence. It does not
enumerate occupation configurations, rescore J22, search seams, add a
fourth fitted mechanism, or run a scientific test suite.
