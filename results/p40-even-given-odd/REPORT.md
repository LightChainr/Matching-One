# The cluster source moves the middle homology sector at fixed matching

## The finite even response is resolved; directional attribution is still open

On the existing million-configuration P40 blocks at N65 and N85, the
absolute black-plus-white cluster source changes the rank-1 sector even
after its first-order change in the matching mean is cancelled. All four
geometry-wise raw responses exceed110 standard errors. Thus this physical
source's **first-order topological marginal is not a pure relative-q
fugacity tangent** at the measured point.

Its H4 direction difference is much smaller. The common raw source's
two-size zero comparison gives nominal p=.21156. An auxiliary,
geometry-adapted compensated source gives p=.04911, driven mostly by N65;
this is a hint in correlated auxiliary views, not a confirmed directional
mechanism. Neither result identifies a continuum energy operator.

This completes the previously missing E_top/source products and the
even-given-odd readout. It reobserves the original configurations once;
there are no new Monte Carlo samples, scientific test suites, server jobs
or GPU operations.

## The source and the fixed-matching experiment

At the archived Bernoulli probability `p=.59274605079`, define

```text
q = rank_black − 1,             E = E_top = q² = 1 − I(rank_black=1),
S = (C_blackNN + C_whiteMatching)/N,
C = Cov(E,S) − Cov(E,q) Cov(q,S)/Var(q).
```

S is the same explicit raw source on every geometry. White connectivity
uses the matching graph. Occupied NN edge count is a control, not E_top.
The source normalization is per site: responses to the extensive cluster
count instead are exactly N times the reported responses, not a fitted
scaling correction.

In the positive measure proportional to
`Bernoulli(p) exp[lambda S + h(lambda)q]`, choose
`h'(0)=−Cov(q,S)/Var(q)`. This holds the matching mean fixed to first
order. Then `d<E>/d lambda=C` and **`dP(rank1)/d lambda=−C`**.
This compensating q source is explicit; it is not a change of Bernoulli p.

N65 uses orientations (8,1)/(7,4); N85 uses (9,2)/(7,6). Within each N,
directions share the original cyclic-site occupations. Each block has
one million old configurations in100 aligned batches. N enters the
engine's random key, giving separate N domains under the usual PRNG
independence assumption. Different analyses of a block remain dependent.

## The raw source has a positive even tangent in all four geometries

For N65, the first and second C responses are respectively
`.00141473398 ±.00001167597` and `.00144263065 ±.00001089886`.
For N85 they are `.00110134755 ±.00000795582` and
`.00109728689 ±.00000996137`. Their z-scores range from110.15 to138.43.
Increasing absolute-cluster fugacity while holding matching fixed therefore
decreases the probability of the intermediate homology sector locally
in lambda, on all four measured geometries.

The three-sector conditional means make the distinction precise. Write
`m_r=E[S|q=r]` and `c=(m_-1+m_+1)/2−m_0`. Since all three sector
probabilities are positive,

```text
C = c [Var(E) − Cov(E,q)²/Var(q)].
```

The raw c values are `.006045591 ±.000049231` and
`.006170530 ±.000047234` at N65, and `.004700097 ±.000033908`
and `.004687375 ±.000042608` at N85. The JSON retains the full
conditional means, probabilities and the algebraic identity residual.
These nonzero curvatures distinguish the first-order finite marginal
from `a+bq`; they do not distinguish every finite-lambda deformation,
high-order response or microscopic mechanism.

## Compensation preserves an even response, with a different source definition

The clock source removes an empirical linear projection on K and K(K−1).
The full source adds occupied NN edges and local Euler count. Coefficients
are fitted separately for each geometry and refitted in every delete-one.
No q or E is included in those source projections.

After full source compensation and then matching compensation, C is
`.00027271023 ±.00000864071` / `.00030158838 ±.00000870608` at N65,
and `.00022792857 ±.00000705241` / `.00023437760 ±.00000822189`
at N85. These auxiliary responses remain28.51–34.64 standard errors from
zero. They are geometry-adapted source diagnostics, not one universal
local field or a proof that the raw effect is independent of every clock
or Euler mechanism.

The operations are sequential. Although S_res has zero empirical
covariance with the original controls Z, subsequently adding hq generally
changes their means by `h'(0) Cov(Z,q)`. This experiment holds q fixed;
it does **not** hold q and every Z fixed simultaneously.

## Directional H4 is the remaining small response

Every direction contrast uses `(first−second)/Delta cos(4theta)` with
the exact rational geometry, and no fitted exponent.

For the primary raw C, the N65 contrast is
`−.00002046240 ±.00001177590` (z=−1.73765), while N85 gives
`+.00000254673 ±.00000863206` (z=.29503).
The nominal joint zero statistic is `3.10647/2`, p=.21156.

For clock-compensated C, the contrasts have z=−1.96226/.05530,
with nominal joint p=.14562. For full-compensated C they are
`−.00002118232 ±.00000885603` and `−.00000404464 ±.00000730597`,
z=−2.39185/−.55361, with nominal joint `6.02744/2`, p=.04911.
The uncompensated JE comparisons have p=.19803/.14377/.05257 for
raw/clock/full respectively. These six views were declared together and
are correlated; the .04911 auxiliary value is not an independent
confirmation or a multiple-comparison-adjusted discovery.

The global even tangent is clear, while its direction difference is not
securely attributed. A shared-stream cancellation allocation is not needed
for this statement: C and its direction difference are marginal quantities,
although the chosen direction pairing determines their sampling covariance.

## One reobservation, paired uncertainty, immutable inputs

The [declaration](../../analysis/p40_even_given_odd_replay.json) was committed
before these missing products were observed. The original production was
already public, so this is a retrospective analysis, not a prereveal test.
The pinned backend blob `b3d2047e35a6840cc7236e4a6088a291c734313b`
is identical at the archived source commit and its recorded original
engine commit. It supplies the exact original geometry, classifier and
counter mapping; no weaker backend is substituted.

The replay adds E*S and four E*control sums. Existing batch first moments
agree to approximately1.3e−12 after normalization, and E counts agree
exactly with the old q² Gram entries. This comparison occurs inside the
single analysis; it is not a separate repeated test campaign.

The scorer uses sample covariances with denominator n−1. Each delete-one
removes the same batch from both directions and refits all source
coefficients. The [JSON](latest.json) retains all100 vectors and the full
49×49 output covariance for each N. Missing cross-direction E Gram cells
are not filled with zeros: the marginal estimands use aligned batch
jackknife covariance directly. Nominal Gaussian statistics do not become
exact finite-sample certificates.

The [run receipt](run.json) records the two deterministic reobservations,
18.29 and23.44 CPU seconds, plus compiler, backend and output hashes.
The saved-moment calculation took about.14 seconds on local ARM64 Python.
Two parent sizes are given as numerical responses rather than an unsupported
scaling plot.

```bash
git fetch origin analysis/p40-production-motif-projection-20260830
python3 scripts/replay_p40_even_given_odd.py --output-dir results/p40-even-given-odd-reproduction
python3 scripts/analyze_p40_even_given_odd.py --output-dir results/p40-even-given-odd-reproduction
```

Both commands refuse to overwrite their existing outputs.

## Next: transport this source to the actual norm-4 question

The pure-q prediction is now explicit, geometry by geometry:
`JE^(q)=Cov(E,q) Jq/Var(q)`. Its difference from the observed response
is C. The next mechanism question is whether the **original norm-4
projection of this difference** is nonzero, using the source definition,
root/thermal coordinate and slope normalization of that comparison.
Form the local ratio before applying the cross-size projection.

The present N65/N85 are two parents, not a dyadic chain. The original
three-size chains are65→130→260 and85→170→340; the existing N65/N13020k
source matrices can supply a lower-statistics fixed-p follow-up, but do
not alone implement that original root-normalized test. Thus the next
missing output is a specific physical source-transport residual, not
another first E_top source test, another generic certificate, or a
new exponent fitted to these two parent sizes.
