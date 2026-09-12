# Completed analysis for existing #613/#276/#582 and #337/#636

2026-09-12. This is a handoff for existing issues, not a request for new production.
No remote issue lifecycle or PR state was changed by this local delivery.

## For #613 / #276 / #582: completed thin-geometry boundary

The exact finite row inequalities P0 <= (1-p^w)^m and
P2 <= [1-(1-p)^w]^m give rank-one concentration on fixed-width long tori and
on w=o(log m). The threshold CDF F=(1+M)/2 tends to 1/2 throughout (0,1),
while its probability law tends to (delta_0+delta_1)/2. In particular,
'w and m both grow' is not a sufficient replacement for #613's geometric
hypothesis: w_j=j, m_j=ceil(exp(j^2)) is a counterexample. The stated
ell/log N -> infinity theorem is NOT contradicted.

For every fixed w>=2 the two births have an independent limiting endpoint
law, with survival exp(-x^w) and exp(-c_w y^w), where
c_w=[z^0](z^-1+1+z)^w. The proof classifies minimal motifs and controls
nonminimal circuits and dependent rare-pattern counts. Joint moments and
permutation birth-rank clocks are covered. No growing-width Weibull uniformity
is asserted. Do not commission another small census or a thin-torus production
merely to rediscover this limit.

For widths 2/3/4 the spectral crossing is an internal balance of two rare
sectors, not concentration of the entire threshold law. The conditional
rank-2 probability is logistic on p=q+x/m, while F remains flat. The probability
window around the median is exponentially small. Its inverse sensitivity is
not improved by the exponentially small deterministic root bias.

Files: notes/thin-torus-two-birth-limits-20260912.md;
notes/thin-torus-central-layer-source-20260912.md; the two thin_torus scripts
and their generated results.

## For #337 / #636: completed source-visible spectrum

Using the preceding all-parameter P0/P2 trace identities and #708's original
509-state certificate, the generic scalar orders are M:16, E=P0+P2:28,
twice the tilted-mean numerator at exp(s)=2:28, twice Z at exp(s)=2:29.
Full integer Hankel minors and independent modular residues prove the lower
bounds; existing exact trace identities prove the upper bounds for all lengths.
The normalized tilted mean is a ratio, NOT claimed to satisfy an order-28
constant-coefficient linear recurrence.

A bounded intrinsic source restores common modes canceled at s=0. Its balance
root has exact derivative -E/M', and, under the proved local spectral contract,
shift -2s/(m*h')+O(m^-2+eta^m/m). For the thin sequences, Z(p,s)->1 everywhere
inside (0,1) at bounded s, so restored source symmetry alone is not a critical
point identifier. The source algebra remains correct; this is a scope correction
for interpreting it.

Files: notes/width4-topological-source-visible-modes-20260912.md;
scripts/width4_topological_source_spectrum.py; full source-spectrum JSON.

## Allocation and limits

No new task is needed for these completed analyses. No GPU, wider state table,
new Monte Carlo, or large retrieval is justified by this delivery alone.

The remaining cross-width question needs estimates uniform in width on the
VISIBLE mode weights and gaps, plus geometry and readout conventions. Neither
a fixed-width theorem nor a generic source-order table provides that uniformity.
A future computation should target one stated separating inequality, not a
width ladder for its own sake. #275's same-source original-U prediction columns
remain a distinct missing physical input; no candidates are promoted here.

The snapshot one-nonzero-event cost formula is restricted to iid rank
observations at one p. It is not a budget for Newman-Ziff histogram reconstruction,
Rao-Blackwellization, importance sampling, or exact transfer algebra.

The preceding parametric definition is an explicit dependency and is not
asserted merged. Existing source and production artifacts are untouched.
