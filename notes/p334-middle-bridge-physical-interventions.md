# P334: marked middle bridges cut every higher-order birth mechanism in two real prefixes

## Main result

Permanently blocking A's single marked middle site198 and B's six marked
middle sites makes their **entire physical birth clocks identical**, at every
one of the174 subset cardinalities. More strongly, the common event reduces
exactly to the original15 direct gates and the same six-edge pair-trigger
forest. Thus the sites initially located from third-order triggering form an
explicit hitting set for **every minimal trigger of order at least three** in
these two real configurations, not merely their first third-order difference.

The intervention is prescribed on counters43042508631(A) and43042514803(B),
both N425/second/k0=252/age10/ell(12,-19). Their original clocks and all
unaffected factors are reused from87b6ca5b; the site sets are taken unchanged
from the marked-middle decomposition fd96cd95. No new sample is generated.

## Preserve the original insertion clock

There remain d=173 insertion labels in every comparison. A blocked site stays
permanently vacant but still occupies its position in the uniform random
insertion order as an inert dummy. If a factor contains q blocked sites, remove
those nodes for its active-site reliability calculation, then multiply its
safe polynomial by `(1+z)^q`. Reuse all other factor polynomials and off-core
dummy factors. The survival denominator remains `C(173,k)`, never `C(173-q,k)`.

The seven single-site blockades and B's joint blockade require only nine
affected-factor DP evaluations. A's singleton is already its joint blockade
and is reused. Neither original baseline was re-solved. The initial factor
evaluation/readout pass took0.030s; later pivotal and interaction readouts use
the saved factor polynomials without another DP.

## Physical clock response

`S40` means P(T>40). All values below are exact rational outputs displayed as
decimals. The single/pair state remains unchanged because these middle sites
are not direct gates or endpoints of any minimal pair trigger.

| Prefix / blocked sites | Remaining g3 | E[T] | Mean increase | S40 | S40 / baseline |
| --- | ---: | ---: | ---: | ---: | ---: |
| A baseline | 5 | 10.48090361 | — | .01078726 | 1 |
| A:198 (also all-middle) | 0 | **10.51477310** | .03386949 | **.01140416** | 1.057188 |
| B baseline | 19 | 10.39836128 | — | .00935215 | 1 |
| B:24 | 18 | 10.40673211 | .00837082 | .00950763 | 1.016625 |
| B:25 | 17 | 10.41127912 | .01291784 | .00957025 | 1.023321 |
| B:184 | 13 | 10.43198834 | .03362705 | .00991585 | 1.060275 |
| B:340 | 14 | 10.42657940 | .02821811 | .00982463 | 1.050521 |
| B:94 | 17 | 10.41173185 | .01337057 | .00959121 | 1.025563 |
| B:361 | 16 | 10.41588858 | .01752730 | .00964958 | 1.031804 |
| B:all six | 0 | **10.51477310** | .11641182 | **.01140416** | 1.219416 |

The original A−B mean gap is0.08254232606. The joint-blocked mean gap and
tail gap are exactly zero, and all174 full safe coefficients agree. The
original B survival lies strictly below A for k=3..154; this entire ordering
disappears after the specified blockades. B's tail40 increases21.94%, so the
middle profile matters beyond the tiny third-step probability discrepancy.

This is an exact comparison under two defined configuration-specific
counterfactual interventions. It is not a unique causal attribution of the
original mean gap to a transferable scalar g3, nor a claim that deleting any
equal-sized set would give the same effect.

## Why the two complete events collapse

The archived minimal-pair graphs are isomorphic, not merely equal in their
first moments. Each is a disjoint union of P3 and a five-vertex double-star
tree with degrees `(3,2,1,1,1)`:

```
A P3:      130--129--286
B P3:      202--204--360

A double star: 39--40--354--355, with additional leaf196 at40
B double star: 180--27--71--72, with additional leaf338 at27
```

These trees have independence polynomials `1+3z+z²` and
`1+5z+6z²+2z³`. There are158 non-direct sites, of which8 meet pair edges and
150 are irrelevant to the pair event. Its safe polynomial is

\[
F_{\rm pair}(z)=(1+z)^{150}(1+3z+z^2)(1+5z+6z^2+2z^3)
=(1+z)^{151}(1+7z+15z^2+10z^3+2z^4).
\]

The additional `(1+z)` in the second form is an algebraic factor, not an
additional identified Boolean-irrelevant site. Both joint-blocked physical
polynomials equal this F_pair coefficient by coefficient. They start

```
1, 158, 12397, 644025, 24920287, 766092007, ...
```

This proves equality of the actual finite events, not just equality in law:
pair triggering is a subset of physical triggering on the same labeled
variables. At each cardinality the pair-safe sets contain the physical-safe
sets. Equal counts therefore force equality of those sets. It follows that
every original inclusion-minimal trigger of order≥3 intersects the prescribed
middle set; otherwise it would survive the blockade as a non-pair trigger.
This hitting-set statement is specific to these two prefixes and does not
assert that third-order middles always hit every higher-order trigger.

## Exact link to the final insertion mark

For one site v, decompose the original safe polynomial as
`F=F0+z F1`, where F0 counts safe subsets with v absent and F1 with v present.
Its inert-site version is `(1+z)F0`, hence

\[
\Delta F=z(F_0-F_1)=zD_v.
\]

The coefficient of z^(k−1) in D_v counts safe preceding subsets on which v
would trigger birth. Uniform insertion order then gives

\[
P(T=k,V_{\rm final}=v)=\frac{[z^{k-1}]D_v}{k\binom{173}{k}},
\quad \Delta S(k)=kP(T=k,V_{\rm final}=v),
\]

and, since the15 untouched direct gates still guarantee eventual birth,

\[
\Delta E[T]=\pi_v E[T\mid V_{\rm final}=v].
\]

Thus the single-site blockade polynomials already supply the baseline final
insertion law without force-on DP or a new random continuation:

| Site | Baseline probability it is the final birth site | E[T given it is final] |
| --- | ---: | ---: |
| A:198 | .001417288 | 23.89739 |
| B:24 | .000347484 | 24.08983 |
| B:25 | .000560160 | 23.06099 |
| B:184 | .001468793 | 22.89434 |
| B:340 | .001232865 | 22.88825 |
| B:94 | .000566825 | 23.58851 |
| B:361 | .000761290 | 23.02314 |

These final-site events are comparatively late (conditional means≈23 rather
than≈10.4 overall). The six B sites together are final in0.00493742 of baseline
permutations. A structural middle is not automatically the last site, but its
complete knockout response now identifies that final-site law exactly.

## Joint interference starts at fourth order and changes sign

Let R(z) be the all-middle safe-polynomial increase minus the sum of the
single-middle increases, each measured against the same baseline. In B,

```
[z^0..z^5] R = 0, 0, 0, 0, -7, -1021.
```

All triples have exactly one marked middle, so third-order blockade effects
are additive. At higher orders they are not: R has negative coefficients for
k=4..27 and positive coefficients for k=28..155. The integrated joint effect
exceeds the sum of singles by0.00238012293 mean steps and0.000105745529 in
tail40. Early overlap and late complementary pathway effects therefore cannot
be represented by one time-independent additive site weight. The observation
is the coefficient sign change; no particular single diagram is inferred from
its sign alone.

## Saved output and scope

`results/p334-middle-bridge-interventions/physical_clock_blockades.json` saves
all eight intervention clocks, nine affected factor polynomials, final-site
joint laws, exact interaction coefficients and the common pair-core identity.
`scripts/p334_middle_bridge_interventions.py` reproduces only these prescribed
cases and reuses saved affected factors on subsequent readout passes.

Same old production source/dependency group; no newMC, server connection,
full-site sweep, baseline re-solve, or repository-wide test campaign. This
advances the interpretation from a missing third-order statistic to an exact
finite intervention on the complete physical birth mechanism.
