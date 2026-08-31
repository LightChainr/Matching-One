# Fixed one-hole decision at the closed-source endpoint

This decision is specified before reading the new defect counts. It follows
the closed-source result0d19179f and the overview's one-hole question in
b8d043fc, without adding a source or fitting an interpolation curve.

Keep the ordinary square graph and bulk source S*=C+F+Bvac fixed. The
Bernoulli chart is pA=s+(1-s)p, pB=p, with thermal p differentiation holding
s,t fixed. Every U(s,t) is evaluated at its own pooled q=0 root and retains
the parent area factor N^(13/8)/2 and normalized cos4 projector.

The fixed parent pair is N50 `(5,5),(1,7)`, delta_cos4=-1152/625. At s=1
it maps to the already completely enumerated N25 `(5,0),(4,3)` pair with
complemented occupations. The endpoint root and baseline coefficients
are taken from that exact calculation, not re-estimated from a new sample.

## One new finite population

At the origin A site fix a vacancy, keep the other24 A sites occupied,
and enumerate all2^25 occupations of B separately in each parent geometry.
Save exact per-K_B sums for1,q,E,S*,qS*,ES*. Translation relates all25
possible A defects, so one defect per geometry supplies the full endpoint
derivative. The free Bernstein degree is25, not50. No baseline replay,
Monte Carlo, additional defect class or parameter sweep is requested.

For each observable O=1,q,E, let H0_O and Hd_O be unnormalized Bernoulli
moments with source exp(tS*) in the saturated and one-defect ensembles.
For x=s-1 the full numerator is

`H_O=H0_O - x*25*(1-p)*(Hd_O-H0_O) + O(x^2)`.

Normalize separately in each geometry by H_1 before pooling/projecting.
The source derivative of that denominator, including its mixed x,t term,
is part of the observable. Differentiate the factor1-p as well.

## Predeclared model comparison

The exact saturated identity alone does not predict the interior.
A specific stronger extension is **source-independent geometric gain**:

`U50(s,t)=g(s)*U25(t)+O((s-1)^2)` near t=0,

with g(1)=2^(13/8) and unknown g'(1), but no t-dependent fitted gain.
Its necessary endpoint condition is

`R = U*U_st - U_s*U_t = 0`.

This is the primary decision. A rational enclosure of R/A50^2 excluding
zero rejects this local gain extension; an enclosure containing zero
leaves it unresolved. The unknown g'(1) is eliminated algebraically, not
estimated to rescue the result. This is a finite equality question, not
a statistical significance threshold or continuum model identification.

Also report the single mixed coefficient Xi=U_t,epsilon=-U_st requested
by b8d043fc, where epsilon=1-s. Nonzero Xi excludes the stated interior
thermal-only extension; zero does not establish full profile closure.
Report U_s and U_t as ingredients, not extra independently selected tests.
All derivatives include the pooled-root displacement and slope denominator.

One exact score is performed after code and count provenance are committed.
No result changes the previous F4 random-block unresolved stop, the P154
lag1 decision, or the P334 prospective decisions. No new source/size/top-up
is specified as a rescue. The geometry and its different Smith classes
remain part of this finite result's scope.
