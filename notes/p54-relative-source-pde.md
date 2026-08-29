# Exact relative-Q source closure and topology inversion

Issue #114 gives a finite relative-fugacity partition function

\[
G(s,u)=\mathbb E_u[e^{s q}],\qquad q\in\{-1,0,1\},\qquad s=\log Q,
\]

whose first `s` derivative at zero is the matching observable.  The three-point
support is stronger than a moment identity: configurationwise, `q^3=q`.
Therefore

\[
\partial_s^3G=\partial_sG,
\qquad
G=P_0(u)+P_+(u)e^s+P_-(u)e^{-s}.
\]

Any thermal or score source `u` is a spectator in this algebra.  Every commuting
mixed derivative inherits the same linear PDE,

\[
\partial_u^\alpha\partial_s^3G
=\partial_u^\alpha\partial_sG.
\]

Equivalently, one may first use the unnormalized numerator
`sum w_0 exp(u t+s q)` and then divide by its `s`-independent partition
function.  That normalization cannot alter the `s`-PDE.

For `F=log G`, division by `G` gives the connected closure

\[
F_{sss}=F_s-3F_sF_{ss}-F_s^3.
\]

One further `u` derivative gives

\[
F_{sssu}=(1-3F_{ss}-3F_s^2)F_{su}-3F_sF_{ssu}.
\]

Thus all mixed-ledger entries with at least three **Q-source** insertions are
derived rows; only the columns with zero, one, or two Q insertions are independent.

## Exact inversion

At `s=0`, write `mu=F_s=E[q]` and `v=F_ss=Var(q)`.  Since
`E[q^2]=v+mu^2`, the complete topology distribution is

\[
P_+=\frac{v+\mu^2+\mu}{2},\qquad
P_-=\frac{v+\mu^2-\mu}{2},\qquad
P_0=1-v-\mu^2.
\]

The first two connected Q-source derivatives therefore contain the entire
finite three-sector distribution, not merely its first two summaries.  Higher
Q derivatives add no independent sector information.

## N=5 exact oracle

For the Gaussian `(2,1)` torus at `p=1/2`, exhaustive enumeration gives

- `(P_-,P_0,P_+)=(1/2,5/16,3/16)`;
- `mu=-5/16`, `v=151/256`;
- `F_sss=555/2048`, both directly and from the closure;
- `partial_s^3 partial_u^b G = partial_s partial_u^b G` for `b=0,...,4`,
  checked first on the unnormalized numerator with `u` coupled to the exact
  Bernoulli score multiple `2 occupied-N` (the normalized identity follows by
  the `s`-independent division above);
- the once-thermal-differentiated logarithmic closure holds exactly.

The machine-readable joint distribution and checks are in
`results/exact-relative-source-pde/latest.json`.

## Boundary: this does not solve thermal kappa3

The closed index is the number of **relative-Q** insertions.  Issue #54 asks for
the third derivative in the thermal direction, schematically `F_suuu`, whereas
the identity closes `F_sss`, `F_sssu`, and higher rows with three Q indices.
Nothing here reduces `F_suuu` to `F_su`, fixes its continuum value, or proves
`kappa3=-5/3`.  The gain is an exact compression of the mixed Q/thermal ledger:
calculate the zero-, one-, and two-Q columns only, then generate every higher-Q
column algebraically.
