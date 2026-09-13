# Sewing with memory: exact site weights, unbiased cuts, and a heat-kernel conjecture

2026-09-13. Continuation of #739, responding to the returns of #740 and #741.

**Status.** Sections 1–5 are finite identities, an elementary infinite-cylinder
limit, or algebraic reanalysis of supplied numbers. Section 6 is explicitly a
research conjecture and its conditional consequences. It is not a proof of the
fixed-subcritical site prefactor, parameter analyticity, or near-critical
universality. No new width-12 calculation was performed.

## 1. What the two returns do and do not establish

The #741 result at commit `3745b13b8e1126017a1e88567a6af44881b8a1ff` supplies
six component densities, at widths 4,8,12, for NN site p=1/4 and matching site
p=1/8. The requested equal-spacing contrast is retained:

| graph and p | beta_eff(4,8,12), recomputed from the supplied logs |
|---|---:|
| NN, 1/4 | 0.792584457189 |
| matching NN+NNN, 1/8 | 0.508452576643 |

These are effective finite-width contrasts, not determinations of an
asymptotic exponent. The width-12 residual certificates are reported by the
external computation; this delivery reads them but does not regenerate them.
The rounded logarithms and those imported error bounds are kept explicitly
in the new script. Extra decimal-input rounding is allowed when propagating
errors. The original density values are not changed.

Three interpretation corrections matter.

For distinct x<y<z, let c=(z-y,x-z,y-x). Then

    sum c_i = 0,      sum c_i*w_i = 0,
    beta_eff = sum c_i*log(nu_i) / [-sum c_i*log(w_i)].             (1)

Both a constant amplitude and a linear mass cancel. Thus the valid unequal
window (2,4,8) gives

    beta_eff = (2 log nu_2 - 3 log nu_4 + log nu_8)/log 2,

namely 0.718245786915 for NN and 0.533377433803 for matching. Applying the
weights (1,-2,1) to (2,4,8), as the return did, retains -2*kappa. Its large
negative numbers are not a prefactor diagnostic and must not be used as
evidence that the model fails at those widths.

Also, even an EXACT nu_w=A*w^(-beta)*exp(-kappa*w) has

    -(log nu_v-log nu_u)/(v-u)
      = kappa + beta*log(v/u)/(v-u).                             (2)

Drift of this adjacent difference is therefore not, by itself, evidence for
additional corrections. And using an identical algorithm on two graphs does
not force their exponents to be identical; that is a structural hypothesis,
not an algorithmic theorem.

A post-return sensitivity calculation, with beta fixed at 1/2 and
log nu=log A-kappa*w-(log w)/2+c1/w, interpolates the three values exactly:

    c1_NN       = 1.010055636135...,
    c1_matching = 0.029179857196... .                            (3)

Three free coefficients for three values leave **zero residual degrees of
freedom**. This is not supporting evidence for beta=1/2. The four available
widths 2,4,8,12 also permit the contrast with weights (-4,15,-20,9), which
annihilates 1,w,1/w. It gives approximately 0.8782607830 and 0.4797263407.
Those are another sensitivity view, not more independent observations or a
reason to commission another width automatically.

The #740 return at commit `787d5d5d40539f115dda23b68be6ff03ee3c82d9` usefully
keeps the site-component sewing unidentified, but its regularity inference
must be corrected. CIV Theorem A gives analyticity in direction on S^(d-1),
and Theorem B concerns the Wulff boundary. Their printed page 13 additionally
mentions inverse-temperature analyticity while deferring its discussion to a
future paper. DM Theorem 4.9 gives directional convexity at fixed p; its model
is bond FK. These statements do not supply the required site-p analyticity
proof. The preceding possible exceptional d set was about nondifferentiability
of kappa_4 or kappa_8 at their inverse images, not zeros of a(d)+b(d)-1.
Both corrections were posted to the existing issues; no result file was erased.

## 2. Complete-component activity has a local three-column weight

Work on the actual infinite cylinder Gamma_w=(Z/wZ) x Z, w>=3. G is either
NN or matching NN+NNN, with distinct lifted edges retained. A finite nonempty
connected vertex set C has probability

    q_p(C) = p^|C| (1-p)^|partial_G C|                            (4)

of being exactly one complete occupied component. The external vertex
boundary consists of DISTINCT sites, not open-to-closed incidences. No
condition is imposed beyond this boundary. In particular the exterior need
not be connected, and entire distant guard rows are not forced to be vacant.

Let S_i={y : (i,y) in C}, i modulo w. For a set S of integer rows write

    E_4(S)=S,
    E_8(S)=S union(S-1)union(S+1).

Define

    B_i = [(S_i-1) union(S_i+1)
           union E_G(S_(i-1)) union E_G(S_(i+1))] minus S_i.     (5)

Then, for EVERY finite set C, whether or not it is connected or winding,

    |C| = sum_i |S_i|,       |partial_G C| = sum_i |B_i|,
    q_p(C) = product_i phi_p(S_(i-1),S_i,S_(i+1)),
    phi_p(A,B,C)=p^|B| (1-p)^|B_i(A,B,C)|.                     (6)

**Proof.** An external neighbour in column i can be adjacent to a site in
column i or one of its two neighbours only. Formula (5) is precisely that
union with the occupied central sites removed. Each external site belongs to
one physical column, so the sum counts it once. Product independence then
gives (4) and (6). Diagonal steps only change E_4 into E_8. End of proof.

This resolves a *weight-bookkeeping* obstruction: the exact local factor is
available. It does not make successive irreducible objects independent. In
column language it is a three-column interaction, equivalently a transfer
on pairs (S_(i-1),S_i). The constraint that C is one connected winding
component remains a global predicate, or must be retained by connectivity
and lift information in an enlarged state. There is no claim that summing
unconstrained column words gives the desired component density.

A six-site example on w=4 consists of a complete row at y=0 and teeth at
(0,1),(2,1). It has 8 distinct external boundary sites but 12 occupied-to-void
incidences. Its true activity is p^6*(1-p)^8, not p^6*(1-p)^12. The error is
already present without any continuum limit or multiple winding paths.

### 2.1 Weight allocation is a gauge choice, not a new physical amplitude

At bounded height, write a pair-state transfer as

    K_(A,B),(B,C)=phi_p(A,B,C).

For any positive function d(A,B), replace it by

    K'_(A,B),(B,C)=d(A,B)*K_(A,B),(B,C)/d(B,C).                  (7)

This is a diagonal similarity transformation. Every cyclic product is
unchanged by telescoping, including a product restricted by the same global
connectivity predicate. Open endpoints acquire the corresponding d factors.
Thus moving boundary weights between adjacent pieces can change an open
endpoint amplitude without changing the closed object. It cannot justify
identifying an open two-point amplitude with a closed component amplitude.
This elementary gauge identity is a tool, not an asserted CIV intertwiner.

## 3. A canonical exact height expansion for the desired cylinder density

At each p<1 an entire empty row has positive probability. Independent such
rows occur arbitrarily far in both directions, so every component on a
fixed-width cylinder is finite. Anchor a winding component by its minimum
row. Let nu_(w,<=H) be the expected number anchored at row 0 with span at most
H rows. Then

    nu_(w,<=H) =
      sum_{C connected,winding; min_y C=0,max_y C<H} q_p(C),
    nu_w = lim_{H->infinity} nu_(w,<=H).                        (8)

All terms are nonnegative. The anchor rate equals the retirement rate used
by the supplied one-frontier engine: stationarity transports one mark per
finite component from its bottom to its retirement row. No winding path is
counted multiple times.

For fixed w a coarse, explicit bound is

    0 <= nu_w-nu_(w,<=H)
      <= w [1-(1-p)^w]^H.                                      (9)

Indeed a component anchored at row 0 and reaching row H requires every row
1,...,H to be nonempty; there are at most w anchors at row 0. Empty rows
separate both adjacency types. This bound is valid but becomes very poor as
w increases; it is NOT a uniform Ornstein–Zernike remainder.

Now define the unanchored strip activity

    Xi_(w,H) = sum_{C connected,winding; C subset rows[0,H)} q_p(C),
    Xi_(w,0)=0.                                                (10)

Its external boundary is still the boundary in the INFINITE cylinder, not a
free-boundary graph of H rows. A shape of vertical span h can be placed at
H-h+1 heights in this strip. Therefore

    Xi_(w,H) = sum_{h=1}^H (H-h+1)*nu_(w,span=h),
    Xi_(w,H)-Xi_(w,H-1) = nu_(w,<=H),
    Xi_(w,H)-2 Xi_(w,H-1)+Xi_(w,H-2) = nu_(w,span=H).           (11)

These are exact finite differences. A strip calculation with the right
external boundary weights can therefore recover the density by a height
difference. It need not guess a unique seam cut or a factor w/n.

### 3.1 A nontrivial finite-height sewing is explicitly closed

At H=1 the only connected winding set is the full horizontal row:

    Xi_(w,1) = [p(1-p)^2]^w.                                  (12)

At H=2 use nonempty column symbols 1,2,3. On NN, adjacent symbols must overlap;
on matching, all adjacent nonempty symbols communicate. In either case
these local compatibility rules are equivalent to the union being connected
and winding. For NN an incompatible interface has no occupied crossing, so
winding is impossible. If all interfaces overlap, either the word is a
single constant row or a double-occupied column joins all strands. For
matching, any two occupied sites in consecutive two-row columns are adjacent.

There are 7 allowed pair states on NN and 9 on matching. Use the weight in
(6) on pair-state transitions. Then

    Xi_(w,2)=trace K_(H=2)(p)^w,
    nu_(w,<=2)=trace K_(H=2)(p)^w - [p(1-p)^2]^w.               (13)

The new script verifies this identity against actual anchored component
activity at w=3,4,5 and p=1/4,1/2 for BOTH graphs. It is an exact finite-height
representation of actual site components, not the unknown all-height renewal
sewing. Fixed H has finitely many transverse states; its leading large-w
behaviour need not have the all-height w^(-1/2) factor. Exchanging the two
limits is a substantive step.

## 4. Multiple seam marks are removable exactly, not a no-go theorem

Fix the seam between columns w-1 and 0. For each winding component C let
M(C) be any nonempty finite set of specified seam marks, and c(C)=|M(C)|.
Marks may be seam edges, or distinct seam rows, but the convention cannot be
changed halfway. A raw seam edge is not automatically an OZ regeneration cut.

For every single configuration,

    sum_{e in M(C)} 1/c(C) = 1.                                (14)

Let component-Palm mean selecting a component proportionally to its activity
per anchor row, and mark-Palm mean selecting a mark proportionally to its
activity. If nu is component intensity and mu is marked intensity, then

    mu = nu * E_component[c],
    dP_mark(C) = c(C)/E_component[c] * dP_component(C),
    nu = mu * E_mark[1/c].                                    (15)

Importantly,

    E_mark[1/c]=1/E_component[c]

but in general E_component[1/c] is different. Multiplying mu by the latter
produces a bias; Jensen's inequality puts that bias above the correct value.
Formula (15) follows from sums in a finite height window, or from stationary
mass transport in the infinite cylinder. It requires the same model, support,
marking rule and sampling law on both sides.

There is also an exact generating-function form. If

    Z_H(t)=sum_C q_p(C)*t^{c(C)},    0<=t<=1,

over the anchored finite-height class, then

    nu_(w,<=H)=integral_0^1 Z_H'(t) dt.                         (16)

The inverse-mark factor is the identity 1/c=integral_0^1 t^(c-1)dt. Multiple
cuts are therefore not an impossibility result for sewing. They demand a
specified mark law and a correct unrooting operator. Whether an OZ chain's
marks agree with the raw seam marks remains to be proved.

### 4.1 A true complete-component small control

For NN, w=4, span at most 3, p=1/2, the exact anchored density and marked
intensity for seam edges are

    nu = 9087/1048576,
    mu = 5601/524288,
    E_component[c] = 11202/9087,
    E_mark[1/c] = 9087/11202.                                  (17)

The product in (15) recovers nu exactly. Using E_component[1/c] instead
would overestimate it by about 9.53%. This is a truncated full-component
activity with the exact external void weight, not a free strip with cut
components retained. The full infinite-height density is NOT (17).

Finite numbers such as 1.5 or 2.1 at widths at most four prove neither tightness
nor growth of c in the asymptotic regime, and do not by themselves determine
an amplitude correction to a differently marked renewal model.

## 5. A modest centre-order consequence, distinct from regularity

Using the preceding definitions of a(d)=kappa_4^(-1)(d) and
c(d)=1-b(d)=kappa_8^(-1)(d), and their continuous strictly decreasing masses,
monotone graph inclusion gives kappa_8(p)<=kappa_4(p) whenever both are
subcritical. If a(d)>=pc_8 then c(d)<a(d) immediately. Otherwise comparing
masses at a(d) gives c(d)<=a(d). Thus

    a(d)+b(d)>=1.                                             (18)

This uses the previous mass-inversion theorem as an input and does not prove
strict inequality at every d. Strictness is a plausible graph-enhancement
question, not the definition of an exceptional Gumbel parameter. In particular
c=a only makes the black and white *occupation probabilities* equal; it does
not make their original-label thresholds a and 1-a coincide unless a=1/2.
No p-analyticity theorem is obtained by this observation.

## 6. A research conjecture: heat-kernel sewing rather than independent pieces

The following is deliberately recorded as a conjectural continuation, not an
accepted result. It remains within the same density problem and does not
request a fourth width or another source programme.

### 6.1 Why this formulation is preferable

The actual site weights have the local memory (5). The correct candidate is
therefore a cyclic Gibbs/Markov-renewal description that retains overlap state,
not necessarily an independent sequence of geometrical blobs. Multiple raw
seam crossings need not be excluded; their weights must be unrooted as in (14).

**Conjecture HK (fixed-subcritical closed-component scaling).** For each
adjacency G and fixed 0<p<pc_G, after the microscopic overlap and marking
weights are included, there are kappa_G(p), D_G(p)>0 and zeta_G(p)>0 such that
for H/sqrt(D_G(p)w)->h in (0,infinity),

    Xi_(w,H)(p)
      ~ zeta_G(p) exp[-w kappa_G(p)]
         sum_{n>=1} exp[-pi^2 n^2 D_G(p) w/(2H^2)].             (19)

This is the Dirichlet heat trace of ONE diffusive transverse mode. D is the
long-time transverse diffusion coefficient, including inter-piece correlations,
not merely a one-step variance. Zeta contains the still-unidentified
microscopic sewing/mark normalization. These quantities are not fitted or
computed for the real site process in this delivery.

The equivalent anchored version is safer as a target, because relative
asymptotics of Xi alone do not justify taking a finite difference. Require
in addition enough uniform remainder control, or directly conjecture

    nu_(w,<=H)(p) / nu_w(p) -> F_range(h),
    nu_w(p) ~ zeta_G(p) exp[-w kappa_G(p)]/sqrt(2*pi D_G(p)w).   (20)

A claim of (19) without finite-difference control is NOT a proof of (20).

The limiting cutoff function is the range distribution of a standard
Brownian bridge, with two equivalent series:

    F_range(h)
      = sqrt(2*pi)*pi^2/h^3 * sum_{n>=1} n^2 exp[-pi^2 n^2/(2h^2)]
      = 1 + 2 sum_{n>=1}(1-4n^2h^2) exp[-2n^2h^2].             (21)

The first converges well at small h, the second at large h. This formula is
an exact Brownian identity; its application to SITE components is conjectural.

**Derivation of the Brownian target.** For a unit-time standard bridge with
range R, integrating the free bridge kernel over starting positions whose
translated bridge fits in (0,h) gives

    trace exp(Delta_(0,h)/2)
      = E[(h-R)_+]/sqrt(2*pi)
      = sum_{n>=1}exp[-pi^2 n^2/(2h^2)].                       (22)

Differentiate in h to get the first series in (21). Poisson summation gives
the second. The discrete identities (11) are exactly the height-counting
analogue of this derivative. This supplies a concrete geometric reason for
the *anchored* w^(-1/2) factor, not just a matching of powers.

### 6.2 What could falsify this conjecture

A second soft transverse mode, slowly decaying interaction between irreducible
pieces, a nontrivial w-dependent unrooting factor, or a different limit shape
of complete winding components can invalidate (19)–(20). Raw seam multiplicity
larger than one does not do so by itself. Conditional support at fixed H also
does not test the Brownian H/sqrt(w) regime.

A specific sufficient marking hypothesis worth examining is tightness with
an exponential moment of c under the component-Palm law, uniformly at fixed
p as w grows, plus convergence of its relevant joint law with the sewing
state. This would preserve the power when converting between mark and component
intensities. It would not identify the constant without that joint law.
This marking hypothesis is **not proved here** and is not implied by a few
small free-boundary means.

Targets from (21) are approximately F_range(1)=0.178, F_range(2)=0.98994.
The result JSON contains higher precision evaluations of the explicit series.
They are not measurements of a percolation cluster.

### 6.3 A geometry prediction that does not fit a mass or an amplitude

If the range convergence in (20) also holds with its first two moments
(the needed uniform integrability is part of this conjecture, not automatic),
let L(C)=max_y C-min_y C+1 for the complete component selected under the
component-Palm law. The Brownian bridge range has

    E R=sqrt(pi/2),       E R^2=pi^2/6.

The first identity follows by reflection from the maximum of a bridge and
symmetry of its minimum; the second follows by integrating (21), using an
absolutely convergent second-moment series. Consequently HK predicts

    E_component L / sqrt(Dw) -> sqrt(pi/2),
    E_component L^2/(Dw) -> pi^2/6,
    Var_component(L)/(E_component L)^2 -> pi/3 - 1
                                            = 0.0471975511966... .       (23)

This last ratio is independent of kappa, D and zeta. It tests the proposed
GEOMETRY of a complete winding component, not merely another fit of density
versus width. A seam-marked sample must first be debiased by (15); otherwise
it probes a different distribution. The existing free-strip seam means are
not this statistic. No new acquisition is commissioned here, and (25) is
not scored against fixed tiny widths as an asymptotic theorem.

### 6.4 A more distant, explicitly conditional critical crossover

A possible extension on approaching pc is

    nu_w^G(p) ~ w^(-1) F(w*kappa_G(p)),                         (24)

after the appropriate lattice metric is fixed. If D_G(p)*kappa_G(p)->D_*>0
and zeta_G(p)->zeta_*>0, (20) would match a large-x tail

    F(x) ~ zeta_* sqrt[x/(2*pi*D_*)] exp(-x).                   (25)

A finite nonzero F(0) would describe a critical cylinder density. Neither the
existence of this scaling function nor its critical value, common metric,
or equality between NN and matching functions is proved here. Fixed-p OZ
estimates are not uniform critical estimates. Moreover the same-label black/
white winding counts obey the exact topology constraint from the earlier
handoff, so two independent Poisson processes cannot simply be continued into
this common window.

Equations (19)–(25) are stored so they can guide a later proof or a specific
comparison; they are not to be copied into the claim ledger as data or theorems.

## 7. Executed scope

The standalone script implements physical lifted BFS and independent gain
union/find, the literal external-neighbour set and an independent column-union
formula, component-Palm/mark-Palm arithmetic, and the two-row trace.

- 11,938 nonempty shape masks across both adjacencies: all boundary sets agree
  and both winding detectors agree. Connected anchored winding shapes are
  then weighted with (4).
- 17,408 actual extended-strip configurations with random guard rows:
  direct complete-component expectations agree with the anchored activity
  sum at p=1/3. Unrelated guard sites are integrated out, not forced vacant.
- Twelve exact two-row trace comparisons at widths 3,4,5 and p=1/4,1/2.
- The returned #741 logs are rescored, not regenerated by a new large solve.
  All displayed long decimals in the JSON beyond input precision are arithmetic
  outputs, not extra precision of the supplied densities.

Run from the repository root:

    python -m unittest discover -s tests -p 'test_cluster_sewing_identity.py' -v
    python scripts/cluster_sewing_identity.py --output /tmp/sewing-new.json

The output path must not exist. No Monte Carlo or full repository CI is included.
No theorem in Section 6 is inferred from these finite checks.

## Sources and integration references

- #740 return: `787d5d5d40539f115dda23b68be6ff03ee3c82d9`, issue comment
  `5651941180`; keep its finite seam-count observations, correct the
  p-regularity and marking inferences.
- #741 return: `3745b13b8e1126017a1e88567a6af44881b8a1ff`,
  `results/geometric-consistency/winding-prefactor-contrast.json`.
- Campanino–Ioffe–Velenik, *Fluctuation theory of connectivities for
  subcritical random cluster models*, Ann. Probab. 36 (2008), arXiv:math/0610100v2.
  Primary PDF Theorems A/B and printed page 13 read; pages 10 and 13 rendered.
  Directional regularity and the deferred temperature remark are distinguished.
- D'Alimonte–Manolescu, arXiv:2510.13648v3, primary HTML sections 1.1–1.2
  and Theorem 4.9 read. Bond-FK model and direction-versus-p distinction checked.

These are bounded source readings, not a completed literature-priority search.
The exact bookkeeping identities here are elementary derivations; no novelty
claim for mass transport, Gibbs transfers, Palm debiasing or heat kernels is made.
