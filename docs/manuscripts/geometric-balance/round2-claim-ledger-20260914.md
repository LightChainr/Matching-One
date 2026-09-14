# 2026-09-14 continuation claim ledger

Status ledger for stacked Draft PR #771.  The purpose is to keep exact topology, author-level probability proofs, conditional consequences and research conjectures visibly separate while the branch is still moving.

This ledger covers the main claims produced in the geometric-balance continuation.  Other contemporaneous branch files should be reviewed on their own terms if they are not named here.

## Status vocabulary

- **EXACT**: finite deterministic/probability identity, or elementary finite-state/convex consequence with no asymptotic model input beyond stated definitions.
- **AUTHOR-PROOF**: a complete proof is written on the branch, but it uses published inputs and/or author-level results from #739 and still needs independent proof review before publication.
- **CONDITIONAL**: deduction is complete once one explicitly named hypothesis/theorem interface is granted.
- **CONJECTURE / PROGRAMME**: falsifiable research mechanism or proposed proof route, not promoted as a theorem.
- **LITERATURE BOUNDARY**: audited statement about what an external theorem does and does not imply for square-site Matching-One.

## A. Exact finite topology and persistent birth structure

| Claim | Status | Main file | Consequence |
|---|---|---|---|
| Digital-Alexander birth reflection `T1_8(1-U)=1-T2_4(U)`, `T2_8(1-U)=1-T1_4(U)` | EXACT | `structural-consequences-20260914.md` | Entire 4/8 birth process is paired samplewise, not only in law. |
| In rank one, complementary NN/matching homology lines coincide | EXACT | `projective-homology-gas-20260914.md` | Same projective slope on the two colours. |
| Rank-one slope is born at `T1` and remains frozen until `T2` | EXACT | `persistent-slope-marked-birth-20260914.md` | Direction is a persistent first-birth mark; it cannot rotate during the plateau. |
| `(G,L,C)_8 =d (G,L,-C)_4` with `G=T2-T1`, `C=(T1+T2-1)/2` | EXACT | `persistent-slope-marked-birth-20260914.md` | NN/matching share the full plateau-width/slope law; only the centre coordinate is complement-odd. |
| Same-parameter count support is only `(0,1)`, `(1,0)`, `(K,K)` | EXACT | `structural-consequences-20260914.md` | Generic bivariate Poisson/copula models are inadmissible. |
| Rank-one state is `(slope,K)` and different slope species have global hard-core exclusion | EXACT | `projective-poisson-hardcore-crossover-20260914.md` | Multi-direction models must enforce torus intersection topology. |
| Neutral count fluctuations are one-dimensional; `W4-W8` is bounded | EXACT | `neutral-gas-limit-collapse-20260914.md` | Every growing-scale joint fluctuation limit lies on the diagonal; extensive LDP rate is `+infinity` off diagonal. |
| Finite reflection dominance `1-T2 <=st T1` | EXACT | `finite-reflection-dominance-20260914.md` | `P2(1-p)<=P0(p)`, `M(p)+M(1-p)<=0`, finite root `>=1/2`. |

Finite controls: `scripts/persistent_alexander_birth_reflection.py` exhausts all `L=3` configurations and all `9!` strict orders; at `L=4` it exhausts configurations and checks 20,000 fixed-seed orders, with zero violations of the recorded identities.

## B. Exact common-window charge coordinates

| Claim | Status | Main file | Consequence |
|---|---|---|---|
| `chi=P0+P2`, `theta=log(P2/P0)`, `M=chi tanh(theta/2)` | EXACT | `charge-neutral-crossover-coordinates-20260914.md` | Root location and charged-sector rarity are separated coordinates. |
| Matching root is exactly `theta=0` | EXACT | same | Common-window modelling should evolve a charge fugacity, not two independent colour means. |
| Joint PGF `G(s,t)=chi[e^{theta/2}s+e^{-theta/2}t]/(2cosh(theta/2))+(1-chi)H(st)` | EXACT | same | Minimal same-parameter state is one scalar susceptibility, one scalar fugacity, and one positive-integer neutral count law. |
| At the root, `M'=chi theta'/2`, `F'=chi theta'/4` | EXACT | same | Algebraic mechanism for balance without concentration. |
| `chi_8(p)=chi_4(1-p)`, `theta_8(p)=-theta_4(1-p)`, neutral law transported unchanged | EXACT | `charge-fugacity-reflection-20260914.md` | Exact 4/8 involution in crossover coordinates. |
| Graph inclusion gives `theta_8(p)>=theta_4(p)`, hence `theta_4(p)+theta_4(1-p)<=0` | EXACT | same | Strong charge-field reflection constraint on any scaling ansatz. |

## C. First-exit / Wulff geometry

| Claim | Status | Main file | Review dependency |
|---|---|---|---|
| `B_S(t)<1` certifies finite exponential susceptibility | AUTHOR-PROOF | `vector-first-exit-domain-20260914.md` | Site-BK first-exit skeleton. |
| True exponential-moment domain is the interior of the polar/Wulff body; large finite boxes exhaust compact interiors | AUTHOR-PROOF | same | Standard uniform directional subcritical exponential estimate. |
| Certified bodies may be convex-hulled; support function, not radial intercept, equals `tau` | EXACT / AUTHOR-PROOF | same | Convexity is exact; domain equality uses the preceding input. |
| Torus winding upper bound `P(r>0)<=N C_T sum_{lambda!=0} exp[-h_T(lambda)]` | AUTHOR-PROOF | `first-exit-torus-winding-upper-20260914.md` | Canonical minimum homology witness is a quotient-vertex-simple cycle; local first-exit witnesses are disjoint SITE variables. |
| With `T=(1-eps)K_p`, full-period bound uses complete period cost `(1-eps)tau_p(lambda)` | AUTHOR-PROOF | same | Depends on Wulff-domain exhaustion. |
| Correlation-norm successive minima organize rank/slope suppression | AUTHOR-PROOF | `correlation-norm-successive-minima-20260914.md` | Deterministic norm-ball packing + torus upper bound. |

Important correction preserved in the branch: a direct global comparison of two quotient cycle arcs with two independent full-plane connection events is **not** used; periodic lifts can reuse the same Bernoulli variable.  The local first-exit construction is the repaired route.

## D. Arbitrary-shape fixed-p homological free energy

For fixed subcritical `p`, put

`rho_n = min_{lambda in Lambda_n\{0}} tau_p(lambda)`.

If `rho_n->infinity` and `log N_n/rho_n->alpha`, then

`(1/rho_n) log P_p(r>0) -> -max(1-alpha,0)`.

Equivalently,

`log P_p(r>0)=-(rho_n-log N_n)_+ + o(rho_n)`.

Status: **AUTHOR-PROOF** in `general-period-homological-free-energy-20260914.md`.

Upper side: first-exit theta bound + correlation-norm lattice packing.  Lower side: finite angular net of fixed local connection seeds + deterministic correction + Harris association + finite-group translation packing.  No fixed direction or OZ prefactor is required.

This promotes the earlier heuristic `connection energy - log opportunities` to a two-sided fixed-parameter theorem at exponential scale.

## E. Directional centres and mass monotonicity

| Claim | Status | Main file |
|---|---|---|
| Fixed primitive direction: `(1/n)log P(r>0)=-max{tau_p(u)-d,0}` for `Lambda=<nu,mv>`, `log m/n->d` | AUTHOR-PROOF | `fixed-direction-exponential-centres-20260914.md` |
| Varying shortest directions `u_n/|u_n|->e`, `log(N/|u_n|)/|u_n|->d`: rate is `-max{tau_p(e)-d,0}` | AUTHOR-PROOF | `varying-direction-exponential-centres-20260914.md` |
| Any nonparallel period in the exponential shortest-period regime is at least transverse height `h=N/ell`, hence direction competition is superexponentially suppressed | EXACT geometry + AUTHOR-PROOF probability | `exponential-homology-class-selection-20260914.md` |
| First-birth slope equals the shortest-period projective line with probability `->1` in that regime | AUTHOR-PROOF | same |
| `p -> tau_{G,p}(e)` strictly decreases for every direction | AUTHOR-PROOF | `varying-direction-exponential-centres-20260914.md` |
| Quantitative FK comparison `tau_p(e)-tau_q(e) >= (q-p)tau_q(e)/rho_FK` uniformly in direction | AUTHOR-PROOF | `direction-uniform-mass-slope-20260914.md` |
| Lower/upper centres satisfy `tau_4,a(e)=d`, `tau_8,1-b(e)=d` | AUTHOR-PROOF | `varying-direction-exponential-centres-20260914.md` |

Issue impact: the fixed/varying direction centre problem in #765 is reduced to proof review rather than a missing direction case.

## F. Strict NN/matching directional separation

| Claim | Status | Main file |
|---|---|---|
| Matching facial enhancement beats a positive NN `p` sprinkling for two-terminal events | AUTHOR-PROOF | `matching-enhancement-mass-gap-20260914.md` |
| `tau_8,p(e)<tau_4,p(e)` for every direction and `0<p<pc(G8)` | AUTHOR-PROOF | `varying-direction-exponential-centres-20260914.md` + enhancement |
| On compact `p` intervals, dual bodies have a positive Euclidean buffer `K_8,p + eta B_2 subset K_4,p` | AUTHOR-PROOF / convex consequence | `strict-wulff-body-separation-20260914.md` |

Primary review hotspot: finite-to-one pivotal bookkeeping / endpoint relocation in the two-terminal adaptation of Grimmett--Li using the corrected Balister--Bollobas--Riordan square-lattice rerouting theorem.

## G. No-prefactor Poisson--Gumbel windows

| Claim | Status | Main file |
|---|---|---|
| Axial complete-component Poisson process and median-centred Gumbel law at regular `d` | parent AUTHOR-PROOF | `poisson-birth-windows.md` |
| Same proof extends to any **fixed primitive integer direction** after exact `SL_2(Z)` straightening to a fixed anisotropic finite-range SITE graph | AUTHOR-PROOF | `fixed-direction-poisson-gumbel-20260914.md` |
| Varying-direction centre theorem does **not** automatically imply varying-direction Gumbel theorem | EXACT scope boundary | same |

The obstruction for moving directions is local rather than energetic: after straightening, the transformed interaction range can grow with the arithmetic of `u_n`, so guard-slab, dependency-neighbourhood and Palm-volume constants need new uniform control.

The at-most-countable exceptional set is retained: no direct whole-subcritical square-SITE theorem proving `p`-analyticity of the mass was identified.

## H. Alternating black/white components (#764)

| Claim | Status | Main file |
|---|---|---|
| Alternating essential regions + component correspondence | AUTHOR-PROOF from existing digital-Alexander CW construction | structural notes / `poisson-tessellation-consequence-20260914.md` |
| `nu E_W L -> 1` needs stationarity/Palm mean spacing, not Poisson | AUTHOR-PROOF | same |
| Poisson anchor process upgrades normalized white span to `Exp(1)` and uniform-location span to `Gamma(2,1)` | AUTHOR-PROOF conditional on parent process theorem | `poisson-tessellation-consequence-20260914.md` |
| Infinite-volume boundary density `beta(q)=(1-q)theta(q)/q` | EXACT | `supercritical-white-slab-bulk-20260914.md` |
| Mean rewards `(nu E L, nu E N/w, nu E B/w)->(1,theta_8(q),(p/q)theta_8(q))` | AUTHOR-PROOF | `white-component-mean-reward-20260914.md` |
| Samplewise `N/(wL)->theta`, `B/(wL)->(p/q)theta` and conditional slab CLT | CONJECTURE | `supercritical-white-slab-bulk-20260914.md` |

The mean boundary/volume ratio `E B/E N -> p/(1-p)` is now theorem-level on the branch even though the stronger samplewise slab LLN remains open.

## I. Fixed-width charge free energy: why the root survives a wide plateau

For fixed circumference `w`, let `Q^G_{w,m}(p)` be the open-strip probability of no horizontal essential component.  Concatenation and an empty separator row give a free energy

`I^0_{G,w}(p)=lim_{m->infinity}-(1/m)log Q^G_{w,m}(p)`.

The torus rank-zero probability has the same rate.  Digital Alexander gives

`theta_{w,m}(p)/m -> Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p)`.

Status: **EXACT free-energy existence + AUTHOR-PROOF transfer/Perron regularity** in `fixed-width-charge-free-energy-20260914.md`.

Consequences:

- unique fixed-width charge-coexistence point `p_w^ch` from `Theta_w=0`;
- finite-torus roots `p^*_{w,m}->p_w^ch` as `m->infinity`;
- a diagonal use of the parent root-consistency theorem gives `p_w^ch->pc` as `w->infinity`.

Interpretation: the matching root balances two rare **topological void-sector free energies**.  It is not balancing the black/white complete-component intensities, which are equal in the infinite-cylinder alternation identity.  This is the spectral explanation of balance without birth-law concentration.

Review hotspot: primitive Perron block / strict monotonicity for the fixed-width safe-homology transfer should be audited if this is promoted as a production estimator.

## J. Projective slope / modular positive control

| Claim | Status | Main file |
|---|---|---|
| Embedded even-spin projective slope harmonic is a literal lattice observable | EXACT | `projective-slope-harmonic-control-20260914.md` |
| 4/8 complement transports the full harmonic exactly | EXACT | same |
| Pinson--Arguin primitive-sector baseline yields a parameter-free critical continuum harmonic | EXACT given merged PR #213 formula | same |
| Modular covariance `A_s(gamma tau)=(|c tau+d|/(c tau+d))^s A_s(tau)` | EXACT sector-reindexing consequence | `projective-slope-modular-covariance-20260914.md` |
| Elliptic-point selection rules, e.g. `A_4(rho_hex)=0` | EXACT continuum symmetry consequence | same |

This is explicitly a positive-control channel for #585/#589; it is **not** original-U and must not be inserted as its surrogate.

## K. Sewing, prefactor and locality

| Claim | Status | Main file |
|---|---|---|
| Finite-state cyclic matrix kernel with one simple Perron band gives `e^{-kappa w}/sqrt(2 pi D w)` and unit logdet residue | EXACT matrix theorem | `matrix-sewing-unit-residue-20260914.md` |
| Finite local memory alone cannot create a continuously varying extra residue | EXACT within that matrix class | `sewing-amplitude-diagnostic-20260914.md` |
| Pure periodization preserves zero-Fourier Perron mass exactly; nonzero `gamma_w-kappa` must come from wrap-sensitive decorations/irreducible pieces | EXACT mechanism statement / CONDITIONAL SITE conclusion | `periodic-mass-locality-mechanism-20260914.md` |
| Exponential decoration locality would imply `gamma_w-kappa=O(e^{-cw})` | CONDITIONAL | same |

This reduces #760 to a concrete decoration-locality lemma rather than another width fit.

## L. Loop/branch morphology and #758/#762

| Claim | Status | Main file |
|---|---|---|
| For candidate `I(A)=min_r[tau(1,2r)-kappa+kappa(A-r)]`, strict convexity gives fixed bulge saturation and an exactly linear tail | CONDITIONAL on the candidate variational formula | structural notes / `loop-branch-linear-tail-tests-20260914.md` |
| Small-`A` overlap forces `D^{-1}=partial_yy tau(1,0)=kappa+kappa_angle''(0)` | CONDITIONAL consistency relation | same |
| Since `r_*<1/2`, #762 tests at `A=1,2` are necessarily in the linear-branch regime if the candidate is correct | CONDITIONAL | `loop-branch-linear-tail-tests-20260914.md` |
| Under near-critical normalized Wulff isotropy, `r_*=1/(2sqrt3)` and the whole normalized rate is explicit | CONDITIONAL / UNIVERSALITY HYPOTHESIS | `near-critical-isotropic-loop-branch-rate-20260914.md` |

The explicit isotropic target is

`I(A)/kappa -> sqrt(1+4A^2)-1` below `1/(2sqrt3)`, and `A+sqrt(3)/2-1` above it.  In particular `I(1)/kappa->sqrt(3)/2` and `I(2)-I(1)->kappa`.

Square-site rotational restoration is **not** claimed as a theorem.

## M. Near-critical literature boundary and next theorem engine (#740/#767)

`near-critical-oz-literature-boundary-20260914.md` is **LITERATURE BOUNDARY**, not a SITE theorem.

Audited result: D'Alimonte--Manolescu v3 proves for square-lattice random-cluster (q=1 = Bernoulli **bond**) a uniform near-critical two-point comparability

`P(0<->r e) asymp pi_1(xi)^2 sqrt(xi/r) exp(-r/xi)`.

What transfers safely as architecture:

- renewal blocks at correlation-length scale;
- Gaussian sewing / Brownian transverse scale `D asymp xi`;
- strict Wulff geometry;
- explicit endpoint arm insertions.

What does **not** transfer as a theorem:

- square-site uniform near-critical OZ;
- exact relative-error amplitude;
- complete-component insertion residue;
- square-site critical arm exponents.

`near-critical-loop-insertion-diagnostic-20260914.md` is **CONJECTURE / DIAGNOSTIC**.  It separates:

- pure cyclic closure: `nu_w ~ xi^{-1}s^{-1/2}e^{-s}`;
- two-endpoint-dressed closure: `nu_w ~ pi_1(xi)^2 s^{-1/2}e^{-s}`;

with `s=w/xi`.  The statistic `J_emp=xi sqrt(s)e^s nu_w` is designed to distinguish insertion semantics near criticality; fixed-`p` data cannot.

`site-near-critical-renewal-programme-20260914.md` is **PROGRAMME**.  It isolates five square-SITE gates: sub-characteristic 4/8 RSW, one-arm stability (or an explicit near-critical arm insertion), uniform killed renewal/barriers, nondegenerate transverse variance, then complete-component cyclic insertion.  This is intentionally smaller than proving conformal universality.

## N. Common-window crossover (#767)

What is now exact:

- support and projective slope constraints;
- charge coordinates `(chi,theta,H)`;
- 4/8 reflection and graph-inclusion inequalities;
- separated Poisson-window boundary conditions.

What remains genuinely open:

- a square-site near-critical scaling law for `(chi,theta,H)` in the merging regime `rho->infinity`, `log rho=o(w)`;
- the complete-component insertion factor at correlation-length scale;
- a proof-level square-site near-critical renewal theorem or an explicit universality input strong enough to replace it.

`projective-poisson-hardcore-crossover-20260914.md` is a **CONJECTURAL closure** only for geometries in which several nonparallel period classes genuinely remain competitive.  It should **not** be used in the exponential shortest-period regime, where `exponential-homology-class-selection-20260914.md` proves deterministic slope selection.

## O. Corrections / superseded shortcuts

1. **Quotient two-arc shortcut rejected.**  A torus simple cycle cannot be globally treated as two independent planar connections because periodic lifts can reuse Bernoulli variables.  Replaced by local first-exit BK witnesses.
2. **Same-p two-subcritical-free-energy phase diagram rejected.**  Because `pc(G8)=1-pc(G4)`, black NN at `p` and complementary matching at `1-p` are not both subcritical away from criticality.  `two-free-energy-rank-phase-diagram-20260914.md` now contains the corrected one-sided lower/upper formulation.
3. **Fixed-p OZ amplitude is not extrapolated to criticality.**  Near-critical bond-FK already shows a moving critical-arm insertion; square-site complete-component insertion remains separate.
4. **Direction centre does not imply direction window.**  Varying-direction centre is coordinate-free; varying-direction Gumbel still needs uniform local component controls.

## P. Issue-level status after this continuation

- **#765 directional centres:** main mathematical direction problem is author-level closed; review/acceptance remains.
- **#766 vector first-exit/Wulff:** main theoretical certificate and torus application are author-level closed; numerical certification can now be built on them.
- **#763 no-prefactor Gumbel:** axial parent proof retained; fixed primitive directions extended; varying directions remain a local-uniformity problem.
- **#764 black/white gap law:** reciprocal mean, Poisson gap laws and mean `(L,N,B)` rewards substantially closed; samplewise white-slab LLN/CLT remains.
- **#767 common near-critical crossover:** exact finite state space/coordinates are closed; square-site merging scaling law remains open.
- **#760 periodic mass locality:** reduced to decoration locality / Perron perturbation.
- **#758/#762 morphology:** candidate mechanism now has sharp linear-tail and conditional isotropic no-parameter tests, but the LDP mechanism itself remains conjectural.
- **#740 prefactor/sewing:** finite-memory sewing theorem and near-critical insertion alternatives are separated; actual SITE complete-component residue remains open.

## Q. Recommended review order

1. Persistent Alexander / projective slope exact identities and finite controls.
2. Vector first-exit theorem and torus cycle application.
3. General-period fixed-p free energy and direction-uniform FK slope inequality.
4. Matching enhancement pivotal conversion and all-direction strict mass gap.
5. Alternating black/white Palm mean results.
6. Fixed-direction Poisson--Gumbel extension.
7. Fixed-width charge free energy / Perron interpretation.
8. Only then review conditional near-critical and morphology mechanisms.

This ordering maximizes the amount of downstream material validated by each proof audit and keeps conjectural near-critical work from blocking the geometric/probability core.