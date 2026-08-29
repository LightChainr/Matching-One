# Scientific Frontier Map

**Updated:** 2026-08-29

This map is the repository's scientific coordinate system. It is not a permission system: exact work, reanalysis, pilots, production and independent theory may proceed in parallel. Priority moves attention; it does not lock, close, veto or demote a task. `STATUS` owns claim language and lifecycle boundaries, `ROADMAP` ranks attention, and `NEXT-TARGETS` specifies decision outputs.

## The five axes

Every comparison must name its state, source, observer, geometry and acquisition. Shared exponents, harmonics or lattice sizes do not make two observers equivalent.

| Axis | Live coordinates | Boundary that must remain visible |
|---|---|---|
| **state** | finite ambient homology; rank-one lifetime and birth filtration; continuum Q4/Jordan or other module; low-dimensional transfer state; microscopic threshold structure | an exact finite state is not yet a continuum-field identification |
| **source** | Bernoulli thermal; rank-birth geometric tilt; Potts-Q measure; confluent projector; explicit field derivative; relative fugacity; deck charge; boundary source | measure, projector and explicit-field Q derivatives are separate terms |
| **observer** | unmarked endpoint; rank-plane projector; `I01/I12` birth gate; matching-even/odd; line H4; landing H4; local pivotal; charged; boundary; cross-microscopic | line, landing, ordinary and charged rows are not interchangeable |
| **geometry** | Gaussian covers; annulus; modulus/Hecke children; boundary cross-ratio; triangular and square controls; finite transfer cylinder | a changed Smith class, modulus, path or marked context is part of the experiment |
| **acquisition** | existing-data reuse; exact construction; small targeted pilot; new production | correlated views of one stream remain one dependency group |

The typed Potts-Q interface is

```text
measure covariance
+ finite confluent-projector derivative
+ explicit bare-field derivative
+ boundary prefactor/gauge term where applicable.
```

## The main mechanism chain

The highest-information line is now a connected chain rather than a collection of H4 fits:

```text
Digital Alexander rank balance and index-2..7 quotient frontier   [main]
    -> essential H1 births K1,K2 and P0/P1/P2 reconstruction      [main]
    -> complement coordinates C (clock translation), W (lifetime) [branch_only]
    -> full finite K2 curve: positive area, three local node crossings [open_pr]
    -> local gates I01,I12 and parity channels S,D                 [branch_only]
    -> tiny nonzero Cov(A_top,J_D4) lattice coupling               [branch_only]
    -> P205 ordinary H4/H4 K1/K2 prism                             [branch_only + PR #267 reuse]
    -> common-field large-N coupling and Q4/Jordan fingerprint     [next]
```

### 1. Digital topology and essential births

On the declared regular square-cell class, Digital Alexander duality gives configurationwise

```text
r_black + r_white = 2,
q = r_black - 1 = A_top.
```

`main` now also contains the filtration oracle, the essential-birth histogram reconstruction and exhaustive HNF checks through quotient index 7. Across 40 quotient representatives and 49,878 filtration paths there are no birth, reflection, rank-sum, reconstruction or line failures. That is a strong finite frontier, not an unrestricted theorem for every degenerate quotient.

For a monotone site ordering,

```text
K1 = K_minus = first essential ambient-H1 birth,
K2 = K_plus  = second essential ambient-H1 birth,
P0 = 1-F1,  P1 = F1-F2,  P2 = F2,
M  = P2-P0 = -1+F1+F2.
```

The positive beta-density derivatives make the physical finite-volume matching root unique and simple. The neutral area is the mean rank-one lifetime. Existing endpoint histograms recover these marginals but do not contain the plateau line `ell`, Smith index `iota`, landing mark or same-sample global/local products.

### 2. Clock/lifetime coordinates and local birth parity

The paired activations have exact normalized coordinates

```text
C = (K1+K2)/(2N+2): complement-odd clock translation,
W = (K2-K1)/(N+1):  complement-even rank-one lifetime.
```

Across the ten archived direction pairs, the exact `Delta cos(4 theta)` normalization gives negative `C` and positive `W`. Those are correlated coordinates reconstructed from the same archives, not independent confirmations or fitted scaling exponents.

The Draft full-curve reuse sharpens the earlier three negative `K2` values at the pooled root. Its exact Bernstein-area identity gives the same positive integrated `K2` direction at all ten sizes. N=265/325/425 each have a delete-one-stable upper node near `.59462/.59635/.59694`, so the negative point value is a local lobe before a zero crossing rather than a reversal of the second activation as a whole. This finite node topology still cannot identify `ell`, `iota` or a continuum carrier.

The local insertion refines the two births:

```text
I01: rank 0 -> at least 1,
I12: rank at most 1 -> 2,
S = I01+I12 = Delta r and integrates to M',
D = I12-I01 = -partial_p P1.
```

A direct `0 -> 2` birth enters both gates, hence contributes `S=2` and `D=0`. Every nonzero `D` event has a canonical rank-one plateau line `ell`; `iota`, the projective line character `chi4(ell)` and the local landing H4 mark are additional typed fields. In particular,

```text
J_D4 = sum_v chi4(ell_v) (I12-I01)_v
```

is a complete matching-odd line-H4 observer. It is not the same object as landing H4.

### 3. A finite lattice bridge exists

With `A_top=q=r-1`, the exact geometric tilt obeys

```text
partial_h <A_top>|h=0 = Cov(A_top, Re[e^(-4 i phi) J_D4]).
```

Tiny axis and Gaussian enumerations give nonzero real/complex couplings with exact complement sign controls. This closes a finite-lattice symmetry obstruction: the global rank projector can couple to the matching-odd line-birth source. It does not establish survival at large N, Q4 dominance or a Jordan logarithmic partner.

### 4. P205 resolves the ordinary character competition

The frozen P205 equal-area, Smith-changing prism at `N=25,50,125` selects one ordinary H4 amplitude over the frozen H8 and H12 alternatives. The retrospective activation-resolved reuse then scores the complete six-coordinate covariance:

```text
K1 line = H4,
K2 line = H4,
chi-square = 2.585155 / 4 df,  p = 0.629455,
K2 share of signed fitted H4 amplitude = 36.732%.
```

Both activations reinforce at every P205 size. The runner-up `H4/H12` is worse by `Delta chi-square=8.918835`; all other pairs are farther away. The result is decisive inside the frozen small-N ordinary H4/H8/H12 character family, but the component split is retrospective, the raw archives are `branch_only`, and neither a continuum operator nor an asymptotic theorem has been identified.

This changes the information ranking of `rho`/C3 geometries. Another ordinary `A_top` H4/H8/H12 selector would largely repeat P205. A C3 acquisition gains substantially more information when it names a distinct typed non-`A_top` module row, a charged row or another observable whose complex character answers a different question.

### 5. The remaining ordinary module question is Q4/Jordan

Two independent selection steps narrow the regular ordinary carrier:

- vacuum/KdV spin 4 lies on the Alexander-even rank line and has exact zero `A_top` projection;
- the regular unlabeled four-leg `[2]` endpoint has an exact ordinary linear selection zero under its declared generic-Q hypothesis.

The transverse thermal `Q4 epsilon`/Jordan family is therefore the first listed candidate allowed by both selectors. “Allowed” is not “identified”: the missing evidence is a nonzero large-N lattice overlap plus a module fingerprint. The sharp fingerprint is the thermal-Q4 positive-mode Ward vector `4:-6:3`; an ordinary primary gives zero, and Jordan residues add a second gate. Modulus/phase dependence and the log-partner response must accompany any scaling power.

## Rank-1 acquisition contract

The next common-field stream should retain, in the same orientation-paired replicas and delete-one batches:

```text
q and q^2;
I01, I12 and direct 0->2;
ell, iota and Re/Im chi4(ell);
line H4 and landing H4 separately;
N*S_H4, N*D_H4 and the unmarked S birth mass;
q*(N*S_H4), q*(N*D_H4) and the required qJ cross moments;
K1,K2,ell,iota when a permutation/threshold stream is used.
```

This one stream can estimate `gamma_D4=Cov(A_top,J_D4)/M'`, preserve the even `S` sign control, distinguish line from landing mechanisms, and connect the P205 activation split to the Q4/Jordan modulus/Ward fingerprint. Endpoint histograms or old marked-pivotal batches cannot reconstruct the missing same-sample `qJ` covariance.

## Parallel programs

### Ordinary versus charged sectors

The ordinary P205 prism selects H4, while the branch-only N325 deck-charged likelihood weights H8/H12/H4 at approximately `71/21/8`. The charged H8 row is a likelihood leader, not a discovered spin-8 field; H12 and H4 are not both rejected at strict 5%. This is sector tension, not contradictory replications. The next charged observations are a typed complex C3 row or the two primitive Z5 cubic fusion channels, with full complex covariance.

### Low-dimensional state and the targeted context rectangle

Norm-4 rejects the frozen analytic q2 and common-generator forms; Jordan remains borderline, and one post-reveal conjugation-even mode closes the visible residual. In the annulus, the simple odd-shell scalar law fails while a two-state recurrence trained on N325/N425 passes held-out N365. Separately, an unmarked final CRT endpoint is factor-order blind, whereas an intermediate rank filtration can retain activation time.

These results motivate a **targeted pilot**, not another retrospective fit: acquire a common source and common readout basis over at least two Gaussian-cover and two annulus contexts, with shared replicas and held-out prediction. Compare one common generator against a context/morphism-enriched realization. Existing disjoint archives do not form the rectangle, and an endpoint-only row cannot test chronological memory.

### Potts-Q tangent and the P263 variance gate

The exact boundary tangent supplies a function-valued positive control, but the lattice estimator still needs an efficiency gate. P263 constructs an exactly unbiased stopped-transcript covariance estimator and proves cancellation of irrelevant spectator noise. Its frozen 20k pilot reveals only about 20%--23% of edges, yet the comparison to the earlier global-score smoke used different random streams and shows no demonstrated variance advantage.

The next P263 output is therefore paired: on each outer transcript record both the stopped estimator and its matched global score, vary completion count `K` and scale, and report their covariance and variance ratio. This is an estimator-mechanism gate, not a rescore of the boundary shape.

Before any higher Q derivative is called a new state, subtract the topology-forced Pascal jet

```text
D^n H|Q=1 = sum_(k<n) binom(n,k) D^k P|Q=1.
```

Only a typed residual after this subtraction can be assigned to an explicit field derivative, singular projector residue or finite-lattice failure of the continuum identity.

### Cross-microscopic logarithmic control

The triangular log-pair cross-cutoff result is positive but its `kappa_proxy` is gauge dependent. The next statistics are a top field with two two-point-normalized spin insertions and the two-modulus weight-12 H4/E6 double ratio. These close the declared gauge before comparison across triangular, square-bond and square-site systems.

### Relative source and threshold origin

The scalar relative source closes on `1,q,q^2`, its three-sector algebra is semisimple, and `HH^2=0`. Higher powers and unsubtracted Pascal Q-jets do not add new representation content. The theory frontier is a connectivity/defect radical, marked junction or singular projector extension.

The numerical origin of square-site `p_c` remains a separate program. Exact block semantics, correlated-hyperedge closures and symbolic obstructions can run in parallel; success or failure there neither identifies nor refutes the continuum carrier.

## Lifecycle boundaries

The mechanism chain deliberately mixes lifecycle states:

| Layer | Current status | Scientific use |
|---|---|---|
| Digital Alexander theorem, filtration, essential births, index-7 frontier | `main_integrated` | exact finite-state backbone within declared scope |
| `C/W`, `I01/I12`, `S/D`, tiny `Cov(A_top,J_D4)` | `branch_only` | exact coordinates and finite controls, not large-N field identity |
| P205 total prism | `branch_only` | frozen ordinary character selection at small N |
| ten-size K1/K2 root and full-curve reuse | `open_pr` #267 | retrospective component, area and stable-node decomposition with shared-stream covariance |
| P205 K1/K2 prism reuse | `open_pr` #267 | retrospective component-character decomposition with branch-only source data |
| Q4/Jordan overlap and fingerprint | `hypothesis` plus branch-only exact predictions | next mechanism discriminator |

Use `docs/STATUS.md` and `analysis/research_ledger.yaml` for exact commit/path pointers. Citation here never promotes an unmerged result to `main`.

## Update rule

Integrate useful science into the map before forcing issue taxonomy to catch up. Record what is supported, what is in tension, and exactly which parameterization was excluded. Change attention when a sharper observable appears; do not lock, close or rhetorically erase a distinct attempt.
