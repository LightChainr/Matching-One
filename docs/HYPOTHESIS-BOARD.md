# Hypothesis and Decision Board

**Frontier date:** 2026-08-29

This board turns the repository into a mechanism-discovery surface. It does not grant permission, lock tasks, close alternatives or promote branch results to `main`. Priority means expected information gained per unit of shared acquisition.

## What changed the decision surface

Four newly recovered results make the previous “collect a larger global `gamma_D4`” plan obsolete:

1. The marked-birth pilot already acquired a common field at q2 N65/N130 and P50 N145. Its mean line source has a noisy but striking q2 transfer `J_D4(130)/J_D4(65)=-0.33074-0.00823i`, while the connected `A_top` response scales approximately with the thermal birth mass.
2. Gate algebra gives exact contact identities for `Cov(q,J_D)` and `Cov(q,J_S)`. The all-order extension shows that every connected insertion using only `f(q)`, `q in {-1,0,1}`, and a root-independent mark reduces to contact moments. No nonlinear lookup table in the global rank charge creates an independent field matrix element.
3. The nine-geometry P275 production rejects ordinary Q4, Q4-Jordan, pure H4, affine-log H4 and zero as frozen models for the global-line `Gamma`; its order-one response is a projective homology-line polarization, not an `N^-13/8` local coupling.
4. Main now contains the exact rank-four operator-mixing audit. Matching parity fixes the four channel assignments and powers but does not supply the two Taylor ratios needed for a two-amplitude closure.

The global-line response remains useful as a projective/contact control. It is not deleted or called useless; independent observers, separated local marks and mean-source transfer now receive more attention.

## Three-object synthesis

```text
rank-one lifetime fiber
    K1/K2, C/W, I01/I12, ell, landing marks
                 |
                 | acted on by geometry/context
                 v
low-dimensional context action
    Gaussian A, annulus U, mixed orders AU and UA
                 |
                 | scalar q quotient forgets
                 v
connectivity / defect radical
    partial partitions, junctions, singular projectors, charged seams
```

This is a conjectural organizing skeleton, not a theorem. It is valuable because each piece has a small, decisive observation that can break it.

## Active hypotheses

| ID | Bold mechanism | Current support and boundary | Smallest decisive observation | Sharp falsifier | Attention |
|---|---|---|---|---|---:|
| H-LIFETIME-POLARIZATION | The local negative `K2` lobes are a shape/polarization effect inside a globally positive second-activation response, not reversal of the second birth mechanism | PR #267 has positive integrated `K2` area at 10/10 sizes and stable local nodes at N265/325/425; line-resolved attribution is unmeasured | Existing-data translation/shape decomposition, then condition N325/N425 second births on `ell`, `W` and landing H4 around the node | The node remains after conditioning while both line and landing responses have no corresponding phase/shape change | 1 |
| H-LOCAL-H4-CARRIER | A high-dimension H4 signal, if present in the birth stream, lives in the mean `J_D4` or in variation within a fixed-q sector, not in connected `f(q)` contact response | Mean q2 ratio is a preregistered lead; exact contact/no-go algebra and P275 rule out the old global-line field-identity interpretation | Frozen third-point mean-`J_D4` q2 chain, or a two-radius landing-H4 UV annihilator / independent global observer with a two-root field Gram | Mean transfer misses the frozen phase/magnitude law and separated-local or fixed-q observers close to zero | 2 |
| H-TWO-CONTEXT-ACTION | Gaussian and annulus data use a small common state, but geometry may act through non-aligned or noncommuting context maps | Gaussian low-rank family survives; annulus odd sector needs rank two; separate spectra cannot align latent bases | Shared 2x2 source/readout rectangle measuring both `c^T A U b` and `c^T U A b`, plus a commensurate clock cell | A shared commuting generator predicts held-out mixed cells, or the rectangle requires rank greater than the proposed state | 3 |
| H-ODD-RADIAL-PLANE | The annulus matching-odd sector is a rotating rank-two plane while the plus control lies on a simpler boundary | Branch-only sector score rejects bounded-window R1 for `A_minus` (`p=.0192`), all named R2 classes close, and separated J2 gains weak support; `A_plus` does not require R2 | One source/readout chosen to separate shared-J2, parity-specific-J2 and gap-one classes, not more replicas of the same shell row | A held-out parity-sensitive cell is predicted by one common scalar/rank-one generator | 3 |
| H-SEAM-TYPED-SPIN4 | A spatially identifiable angular row and a colour seam can jointly separate ordinary singlet from charged spin-4 candidates before an exponent fit | Exact Q=4 preflight gives colour ratios 1/0; open PR #291 proves one C4 orbit aliases scalar/spin 4 | Add calibrated phase diversity or typed internal data, report the necessary angular-rank gate separately, then emit constituent numerators and six raw seam `Z/N` values | Angular rank remains one after calibration, or the independent colour ratios fail 1/0 | 2 |
| H-CHARGED-QUADRATIC | The charged H8 preference may be a quadratic `E4^2` response rather than a new spin-8 primary | Charged norm-5 ranks H8/H12/H4 about 71/21/8, while ordinary P205 selects H4; this is a sector tension, not a contradiction | Complex C3 parent/children including zero order, followed by primitive Z5 `113/122` fusion in the same transported basis | C3 DFT is not r=2, the parent zero is not quadratic, or the two cubic rows require incompatible phases | 4 |
| H-CONNECTIVITY-RADICAL | Q-tangent singularities and threshold nonmonotonicity may both live in connectivity information discarded by the semisimple `1,q,q^2` quotient | Scalar relative source has rank-three closure and `HH^2=0`; the common-radical interpretation remains conjectural | Build a tiny partition/junction algebra, quotient the scalar sector and project both defect vectors into its radical/Jordan blocks | The extension remains semisimple or the two defects lie in unrelated simple components | 5 |
| H-CROSS-MICRO-JORDAN | Triangular log-pair and square spin-4 response may be different microscopic realizations of one descendant extension | Triangular cross-cutoff control is positive; raw `kappa_proxy` and raw amplitudes remain gauge-dependent; square global `Gamma` is not a valid field row | A normalization-free top-field statistic and two-modulus H4/E6 polynomial null using a separated-local or typed square observer | Both systems pass internally but disagree after all allowed shear and metric gauges are cancelled | 6 |
| H-IMPROVED-ACTION | An exactly critical oblique family can tune ordinary H4 through zero; a quadratic charged response would vanish to second order there | Main proves only the microscopic stencil feasibility gate; no exact-critical tunable family yet exists | Exact construction first, then three common-frame coupling values bracketing the zero | The exact microscopic proxy crosses zero without a large-scale ordinary sign change | opportunity |
| H-FACTOR-FOUR | The factor four between the main Q-velocity gap and the boundary high-branch log coefficient is a Ward-descendant collision rule | `sqrt(3)/(4pi)` versus `sqrt(3)/pi` is an exact numerical coincidence across separate derivations | Symbolically apply the full Q4 Ward operator to the colliding generic-Q pair, retaining projector and explicit-field terms | A different residue, colour tensor or Q parametrization removes the factor four | exact sprint |

## Minimum decision portfolio

The next work should share acquisitions instead of opening one stream per hypothesis.

| Work package | Reuse first | Smallest new acquisition | Decisions covered | Cost |
|---|---|---|---|---:|
| Activation shape and local carrier | Existing ten K1/K2 archives; existing marked-birth pilot | Only after reuse: mean-JD4 third q2 point or two-radius landing marks | H-LIFETIME-POLARIZATION, H-LOCAL-H4-CARRIER | none -> moderate |
| Typed exact spin-4 | Existing Q=4 seam preflight and PR #291 alias gate | Calibrated multi-phase or typed-edge kernel; spatial angular rank and colour-seam character scored separately | H-SEAM-TYPED-SPIN4; supplies a positive Q-tangent control | low |
| Mixed context rectangle | Existing Gaussian summaries, annulus R2 fits and P200 static residuals | Both mixed orders in a shared source/readout basis with intermediate state | H-TWO-CONTEXT-ACTION, H-ODD-RADIAL-PLANE | moderate |
| Charged complex geometry | Existing N325 charged engine | One parent plus three C3 children; reuse for Z5 cubic rows | H-CHARGED-QUADRATIC | moderate |
| Exact theory sprint | Existing Q velocities, Ward row, scalar closure and threshold semantics | Symbolic residue, tiny connectivity radical and exact-critical oblique construction | H-CONNECTIVITY-RADICAL, H-FACTOR-FOUR, H-IMPROVED-ACTION | low |
| Gauge-free cross-microscopic control | Existing triangular Phase B sufficient statistics | Normalized top statistic plus typed square row at two moduli | H-CROSS-MICRO-JORDAN | moderate |

## Reinterpretation rules

- Exact contact closure reclassifies `gamma_D4`; it does not negate the observed finite response or the utility of the marked stream.
- P275 rejects the declared global-line field models; it does not reject ordinary global H4, local Q4, charged four-leg or separated landing observers.
- Annulus R1 failure is sector-specific. It does not imply that every plus/minus row needs the same rank-two generator.
- P200 identifies a static mixed-factor interaction. Without `h0,h2,h5,h25` or both path orders it does not identify chronology or morphism memory.
- Main's rank-four mixing audit rejects an assumed amplitude closure, not the exact parity zero pattern or the four powers.
- A lower-ranked route stays open and can immediately rise when it offers a sharper observer, exact map or shared acquisition.

## Team handoff

Every new result should state `state/source/observer/geometry/acquisition`, lifecycle status, immutable source commit, dependency group, exact model pair separated, and the next observation that would change the conclusion. “More samples” is not a next target unless the current uncertainty, rather than model aliasing or missing semantics, is the decision bottleneck.
