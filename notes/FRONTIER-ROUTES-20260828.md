# Matching One research frontier — post-N185/N265 map

Status: exploratory roadmap, 2026-08-28. No existing preregistration order is changed by this note.

## Where the project actually stands

The signal-discovery phase is over. The main matching-odd result now has several independent supports: same-N orientation data, fixed-coordinate Gaussian multiplier tests, local residual-to-root closure, and the genuinely new N=185/265 prospective target block.

The N=185/265 matching-odd `DeltaM` result remains compatible with the frozen `x=21/4` H4-like law and strongly disfavors zero and the larger frozen `x=17/4` adversary.

The initially reported matching-even sign reversal in #108 must **not** be treated as a physical falsification: the frozen even amplitude was fitted from P31 `either/even`, while the threshold-rank target is `cross/even`. P31 itself gives the exact opposite sign between those channels. #132 / PR #134 applies only the exact cross/either channel map and gives a corrected no-refit score near `chi2=0.57/2`.

P48 `P4[S']` is different: its pure `N^-5/4` law genuinely fails on the new N=185/265 geometries, while both predeclared q=2 and Jordan-log corrections remain viable. This is now a real mechanism-discrimination problem.

The central project risk is therefore **non-identifiability of subleading mechanisms**, not lack of numerical precision.

## Immediate order

1. Merge/audit the cross/either protocol repair (#132 / PR #134). Audit channel semantics in every scorer before interpreting even-sector signs.
2. Run norm-5 Gaussian spectroscopy (#57). This is now the highest-value new-compute experiment: it discriminates H4 from H12 and also separates live radial competitors through parameter-free multiplier ratios.
3. Run the prospective 145->290 full-curve slope/root test using the already-frozen finite-size correction model on canonical main.
4. In parallel, extract more information from existing curves rather than simulating more of the same sizes.
5. Let an explicit information optimizer (#102) choose the next expensive geometry after #57 rather than selecting a larger N manually.

## P0 analysis — mostly zero additional simulation

### Prequential evidence ledger — #95

Build a chronology-locked predictive evidence ledger. Score only endpoints that a model predicted before reveal; use joint target+source covariance and predictive log score, not local chi-square alone. Do not double-count root/slope/DeltaM derived from one histogram as independent evidence.

### Pivotal-measure / Russo bridge — #100

For increasing wrapping events Russo's formula makes the slope a total pivotal mass. For the matching function,

`M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p)`.

This gives a probability-theory explanation of the observed `L^(3/4)` denominator and connects directly to rigorous near-critical pivotal-measure theory. The high-value extension is an orientation-resolved four-arm/pivotal H4 moment.

### Intrinsic quantile-center spectroscopy — #101

Use the already-solved intrinsic levels `Mbar(p_-^u)=-u`, `Mbar(p_+^u)=+u` themselves. For

`t(delta)=a1 delta + a2 delta^2 + ...`,

`w_u=(p_+-p_-)/2 ~ L^-3/4`, while differences of midpoint coordinates

`c_u=(p_++p_-)/2`

at two fixed levels isolate the even nonlinear thermal coordinate with leading `L^-3/2=N^-3/4`. This can separate bare-coordinate nonlinearity from genuine irrelevant/Jordan corrections.

### Multi-u thermal response — #119

Use the frozen `u={0,.025,.05}` vector as one correlated functional observable. q=2, Jordan and coordinate-nonlinearity mechanisms need not have the same near-critical u-shape even when their radial laws are confounded.

### Joint operator-mixing matrix — #125

Stop fitting four derivative/parity channels independently. Model `[P4[S],P4[D],P4[S'],P4[D']]` as a correlated vector produced by the smallest justified RG basis. One physical correction should move several channels coherently. Use new N/multipliers for validation; do not expand the basis after every residual.

### Information-optimal Gaussian design — #102

Enumerate primitive Gaussian tori/multipliers and rank candidate runs by expected model KL/log-score separation per CPU-second using exact harmonic arithmetic and measured variance/throughput. This is the natural successor to the synthetic design red-team.

## P1 theory and exact controls

### FK/Potts torus-sector derivation — #114

Try to derive the matching scaling observable from `Q` derivatives of the random-cluster/Potts torus partition function and homology sectors at `Q->1`. A successful derivation would identify the continuum channel from the observable itself rather than from exponent matching.

### Four-arm anisotropy — #121

Pivotality is a four-arm event with thermal dimension `5/4`. Derive or measure an orientation-resolved four-arm H4 correction. If its correction really adds four units of conformal level, `5/4+4=21/4` becomes a geometric arm-event statement.

### Torus-modulus spectroscopy — #103

Current Gaussian work varies microscopic orientation at fixed `tau=i`. Vary continuum modulus with the general period-matrix engine. Shape dependence is an orthogonal fingerprint of a conformal family; seek parameter-free amplitude ratios across moduli rather than another free exponent.

### Exactly-critical tunable anisotropy — #106

Use canonical critical isoradial bond percolation/star-triangle families as a controlled anisotropy laboratory. Search for an improved discretization where the leading H4 amplitude crosses zero, exposing the subleading spectrum.

### Configuration-level Euler/Betti identity — #111

Lift the matching identity before expectation to Euler-Poincare/Betti/homology statistics on each torus configuration. This may yield exact controls, explain configuration-identical wrapping differences, and connect matching to persistent topology.

### Universal amplitude ratios — #118

Derive truly dimensionless combinations that cancel thermal metric and lattice couplings before comparing square-site, self-matching, square-bond and isoradial controls.

### Exact self-matching Beta-family falsification — #115

N=10 has the exact threshold CDF `Beta(3,3)`. The next `(5,1), N=26` quotient is only `2^26` configurations. Compute the geometry-only shortest wrapping occupation `s`, freeze `2 I_p(s,s)-1`, then enumerate. Do not choose the beta parameter after seeing the target polynomial.

### Finite matching Galois complexity — #104

Axis L=2..5 finite matching polynomials are already irreducible of degree `L^2`. Use discriminants and modular factor-cycle certificates to determine whether their Galois groups are full/large, further excluding a simple finite-cell radicals mechanism.

## P2 high-risk routes

- discrete-holomorphic/parafermionic defect — #109;
- correlated-hyperedge self-dual embedding beyond the independent-bond no-go — #123;
- symmetry-resolved transfer-matrix eigenoperator spectroscopy — #120;
- theorem-grade finite-size bound feasibility — #112;
- local rather than global matching-complex-zero scaling — #113;
- full standardized threshold-distribution collapse — #122.

Each is gated: exact/theory/control evidence first, no large simulation until it produces a distinct falsifiable target.

## Compute methodology we were underusing

For any future billion-replica campaign, consider a **predeclared optional-stopping-safe** likelihood/e-value design (#126). H4-vs-H12 may become decisive far earlier than the sample count needed to show the tiny H4 child amplitude is nonzero at 3 sigma.

## What we should stop doing

- broad PSLQ searches on the decimal threshold;
- treating every residual as permission to add another free power;
- increasing transfer-matrix width solely for extra digits;
- N=1105 production before norm-5/new-geometry discriminators justify it;
- GPU use without an information target and CPU/GPU exact agreement;
- comparing predictions/targets across wrapping channels without an explicit channel map.

## Decision tree

1. If norm-5 supports H4 and the live radial `x=21/4` law, concentrate theory on pivotal/FK/four-arm bridges (#100/#114/#121) and independent controls (#103/#106/#111).
2. If H4 wins but `13/8` loses, retain the angular mechanism and reopen only the radial/operator identity.
3. If H12/another odd harmonic wins, preserve the Gaussian-semigroup framework but abandon the thermal level-4 identification.
4. If q=2 versus Jordan remains confounded, use intrinsic quantile/multi-u/modulus information rather than adding more same-N replicas.
5. If several mechanisms remain viable, let #95 and #102 choose the next discriminating evidence block rather than expanding the model set manually.
