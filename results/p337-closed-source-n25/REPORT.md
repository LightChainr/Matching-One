# Exact closed source Sstar reaches the N25 original global U

**V_Sstar=0.12616536341417; V_Bvac=0.33291070842057.**
The exact rational enclosure of V_Sstar/A excludes zero.
Decision: `closed_Sstar_common_thermal_alias_rejected`.

## Fixed mechanism, not a source scan

Sstar=Ctot+F4+Bvac is selected by T(Ctot)=Ctot+F4, T(F4)=Bvac, T(Bvac)=0.
It is fixed under this map, and bare Ctot reaches it after two steps.
The only companion source is Bvac, the next endpoint of F4. Both are bulk
sources exp(t*S); no amplitude or density normalization was fitted.

One full pass per geometry enumerated all 2^25 configurations of (5,0) and
(4,3). The inherited geometry, component updates and traversal were unchanged.
At each leaf Bvac=50−4K+occupied_NN_edges and Sstar=Ctot+F4+Bvac.
Only their sufficient statistics plus q/E were saved; old C/F4 source
responses were not recomputed. There were no random samples or cloud jobs.

## Root-complete finite response

| Quantity | Numerical evaluation of exact coefficients |
|---|---:|
| Fresh exact pooled root | 0.5926655393282267 |
| Native U25 | 0.88046615696337 |
| V25 for Sstar | 0.12616536341417 |
| V25 for Bvac | 0.33291070842057 |
| N50 F4 endpoint predicted from Bvac | 1.0268369968409 |

Let Q=mean(q), Y=P4(E), D=Q_p, A=25^(13/8)/2 and j_O=Cov(O,S).
Every source uses the exact original-U derivative

`V_S/A = jY_p/D − Y_pp*jQ/D² − Y_p*jQ_p/D² + Y_p*Q_pp*jQ/D³`.

Root motion, slope motion and per-geometry covariance centering are all included.
The exact integer counts contain their configuration multiplicities; no extra
binomial factor is applied. The root uses 128 rational bisections in [11/20,13/20].
`latest.json` preserves outward rational enclosures for the root, V/A, and each
of its four terms. These are computational bounds, not confidence intervals.

## What changes

A nonzero Sstar response excludes common-thermal invisibility of this uniquely
closed source **for this finite global observable**. Closure of the source
under decimation does not force it to disappear after root normalization.
The companion predicts the next F4 endpoint through
`V50_endpoint_F4=2^(13/8)*V25_Bvac`; this is a transported derivative, not
an independently simulated parent endpoint.

The axis Z5xZ5 and tilted Z25 quotients have different Smith classes. Nothing
here identifies a continuum field or an asymptotic exponent. The calculation
does not change the independent larger-N F4 experiment or revive a lag-one
source. If the reported enclosure contains zero, the finite alias is unresolved.

## Provenance and reproduction

Contract freeze: `5598812612d176e30c0e9ee50d2fd78f382db632`.
Enumerator parent: `2bfe9b90:scripts/exact_decimation_plaquette_u.cpp`.
Corrected scorer parent: `c76b038b:scripts/analyze_decimation_plaquette_u.py`.
Full commit/blob IDs, source hashes and enumeration receipts are saved alongside
the exact profiles. Run once with
`python scripts/p337_closed_source_score.py --output-dir NEW_DIRECTORY`.
For future inspection only, `--counts-dir results/p337-closed-source-n25`
consumes these saved new-source profiles without enumerating again.
