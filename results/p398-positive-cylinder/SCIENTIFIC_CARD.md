# P398: positive-cylinder rooted/landing response

- **Lifecycle:** branch-only / exact finite positive-measure calculation / completed / no new random samples.
- **Observer/source/geometry:** existing AP and landing cyclic-charge-one functions; unique stationary past-frontier distribution; Q=1 square-bond width-4 cylinder, p=1/2, after-horizontal-layer slice.
- **Result:** the entire two-dimensional complex charge-one sector propagates with distinct eigenvalues `(3±sqrt(5))/64`. `det C(d)=(73216/1940449)1024^(-d)>0` for every finite d>=1. No exact common ray, no Jordan block in this sector.
- **New phenomenon:** ordinary fast/slow mode separation makes the response approximately rank one at large distance; the fast/slow decay ratio is only 0.145898. The ordered cross asymmetry at d=1 is exactly `(1-i)/1393`, consistent with a past-connected frontier readout rather than a time-reversal-symmetric local observable.
- **Dependency:** one deterministic 14-state transfer construction, inherited from afc619c; all distances are coordinates of the same exact model, not independent evidence. No P250/P337 production covariance is reused or counted.
- **Not proved:** site-Matching coupling, continuum/Jordan identity, universality, field multiplicity, or intrinsic temporal memory. Equal-time positive definiteness is only a baseline.
- **Next observation that could change the interpretation:** a width or microscopic-readout comparison of these two explicit eigenmixtures, not another fixed-radical jet.
- **Artifacts:** `latest.json`, `analysis/p398_positive_cylinder_protocol.json`, `scripts/p398_positive_cylinder.py`, `notes/p398-positive-cylinder-propagation.md`.
- **Checks:** three focused arithmetic tests, including the direct 14×256 physical bond-mask transfer sum; no broad suite.
