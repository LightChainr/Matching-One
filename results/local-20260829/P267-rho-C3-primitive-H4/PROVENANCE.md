# P267 rho-child primitive-H4 pilot provenance

- Prereveal design commit: `8eddecf4aa053132a434ab6d4454ad493cf1ee72`.
- Authorized manifest commit: `afdfb6b0e30abbf91fa8021ec608db7e7e334832`.
- Execution: local macOS, Python 3.13, eight worker processes.
- Acquisition: 100,000 common-field replicas, 100 batches, seed `267156112`.
- Counter field: the same deterministic 224-bit bond vector drives all three N112 children for every replica.
- Wall time: 6.40 seconds; aggregate worker CPU time: 47.88 seconds.
- Exact preflight: 65,536 masks on each tiny child, zero homology invariant failures, digest `54f6d2cc80670db23d86b73d3d4d1f934bde1d14ed80a0c813b02a311b268f74`.
- Numerical baseline: primitive-line cutoff 9; maximum change from cutoff 7 is below `1.05e-27`.

The batch CSV is the full acquisition record. The run JSON contains the raw
complex means, exact geometry metadata and common-field 6x6 covariance. The
score JSON retains that covariance and its transformed C3 covariance.
