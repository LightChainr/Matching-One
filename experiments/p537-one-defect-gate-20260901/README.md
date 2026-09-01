# P537 physical one-defect diagonal-edge gate

This artifact implements the existential stop rule in Issue #537.  It uses
the complete pooled-root coefficients from both N25 geometries, but searches
for one literal physical edge in the axis geometry.  Finding one is enough to
reject a blanket factorization of the full graph into an independent
landing/rank defect and an independent source/Bell defect.

The producer fixes `x=0`, `z=East(x)`, enumerates the other 23 Bernoulli sites
in increasing integer-mask order, and scans source mark `y` in vertex order.
It stops at the first edge for which all three observables change:

```text
rank_index=q+1, source Bell(x4+y4), and g16(source Bell).
```

Requiring `Delta g16 != 0` is stricter than the issue's topological Bell
condition.  It prevents a nonzero `-beta*b` allocation on a kernel-zero Bell
transition from deciding the gate by itself.

The witness records a common first-occurrence canonical component map on
`x4+y4+z4` ports (`outer_C`), the global NN-black partition of the four
thermal ports (`outer_B`), the off-z matching-white partition (`outer_W`),
rank and Bell before/after, off-port contact, the physical transition ID and
the exact integer sufficient statistics.  The scorer independently derives
both Bell keys by restricting the joint map.

Run the complete verifier from the repository root:

```bash
python3 experiments/p537-one-defect-gate-20260901/verify.py
```

The verifier compiles the C++17 producer in a temporary directory, reproduces
`witness.json` byte for byte, independently runs the pooled-root scorer and
compares `result.json` byte for byte, then runs the Python topology oracle.

The compact `global-diag1.csv` contains exactly the 48 `__GLOBAL__/diag1`
rows extracted from PR #544 commit `e8e9c7cf`'s complete N25 Schur aggregate.
It is sufficient to freeze `mu_a` and `beta_diag1`; no landing-cell root or
counterterm is re-estimated.

This witness has `arm_mask=3`, is adjacent to the source, and has explicit
global source/landing contact.  It decides the full-graph stop rule through a
contact/collision channel.  It does not claim that the ordinary separated
four-arm sector contains a diagonal edge, nor does it give an asymptotic
lower bound.
