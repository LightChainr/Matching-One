# Two named Q paths already have opposite original-U tangents

Decision: **named_Q_paths_have_opposite_original_U_tangents** on the same N25 iid baseline.

| Root-adjusted derivative at Q=1 | Response | Exact sign |
|---|---:|---|
| tied_Q_path | +0.0630826817070846 | positive |
| local_edge_control | +0.332910708420572 | positive |
| rank_projected_site_RC_Q_path | -0.269828026713487 | negative |

The identity is `V_Sstar/2 = V_(C_B-r/2) + V_B`, with
`V_B=V_Bvac` on the thermal quotient. The first term is the Q tangent of
the rank-projected ordinary site-RC law; the left side is the fixed tied-edge
family. Both use the same rank projection and original q/E observer.
The factor1/2 converts the pre-existing t derivative to log-Q, since Q=exp(2t).
At Q=1 the log-Q and Q derivatives coincide. Each path follows its perturbed
pooled root, retaining its individual geometry normalizers and slope motion.

The local edge control reverses the response sign; it cannot be removed by
a common density reparameterization. This excludes equality of these two
named finite-observer tangents, not Potts universality or a continuum field.
It does not compare the tied law with *unprojected* ordinary site-RC.

Only the saved rational source enclosures were transformed. No profile,
root or scorer was rerun. Both source values were known before this reduction;
this is an exact algebraic consequence, not a prospective independent test.
Full reduced-response bounds, common-root provenance and the fixed matrix
are stored in latest.json. The area factor is displayed numerically only.

The ordinary endpoint's regular Q-activation exclusion and the remaining
trace/confluent interface are explained in
[the mechanism note](../../notes/weak-q-paths-and-regular-selection.md).
Neither this N25 sign comparison nor its source control establishes a
cross-size logarithmic velocity. P154/P334/F4 production decisions stand.

Reproduce: python scripts/analyze_weak_q_path_comparison.py --output-dir NEW_DIRECTORY.
