# Fixed two-coupling endpoint profile test

**finite_root_two_coupling_profile_closure_rejected**. The positive rational lower bound rejects common thermal-plus-S profile closure at this finite root, including p-dependent coordinate changes. At least one additional finite response direction is needed beyond the two endpoint tangents.

| Profile | T | C | H |
|---|---:|---:|---:|
| q_first | 5.8824957165133 | 1.1168901349429 | -0.57701112915732 |
| E_first | -0.056485379803595 | 0.13098621519115 | 0.0043952012039549 |
| q_second | 5.8995287019857 | 1.1343383460836 | -0.59003162716463 |
| E_second | 0.045809178985119 | 0.16334206048399 | -0.01974510457728 |

All four predeclared 3x3 minors enter D3; D3 ≈ 0.00043915423800966.
Its outward rational enclosure is [4391542380096601812723550037505088903/10000000000000000000000000000000000000000, 4391542380096601812723550037505093099/10000000000000000000000000000000000000000].
The positive lower-bound witness, if any, is `1/10000`.
The auxiliary sum of squared (T,C) minors is ≈ 3.0611227187792;
endpoint rank two certified: True. Full bounds and the equivalent root-comoving columns
are in `latest.json`. D3 magnitude depends on units; its zero criterion does not.

The four pinned tables pass every K/count check. Parent endpoint coefficients use the
fixed child complement/reversal. Exact pooled-root endpoint signs are negative/positive
and the slope is strictly positive. All arithmetic before serialization uses Fraction
intervals; reported bounds are rounded outward on the 1e-40 rational grid. The primary
decision uses the unrounded rational lower bound, with a short exact positive witness.

This tests one common map for all four profiles, not scalar U alone or separate maps
per profile. An effective source shift can survive an earlier scalar-gain rejection:
U(epsilon,t)=u(t+c*epsilon) gives U*U_t_epsilon-U_epsilon*U_t=c*(u*u''-u'^2),
which need not vanish. Thus that earlier rejection does not settle this tangent test.

The result is conditional on the published finite graph counts and their different
quotient Smith classes. It is not a continuum field count, an independent confirmation,
or a revision of previous P154/P334/F4 stop decisions. No samples, enumeration or cloud
jobs were added. Freeze commit: `76a070d4aa95866f129572940f591b197cda064d`.
