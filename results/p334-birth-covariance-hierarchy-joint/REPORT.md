# Covariance response hierarchy: shared-batch shares and estimator difference

## N325

| Coordinate | Estimate | Shared-batch SE |
|---|---:|---:|
| covariance_hierarchy.all.plus->S.cov_xy.total | 4.457782782e-07 | 3.56657e-08 |
| covariance_hierarchy.all.plus->S.cov_xy.between_prefixes | 4.36457711e-07 | 3.29074e-08 |
| covariance_hierarchy.all.plus->S.cov_xy.within_prefix | 9.320567234e-09 | 6.40755e-09 |
| hierarchy_joint.plus->S.between_prefixes.signed_response_share | 0.9790914728 | 0.0134141 |
| hierarchy_joint.plus->S.within_prefix.signed_response_share | 0.02090852715 | 0.0134141 |
| rankcell_transport.all.plus->S.cov_xy.within_rankcell_prefixes | 3.698558531e-08 | 9.35027e-09 |
| rankcell_transport.all.plus->S.cov_xy.between_rankcells | 3.994721257e-07 | 2.93879e-08 |
| rankcell_transport.all.plus->S.cov_xy.within_rankcell_total | 4.630615255e-08 | 1.167e-08 |
| hierarchy_joint.plus->S.between_rankcells.signed_response_share | 0.896122905 | 0.0217603 |
| connected_clock.plus->S.delta_intrinsic_rank_cov_tau12 | 4.591451791e-07 | 7.51861e-08 |
| hierarchy_joint.plus->S.exactscore_total_minus_matchedmask_intrinsic | -1.336690092e-08 | 6.15837e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.total | -1.309059333e-06 | 9.20886e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.between_prefixes | -1.259345434e-06 | 8.29938e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.within_prefix | -4.971389877e-08 | 2.31287e-08 |
| hierarchy_joint.minus->D.between_prefixes.signed_response_share | 0.9620231891 | 0.0164508 |
| hierarchy_joint.minus->D.within_prefix.signed_response_share | 0.03797681091 | 0.0164508 |
| rankcell_transport.all.minus->D.cov_xy.within_rankcell_prefixes | -1.334367975e-07 | 2.75886e-08 |
| rankcell_transport.all.minus->D.cov_xy.between_rankcells | -1.125908637e-06 | 7.22159e-08 |
| rankcell_transport.all.minus->D.cov_xy.within_rankcell_total | -1.831506963e-07 | 3.94035e-08 |
| hierarchy_joint.minus->D.between_rankcells.signed_response_share | 0.8600898434 | 0.0245803 |
| connected_clock.minus->D.delta_intrinsic_rank_cov_tau12 | -1.194918087e-06 | 3.25239e-07 |
| hierarchy_joint.minus->D.exactscore_total_minus_matchedmask_intrinsic | -1.141412459e-07 | 2.9828e-07 |

## N425

| Coordinate | Estimate | Shared-batch SE |
|---|---:|---:|
| covariance_hierarchy.all.plus->S.cov_xy.total | 4.219105362e-07 | 3.39756e-08 |
| covariance_hierarchy.all.plus->S.cov_xy.between_prefixes | 4.095189781e-07 | 3.04818e-08 |
| covariance_hierarchy.all.plus->S.cov_xy.within_prefix | 1.239155804e-08 | 8.22064e-09 |
| hierarchy_joint.plus->S.between_prefixes.signed_response_share | 0.9706298919 | 0.0183007 |
| hierarchy_joint.plus->S.within_prefix.signed_response_share | 0.02937010805 | 0.0183007 |
| rankcell_transport.all.plus->S.cov_xy.within_rankcell_prefixes | 3.191765512e-08 | 8.33132e-09 |
| rankcell_transport.all.plus->S.cov_xy.between_rankcells | 3.77601323e-07 | 2.69348e-08 |
| rankcell_transport.all.plus->S.cov_xy.within_rankcell_total | 4.430921315e-08 | 1.36856e-08 |
| hierarchy_joint.plus->S.between_rankcells.signed_response_share | 0.8949796003 | 0.0275921 |
| connected_clock.plus->S.delta_intrinsic_rank_cov_tau12 | 5.347939641e-07 | 6.33625e-08 |
| hierarchy_joint.plus->S.exactscore_total_minus_matchedmask_intrinsic | -1.12883428e-07 | 5.79992e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.total | -9.733073551e-07 | 5.82399e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.between_prefixes | -9.343339758e-07 | 5.40117e-08 |
| covariance_hierarchy.all.minus->D.cov_xy.within_prefix | -3.897337926e-08 | 1.88738e-08 |
| hierarchy_joint.minus->D.between_prefixes.signed_response_share | 0.9599577882 | 0.0186969 |
| hierarchy_joint.minus->D.within_prefix.signed_response_share | 0.04004221181 | 0.0186969 |
| rankcell_transport.all.minus->D.cov_xy.within_rankcell_prefixes | -1.064304836e-07 | 2.12545e-08 |
| rankcell_transport.all.minus->D.cov_xy.between_rankcells | -8.279034922e-07 | 5.29749e-08 |
| rankcell_transport.all.minus->D.cov_xy.within_rankcell_total | -1.454038628e-07 | 2.32797e-08 |
| hierarchy_joint.minus->D.between_rankcells.signed_response_share | 0.8506084824 | 0.0217158 |
| connected_clock.minus->D.delta_intrinsic_rank_cov_tau12 | -1.088798783e-06 | 1.61897e-07 |
| hierarchy_joint.minus->D.exactscore_total_minus_matchedmask_intrinsic | 1.154914281e-07 | 1.62115e-07 |

Signed response shares are not baseline variance fractions or probabilities. Prefix law is unchanged; between-prefix transport is the response of relationships among conditional means. Exact-score and matched-mask values reuse the same e32a8593/959a7fa2 block and estimate the same physical target with different estimators/products. No fit, inverse, cell reanalysis, new sampling or independent evidence.
