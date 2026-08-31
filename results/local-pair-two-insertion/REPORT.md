# Single-insertion regularity does not extend to the fixed two-copy local tensor

The exact fixed contraction is

`G(Q)=Q(Q−3)(3Q²−9Q+8)/[8(Q−2)(Q−1)]`.

It has a nonremovable Q1 simple pole with residue **1/2**, although every
one-insertion Bell4 closure is finite. The 17x17 two-hole witness has four
disjoint occupied NN paths, K=52, four components, and rank0. Its physical
reflected gluing equals this Frobenius closure. No occupation ensemble was
enumerated.

After summing the two holes' vacant/occupied states, the conditional
partition, apart from common exterior factors, is

`Q^4+(v_x+v_y+v_x*v_y)Q+epsilon_x epsilon_y G(Q)`.

Both first-insertion terms vanish. The connected mixed log-partition
susceptibility has residue **1/[2(1+v_x)(1+v_y)]**. A common partition
normalizer therefore does not remove this pole. Separate single-site
quadratic counterterms cannot change this mixed derivative.

This excludes an unrenormalized finite-strength Q1 continuation in every
physical exterior for the specified tensor. It does **not** negate the
completed finite linear original-U response, or prove divergence after all
exterior configurations are summed. Rescaling each insertion by
sqrt(Q−1) makes this two-copy limit 1/2 but sends its finite one-insertion
response to zero; that is a different normalization and mechanism.

The 15 exact colour equality patterns, physical coordinates, rational
formulas and provenance are in latest.json; execution and hashes are in
run.json. This is exact algebra plus one graph construction, not Monte
Carlo, a new lattice enumeration, or a scientific test suite.
