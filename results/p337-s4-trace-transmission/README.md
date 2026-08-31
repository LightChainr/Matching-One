# Fixed Q4 S4 seam trace: complete missing per-component counts

These two raw histograms complete the producer side of frozen contract
`55fdba789a576d8d4c507372b7834f92cf506c80`. They are **not a response score**.
The code was committed as `0b2ad6ff9675185789256b0a758f8eb5bc7d02ed`
before the single full enumeration of each fixed N25 geometry.

Columns are `k,g,q,bad2,n_bad3,count`. The old k/g/q quantities retain
the white-Alexander/occupied-NN definition and source convention.
`bad2` is an any-component Boolean. `n_bad3` counts the distinct occupied
components with nontrivial P1-cycle gain modulo3. They are not the global
ambient rank and are not computed by reducing Cartesian displacement.

The black rollback DSU alone adds Z/6 vertex potentials and two root cycle
flags. Union joins OR the flags and update the two global component counts;
rollback restores both flags and counters. The geometry representatives,
white graph, binary subset traversal, edge counts and face counts are those
of `scripts/p337_closed_source_finite_exact.cpp`. No full homology module
or extra source has been introduced.

The prescribed insertion can now be recovered exactly as

`beta = 1/6 + (1/2)*(1-bad2) - (2/3)*4^(-n_bad3)`.

The root coordinator will perform the single original-U score at the
already saved m2 root, retaining separately normalized geometries and
root/slope motion. This producer did not examine that new response,
rescore the old four-point baseline, or search for a new root.

Axis and tilted each cover33,554,432 configurations; runtime was1.835s
and1.61768s, respectively, as simultaneous local single-thread processes.
Full source/output hashes and commands are in [run.json](run.json).
No random sampling, cloud task or test suite was run.
