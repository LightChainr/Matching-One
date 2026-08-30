# P250 multi-separation charged-cubic smoke

This package separates three questions that the original compact three-point
score combined:

1. is the charged two-point denominator resolved at the requested separation;
2. is the local-variance-normalized cubic vector detectably nonzero;
3. conditional on both gates, does the separation-normalized cubic satisfy the
   common-phase closure?

`existing_archive_reveal.*` reclassifies the one-million-replica result at
`be80f25`.  `response_4k.*` is a capped model-development stream at separations
1, 2 and 3; `score_4k.*` applies the support-first decision order.
`exact_mapping_gate.json` verifies every translated root label used by the new
observable.  The scientific interpretation and next selector are in
`notes/p250-separation-normalized-near-zero-reveal.md`.

This is not a production acquisition.  Its sample cap is 5,000 replicas.
