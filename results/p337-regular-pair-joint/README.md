# Canonical Kreg joint source: raw exact N25 cross moments

Both frozen traversals completed once, locally in parallel: axis `(5,0)` **3.37178 s**, tilted `(4,3)` **3.19402 s**, each covering all **33,554,432 configurations**, both exit zero. The executable was compiled once. No benchmark, Monte Carlo, old-source scoring, root search or cloud run occurred. This task has not scored J2.

Each CSV contains the prescribed 13 columns: K/count/q/E plus the three signed source cross moments for total, adjacent and nonadjacent pairs. `b16=sum_y g16(origin,y)`; divide each of the nine source columns by **400 exactly once**. No further factor of two or N is present. Origin-occupied configurations remain in the population with zero source. Adjacent marks share their physical edge-node ID; the original q/E remain unchanged.

Pre-data provenance: contract `4ce4dfe894c9fe96f268c61cf21eb6585dba5418`, scorer `5da4749245450048625a2da43e8f73da1ee9275c`, producer `30891e0468dc874084be2eded58426087ad04c45` (public `8771d6ec`), accepted response proof `93651d61` (public `8212a6b3`) and adjacent-edge proof `4bdf275f` (public `6d1e453e`). `run.json` records their complete SHAs, actual commands, compiler, times and producer/kernel/CSV hashes.

These are new joint contractions on the same complete N25 populations as earlier work, not independent stochastic confirmation. Total/adjacent/nonadjacent are one exact linear decomposition. The parent task alone will perform the locked moving-root/slope readout.
