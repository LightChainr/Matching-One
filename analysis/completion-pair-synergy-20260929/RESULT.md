# Completion-pair synergy: exact L4 successor calculation

Uniform-permutation count clock. Exhaustive occupied subsets and integer prefix DP, not 16! paths.

## square

19932 rank-one configurations; 72 (k,D,nu2) cells.
34 cells have different pair counts; 34 different successor laws.

First weak-history witness: k=6, a=4, D=[0, 1], nu=0.
Identical one-step survival: 1.
Two-step early-minus-late survival: -4/4131.
Mean synergy-edge early-minus-late: 20/459.

Wall time: 0.540s.

## triangular

23502 rank-one configurations; 90 (k,D,nu2) cells.
48 cells have different pair counts; 48 different successor laws.

First weak-history witness: k=5, a=4, D=[0, 1], nu=0.
Identical one-step survival: 1.
Two-step early-minus-late survival: 2/121.
Mean synergy-edge early-minus-late: -10/11.

Wall time: 0.693s.

Complete per-cell integer birth masses, edge sums and successor counts are in result.json.
Existence witnesses are selected lexicographically from an exact finite census, not significance-ranked samples.
