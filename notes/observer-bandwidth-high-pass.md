# Six-level finite-noise high-pass control

For `0<rho<1`, define

`H_{rho,d} = product_{m=0}^d [I-rho^(-m)T_rho]`.

On Walsh degree `j`, its exact multiplier is

`h_j = product_{m=0}^d (1-rho^(j-m))`.

The committed `rho=1/2,d=4` oracle expands this polynomial as six noise
levels `T_1,T_(1/2),...,T_(1/32)`. It annihilates every raw monomial of degree
at most four on a five-site cube (992 exact output-point checks), while the
centered degree-five product survives at all 32 points with the predicted
positive multiplier. An independent `L=3` square-torus Euler control has no
Walsh energy above degree four and therefore has exactly zero filtered L2
energy.

The alternating representation has a large reported coefficient L1 norm.
Exact population contraction therefore does not promise a low-variance Monte
Carlo estimator; common-noise coupling and variance per wall time still need a
separate pilot.
