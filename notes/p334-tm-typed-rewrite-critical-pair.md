# Finite square types and the last TM critical pair

The curvature-corrected Rayleigh polynomial has a finite local alphabet. List
the vertices of an ordered missing-site square as `(00,10,01,11)` and record
their ambient ranks. Above a rank-one fixed-line base, monotonicity permits
exactly:

- `D=1222` (both singles coexit);
- `M_left=1212`;
- `M_right=1122`;
- `Y=1112` (the double insertion is the first exit);
- `F=1111` (the entire square is flat in the fixed-line sector).

Alexander complement reverses the square and sends each rank `r` to `2-r`,
giving respectively `0001`, `0101`, `0011`, `0111`, and `1111`. These are the
complete exit/birth local types; there is no omitted sixth case.

Across the existing atlas their ordered counts are:

- `D`: 361,632;
- `M_left`: 517,168;
- `M_right`: 517,168;
- `Y`: 310,920;
- `F`: 1,708,048.

## The finite aggregate rewrite system

The unique negative product in

`M^2+4Y(T-D)-4DF`

is `D x F`. Give its four orientation replicas lexicographic labels. There
are two positive token reservoirs:

1. `R_M`: ordered `M x M` tokens;
2. `R_Y`: four replicas of `Y x nonD` tokens.

Apply `R_M` first, then `R_Y` to any residual. The number of unmatched hard
tokens strictly decreases, so this count-level rewrite terminates. It closes
all 984 rows:

- 968 use `R_M` only;
- 16 use both rules;
- zero hard tokens remain.

The 16 synergy-rescue rows have only four exact `(T,D,M,Y,F)` signatures,
four quotient realizations each. At most `133/2880` of the available synergy
pool is required after the mixed pool is exhausted.

## Extreme-ray classification

Each of the nine bounded Pareto rays now has an explicit minimal witness and
a deterministic rule. The three mixed-deficient rays use `R_M` then `R_Y`,
requiring respectively `133/2880`, `539/14064`, and `9/277` of their synergy
pool. The other six rays close by `R_M`. Four rays form the exact convex lower
hull: one synergy-rescue endpoint followed by three mixed-only rays, including
the synergy-free endpoint.

Thus every bounded cone regime is represented by the same two-rule grammar;
there is no third hidden mechanism.

## The one unclosed general critical pair

The rewrite above is machine-verifiable but aggregate: lexicographic token
labels do not construct images from the underlying site configurations. The
single unresolved topology rule is still the hard pair `D x F`.

Given a coexit square and a flat square, one must cross-switch their two
ordered missing-site pairs, or pass through the Alexander-dual birth square,
to create either two mixed squares or a synergy square plus a non-coexit
square, without collisions between images.

The exact unresolved capacity is

`K=max(0,4DF-M^2)-4Y(T-D)`.

All bounded rows have `K<=0`. Nevertheless the known `N=6` opposite-pair
table has `D,F>0` and `M=Y=0` at that displacement and even after grouping by
torsion order. Therefore the missing realization must exchange mass across
relative-displacement classes. A purely local square rewrite cannot close the
critical pair.

Conditional theorem: a globally injective topology realization of `R_M` and
`R_Y` for every `D x F` token proves aggregate TM. Translation regularity then
produces the explicit one-mark Hall injection and every proper Hall cut.
