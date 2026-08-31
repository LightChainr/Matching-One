# P398: T4's plus-tail repair is current-enabled, not generic reversible geometry

The main T4 improvement **does not survive current deletion**. On plus,
adding T4 to the fixed triplet-incidence span improves the original t=4
error from +8.0624% to +2.3768%; in the reversible process the corresponding
error moves from -5.7682% to -5.7891%, a slight worsening. In contrast,
the earlier triplet-incidence extension makes large improvements in the
reversible model. These are distinct roles for already named geometries.

This does not contradict the preceding result that the *full-process*
fast-to-slow inversion survives current deletion. The existence of the
inversion is not the same question as which omitted geometry repairs a
particular reduced approximation.

## Exactly one geometry x current comparison

The 5+5, 7+7, and 8+8 spans are unchanged from `2385062d`/`30eef34a`:

- 5+5: A,T2,T3,S11,B2 and their Kreweras partners;
- 7+7: additionally the already defined triplet boundary-incidence classes;
- 8+8: additionally T4 and its Kreweras partner.

We reused the archived exact independent columns, without reselecting rank
or features. Let P be each stationary-L2 orthonormal basis and M=-G. In
stationary-whitened coordinates,

\[
B=P^*MP,\qquad B_{\rm rev}=(B+B^*)/2=P^*(-S)P.
\]

Thus the projector is identical with or without current. Original G scores
are reused from their archived results; only the small Hermitian projected
matrices are newly diagonalized. Full G and full S references are reused
from `520a9d21`, with the same sources and distance grid. There is no new
full-sector spectrum solve, fitted parameter, width, or Monte Carlo.

These are observable Galerkin projections, **not Markov state chains**.
Even though full S is a reversible generator, a matrix in a complex named
observable basis need not be a transition-rate matrix.

## Crossing and slow masses

| Dynamics / dimension per ray | Crossing | Minus slow mass | Plus slow mass |
|---|---:|---:|---:|
| Full G | .2656573200 | 2.8196586326 | 1.9557501384 |
| G, 5 | .2650793593 | 2.8780045603 | 1.9657576829 |
| G, 7 | .2654226141 | 2.8428902867 | 1.9313697560 |
| G, 8 | .2656408870 | 2.8404109327 | 1.9479283947 |
| Full S | .2722634760 | 2.5407959794 | 1.8363180504 |
| S, 5 | .2698038155 | 2.6874758524 | 1.8865264703 |
| S, 7 | .2754367737 | 2.6040503689 | 1.8535863593 |
| S, 8 | .2752921923 | 2.6024649989 | 1.8534155439 |

The plus G pole in the 7-dimensional projection is artificially too slow.
T4 moves it **up** toward the full pole, reducing the relative mass error
from -1.24660% to -.39994%. Under S, T4 only moves the pole slightly **down**,
from +.94038% to +.93107% error. The latter is an ordinary Hermitian Ritz
improvement in a nested span; it has almost none of the original repair's
magnitude. The current-enabled correction cannot be explained as just a
better symmetric slow-mode approximation.

The crossing improvement partly persists, but the S crossing remains
1.11242% high, whereas the G crossing is only .00619% low. A ratio crossing
is not a proxy for accurate individual tails.

## Original fixed t=2 and t=4 errors

All numbers are percentage errors against **that process's own full
propagator**, not against a shared G reference.

| Dynamics / dimension | Minus t=2 | Minus t=4 | Plus t=2 | Plus t=4 |
|---|---:|---:|---:|---:|
| G, 5 | -2.84954 | -11.83004 | -.18179 | -2.09366 |
| G, 7 | -.66136 | -4.01145 | +2.92025 | +8.06241 |
| G, 8 | -.47399 | -3.39155 | +.81666 | +2.37679 |
| S, 5 | -13.96362 | -35.03716 | -5.32225 | -14.22295 |
| S, 7 | -5.38792 | -15.83919 | -2.55507 | -5.76819 |
| S, 8 | -5.26058 | -15.47353 | -2.58069 | -5.78914 |

The **5->7 triplet-incidence step** repairs reversible t=4 errors by 19.198
percentage points on minus and 8.455 on plus. The **7->8 T4 step** adds only
.366 points on reversible minus and worsens reversible plus by .02095
points. T4 is therefore not a general replacement for the missing symmetric
geometry, especially the sizable reversible minus tail.

The tiny reversible plus-tail worsening is consistent with its reduced
leading residue, .4622953120 -> .4618736107, offsetting the tiny mass
improvement. It is not a contradiction to the variational eigenvalue bound:
the bound does not fix source overlaps or enforce monotone error of a
particular correlation at every distance.

## The requested 2x2 difference

For each metric use e_(n,dyn)=prediction/full_dyn-1. Define the signed
interaction as

\[
\Delta=(e_{8,G}-e_{7,G})-(e_{8,S}-e_{7,S}).
\]

The score retains this signed quantity and all four cells. For a more
readable account of repair, the table gives absolute error reduction
|e7|-|e8| (positive is improvement), in **percentage points**:

| Metric | T4 repair in G | T4 repair in S | Repair difference G-S |
|---|---:|---:|---:|
| Crossing | .08216 | .05310 | .02906 |
| Minus mass | .08793 | .06240 | .02553 |
| Minus t=2 | .18737 | .12734 | .06003 |
| Minus t=4 | .61990 | .36566 | .25424 |
| Plus mass | .84666 | .00930 | .83736 |
| Plus t=2 | 2.10359 | -.02562 | 2.12921 |
| Plus t=4 | 5.68562 | -.02095 | **5.70657** |

The signed plus t=4 interaction is -5.66467 percentage points. These
nonlinear output differences establish a *current-enabled geometry effect*
in this controlled finite process; they are not a literal percentage
decomposition of microscopic reversible and irreversible mechanisms.
T4 is still a geometric readout, not itself a pure antisymmetric operator.

Complete inputs, compressed matrices, scores and hashes are in
`results/p398-width8-geometry-current-interaction/latest.json`.
No new candidate was selected, and no independent statistical, continuum,
Jordan, or morphism-history evidence is being claimed.
