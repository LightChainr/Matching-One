# The canonical joint interaction reaches original U beyond adjacent contacts

**Both frozen zero-transmission closures are excluded.** For the unchanged
canonical regular local tensor, the complete joint first-Q response on
the original N25 pair is strictly negative. Its nonadjacent contribution
is also strictly negative and larger in magnitude than the adjacent
contribution. This closes the finite joint-interaction-to-U interface;
it does not identify a continuum field.

## 1. The completed numerical decision

At every original vacant vertex insert `epsilon*Kreg/N`, with
`Kreg=average_C4 i(I-P1)i^dagger` and `c(Q)=1` identically. Keep the
original occupied summand, rank, q/E, separately normalized geometries,
pooled moving root and thermal-slope denominator. The frozen quantity is

```
J2 = d_logQ d_epsilon^2 U | Q1,epsilon=0.
```

| Fixed pair contribution | J2 in original-U units | J2/A25 | Exact interval sign |
|---|---:|---:|---|
| Total | **-0.0055194314248394015** | -0.00005905706006949678 | strictly negative |
| Adjacent NN pairs | -0.0017510744544027990 | -0.000018736225034779993 | strictly negative |
| All other distinct pairs | **-0.0037683569704366022** | -0.00004032083503471679 | strictly negative |

Here `A25=25^(13/8)/2`. The strict rational upper bounds include

```
J2_total/A25 <= -147642650173741949172258069338667607
                 /2500000000000000000000000000000000000000 < 0,

J2_nonadjacent/A25 <= -80641670069433570128095233537434551
                       /2000000000000000000000000000000000000000 < 0.
```

The bounds are finite exact-population enclosures, not confidence
intervals. [The complete score](../results/p337-regular-pair-joint/score/score.json)
retains the bounds, source polynomials and all response components.
The adjacent/nonadjacent split was fixed before collection and sums
linearly to the total; it is not two independent evidence blocks.

The primary result rejects an additive first-Q effective log weight
linear in the specified epsilon as a global closure of this model.
The secondary result rejects the claim that its first-Q joint U
response is carried only by NN contacts. Nonadjacent here includes
plaquette diagonals on N25, not just macroscopically separated sites.

## 2. The actual two-insertion source, not a product of marks

For original occupation A, let g_xy be the first-Q derivative of the
relative joint colour contraction at two vacant marks. The new
[joint-response derivation](regular-pair-joint-original-u-response.md)
proves

```
s2(A) = (2/N^2) sum_{x<y} g_xy(A),
J2 = L[s2],
E[f s2] = E[f (1/N) sum_{y!=0} g_0y]   for translation-invariant f.
```

The last identity supplies all K/q/E moments and their thermal
derivatives without visiting all N(N-1) ordered pairs per configuration.
In the producer `b16=sum_y 16*g_0y`; every source column is divided by
`16*N=400` exactly once. The second derivative's factor two is already
included. There is no same-site term in a tensor linear at each vertex.

Every nonempty Kreg insertion is pointwise zero at Q1. Hence products
of separate one-site responses start at order `(Q-1)^2`. Their covariance
cannot replace the first-Q joint contraction. The surviving s2 response
still requires the original covariance centering, thermal derivatives,
root shift and slope response.

## 3. How the negative response is transmitted

The total response decomposes in the original activity coordinate as:

| Required original-U term | Contribution |
|---|---:|
| Direct centered thermal/source response | -0.005529496208369856 |
| Pooled-root motion | -0.0014591233330132293 |
| Source change of thermal slope | +0.000007758539927132672 |
| Root change of thermal slope | +0.0014614295766165514 |
| **Total** | **-0.0055194314248394015** |

The full mixed probability-root shift is
`d_logQ d_epsilon^2 p0=+0.00037853828485651174`.
The two root terms nearly cancel, but neither was dropped. Both fixed
pair classes transmit with the same negative sign; this is not a total
created by cancellation between adjacent and nonadjacent classes.

At this same root the integrated source means are **positive** in both
geometries: `E[s2]=0.005796796990976299` and `0.005636471996260121`.
Their nonadjacent means are positive as well. Thus even a positive
occupation-averaged joint susceptibility does not determine the sign
of the shape-sensitive, root/slope-normalized U response. The example
is inside one fixed canonical model, without adjusting a counterterm.

## 4. Two additional exact mechanism distinctions

**A contact still needs an external bypass.** The
[adjacent-kernel proof](regular-pair-adjacent-joint-kernel.md) shows that
adjacent vacant vertices share one physical edge-node. Counting that
one component once makes the existing Bell8 table and Q^(-|pi|)
normalization apply unchanged. That shared edge alone has g=0. One
additional shared outside component gives a nonnegative contrast
product; multiple bypass components can give negative weights. An
explicit rectangular-torus occupation realizes g=-2. These statements
explain why neither ignoring contacts nor assigning them a universal
positive sign is justified.

**The joint directions are not copies of temperature or the old source.**
Three elementary K=2 occupations have the same one-insertion activation
`a=(N-2)/N`, but different joint sources:

| Occupied sites, all others vacant | s_adjacent | s_nonadjacent |
|---|---:|---:|
| NN domino | 1/N^2 | 0 |
| Opposite plaquette corners | 0 | 1/(2N^2) |
| Straight distance-two pair | 0 | 0 |

They prove linear independence of these two joint directions modulo
`span{1,K,a}`. In fact no function of (K,a) alone distinguishes the
three occupations. The numerical result now shows that the original-U
functional is nonzero on both of those specified directions.

The [distinct-site coordinate argument](regular-pair-offdiagonal-coordinate-invariance.md)
adds a useful scope distinction. Off-diagonal local-parameter Hessians
are invariant under independent unit-Jacobian local reparameterizations.
Their nonzero sum excludes a first-Q action separable in the individual
site parameters, even if each one-site dependence is nonlinear. A
uniform nonlinear relabeling of epsilon adds a diagonal path-acceleration
term proportional to the old first response W; that would be a different
full path derivative, not removal of the distinct-site Hessian measured
by the declared linear tensor convention.

## 5. Execution, dependencies and the completed boundary

The contract was frozen at `4ce4dfe894c9fe96f268c61cf21eb6585dba5418`,
the scorer at `5da4749245450048625a2da43e8f73da1ee9275c`. The full response
proof `93651d61` is public as `8212a6b3`; the adjacent proof `4bdf275f`
as `6d1e453e`; producer `30891e04` as `8771d6ec`. The coordinate argument
was committed at `401db222` before production. Accepted proofs and code
preceded a single explicit production GO.

The two complete N25 traversals took 3.37178 and 3.19402 seconds on two
local workers. Each covers exactly 2^25 configurations, collecting only
the new joint crossmoments. The raw commit `f67d7646` is public as
`a2307d67e48da4bbd18ff7fd47f731eecdaca560`. The scorer ran once at that
head in 0.615483625 seconds, importing the old root, denominator and U/A
unchanged. [The raw receipt](../results/p337-regular-pair-joint/run.json)
records commands and hashes. No old-source score, spatial MC rerun,
parameter scan, root search, production restart or cloud operation was
performed. The local managed research Python environment was reused.

**Science card:** original moving-root U / canonical C4 regular local
colour interaction / distinct-site first-Q Hessian / N25 `(5,0),(4,3)`.
These are new joint contractions on the same complete exact populations
as earlier N25 work, not independent stochastic confirmation. Total and
its fixed parts share all dependencies. The independently sampled L32/L64
[spatial result](regular-pair-spatial-transmission-result.md) remains a
different observer and is not pooled into this score.

The fixed joint-global and contact-only decisions are complete. This
does not assign H4/Jordan identity, establish an asymptotic exponent, or
show that the observed joint term generates the original anomaly. The
next field claim must supply a named projection and a discriminating
prediction. No further N25 source mixture or completion coefficient is
added to improve this result; stopped P154/P334/F4 production stays stopped.
