# Mesoscopic coexistence: a dilute black gas against a pure full background

The combined theorem and root consequences are summarized in
[the mesoscopic two-cloud synthesis](closed-source-mesoscopic-two-cloud-synthesis.md).

Let k denote the geometric dilation, so the equal-area pair has
`N=25 k^2`; keep m for the coupling. Assume

```
m -> infinity,       N/m^4 -> 0,       ell_N -> infinity,
h=(1-m^-2)^(-1) exp(s/N),              |s|<=S.                (1)
```

The argument below does not require N/m^2 to stay bounded. It therefore
covers both the special path k=m and the genuinely mesoscopic paths on
which N/m^2 diverges. Under the growing-cutoff hypotheses stated below,
the short-contour ensemble splits exactly into rank0/vacant background
and rank2/full background. Uniformly for bounded s,

```
Z0_good = (1-m^-2)^(-N) [1+O_S(m^-2+N/m^4)],

Z2_good/h^N = 1 + N h^-1 m^-4
                 +O_S(N/m^6+(N/m^4)^2),                     (2)

Z2_good/Z0_good = exp(s)[1+O_S(m^-2+N/m^4)].                 (3)
```

Thus rank0 may contain a mesoscopically large number `~N/m^2` of black
particles, while rank2 is conditionally pure: its probability of any
hole is `N/m^4[1+o(1)]`. The common bulk factor in (2) is what permits
coexistence beyond the old bounded-N/m^2 Poisson regime.

This note starts at `410015f5505dc2d8ca0e9ac904f656a4adc9fe86` and
uses the same occupied-NN / matching-white resolved-corner convention as
`closed-source-poisson-double-scaling.md`. It does not enumerate a new
configuration, evaluate a coupling grid, or reuse a fixed-cutoff constant
as though it were uniform in the growing cutoff.

## 1. The growing good event and its exact two backgrounds

Use the resolved mixed-edge curves. Let

```
A_N(m)=108N/[m^2(1-3/m)]
```

and, for an arbitrary fixed M, take the declared cutoff

```
R_N=max(16, ceil{
 [A_N(m)+(M+2)log N+2log m+10]/log(m/3)}).                   (4)
```

The geometric/scaling input is that, eventually,

```
R_N<ell_N,       R_N/m^2 ->0,       R_N/N ->0.               (5)
```

The already proved all-activity contour inequality then makes the
relative long-contour contribution `O(N^-M)`, uniformly for |s|<=S.
M can be fixed as large as required before the limit. This paragraph
uses the finite inequality with its actual `exp(A_N)` factor; it does
not claim that the old numerical choice R=16 is uniform here.

On the good event every resolved curve has length below ell_N and is
therefore contractible. Lift an outermost curve to the square lattice.
It bounds a disk, and its projected disk is embedded because it is a
component of the covering of a contractible region. Disjoint outermost
curves have disjoint disks. After removing those disks, the torus has
one connected genus-one exterior. The occupation colour on that exterior
is constant.

- If the exterior is vacant, every occupied NN component lies in a
  lifted disk, so its ambient image is zero. The configuration has r=0.
- If the exterior is occupied, the black exterior is one NN component
  carrying both torus directions. Every other black component is an
  island in a white hole. The total image has r=2.

This proves an exact good-family dichotomy, not a dominant-background
heuristic. The corner resolution is important: it preserves black NN
components while letting diagonally touching white sites belong to the
matching-white region. No step assumes that a white hole is NN-connected.
Nested curves are retained in the patch weights below.

## 2. Uniform disk and nested-island bound

A lifted simple dual curve gamma of length l<R lies in a coordinate box
with at most

```
Area(D_gamma) <= C0 l^2,       C0=4                              (6)
```

primal sites. This deliberately loose constant follows from the spans of
a length-l lifted walk and is uniform in the quotient and in l.

For an outer white hole in a black/full background, let V be its total
number of vacant sites, let i be the number of black NN island components
strictly inside it, and let the total length of its outer and internal
curves be B. Relative to the full configuration its **exact** activity is

```
h^-V m^(-B+2i).                                               (7)
```

If the outer curve has length l, every internal black component has at
least one resolved boundary curve. Hence, after assigning one cluster
reward m^2 to one boundary of each such component,

```
m^(-B+2i) <= m^-l product_(internal curves eta) m^(2-|eta|). (8)
```

This keeps the actual nested islands; it does not replace them by an
unconditioned gas of occupied sites. Assigning m^2 to every internal
curve only enlarges the right side and avoids any assumption about the
NN connectivity of intervening white regions.

There are at most `4A 3^(n-1)` rooted internal dual walks of length n
in a disk of area A. Dropping disjointness, nesting and cut consistency
therefore bounds the entire internal partition by

```
exp{108 C0 l^2/[m^2(1-3/m)]}.                                (9)
```

Moreover, for m^2>=2,

```
|log h| <= 2/m^2+S/N.                                       (10)
```

Thus the loss h^-V and every internal island together multiply the outer
activity m^-l by at most

```
exp(c_N l^2),
c_N=C0[2/m^2+S/N+108/{m^2(1-3/m)}].                          (11)
```

Equations (5) imply `c_N R_N ->0`. This is the uniform point that a
fixed-R proof did not provide. For all l<R_N the exponent in (11) is
`o(1)l`, so the outer-curve sum retains its geometric ratio
`3 exp(o(1))/m`. All constants below are uniform for |s|<=S and for the
fixed finite collection of equal-area quotient sequences.

## 3. Rank2/full side: exact leading one-hole term

The outer boundary of a white hole has length at least4. A length-four
dual edge-cycle encloses exactly one vacant primal site; it has no room
for an internal curve. There are exactly N such translated cycles. Its
activity is

```
a4=h^-1 m^-4.                                                (12)
```

The configuration with one site removed from full occupancy remains a
connected rank2 black configuration. One may reroute any fundamental NN
walk around the missing site through its injective local square stencil;
the detour preserves its deck displacement. It has Bmix=4,C_B=1,r=2,
confirming (12) directly.

The square dual graph has no odd closed walk. Every other outer curve
therefore has length at least6. By (11) and the usual rooted-walk count,

```
W_ge6 <= (4N/3) sum_(l>=6)
              [3 exp(c_N R_N)/m]^l = O(N/m^6).                (13)
```

Drop compatibility between outer holes. The full state and the N
single-hole states give the lower bound `1+Na4`; the relaxed polymer
product gives the upper bound `exp(Na4+W_ge6)`. Since N/m^4->0,

```
Z2_good/h^N
 =1+N h^-1m^-4+O(N/m^6+(N/m^4)^2).                          (14)
```

In particular the full configuration has conditional rank2 probability

```
1-N/m^4+o(N/m^4).                                           (15)
```

The estimate includes diagonal matching-white contacts and every nested
black island through (8)-(11). A diagonal pair of vacant sites is not
silently treated as two NN-connected holes; it simply appears through
its resolved outer/internal curve family and is covered by (13).

For two equal-area geometries with injective local stencil, the leading
single-hole term is the same N times the same local activity. It cancels
from their full-side shape ratio:

```
[Z2_good,f/h^N]/[Z2_good,s/h^N]
 =1+O(N/m^6+(N/m^4)^2).                                    (16)
```

This is a parameter-free shape prediction of the pure full side, not an
assertion that the rank0 gases already agree to the same smaller order.

## 4. Rank0 side: a dilute gas with common bulk pressure

Now let the genus-one exterior be vacant. An outermost black component
with boundary l receives one cluster factor m^2, so its outer activity
is `m^(2-l)` times its h factors. The same disk/internal-curve argument
as (8)-(11) applies with black and white interchanged while retaining
the prescribed corner resolution.

The length-four patches are precisely isolated occupied sites, of
activity

```
x=h/m^2=exp(s/N)/(m^2-1).                                   (17)
```

Mutually compatible singleton patches form the NN independent-set gas.
Its partition differs from `(1+x)^N` by relative logarithm `O(Nx^2)`.
Every nonsingle outer black patch has l>=6, and its extra outer factor
m^2 changes (13) to a total `O(N/m^4)`. Dropping compatibility gives

```
log Z0_good = N log(1+x)+O(N/m^4).                           (18)
```

This coarse bound deliberately does not claim that every adjacent tree
is a new interaction: some tree weights coincide with their independent
particle weights. It is sufficient and uniform without a separate
large-component catalogue.

Put

```
G_N=(1-m^-2)^(-N).
```

From (17), exactly

```
N log(1+x)
 = log G_N + N log[1+(exp(s/N)-1)/m^2]
 = log G_N + s/m^2+O_S(1/(Nm^2)).                            (19)
```

Combining (18)-(19),

```
Z0_good=G_N[1+O_S(m^-2+N/m^4)].                             (20)
```

The density tends to zero like m^-2. Its total black count need not stay
bounded: the mean scale N/m^2 may diverge. This is the mesoscopic dilute
gas that replaces the finite Poisson cloud of the bounded-N/m^2 theorem.

The full state has weight

```
h^N=G_N exp(s).                                              (21)
```

Equations (14),(20),(21) prove (2)-(3). The factor G_N can diverge like
`exp(N/m^2)`; it cancels between the two phases rather than being bounded
away by a fixed-volume argument.

## 5. Root and shape consequences

The growing-cutoff bound makes every rank1 configuration bad, because
rank1 has an essential boundary of length at least ell_N>R_N. With M
chosen in (4), the relative bad contribution is `O(N^-M)`. Define

```
eta_N=m^-2+N/m^4+N^-M.
```

Uniformly for bounded s, each geometry therefore has

```
q_g(s)=tanh(s/2)+O_S(eta_N),
P1_g=O(N^-M),
P0_g=1/[1+exp(s)]+O_S(eta_N),
P2_g=exp(s)/[1+exp(s)]+O_S(eta_N).                            (22)
```

The pooled root is consequently

```
s_root=O(eta_N),
N log h_root=-N log(1-m^-2)+O(eta_N),                         (23)

logit p_root
 =log m-log(1-m^-2)+O(eta_N/N).                              (24)
```

For the special path k=m, so N=25m^2, this becomes

```
s_root=O(m^-2+N^-M),
logit p_root=log m+m^-2+O(m^-4+N^(-M-1)).                    (25)
```

The coefficient at order m^-4 in (25) is not fixed by this bound.

For an equal-area geometry pair the restricted sector cross-ratio obeys

```
Xi_good=Z2_f Z0_s/(Z0_f Z2_s)=1+O(eta_N),                    (26)
```

while the pure full-side ratio has the sharper cancellation (16). These
are the new shape predictions. They say that a geometry-dependent leading
full-hole amplitude is absent: any observed full-side difference above
`O(N/m^6+(N/m^4)^2)` would contradict the short-curve mechanism or its
scope assumptions.

Since E=q^2 equals1 on both good sectors, original angular E differences
still require rank1/bad transmission. Equations (22)-(26) do not assign
their subleading sign or a continuum exponent. They do show that a large
rank0 particle count is compatible with an asymptotically pure rank2
phase and a common root: the two statements concern conditional structure,
not equality of defect counts.

## Scope

The proof requires (5), an injective square stencil, the existing resolved
corner convention, bounded s, and a finite collection of equal-area
geometries. It is a finite-curve argument followed by a joint limit; it
does not settle fixed m, interchange infinite-volume and strong-coupling
limits, or prove a Potts/CFT interpretation. If the cutoff cannot satisfy
R_N<ell_N and R_N/m^2->0, the internal factor (11) is the exact obstruction:
the present contour sum then no longer yields the claimed uniform bounds.
