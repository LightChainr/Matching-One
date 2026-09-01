# The mesoscopic closed-source root is centered at `h=1+m^-2`

## Result

Put

```text
a=m^-2,       lambda=N a,       tau=N a^2=N/m^4.
```

On the growing N25 axis/tilted pair, the first useful mesoscopic
approximation is not empty versus full. It is a dilute occupied-site gas
versus a dilute vacancy gas on the full background. Their restricted
partition functions are

```text
R0(h)=(1+a h)^N,
R2(h)=(h+a^2)^N.                                                (1)
```

Their equality has the exact solution

```text
h_bar=(1-a^2)/(1-a)=1+a.                                      (2)
```

This gives a new universal finite-coupling fingerprint for the named
closed source. In the older one-cloud chart

```text
h_c=(1-a)^-1,       h=h_c exp(s/N),                            (3)
```

the same center is

```text
h_bar/h_c=1-a^2,
s_bar=N log(1-a^2)=-N/m^4+O(N/m^8).                           (4)
```

Consequently, when `tau ->0` and the local-gas remainder is `o(tau)`,
the actual pooled root obeys

```text
s_root=-tau+o(tau),
h_root/h_c=exp[-m^-4+o(m^-4)],
log h_root=m^-2-(1/2)m^-4+o(m^-4).                             (5)
```

The coefficient `-1` in `s_root/tau`, equivalently `-1/2` in the second
term of `log h_root`, is universal within this fixed `Sstar` completion.
It comes from the exactly N single-vacancy configurations. It is not a
new field amplitude and need not persist after changing the source law.

The power gate `N/m^4 ->0` is therefore correctly sharp if the conclusion
includes a **pure full configuration conditional on rank two**. It is not
the weakest gate for the root or sector sigmoid. Once the vacancy gas in
(1) is retained, the more natural sufficient interaction gate is

```text
N/m^6 ->0,                                                     (6)
```

together with the independently supplied growing-systole contour error.
This permits `tau` to be bounded or even divergent and still gives a
two-sector sigmoid after centering at `h_bar`. Purity of the full state is
then deliberately no longer asserted.

This note only develops consequences of the dilute/full asymptotics and
the activity-uniform contour input. It introduces no configuration,
simulation, coupling scan or fitted exponent.

## 1. Why the two reference gases have exactly these activities

For `Sstar`, a rank-zero occupied graph has

```text
g=B_mix-2 C_B=2K-2 beta_1,
w=h^K m^-g=(a h)^K m^(2 beta_1).                               (7)
```

Every occupied forest therefore has exactly the product-site weight
`(a h)^K`; adjacency by itself does not change its weight. The first
departure from the rank-zero product gas requires an occupied cycle. On
the square local graph its smallest representative is a four-cycle. Its
product activity is order `a^4`, and its cycle reward is `m^2`, so the
first connected pressure correction is order

```text
N a^4 m^2=N/m^6.                                               (8)
```

Essential cycles are not inserted into this local estimate: they are
part of the separately controlled bad-contour/rank-one remainder.

For a rank-two occupied background with v vacancies, divide by the full
weight `h^N`. A single vacancy has four mixed edges and exact relative
activity

```text
z=h^-1 m^-4=h^-1 a^2.                                         (9)
```

N independent vacancy positions give

```text
h^N(1+z)^N=(h+a^2)^N.
```

The first failure of this product gas is an adjacent vacancy pair: its
boundary has length six instead of eight. Relative to two independent
vacancies it receives an `m^2` enhancement, again producing a connected
pressure correction `O(N/m^6)`. Longer connected polymers are bounded by
the same local contour expansion for large m. Thus, uniformly on bounded
root charts and separately in each geometry,

```text
Z0_g=(1+a h)^N exp[r0_g],
Z2_g=(h+a^2)^N exp[r2_g],
max_g |rj_g| <= C N/m^6 + epsilon_bad,N.                        (10)
```

Here `epsilon_bad,N` denotes the additive log-partition effect obtained
after converting normalized bad probability to restricted relative
partition control. Equation (10) is stronger than merely saying that the
two total pressure densities agree. It controls each restricted sector.
On the bounded two-phase chart both reference sectors have a fixed
fraction of the total, so `P(bad)<=epsilon` implies a relative restricted
error `O(epsilon)`. Without this denominator step, a small normalized bad
probability or equal total pressure would not locate the pooled root.

For the N25 sequences the growing-cutoff contour estimate is
superpolynomial, hence smaller than either displayed power error in the
regimes used below. If only `epsilon_bad,N=o(1)` is available, retain it
explicitly in every root bound rather than silently assigning it a power.

## 2. Exact root chart and the limiting sigmoid

Use the natural two-cloud chart

```text
h=h_bar exp(u/N),       h_bar=1+a.                             (11)
```

The reference rank-two versus rank-zero log odds are

```text
L_N(u)=N log[(h+a^2)/(1+a h)].                                 (12)
```

They satisfy `L_N(0)=0` exactly and

```text
L_N'(0)=h_bar/(h_bar+a^2)-a h_bar/(1+a h_bar)
       =1-a+O(a^2) ->1.                                       (13)
```

Therefore `L_N(u)=u+o(1)` uniformly for bounded u. With (10), vanishing
rank-one mass and the same local reference terms in both geometries,

```text
Z_g/[1+a h]^N = 1+exp[L_N(u)]+o(1),
P0_g ->1/(1+e^u),       P2_g ->e^u/(1+e^u),
q_g=P2_g-P0_g ->tanh(u/2).                                    (14)
```

This is the genuine two-term exponential form. It is a statement about
restricted sector partition functions, not an inference from matching
pressures. The pooled root is unique and monotone, so evaluating the
limiting signs at `u=+-eta` brackets it and gives

```text
u_root=O(N/m^6)+O(epsilon_bad,N).                              (15)
```

In the old chart (3), an exact reference expression before taking the
limit is also useful:

```text
R0/h_c^N = exp{N log[1+a(exp(s/N)-1)]},
R2/h_c^N = exp(s) (1+a^2/h)^N.                                (16)
```

For bounded s and `tau->0`, the first exponent is `a s+o(tau)` while the
second extra log is `tau+o(tau)`. Hence

```text
Z_g/h_c^N = 1+exp(s+tau)+o(tau),                              (17)
```

which proves (5), provided the bad/restricted remainders are `o(tau)`.
This last rate condition is essential for resolving the coefficient
`-1`; an undifferentiated `O(tau)` partition estimate locates the scale
but not its universal constant.

More generally define

```text
d_g(u)=r2_g(u)-r0_g(u).
```

A first-order implicit-root expansion gives

```text
u_root=-mean_g d_g(0)/(1-a+O(a^2)) + O(max_g |d_g|^2).          (18)
```

Thus the next term after the universal two-cloud center is determined by
occupied-cycle and adjacent-vacancy polymer counts, plus any residual
geometry-dependent bad sector, through `d_g`. Equal area, equal total
pressure, and the existence of the sigmoid do not determine those counts.
Section 5 now evaluates the first two counts for the square local graph.

## 3. Mesoscopic conditional laws

At the reference root `h_bar=1+a`, the rank-zero product occupation
probability and rank-two vacancy probability are

```text
theta0=a(1+a)/(1+a+a^2)=a-a^3+O(a^4),
theta2=a^2/(1+a+a^2)=a^2-a^3+O(a^5).                          (19)
```

Under the pure-full gate `tau=N a^2 ->0`,

```text
N theta0=lambda+o(1),       N theta0^2=tau+o(tau),
N theta2=tau+o(tau).
```

An elementary Bernoulli-to-Poisson coupling changes at most
`N theta0^2` probability, so conditional on rank zero

```text
K = Poisson(lambda)+o_TV(1).                                  (20)
```

This remains meaningful when `lambda->infinity`: its centered and scaled
characteristic function converges to `exp(-t^2/2)`, giving

```text
(K-lambda)/sqrt(lambda) -> Normal(0,1),
P(K=0 | r=0)=exp[-lambda+o(1)] ->0.                            (21)
```

Conditional on rank two, the vacancy mean is `tau+o(tau)`, so the full
configuration has probability `1-tau+o(tau) ->1`. At the pooled root,
both geometries have

```text
P0=1/2+o(1),       P2=1/2+o(1),       P1=superpolynomially small. (22)
```

The two rank weights are one half even though rank zero is a large
mesoscopic cloud and rank two is a single configuration. The denominator
also remains macroscopic: the two conditional K means differ by
`N-lambda+o(N)`, so `Cov(q,K)/N ->1/2` and the within-geometry sector
variance `kappa->1`. Together with the existing rank-one bound, the
original separately normalized pooled U is superpolynomially small. This
is not obtained by replacing the within-geometry denominator with a
mixture variance.

If only (6) holds and tau does not vanish, equations (14)-(15) still give
the centered two-sector sigmoid. The conditional rank-two phase then has
a vacancy cloud of mean asymptotic to tau; it is intentionally not called
the full configuration. If tau diverges, its vacancy count has its own
Gaussian limit after centering/scaling, while its density still vanishes.

## 4. Which power gate is actually sharp?

There are two different questions.

1. **Pure-full statement.** The ratio of the total one-vacancy weight to
   the full-state weight is exactly

   ```text
   N/(h m^4)=tau[1+o(1)].                                     (23)
   ```

   Therefore conditional concentration on the single full configuration
   forces `tau->0`. The local expansion above also proves sufficiency.
   Thus `N/m^4->0` is the correct and weakest power gate for this purity
   conclusion (within the fixed chart `h->1`).

2. **Root/sigmoid statement.** Single vacancies need not be discarded;
   resumming them changes `h_c` to the exact two-cloud center `1+a`.
   The first interactions not represented in the two product gases occur
   at order `N/m^6`, from a black four-cycle or adjacent vacancy pair.
   Consequently `N/m^6->0`, plus the separate bad-contour control, is the
   natural weaker sufficient gate for the clean centered sigmoid and
   `u_root=o(1)`. Section 5 shows that its first two connected-polymer
   coefficients do not cancel. Hence this gate is locally sharp for the
   named root chart as `N/m^6 ->0`; no global necessity statement is made
   for arbitrary non-small `N/m^6`, where further resummation is required.

The hierarchy is therefore

```text
N/m^4 ->0  => pure full state + Poisson(lambda) rank-zero cloud,
N/m^6 ->0  => two resummed dilute clouds + centered sector sigmoid. (24)
```

The first implication contains the second. Neither is a fixed-m
thermodynamic theorem, and neither assigns the subleading sign of original
U.

## 5. The first connected-pressure coefficients do not cancel

The `N/m^6` terms can be computed without another polymer family. Work at
the two-cloud center `h_bar=1+a`, and let

```text
theta0 = a h_bar/(1+a h_bar),
theta2 = (a^2/h_bar)/(1+a^2/h_bar)=a^2/(h_bar+a^2).             (25)
```

These are the product-gas probabilities of an occupied site in the
rank-zero reference and a vacancy in the rank-two reference.

### Rank zero: one elementary occupied cycle per face

A forest has `beta_1=0` and exactly its reference product weight. The
smallest departure is all four vertices of one elementary square occupied.
The occupied NN graph then has `beta_1=1`, so (7) multiplies this event by
`m^2`. The Mayer increment relative to the forest/product gas is therefore

```text
(m^2-1) theta0^4.                                              (26)
```

There are exactly N elementary faces. The `-1` is the necessary forest
subtraction: it removes the product-gas weight already present in `R0`.
It is order `N/m^8`, while the rewarded term is order `N/m^6`. Two
separated plaquettes exponentiate and disappear from the log; overlapping
cycles first contribute at order `N/m^8`. Consequently

```text
r0 = N(m^2-1)theta0^4 + O(N/m^8) + o_bad(N/m^6)
   = N/m^6 + o(N/m^6).                                        (27)
```

The coefficient is `c0=1`. It is geometry independent once the elementary
square embeds, as it does on both growing N25 sequences.

### Rank two: two adjacent-vacancy polymers per site

Two independent vacancies have relative product weight
`h^-2 m^-8`. If they occupy the endpoints of one NN edge, their union has
six rather than eight mixed boundary edges and exact weight
`h^-2 m^-6`. Thus each NN edge has Mayer increment

```text
(m^2-1) theta2^2.                                              (28)
```

The square torus has exactly `2N` unoriented NN edges. The subtraction by
one again removes the pair already included in the vacancy product gas.
Connected triples and overlapping adjacent pairs begin at `N/m^8`, so

```text
r2 = 2N(m^2-1)theta2^2 + O(N/m^8) + o_bad(N/m^6)
   = 2N/m^6 + o(N/m^6).                                       (29)
```

The coefficient is `c2=2`. Therefore the first restricted log-odds
correction is

```text
d=r2-r0=N/m^6+o(N/m^6).                                       (30)
```

It is positive: the vacancy/full sector receives one more unit of
connected pressure than the occupied/rank-zero sector. Since the slope in
(13) tends to one, the actual pooled root in the natural chart has the
parameter-free displacement

```text
u_root=-N/m^6+o(N/m^6).                                       (31)
```

Equivalently, refining (5),

```text
h_root=(1+m^-2) exp[-m^-6+o(m^-6)],
log h_root=m^-2-(1/2)m^-4-(2/3)m^-6+o(m^-6).                  (32)
```

In the old chart `h=h_c exp(s/N)`,

```text
s_root=-N/m^4-N/m^6+o(N/m^6),                                 (33)
```

because `N log(1-m^-4)` has no `N/m^6` term. Resolving these constants
requires the restricted bad-sector error to be `o(N/m^6)`; a merely
vanishing or `O(N/m^6)` error is insufficient. The growing-systole contour
input supplies that separation here.

The noncancellation `c2-c0=1` makes `N/m^6` the first actual obstruction
to `u_root=o(1)` in the named square-lattice model. The individual numbers
1 and 2 are local incidence counts, not universal across a changed lattice
or completion. Their difference is parameter free only after the model,
port semantics and source have been fixed.

## Scientific card

- **Mechanism changed:** the mesoscopic root balances an occupied-site
  gas against a vacancy gas; it is not empty-versus-full at first
  correction order.
- **New parameter-free predictions:** `h_bar=1+m^-2`,
  `s_root/(N/m^4)->-1` under the pure-full gate, and the first connected
  correction `u_root/(N/m^6)->-1` when errors are rate resolved.
- **First model count:** one elementary face versus two NN edges per site,
  so `c2-c0=2-1=1`; it does not cancel.
- **Boundary:** `N/m^4->0` is sharp for full-state purity; the resummed
  root/sigmoid only needs the more natural sufficient gate `N/m^6->0`.
- **Not claimed:** a fixed-m transition, a continuum field identity, an
  independent data block, or the sign of the vanishing original U.
