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

Thus the next term after the universal two-cloud center is explicitly
model-count dependent: occupied-cycle and adjacent-vacancy polymer
coefficients, plus any residual geometry-dependent bad sector, enter
through `d_g`. Three facts do not determine it: equal area, equal total
pressure, or the existence of the sigmoid.

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
   `u_root=o(1)`. It is not claimed here to be a logically necessary gate
   for every possible observable, since a further exact cancellation of
   the two connected-polymer coefficients would have to be ruled out to
   establish necessity.

The hierarchy is therefore

```text
N/m^4 ->0  => pure full state + Poisson(lambda) rank-zero cloud,
N/m^6 ->0  => two resummed dilute clouds + centered sector sigmoid. (24)
```

The first implication contains the second. Neither is a fixed-m
thermodynamic theorem, and neither assigns the subleading sign of original
U.

## Scientific card

- **Mechanism changed:** the mesoscopic root balances an occupied-site
  gas against a vacancy gas; it is not empty-versus-full at first
  correction order.
- **New parameter-free prediction:** `h_bar=1+m^-2`, equivalently
  `s_root/(N/m^4)->-1` under the pure-full gate and rate-resolved errors.
- **Model-dependent remainder:** the imbalance between occupied-cycle and
  adjacent-vacancy connected polymers in (18).
- **Boundary:** `N/m^4->0` is sharp for full-state purity; the resummed
  root/sigmoid only needs the more natural sufficient gate `N/m^6->0`.
- **Not claimed:** a fixed-m transition, a continuum field identity, an
  independent data block, or the sign of the vanishing original U.
