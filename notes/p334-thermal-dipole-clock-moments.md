# Thermal dipoles separate lifetime mass, plateau center, and clock organization

The complete E response has an exact moment interpretation. Its area is
minus the mean lifetime response, its first moment is minus the response
of **lifetime-weighted birth center**, and its second moment contains both
lifetime-weighted center width and within-window width. These give
parameter-free tests of a rigid center-translation mechanism.

An important distinction follows: `Cov(C_tau,W_tau)` is already determined
by the two endpoint marginal variances. It is not, by itself, new evidence
about their joint dependence. The genuinely additional second-order joint
coordinate is `Cov(tau1,tau2)=Var(C_tau)-Var(W_tau)/4`.

## Exact continuous-priority representation

Work first in **one orientation**, under either frozen common-label tangent
from cfaae36c/ffb70969. Generate N independent U(0,1) priorities, independently
of the tilted label permutation, and attach their sorted values to its
positions. The policy depends on label order/contact geometry, not on these
numeric order statistics. Therefore the usual conditional order-statistic
law remains unchanged when the permutation law is tilted.

Let `tau1=U_(K1) <= tau2=U_(K2)`,
`C_tau=(tau1+tau2)/2`, and `W_tau=tau2-tau1`. Direct births have W_tau=0.
For the frozen even-topology convention,

\[
 E_t(p)=1-P_t(\tau_1\le p<\tau_2)=1-F_{1,t}(p)+F_{2,t}(p).
\]

Integrating the random plateau interval, then taking expectation, proves
for every integer m>=0

\[
 \boxed{\int_0^1p^m E_t(p)\,dp
  ={1\over m+1}
   -\mathbb E_t\!\left[{\tau_2^{m+1}-\tau_1^{m+1}\over m+1}\right].}
\]

Write H_X for `partial_t E_t[X]|0`, allowing a plus or minus source
direction, and `I_m=int p^m partial_t E_t(p)|0 dp`. Since the finite
permutation mixture is differentiable, differentiation commutes with
these bounded integrals. Then

\[
 I_0=-H_{W_\tau},\qquad
 I_1=-H_{C_\tau W_\tau},\qquad
 I_2=-H_{C_\tau^2W_\tau+W_\tau^3/12}.
\]

The last equality uses
`tau2^3-tau1^3=3 C_tau^2 W_tau+W_tau^3/4`. Thus I2 alone does **not** isolate
center width: it also contains the third lifetime moment.

## Discrete K to continuous moments: retain the finite-N corrections

Set `x^overline{r}=x(x+1)...(x+r-1)` and `D_r=(N+1)^overline{r}`.
Conditional on a permutation and its K1,K2,

\[
 \mathbb E[\tau_j^r\mid K_1,K_2]={K_j^{\overline r}\over D_r},\qquad
 \mathbb E[\tau_1^a\tau_2^b\mid K_1,K_2]
 ={K_1^{\overline a}(K_2+a)^{\overline b}\over D_{a+b}}.
\]

For K1<K2, condition on tau2: the ratio tau1/tau2 is an independent
Beta(K1,K2-K1) variable. Its rising-factorial moment cancels the first a
factors of the tau2 moment, proving the second formula. The equality case
K1=K2 follows directly and agrees with the same formula.

In particular
`E[tau1 tau2|K]=K1(K2+1)/((N+1)(N+2))`, not the product of the two
conditional means. With discrete `C=(K1+K2)/2`, `W=K2-K1`, all quantities
on the right of the following table are evaluated per permutation:

| Conditional continuous moment | Exact discrete readout |
|---|---|
| E C_tau, E W_tau | C/(N+1), W/(N+1) |
| E C_tau^2 | (C^2+C-W/4)/D2 |
| E W_tau^2 | W(W+1)/D2 |
| E(C_tau W_tau) | W(C+1/2)/D2 |
| E W_tau^3 | W(W+1)(W+2)/D3 |
| E(C_tau^2 W_tau) | W(C^2+2C+1/2-W/4)/D3 |

The gap moment follows from the Dirichlet spacings:
`E[W_tau^r|K]=W^overline{r}/D_r`, with all positive moments zero when W=0.
For example,

\[
 \operatorname{Cov}(C_\tau,W_\tau\mid K)
 ={W((N+1)/2-C)\over (N+1)^2(N+2)}.
\]

Consequently population continuous covariance is the discrete covariance
divided by `(N+1)^2` **plus the expectation of this conditional term**.
It cannot be obtained by simply rescaling discrete C/W covariance. These
formulas require no resampling of the continuous priorities.

## The dipole and what its connected remainder means

Let `mu_C=E C_tau`, `mu_W=E W_tau`, and
`kappa_CW=Cov(C_tau,W_tau)` in the original full target population. Then

\[
 \boxed{I_1-\mu_C I_0=-\mu_W H_{C_\tau}-\delta\kappa_{CW}.}
\]

For a fixed reference a, including the production p_ref, the exact form is

\[
 I_1-aI_0=-\mu_W H_{C_\tau}-\delta\kappa_{CW}
                    -(\mu_C-a)H_{W_\tau}.
\]

The last term is required when lifetime mass changes; moving the dipole
origin from mu_C to p_ref does not eliminate it.

Since `C_tau W_tau=(tau2^2-tau1^2)/2`,

\[
 \kappa_{CW}=\tfrac12[\operatorname{Var}(\tau_2)
                              -\operatorname{Var}(\tau_1)].
\]

Thus the dipole remainder is a precise **endpoint variance-imbalance
response**. The whole mean E curve is determined by endpoint marginals;
it cannot distinguish different admissible K1/K2 joint laws sharing those
marginals. At second order the missing joint direction is

\[
 \operatorname{Cov}(\tau_1,\tau_2)
       =\operatorname{Var}(C_\tau)-\tfrac14\operatorname{Var}(W_\tau),
\]

whereas `Var(C_tau)+Var(W_tau)/4` and kappa_CW are both fixed by the endpoint
marginals. The saved joint K1,K2 rows make this extra coordinate available;
it should not be inferred from the E curve alone.

## Parameter-free rigid-translation predictions and their precise scope

Define the narrow comparison mechanism by a constant infinitesimal shift
`C_tau -> C_tau+v t`, with W_tau unchanged, of the **whole target output
law**. This is a moment-level rigid-translation ansatz, without boundary
clipping; it is not a claim that every such shift can be implemented by
the frozen next-label policy. It implies

```
H_Ctau = v,
H_Wtau = H_Wtau^2 = 0,
delta Var(C_tau) = delta Var(W_tau) = delta Cov(C_tau,W_tau) = 0,
I0 = 0,     I1 = -mu_W H_Ctau,
I2 = -2 E(C_tau W_tau) H_Ctau.
```

In the second line `H_Wtau^2` means the response of `E[W_tau^2]`, not the
square of H_Wtau. More generally, write
`Q_r=E[(tau2^(r+1)-tau1^(r+1))/(r+1)]`; Q0=mu_W. For m>=1 the ansatz gives
`I_m=-m v Q_(m-1)`, hence the velocity-free hierarchy

\[
 \boxed{\mu_W I_m-m Q_{m-1}I_1=0,\quad m\ge1,\qquad I_0=0.}
\]

For m=2 this is the immediately available test
`mu_W I2-2 E(C_tau W_tau) I1=0`. No center-shift parameter is fitted.
The separate prediction `I1+mu_W H_Ctau=0` compares to the independently
named center response in the *same* source, with its common covariance.
Satisfying a few moments would not prove rigid translation.

A path-dependent center shift `C_tau -> C_tau+t v(path)`, W_tau fixed,
is already more general. It instead gives

```
H_Ctau = E v,
delta Var(C_tau) = 2 Cov(C_tau,v),
delta Cov(C_tau,W_tau) = Cov(v,W_tau),
I1 = -E(v W_tau),       I2 = -2 E(v C_tau W_tau).
```

It can fail the rigid predictions while every trajectory's lifetime stays
fixed. In particular, moving only active prefixes and leaving the others
unchanged is path-dependent at the full-population level. A general thermal
remap of both endpoints can change W_tau as well. Rejection of the narrow
constant-shift null therefore does not reject all center motions or all
thermal remaps, and a mean center response alone never specifies the shape.

## What I2 adds: the normalized plateau's center and width

When mu_W>0, normalize the rank-one plateau measure:

\[
 d\nu(p)={P(\tau_1\le p<\tau_2)\over\mu_W}\,dp.
\]

Its mean eta and variance V are

\[
 \eta={E(C_\tau W_\tau)\over\mu_W},\qquad
 V={E(C_\tau^2W_\tau+W_\tau^3/12)\over\mu_W}-\eta^2
   =\operatorname{Var}_{W\text{-weighted}}(C_\tau)
      +{E W_\tau^3\over12\mu_W}.
\]

The exact moment responses are

\[
 \delta\eta=-{I_1-\eta I_0\over\mu_W},\qquad
 \boxed{\delta V=-{I_2-2\eta I_1+(\eta^2-V)I_0\over\mu_W}.}
\]

These separate plateau mass, lifetime-weighted location and width with no
fitted kernel or adjustable scale. The width has two identifiable pieces
only after supplying the joint lifetime-cubic readout above; I2 alone mixes
between-window center spread with within-window width. A rigid translation
has delta eta=v and delta V=0. These identities are available for the new
thermal histograms but are **not asserted here to hold empirically under
any proposed null**.

## Pairing, covariance and scientific handoff

Define all raw moments, population means, variances, covariances and their
responses separately for f and s first, with the original full-prefix
denominator. Then apply `S=(f+s)/2` and
`D=(f-s)/delta_cos4` to each complete identity. For example the S center
term is `(mu_W,f H_C,f+mu_W,s H_C,s)/2`, not the product of paired means.
For any scalar orientation arrays a,b,

```
S(ab) = S(a)S(b) + (delta_cos4^2/4) D(a)D(b),
D(ab) = S(a)D(b) + D(a)S(b).
```

Covariance of paired C/W variables would additionally contain cross-geometry
terms and is not the paired contrast of the within-orientation covariances.
Pool each population's moments before forming products; delete an original
batch and re-form the entire expression for joint uncertainty. Do not
average products formed separately in the20 batches or treat these moment
readouts as independent evidence.

Scientific card: the new exact hierarchy tests whether the E dipole is a
rigid plateau translation, a lifetime-mass change, or a changing normalized
plateau shape. It also separates endpoint variance imbalance from genuine
additional two-birth joint information. Sources are the common-label policy
cfaae36c/ffb70969 and saved histogram4db356e1; the root's new moment readout
9059776d is the intended consumer. This note is theory only, with no MC,
DP, raw replay, numerical verification or claim about current null outcomes.
