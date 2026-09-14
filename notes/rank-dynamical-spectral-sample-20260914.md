# Dynamical topological rank and its pivotal spectral sample

Date: 2026-09-14. Exact finite product-space identities, small-size controls, and a dynamical universality conjecture. Addendum to draft #773 / #769 / #776.

## 1. p-biased Fourier spectrum of the rank observable

Fix Bernoulli site parameter p, q=1-p, and

    X=r-1.

For each site define the centered normalized coordinate

    phi_i(omega)=(omega_i-p)/sqrt(pq).

The product basis `phi_S=prod_{i in S} phi_i` is orthonormal, and

    X = sum_S Xhat(S) phi_S,
    Xhat(S)=E[X phi_S].                                   (1)

At a balance root `E X=0`, Parseval gives

    Var X = sum_{S!=empty} Xhat(S)^2.                    (2)

Define the spectral-sample law

    P(S_spec=S)=Xhat(S)^2/Var X,  S!=empty,              (3)

and its size

    K_spec=|S_spec|.

This is an exact finite distribution attached to the topological rank observable.

## 2. Heat-bath dynamical percolation turns the spectrum into an autocorrelation law

Let every site be independently refreshed from Bernoulli(p) at rate one. The p-biased Fourier mode on S has semigroup eigenvalue `-|S|`. Therefore

    C_X(t)=Cov[X(omega_0),X(omega_t)]
          =sum_{S!=empty} Xhat(S)^2 e^(-t|S|),            (4)

or, normalized,

    C_X(t)/Var X = E[e^(-t K_spec)].                      (5)

Thus the entire stationary dynamical autocorrelation of rank is the Laplace transform of one positive integer-valued spectral sample.

## 3. First spectral moment = dynamic topological clock

Discrete differentiation of the p-biased Fourier expansion gives

    sum_S |S| Xhat(S)^2
      = pq sum_v E[(Delta_v X)^2].                        (6)

Hence

    E K_spec
      = pq N E[(Delta_v X)^2] / Var X                   (7)

on a transitive N-site torus.

Using the local boundary fluxes from the two-threshold calculus,

    E[(Delta_v X)^2]=alpha+4 beta+gamma,                 (8)

where beta is direct rank `0->2` birth.

At a balance root `P0=P2=a`, the conditional occupation difference `g1=E[K|2]-E[K|0]` satisfies

    pq M'_p = a g1,
    Var X=2a.

Therefore the exact ratio

    2 E K_spec / g1 - 1
       = 2 beta/[alpha+2 beta+gamma]
       = omega_2.                                        (9)

The direct-jump fraction can therefore be read either from local rank flips or from the global Fourier/dynamical spectrum. This is a strong independent check for #769.

## 4. All spectral factorial moments are squared multi-pivotal energies

For distinct sites `i_1,...,i_m`, repeated discrete differentiation of (1) gives

    E[(Delta_{i_1}...Delta_{i_m}X)^2]
      = (pq)^(-m)
        sum_{S superset {i_1,...,i_m}} Xhat(S)^2.         (10)

Summing over ordered distinct m-tuples yields

    E[(K_spec)_m]
      = (pq)^m / Var X
        sum_{i_1,...,i_m distinct}
          E[(Delta_{i_1}...Delta_{i_m}X)^2].             (11)

This hierarchy is important because it uses **squares**, not signed sums. It therefore retains the absolute multi-pivotal geometry which may cancel out of ordinary thermal derivatives.

For m=2, (11) is the L2 analogue of #769's pivotal-pair atlas. Higher m gives a canonical extension without inventing a new family of observables.

## 5. Exact matching parity of the spectrum

Let G and Ghat be a primal/matching pair at complementary parameters p and q. Under configuration complement,

    X_G(omega)=-X_Ghat(omega^c),
    phi_i^(q)(1-omega_i)=-phi_i^(p)(omega_i).

Therefore

    Xhat_Ghat,q(S)=(-1)^(|S|+1) Xhat_G,p(S),              (12)

and hence

    |Xhat_Ghat,q(S)|^2=|Xhat_G,p(S)|^2.                  (13)

The **entire spectral-sample distribution is matching-even**.

For a self-matching model at p=1/2, (12) applies within one model and forces

    Xhat(S)=0 for every even |S|.                         (14)

The triangular L=3,4 controls satisfy this exactly: all even Fourier levels vanish to arithmetic zero.

This is the dynamical counterpart of the exact matching-even persistence-gap law.

## 6. Critical scaling and time normalisation

A four-arm pivotal density of order `L^-5/4` gives

    E K_spec = O(L^(3/4)).                               (15)

Thus the topological rank changes on heat-bath time scale `L^-3/4`, the standard dynamical-percolation scale.

However the absolute constant in that time scale can depend on microscopic normalisation. A cleaner observable is the mean-normalized spectral size

    Y_L=K_spec/E K_spec.                                  (16)

Then

    C_X(tau/E K_spec)/Var X
       = E[e^(-tau Y_L)].                                (17)

**Conjecture DS1.** At fixed torus modulus, `Y_L` converges to a lattice-independent law in the percolation universality class. Equivalently, the mean-clock-normalized autocorrelation (17) converges to a universal function.

For triangular site this should be approachable from the Garban--Pete--Schramm dynamical scaling limit once torus rank is shown to be an appropriate continuity observable.

## 7. Small exact controls

Direct p-biased Fourier transforms of all configurations give:

| lattice | L | E K_spec | E K_spec/L^(3/4) | CV^2(K_spec) |
|---|---:|---:|---:|---:|
| square |3|1.62559858|0.71313601|0.38762459|
| square |4|2.00362886|0.70838978|0.53000273|
| triangular |3|1.68072289|0.73731857|0.47785871|
| triangular |4|2.03530475|0.71958890|0.56021794|

The absolute mean amplitude is not claimed universal. More relevant is the mean-normalized Laplace function from (17):

| tau | square L3 | square L4 | triangular L3 | triangular L4 |
|---:|---:|---:|---:|---:|
|0.5|0.631779|0.640893|0.637511|0.642876|
|1|0.423016|0.443515|0.435346|0.448147|
|2|0.207493|0.238291|0.223579|0.245632|
|4|0.056991|0.081198|0.066131|0.087007|

The cross-lattice agreement is suggestive but not an asymptotic result.

The exact direct-jump identity (9) gives at the respective balance points

    omega_2(square L3)    ~0.10397,
    omega_2(square L4)    ~0.07432,
    omega_2(triangular L3)=2/29 ~0.06897,
    omega_2(triangular L4)~0.03886.

This is a clean finite-size measure of the subleading shared-boundary fusion channel.

## 8. Six-arm versus eight-arm placement

If beta is a six-arm shared-boundary event, then

    omega_2 ~ L^-5/3.                                    (18)

If direct simultaneous birth itself requires eight arms,

    omega_2 ~ L^-4.                                      (19)

Because the entire dynamical spectral law is matching-even, either behaviour belongs to an even process channel. It does not by itself decide the matching-odd one-point root correction.

A useful split is therefore:

    dynamic spectral sample / persistence gap  -> absolute/even fusion geometry,
    alpha-gamma / odd rank-curve correction    -> matching-odd selection rule. (20)

## 9. Research interfaces

**#769:** compare the real-space squared pair-influence sum with the second factorial spectral moment from (11). This independently checks the pivotal-pair implementation.

**#776:** self-matching triangular site has the exact even-level vanishing (14). A larger exact/transfer control can test the normalized spectral-size law before any square-site universality claim.

**#768:** spectral data can establish the presence/scaling of six-arm absolute geometry, but matching parity (13) means it cannot alone kill the odd 8-arm hypothesis.

A future dedicated dynamical task should target the torus continuity interface to rigorous dynamical percolation and the universal mean-normalized autocorrelation, not a generic time-series simulation.

## 10. Boundaries

- The p-biased Fourier identities are exact finite product-space facts.
- The L^(3/4) clock and cross-lattice convergence are universality/critical-scaling statements.
- Mean-clock normalization removes one dynamical metric factor but does not prove the whole spectral-size law is universal.
