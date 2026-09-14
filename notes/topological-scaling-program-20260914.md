# Topological scaling program: integrated state after the 2026-09-14 analysis rounds

Date: 2026-09-14. Intended branch: `analysis/topological-scaling-program-20260914`, based on `main@d31fa5fe9584b77595e4a78b9560aa2ff77cf0b3`.

**Scope and epistemic status.** This note integrates the owner's delegated analysis rounds into one research program. It cites unmerged research assets by pinned PR/commit and does not promote them to merged or publication-certified status. Exact finite identities are labelled separately from conditional asymptotic theorems and conjectures. No new Monte Carlo, GPU run, numerical `p_c`, STATUS claim, or original-U identification is introduced here.

The main working assets remain:

- PR #739, `907a9d94ee0c3f20b00188647872a067abffd275`: arbitrary-period geometric-balance manuscript plus complete-component tools;
- PR #737, `7b459c4809baee9dfc2b910091e00a3a8d508d1c`: spatial-source Hessian, rank-source zeros and conditional odds;
- merged #705 / commit `fc19cc748527249f0ce1c69d7cd2a88a4879427b`: exact width-two torus/cylinder bridge;
- open tasks #768, #769, #770 created from this program.

---

## 1. One architecture instead of several unrelated projects

The present program has five layers.

### Layer G — geometry

For honest integer-period square-site tori, the balance root compares rank-0 and rank-2 rare sectors. The current author proof in #739 gives root convergence when the genuine shortest period tends to infinity, with no aspect/shear condition, while the full birth-time law has the sharper geometric boundary `log N / ell -> 0`. The root and the full law are therefore distinct asymptotic objects.

### Layer C — complete winding components

On a fixed-width cylinder, complete winding components have a unique Palm anchor. PR #739's tagged-span resolvent gives an all-height phase-type law, direct two-fugacity activity `p^K(1-p)^B`, exact moments, and exact Palm sampling. The fixed-subcritical regime is naturally a dilute gas of complete winding cuts; complementary white components fill the intervals between cuts.

The earlier analysis rounds derived the deterministic black-gap/white-span sandwich and, conditional on the existing Poisson-window input, the limiting white-gap laws and their moments. These should remain a separate fixed-subcritical package from the near-critical program below.

### Layer P — topological pivotal geometry

For `X=r-1`, finite Bernoulli multilinearity gives

    M'(p)  = sum_v E[Delta_v X],
    M''(p) = sum_{u != v} E[Delta_u Delta_v X].

The first derivative is a weighted topological pivotal measure; the second derivative is a signed two-pivotal interaction. Exact L=3,4 controls show that rank jumps of two occur and that `M''` hides large positive/negative cancellation. The real-space pair kernel and PR #737's spatial-source Hessian are the same object in two representations:

    J_L(d) = E[Delta_0 Delta_d X] = M_{0d},
    S_piv(k) = sum_d J_L(d) exp(i k.d).

Hence real space is the natural place to classify arm fusion, while Fourier space is the natural place to compare source sectors. Issue #769 owns the exact L=5 pair atlas.

### Layer T — matching-odd continuum selection

Configuration-level digital Alexander duality makes the matching function an antisymmetric pair projection rather than a generic critical observable. If `Ghat` is the matching complement,

    r_G(omega) + r_Ghat(omega^c) = 2,
    M_G(L,p) = -M_Ghat(L,1-p),
    M_G^(n)(L,p) = (-1)^(n+1) M_Ghat^(n)(L,1-p).

Therefore the correct correction-to-scaling question is not "what is the smallest correction exponent in percolation?" but "what is the first correction allowed in the matching-odd topological sector?"

The strong candidate is the 8-arm/8-leg field. With the standard polychromatic arm exponents

    alpha_4 = 5/4,
    alpha_8 = 21/4,

one has

    alpha_8 - alpha_4 = 4.

Thus the candidate mechanism is

    M_L(p_c) ~ L^(2-alpha_8) = L^(-13/4),
    M'_L(p_c) ~ L^(2-alpha_4) = L^(3/4),
    p_L^* - p_c ~ L^-4.

If the 8-arm bulk channel is spinless then `h=hbar=21/8`; level-(n,n) descendants naturally produce a `4,6,8,...` root-correction tower. This is not accepted until lower matching-odd channels, especially the 6-arm channel, are excluded. Issue #768 is the selection-rule audit and has an explicit 6-arm kill test.

### Layer N — thermal normal coordinate and disorder/source calibration

Raw `p` need not be the RG thermal normal coordinate. Writing

    t(p)=a1 epsilon + a2 epsilon^2 + ...,
    u=L^(3/4) t(p),

with an odd leading scaling function explains a leading `M'' ~ L^(3/4)` without introducing a new even primary. The finite diagnostic `M''/(2M')` estimates `a2/a1` only if the normal-form hypothesis is correct; it must be calibrated independently.

Issue #770 uses the subcritical mass coordinate

    t_mass(p) = - kappa(p)^(3/4)

and the amplitude-free curvature

    c2_mass = 1/2 [kappa''/kappa' - (1/4) kappa'/kappa]

as the independent test. The same issue records the exact-finite Harris positive control obtained from PR #737's quenched-root Hessian.

---

## 2. New exact result: the rank-source sector is completely solvable

This section is independent of the 8-arm conjecture.

Put

    X = r-1 in {-1,0,1},
    P_j = P(r=j).

The three exact rank projectors are

    1_{r=2} = (X^2+X)/2,
    1_{r=1} = 1-X^2,
    1_{r=0} = (X^2-X)/2.                                  (2.1)

Consequently any rank-conditioned observable can be reconstructed from ordinary unconditioned moments with `X` and `X^2`; no state-space identity beyond the three-valued rank is required.

At a balance root `P_2=P_0=a`,

    E X = 0,
    E X^2 = 2a,

and for any integrable observable A,

    E[A|2] = E[(X^2+X)A] / E[X^2],
    E[A|0] = E[(X^2-X)A] / E[X^2],
    Delta E[A] = 2 E[XA] / E[X^2].                        (2.2)

Since `E X=0`, the numerator in the last formula is also `2 Cov(X,A)`.

### 2.1 Rank-source equation of state

Couple a finite topological source `s` to `X`:

    Z_X(s)=E exp(sX)=P_0 exp(-s)+P_1+P_2 exp(s).

At a balance root this becomes

    Z_X(s)=1+2a(cosh s-1),                                (2.3)
    m_X(s)=d/ds log Z_X(s)
          = 2a sinh s / [1+2a(cosh s-1)].                 (2.4)

All odd source cumulants vanish exactly. The first nonzero ones are

    kappa_2 = 2a,
    kappa_4 = 2a(1-6a),
    kappa_6 = 2a(1-30a+120a^2),
    kappa_8 = 2a(1-126a+1680a^2-5040a^3).                (2.5)

Thus the entire finite rank-source law at balance is controlled by the single scalar `a=P_0=P_2`.

### 2.2 Exact source zeros at balance

Zeros of (2.3) satisfy

    cosh s = 1 - 1/(2a).                                  (2.6)

Hence:

- if `a>1/4`, the nearest zeros are purely imaginary,

      s = +/- i theta + 2 pi i n,
      theta = acos(1-1/(2a));

- if `a=1/4`, `Z_X(s)=cosh^2(s/2)` and the zeros at `(2n+1)i pi` are double;
- if `a<1/4`, the nearest zeros have imaginary part `pi` and nonzero real parts.

This sharpens, at balance, PR #737's general zero-free strip `|Im s|<pi/2`. A double static source zero is still not an operator Jordan diagnostic.

### 2.3 The source-tilted balance root is an exact inverse conditional-odds curve

Use logit thermal coordinate `z=log(p/(1-p))` and unnormalised rank-sector sums `Z_j(z)`. Define

    g(z)=log[Z_2(z)/Z_0(z)].

After adding the topological source `s`, the sourced balance condition is exactly

    g(z_*(s)) + 2s = 0.                                   (2.7)

The derivatives of `g` are conditional occupation cumulant differences:

    g^(n)(z) = kappa_n(K|r=2)-kappa_n(K|r=0) =: g_n.      (2.8)

At `s=0`,

    z_*'  = -2/g_1,
    z_*'' = -4 g_2/g_1^3,
    z_*'''= 8(g_1 g_3-3g_2^2)/g_1^5.                     (2.9)

So a finite topological-source root curve is an exact generator for rank-conditioned thermal cumulant differences. This is complementary to PR #737's conditional-odds integration: the latter reconstructs `g(z)`, while (2.7) uses `s` as the canonical coordinate conjugate to that odds.

### 2.4 Exact unconditioned formula for the spatial-source Hessian

Let `h_i` be a mean-zero spatial logit source and `H=sum_i h_i n_i`. Translation invariance within each rank sector gives `E[H|r=j]=0`. PR #737 gives

    z_*''(h) = - Delta Var(H) / Delta E[K].                (2.10)

Using the rank projectors (2.1), at balance

    Delta Var(H) = E[X H^2]/a,
    Delta E[K]   = E[X K]/a.

Therefore the same root Hessian has the exact unconditioned representation

    z_*''(h) = - E[X H^2] / E[X K].                       (2.11)

This directly links three previously separate objects:

- the ambient-topology source `X` used in the graded Q->1 analysis;
- the rank-conditioned structure factor of PR #737;
- an ordinary full-ensemble mixed moment requiring no explicit rare-sector conditioning.

It does **not** make `X` a continuum field identity; it is a finite algebraic projector relation.

---

## 3. Continuum positive control from Pinson wrapping probabilities

For the critical continuum square torus, Pinson's topology result gives equality of the no-wrapping and cross-topology probabilities. Newman--Ziff quote for the square torus

    R_e = P(any wrapping) = 0.690473724570168677...,

hence

    a_* = P(no wrap)=P(cross)=1-R_e
        = 0.30952627542983132276969363335....              (3.1)

At large square tori, if the finite matching root converges to the critical continuum topology law, the full rank distribution at the root must approach

    (P_0,P_1,P_2)
      -> (a_*, 1-2a_*, a_*).                              (3.2)

Therefore the rank-source continuum family is explicitly

    Z_*(s)=1+2a_*(cosh s-1).                               (3.3)

It predicts

    Var(X)      -> 0.61905255085966264554...,
    E[X^4]/E[X^2]^2 -> 1/(2a_*)
                     = 1.61537174608411717427...,
    excess kurtosis -> -1.38462825391588282573...,         (3.4)

and a nearest source zero

    theta_* = 2.23365382365240413977...,
    theta_*/pi = 0.71099409438078560159....                (3.5)

These are not new independent universal constants; they are exact transforms of the known topological wrapping probability. Their value is as **typed positive controls** for source-response, zero-finding and map-resolved torus machinery.

The included exact finite enumeration gives:

| L | balance root p_L^* | a_L=P0=P2 | theta_L/pi |
|---|---:|---:|---:|
| 3 | 0.58651145511267563565 | 0.32912802663079765386 | 0.673757239449169 |
| 4 | 0.59067211233102829690 | 0.32274409827319622603 | 0.685073198471710 |
| continuum square torus | p_c | 0.30952627542983132277 | 0.710994094380786 |

Two finite sizes do not establish a convergence law. Issue #769 has been asked to report the essentially free L=5 continuation while constructing its pivotal-pair atlas.

### 3.1 Interface to map-resolved torus tomography

For arbitrary torus modulus/twist, the same construction uses the universal cross/no-wrap probability `a_*(tau)`. Thus

    Z_tau(s)=1+2 a_*(tau)(cosh s-1)                         (3.6)

is a known topological-source family with explicit map/connectivity semantics. It is a useful positive control for #585's map-resolved modular-covariance machinery. It is not an original-U candidate.

---

## 4. The near-critical program to test next

The preceding exact finite algebra suggests three higher-level scaling conjectures.

### Conjecture A — near-critical winding intensity

There is a thermal metric `a_t` and a universal intensity crossover `I(lambda)` such that

    w nu_w(p_c + lambda/(a_t w^(3/4))) -> I(lambda).       (4.1)

On the subcritical tail, consistency with `kappa(p)~C(p_c-p)^(4/3)` requires

    -log I(lambda) ~ C_- |lambda|^(4/3), lambda -> -infty. (4.2)

For aspect parameter `r=N/w^2`, the first-winding center solves roughly `r I(lambda)~1`, recovering the #739 full-law boundary `log N / w -> 0` as the condition that the center returns to the ordinary near-critical window.

### Conjecture B — complete-cluster Palm score Ward cancellation

PR #739 gives exactly, at fixed width,

    d log nu_w / dp
      = E_Palm[K/p - B/(1-p)].                             (4.3)

If (4.1) holds and the critical complete winding Palm mass and distinct vacant boundary both scale as `w^(91/48)`, then the two extensive geometric scores must cancel to leave only `w^(3/4)`:

    E[K]/p_c - E[B]/(1-p_c) = O(w^(3/4)).                 (4.4)

In particular,

    E[B]/E[K] -> (1-p_c)/p_c.                             (4.5)

This is a high-value falsifiable prediction because the two-fugacity complete-component operator already exposes K and B derivatives. It requires no continuum field naming.

### Conjecture C — matching-odd 8-arm defect

The first nonzero scalar correction in the matching-odd rank sector is the spinless 8-arm/8-leg family `x=21/4`, while all lower identity/thermal/6-arm/generic irrelevant contributions are pair-even or map-sector-forbidden. This produces the observed candidate root tower `L^-4,L^-6,L^-8,...`. Issue #768 decides whether this survives a real selection-rule audit.

---

## 5. Original-U boundary and how this program should be used

Nothing above solves #275. The original-U contract still requires two named physical candidates to provide same-source raw forward columns, the physical normalizer, moving-root counterterm, allowed amplitude class and nuisance transfer.

The present program adds **calibration gates**, not replacement observables:

1. declare the source type (`thermal`, `B_even`, ambient-topology `X`, spatial field, etc.);
2. satisfy exact matching-pair parity;
3. separate raw-p coordinate curvature from an independently calibrated thermal normal coordinate;
4. pass at least one source-explicit topological/Fourier positive control;
5. only then enter the frozen original-U forward-map comparison.

The rank-source family (Section 3) is particularly useful because its continuum topology semantics are known exactly; failure to reproduce it should block interpretation of a more complicated map-resolved candidate.

---

## 6. Active external tasks and stop rules

- **#768 — matching-odd 8-arm sector audit.** Main kill test: an allowed nonzero 6-arm scalar invalidates the claimed 8-arm leading mechanism.
- **#769 — exact pivotal-pair atlas.** Main target: L=5 exact signed/absolute pair influence, distance/orbit structure, Fourier agreement with spatial Hessian, and the rank-source L=5 positive control.
- **#770 — mass thermal coordinate and Harris control.** Main target: independent `c2_mass`, matching-pair parity, and uncorrelated/long-range disorder scaling.

Stop rules:

- Do not infer a continuum field from one exponent or one Fourier sign.
- Do not infer a multi-arm exponent from three finite sizes.
- Do not call raw `M''` a two-pivotal amplitude without quantifying its signed cancellation.
- Do not replace the physical source of #275 by the rank source `X`.
- If #768 finds a lower allowed matching-odd channel, correct the 8-arm story instead of relabelling it.
- If #770 cannot control `gamma_w-kappa` close enough to `p_c`, keep the thermal-coordinate comparison explicitly blocked on mass locality.

---

## 7. Sources used in this integrated note

Repository:

- PR #739 / `907a9d94`: `docs/manuscripts/geometric-balance/*`, especially `tagged-span-resolvent.md` and `poisson-birth-windows.md`.
- PR #737 / `7b459c48`: `spatial-source-hessian-20260913.md`, `rank-source-zero-free-strip-20260913.md`, `conditional-odds-integration-20260913.md`, `annealed-quenched-root-20260913.md`.
- #705 / `fc19cc74`: exact width-two cylinder-sector bridge.
- PR #746: exact graded split of duality-even `B_even` and ambient topology `X`.

External primary anchors:

- Smirnov--Werner, critical arm exponents, arXiv:math/0109120.
- Garban--Pete--Schramm, pivotal and near-critical measures, arXiv:1008.1378 and related work.
- Pinson / Newman--Ziff torus wrapping probabilities; Newman--Ziff quote `R_e=0.690473724570...` for the square torus.
- Jacobsen 2015, arXiv:1507.03027, open/closed graph-polynomial sectors and the observed fast correction sequence.
- Mertens--Ziff 2016, arXiv:1603.07289, matching function scaling.
- Mertens 2022, arXiv:2109.12102, square-site threshold estimator and correction tower.

The square-lattice use of exact percolation exponents remains a universality input unless separately proved for that lattice/model. The finite algebra and rank-source formulas in Sections 2--3 do not depend on that universality assumption.
