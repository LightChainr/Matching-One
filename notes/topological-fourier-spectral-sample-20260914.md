# Topological Fourier spectral sample: pivotal geometry, source Hessians and noise sensitivity

Date: 2026-09-14. Exact finite algebra and tiny controls; asymptotic spectral-sample statements are explicitly conditional. Addendum to draft #773 / #737 / #769 / #784.

## 1. p-biased Fourier expansion of the rank observable

Let

    X(omega)=r_black(omega)-1 in {-1,0,1}

under product Bernoulli(p), q=1-p. Use the standard p-biased orthonormal basis

    phi_i(omega_i)=(omega_i-p)/sqrt(pq),
    phi_S=prod_(i in S) phi_i.

Write

    X-E X = sum_(S != empty) Xhat_p(S) phi_S.              (1)

Conditioning on all other bits gives exactly

    Xhat_p({i}) = sqrt(pq) E_p[Delta_i X].                 (2)

For i != j,

    Xhat_p({i,j}) = pq E_p[Delta_i Delta_j X].             (3)

No Boolean-valued assumption is used; Delta_i X may be 0,1,2 and the mixed difference may have either sign.

Consequently

    M'(p) = (pq)^(-1/2) sum_i Xhat_p({i}),                 (4)

and PR #737's spatial Hessian entries are

    M_ij = (pq)^(-1) Xhat_p({i,j}), i != j.                (5)

Thus the thermal derivative is the coherent level-1 Fourier sum, while the spatial-source Hessian is the signed level-2 Fourier tensor.

## 2. Exact matching/complement relation at every Fourier level

Configuration complement gives

    X_G(omega) = - X_Ghat(1-omega).

Complement maps Bernoulli(p) to Bernoulli(q) and

    phi_i^p(omega_i) = - phi_i^q(1-omega_i).

Therefore, for every S,

    Xhat_(G,p)(S) = (-1)^(|S|+1) Xhat_(Ghat,q)(S).        (6)

This is the all-level version of the derivative parity already used elsewhere:

- odd Fourier levels are matching-even between the paired models;
- even Fourier levels are matching-odd.

For a self-matching model at p=1/2,

    X(1-omega)=-X(omega)

inside the same model. Hence

    boxed: Xhat(S)=0 for every even |S|.                   (7)

This exact theorem is much stronger than the special statements M''=0 or g_2=0.

## 3. The topological spectral sample

Define a random nonempty subset S_X by

    P(S_X=S) = Xhat_p(S)^2 / Var_p(X).                    (8)

For standard independent resampling noise with correlation rho=1-epsilon,

    Cov[X(omega),X(omega_epsilon)]/Var(X)
      = E[rho^|S_X|].                                     (9)

Discrete differentiation and Parseval give

    sum_S |S| Xhat(S)^2
      = pq sum_i E[(Delta_i X)^2],                        (10)

and

    sum_S |S|(|S|-1) Xhat(S)^2
      = (pq)^2 sum_(i != j) E[(Delta_i Delta_j X)^2].     (11)

Therefore

    E |S_X|
      = [pq/Var(X)] sum_i E[(Delta_i X)^2].               (12)

A direct rank-0 to rank-2 insertion enters (12) with weight four because `(Delta_i X)^2=4`. This differs from an ordinary Boolean crossing influence.

Equation (11) is useful for #769: the three pair observables

    sum J_ij,
    sum |J_ij|,
    sum J_ij^2

with `J_ij=E Delta_i Delta_j X` have distinct meanings. The first is a signed coherent derivative, the second an absolute mixed influence, and the third contributes to Fourier/spectral energy. They should not be assigned the same arm exponent by default.

## 4. Translation-invariant source Hessian = spatial level-2 Fourier transform

On a transitive torus write

    J(d)=E[Delta_0 Delta_d X].

Then

    Xhat_2(d)=pq J(d).                                    (13)

PR #737 diagonalizes the translation-convolutional Hessian by spatial momentum. Its eigenvalue at wavevector k is

    lambda(k)=sum_d J(d) exp(i k.d)
             =(pq)^(-1) sum_d Xhat_2(d) exp(i k.d).       (14)

So real-space pivotal-pair fusion and Fourier source response are exactly the same level-2 tensor in two bases.

## 5. Exact L=3,4 controls

Independent lifted-homology enumeration gives the following spectral-sample size moments at the finite square balance root and at triangular p=1/2.

| model | L | E|S_X| | E|S_X| / L^(3/4) | Var(|S_X|) |
|---|---:|---:|---:|---:|
| square NN |3|1.62559858|0.71313601|1.02432540|
| square NN |4|2.00362886|0.70838978|2.12771111|
| triangular self-matching |3|1.68072289|0.73731857|1.34986936|
| triangular self-matching |4|2.03530475|0.71958890|2.32068344|

Two widths are only controls. The closeness of the scaled first moment is suggestive, not an exponent or amplitude fit.

For triangular self-matching, a full exact Walsh-Hadamard transform was also performed. Every nonconstant even-level numerator is exactly zero.

At L=3, spectral mass by level is approximately

    level 1: 0.712444
    level 3: 0.238328
    level 5: 0.045745
    level 7: 0.003389
    level 9: 0.000094.

At L=4 it is approximately

    1: 0.615471,
    3: 0.278853,
    5: 0.082393,
    7: 0.019589,
    9: 0.003252,
    11: 0.000411,
    13: 0.000030,
    15: 0.000001,

with all even levels exactly absent. The migration away from fixed low levels is qualitatively consistent with a growing critical spectral sample, but these sizes do not certify its scaling law.

## 6. Conditional critical conjecture

For ordinary critical crossing events, the percolation Fourier spectrum and pivotal scaling imply a characteristic spectral size on the L^(3/4) scale. The rank observable is not Boolean and the crossing theorem does not automatically transfer.

The natural conjecture is

    E |S_X| ~ C_spec L^(3/4),                              (15)

with a nondegenerate rescaled spectral-sample law. If true, the dynamical decorrelation/noise scale is

    epsilon_L ~ L^(-3/4).                                 (16)

The equality of this exponent with the thermal-window exponent should be interpreted through the common pivotal measure, not as an identity between thermal perturbation and dynamical resampling.

## 7. Why this helps the 6-arm / 8-arm question

The 8-arm conjecture concerns a matching-odd one-point finite-size correction. The spectral sample measures squared Fourier activity and need not obey the same cancellation.

In particular:

- abundant theta / six-arm direct-jump geometry can contribute strongly to spectral energy;
- it can still cancel out of the matching-odd one-point amplitude;
- self-matching triangular parity removes every even Fourier level exactly, providing a strong control on which cancellation is symmetry-enforced.

Thus #768 should distinguish "a six-arm structure exists in the spectral/pivotal geometry" from "a six-arm matching-odd scalar survives in M(p_c)".

## 8. Literature boundary

Primary methodological anchors: Benjamini--Kalai--Schramm on percolation noise sensitivity; Schramm--Steif on quantitative Fourier-level control; Garban--Pete--Schramm on the critical percolation Fourier spectrum/spectral sample and pivotal scaling. None is treated as already proving (15) for this bounded rank-valued observable.

## 9. Boundaries

- Equations (1)--(14) are exact finite product-measure/Fourier identities.
- The L=3,4 numbers are exact-enumeration/high-precision controls, not scaling estimates.
- Equation (15) is a conjectural transfer of the pivotal/spectral paradigm to X.
- Squared spectral energy, absolute mixed influence and signed source Hessians are distinct observables.
