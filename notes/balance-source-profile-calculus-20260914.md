# Balance-line source profiling: an exact quotient by thermal nuisance

Date: 2026-09-14. Additive continuation of draft PR #773. This note gives finite-system identities for arbitrary exponential sources at a rank balance point. It formalizes the nuisance quotient that was implicit in the spatial-source Hessian work.

## 1. General source and balance root

For the rank sectors `r=0,2`, introduce logit thermal coordinate z and a physical source A:

    Z_j(z,eta) = sum_{r(omega)=j} exp[z K(omega) + eta A(omega)].

Let

    g(z,eta)=log Z_2-log Z_0,

and let `z_*(eta)` be the local balance root `g(z_*(eta),eta)=0`. Write `Delta` for rank-2-conditioned minus rank-0-conditioned expectation at the unperturbed balance root, and put

    D = Delta E[K] >0.                                   (1)

Then

    g_z = D,
    g_eta = Delta E[A],
    g_zz = Delta Var(K),
    g_zeta = Delta Cov(K,A),
    g_etaeta = Delta Var(A).                             (2)

Therefore

    beta_A := Delta E[A]/D,
    z_*'(0) = - beta_A.                                  (3)

The number `beta_A` is the exact first-order thermal component of A as seen by the balance equation.

## 2. Thermal profiling makes the second derivative a residual variance

Define the profiled source

    R_A = A - beta_A K.                                  (4)

By construction

    Delta E[R_A]=0.                                     (5)

Implicit differentiation gives

    z_*''(0)
      = -[Delta Var(A)-2 beta_A Delta Cov(A,K)
           + beta_A^2 Delta Var(K)]/D
      = - Delta Var(R_A)/D.                              (6)

This is the exact finite analogue of profiling a nuisance direction before reading a curvature. No approximation, CFT input or Gaussian assumption is used.

It immediately passes three controls:

- `A=K`: `beta=1`, `R_A=0`, hence a thermal source is exactly absorbed by `z` and `z_*''=0`;
- a zero-mean Fourier site field at a transitive root has `beta=0`, recovering the PR #737 conditional structure-factor formula;
- `A=X=r-1`: it reproduces the topological-source root derivatives obtained from `g(z)+2s=0`.

## 3. Multisource Hessian and the quotient space

For sources `A_a`, define

    beta_a = Delta E[A_a]/D,
    R_a = A_a-beta_a K.                                  (7)

Then the balance-root Hessian is

    partial_a partial_b z_*
      = - Delta Cov(R_a,R_b)/D.                          (8)

Adding any multiple of K to a source leaves `R_a` unchanged. Adding a constant also has no effect because it multiplies both sector partition functions by the same factor. Thus (8) descends exactly to the source quotient

    source space / span{1,K}.                            (9)

This is the correct finite-system statement behind the phrase “thermal nuisance profiling.” It is stronger and cleaner than checking whether a raw source happens to have zero first derivative.

The bilinear form in (8) is signed: `Delta Cov` need not be positive. Indefinite curvature is therefore expected and does not by itself signal a bad source.

## 4. Remove rare-sector conditioning with the rank projector algebra

At balance let

    X=r-1,
    P0=P2=a,
    E X=0,
    E X^2=2a.

The exact projectors are

    1_{r=2}=(X^2+X)/2,
    1_{r=1}=1-X^2,
    1_{r=0}=(X^2-X)/2.                                  (10)

For any observable A,

    Delta E[A] = 2 E[X A]/E[X^2].                        (11)

Hence

    beta_A = E[X A]/E[X K].                              (12)

Because `Delta E[R_A]=0`, the difference of squared sector means drops out and

    Delta Var(R_A)=2 E[X R_A^2]/E[X^2].                  (13)

Combining (6),(11) gives the unconditioned formula

    z_*'(0)  = - E[X A]/E[X K],
    z_*''(0) = - E[X (A-beta_A K)^2]/E[X K].             (14)

Likewise, for two sources

    partial_ab z_*
      = - E[X R_a R_b]/E[X K].                           (15)

Thus every balance-line source response through second order can be computed as an ordinary expectation under the unconditioned balance ensemble with one signed rank weight X. Rare rank-sector conditional sampling is a computational option, not a mathematical necessity.

## 5. Exact odd/even rank ANOVA

At balance the rank sigma-algebra has only three states. The two centered variables

    X,
    Y = X^2-2a                                          (16)

are exactly orthogonal:

    E[X Y]=0.                                           (17)

Together with the constant, `{1,X,Y}` spans every function of rank. Therefore for any square-integrable A,

    E[A | r]
      = E A
        + [E(AX)/E(X^2)] X
        + [E(AY)/E(Y^2)] Y.                              (18)

Consequently the variance explained by conditioning on rank splits exactly into only two nontrivial channels,

    Var(E[A|r])
      = E(AX)^2/E(X^2) + E(AY)^2/E(Y^2).                 (19)

There is no hidden third rank direction. `X` is the matching-odd/topological channel; `Y` is the matching-even rank-shape channel. This gives a precise finite interpretation of the graded split used in #746 and a clean audit target for Rao–Blackwell hierarchies such as #578.

## 6. Information-geometric consequence for K-only conditioning

For A=K, the odd rank-channel squared correlation is

    rho_KX^2 = Cov(K,X)^2/[Var(K) Var(X)].                (20)

Under the standard critical scalings `Delta E K ~ L^(3/4)`, `Var(K)~L^2`, and `Var(X)=O(1)`, this gives

    rho_KX^2 = O(L^(-1/2)).                              (21)

So any linear variance reduction based only on total occupation K becomes asymptotically weak for the topological rank signal. This does not contradict strong finite-size gains. It says that a long-horizon conditional-integration state must retain connectivity/cut geometry if it is to preserve topological information.

## 7. Use as an original-U pre-filter, not an identification

For a proposed microscopic source A intended to feed the original-U contract:

1. compute `beta_A` and remove the exact thermal nuisance `beta_A K`;
2. evaluate the signed quotient Hessian (15) or its represented raw coordinates;
3. only residual directions not annihilated in the quotient can supply source-visible information beyond a root shift.

Passing this filter is not a continuum-field identification. It is an exact finite obstruction against candidates whose purported signal lies entirely in `span{1,K}` or in an already-declared nuisance direction.
