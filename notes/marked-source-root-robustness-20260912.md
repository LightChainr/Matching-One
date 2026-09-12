# Positive marks: robust limiting roots, but not automatically unique finite roots

2026-09-12. Corollary of the arbitrary-period odds theorem. This addresses
actual microscopic local-source and homology-source variants without replacing
#275's frozen source or asserting a continuum identification.

## 1. General bounded-oscillation theorem

Let W_Lambda(omega)>0 be a fixed, p-independent weight on configurations, and

    Omega_Lambda = max_omega log W_Lambda - min_omega log W_Lambda.

Keep the physical normalizer:

    P_j^W(p)=E_p[W 1{r=j}]/E_p[W],    M^W=P_2^W-P_0^W.

Then exactly

    exp(-Omega) P_2/P_0 <= P_2^W/P_0^W <= exp(Omega) P_2/P_0.    (1)

The denominator cancels in this RATIO, not from the probability law. The proof
is the elementary min/max bound on W in each of the two event sums.

**Theorem B.** On any integer-period sequence with ell->infinity, if

    Omega_Lambda = o(N/ell),                                    (2)

then every zero of M^W tends to p_c(NN). There is at least one zero, since
M^W is continuous and its endpoint values are -1 and +1. Uniqueness is NOT
claimed for general W.

Proof. Below p_c the unweighted log odds are <=-kappa N/ell; above p_c they
are >=kappa* N/ell. The oscillation in (2) is smaller than either fixed-p
margin. The unweighted odds are increasing in p, so the same signs hold on
all p outside [p_c-epsilon,p_c+epsilon], not just at its endpoints. Consequently
ALL weighted zeros are trapped there. Also N/ell>=sqrt(3)ell/2->infinity.

This is a statement about a source family. It does not convert an unseen
source into new evidence or authorize replacing an old score's source.

## 2. Local product fields and intrinsic rank sources

For a local log-odds field

    W(omega)=exp(sum_v eta_v omega_v),

Omega=sum_v |eta_v|, and the normalized measure is still a product law with
site probabilities

    p_v(p)=p exp(eta_v)/(1-p+p exp(eta_v)).

Each p_v is strictly increasing in p. Hence M^W is strictly increasing, and
Theorem B supplies a UNIQUE consistent root whenever

    sum_v |eta_v|=o(N/ell).                                     (3)

A fixed number of bounded local insertions, placed ANYWHERE, obey (3). More
generally k_Lambda=o(N/ell) sites with uniformly bounded fields suffice. This
connects the finite local-source controls to a limiting-root statement, not a
prediction of their finite shift amplitude or harmonic type. The field here is
in log odds; it is not silently equated to a fixed additive p+-epsilon contract.

For W=exp(s(r-1)), Omega=2|s| and the rank odds are multiplied by exp(2s).
The source-balanced root is again unique and consistent if |s|=o(N/ell).
For axis rectangles N/ell=m, recovering the earlier |s|=o(m) statement.

For a bounded nonlocal mark f_Lambda and W=exp(s_Lambda f_Lambda), it is enough
that |s_Lambda| osc(f_Lambda)=o(N/ell). No monotonicity or locality of f is
needed to localize all zeros, but uniqueness does not follow.

## 3. An exact finite warning: a positive fixed mark creates three roots

On the axis 4x4 torus choose one occupied row-column cross A (7 sites, rank 2)
and the complementary 3x3 block B (9 sites, rank 0). Set

    W(omega)=1+10^12 [1{omega=A}+1{omega=B}].

It is positive and independent of p. The unnormalized weighted balance is

    M(p)+10^12 p^7(1-p)^7(1-2p).                               (4)

Its normalizer is 1+10^12[p^7(1-p)^9+p^9(1-p)^7]>0.
Independent physical lifted-graph enumeration reproduces the base Bernstein
coefficients, then Fraction arithmetic gives signs

    p=1/1000: negative;  p=1/4: positive;
    p=3/4: negative;     p=999/1000: positive.

So there are at least THREE distinct finite balance roots. This prevents
silently importing monotonicity from the original product law into an arbitrary
weighted ensemble. The example does not refute Theorem B: the latter localizes
all roots asymptotically under (2), not their finite multiplicity.

Files: scripts/marked_balance_controls.py and its saved exact-rational report.
The 65,536-configuration axis census here is only this additional finite
regression, separate from the oblique HNF controls.
