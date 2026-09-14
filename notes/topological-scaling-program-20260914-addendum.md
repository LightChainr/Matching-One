# Topological scaling program — source geometry and near-critical equation-of-state addendum

Date: 2026-09-14. Addendum to `notes/topological-scaling-program-20260914.md` on draft PR #773. The finite identities below are exact. Statements invoking a smooth near-critical scaling function are conditional; square-site use of the critical exponents remains a universality input unless separately proved.

## 1. General balance-line source calculus

Let `z=log(p/(1-p))` and let `A(omega)` be any configuration observable coupled exponentially inside each rank sector,

    Z_j(z,eta)=sum_(r=j) exp[z K(omega)+eta A(omega)],
    g(z,eta)=log Z_2(z,eta)-log Z_0(z,eta).

Let `z_*(eta)` solve the sourced balance equation `g(z_*(eta),eta)=0`. At `eta=0`, write `Delta` for rank-2 minus rank-0 conditional expectation. Exact implicit differentiation gives

    z_*' = - Delta E[A] / Delta E[K],                       (1.1)

and

    z_*''
      = -{ Delta Var(A)
           + 2 Delta Cov(K,A) z_*'
           + Delta Var(K) (z_*')^2 }
        / Delta E[K].                                      (1.2)

At an unsourced balance root `P_2=P_0=a` and `X=r-1`, the projector identity gives

    Delta E[A] = E[X A]/a,
    Delta E[K] = E[X K]/a,

hence

    z_*' = - E[X A]/E[X K].                                (1.3)

This single theorem contains several previous formulas:

- `A=X`: topological chemical potential, `z_*'=-2/Delta E[K]`;
- `A=n_i`: one addressed site, translation gives `z_*'=-1/N` exactly;
- `A=H=sum h_i n_i`, `sum h_i=0`: `z_*'=0`, and (1.2) becomes the spatial-source Hessian of PR #737.

For a vector source `A_a`, the Jacobian is simply

    partial z_*/partial eta_a = -E[X A_a]/E[XK].            (1.4)

Thus first-order balance-root motion only measures overlap with the finite matching-odd rank direction; it is not a continuum field identifier.

## 2. Fisher geometry of the thermal and topological sources

At a homogeneous balance root the two canonical finite scores are

    thermal:     K-EK,
    topological: X.

The full Bernoulli ensemble gives exactly

    Var(K)=N p(1-p),
    Var(X)=2a,
    Cov(K,X)=E[XK]=a Delta E[K]=p(1-p) M_p'(p).             (2.1)

Therefore their squared Fisher correlation is

    rho_KX^2
      = p(1-p) [M_p'(p)]^2 / (2 a N)
      = a [Delta E[K]]^2 / [2N p(1-p)].                    (2.2)

If the standard critical law `M_p'~C L^(3/4)` holds on an `L x L` torus and `a->a_*>0`, then

    rho_KX^2 ~ const L^(-1/2),
    rho_KX   ~ const L^(-1/4).                             (2.3)

So the global thermal density score and the bounded topology score become asymptotically orthogonal in Fisher angle even though a finite topological source shifts the balance root by the critical thermal-window scale.

The source-balanced slope has the exact information-geometric form

    dz_*/ds = -Var(X)/Cov(K,X) = -2/Delta E[K].             (2.4)

Hence under the same scaling,

    L^(3/4) dz_*/ds = O(1).                                (2.5)

The exact L=3,4 controls committed with this addendum are:

| L | `sqrt(L) rho^2` | `L^(3/4) dz_*/ds` |
|---|---:|---:|
| 3 | 1.13262992139208 | -1.54804444358016 |
| 4 | 1.16077682539732 | -1.51656450631935 |

Two sizes are only a structural positive control, not an exponent fit.

This same `L^-1/2` squared-correlation scale is the finite information-geometry counterpart of the Harris disorder-variance exponent `2/nu-d=-1/2` when `nu=4/3`.

## 3. Near-critical rank-source equation of state

For fixed torus modulus/twist `tau`, introduce a correctly normalised thermal coordinate `lambda` and the limiting rank probabilities

    Pi_j(lambda;tau), j=0,1,2.

The proposed two-variable continuum topological partition function is

    Z_tau(lambda,s)
      = Pi_0(lambda;tau)e^(-s)
        + Pi_1(lambda;tau)
        + Pi_2(lambda;tau)e^s.                             (3.1)

Matching/complement exchange, together with an odd thermal normal coordinate, predicts

    Pi_2(lambda;tau)=Pi_0(-lambda;tau),
    Pi_1(lambda;tau)=Pi_1(-lambda;tau),
    Z_tau(lambda,s)=Z_tau(-lambda,-s).                     (3.2)

At `lambda=0`, this reduces to the exactly solvable critical rank-source law in the main note. The horizontal centre of the two source zeros is

    s_b(lambda;tau)=1/2 log[Pi_0/Pi_2],                    (3.3)

an odd function of `lambda`. The finite source-balanced curve is the inverse relation `s=-s_b`.

This separates the leading near-critical crossover from the proposed fast matching-root correction. If

    p_L^*-p_c=O(L^-4),

then the leading thermal coordinate satisfies

    lambda_L^*=O(L^(3/4)L^-4)=O(L^-13/4)->0.               (3.4)

Thus an `L^-4` displacement is too small to be a shifted leading near-critical centre; it must arise from an irrelevant matching-odd correction to (3.2). A nonzero limiting `lambda_L^*` would directly falsify that interpretation.

## 4. Odd log-odds normal form and the conditional-cumulant parity ladder

Define the rank-sector log odds at zero topological source,

    g_L(z)=log[Z_2(z)/Z_0(z)].

Assume a smooth leading near-critical form

    g_L(z)=G(lambda_L(z))+o(1),
    lambda_L(z)=L^(3/4)t(z),
    G(-lambda)=-G(lambda).                                 (4.1)

At the balance point `lambda=0`, sector thermodynamic identities give

    g_L^(n)(z)=Delta kappa_n(K),                            (4.2)

where `kappa_n(K)` is the nth occupation-number cumulant conditioned on the rank sector.

Because `G` is odd while `t(z)` need not be linear, Faà di Bruno immediately predicts the derivative parity ladder

    g_(2m+1) = O(L^(3(2m+1)/4)),
    g_(2m)   = O(L^(3(2m-1)/4)).                           (4.3)

In particular

    g_1 ~ L^(3/4),
    g_2 ~ L^(3/4),
    g_3 ~ L^(9/4),
    g_4 ~ L^(9/4), ...                                    (4.4)

The anomalously small exponent of the even derivative is therefore a consequence of matching oddness plus nonlinear thermal coordinates, not evidence by itself for a new even primary.

Two especially useful ratios follow:

    g_2/g_1 -> t''(z_c)/t'(z_c),                           (4.5)

which is nonuniversal thermal-coordinate curvature, whereas

    g_3/g_1^3 -> G'''(0)/G'(0)^3                           (4.6)

is independent of the thermal metric and is a candidate universal number for each torus modulus/observable convention.

The exact L=3,4 controls are:

| L | `g1/L^(3/4)` | `g2/L^(3/4)` | `g3/L^(9/4)` | `g2/g1` | `g3/g1^3` |
|---|---:|---:|---:|---:|---:|
| 3 | 1.29195257 | -0.08509569 | 0.05542263 | -0.06586595 | 0.02570085 |
| 4 | 1.31877015 | -0.08831509 | 0.06321249 | -0.06696777 | 0.02756104 |

Again, this is not a two-point exponent fit. Issue #769 can extend these quantities to L=5 at essentially zero marginal cost once the rank-by-occupation histogram is available.

Higher odd ratios `g_(2m+1)/g_1^(2m+1)` are natural metric-free targets. Even/odd neighbouring derivative ratios additionally test the thermal-coordinate normal form, but small widths can have large subleading corrections.

## 5. Conditional global-density CLT implied by the same scaling ansatz

There is a further non-obvious consequence. Under homogeneous Bernoulli `p=logistic(z)`, the exact conditional mgf in rank sector `j` is

    E[ exp(tK/L) | r=j,z ]
      = [(1+e^(z+t/L))/(1+e^z)]^(L^2)
        * P_j(z+t/L)/P_j(z).                               (5.1)

After centering by `L^2 p`, the first factor tends to the Gaussian mgf `exp[p(1-p)t^2/2]`. Under a smooth near-critical scaling law, `z -> z+t/L` changes `lambda` only by `O(L^-1/4)`, so the second factor tends to one for every sector with positive limiting probability. Therefore, conditionally on the scaling hypothesis,

    (K-L^2 p)/L | (r=j)
       => Normal(0,p(1-p))                                 (5.2)

for every `j=0,1,2` in the interior near-critical regime.

Topology and total density therefore decouple at the central-limit scale even though their subleading covariance controls thermal response. More precisely,

    E[K|j]-L^2p = O(L^(3/4)),
    Var(K|j)-L^2p(1-p) = O(L^(3/2))                        (5.3)

in a generic near-critical point; at the matching-symmetric centre, parity can suppress the leading even-sector variance difference to the smaller coordinate-curvature scale described in Section 4.

This provides a new interpretation of the Fisher result (2.3): the topology-density correlation vanishes at CLT scale, while the `L^(3/4)` conditional mean separation remains exactly the quantity that moves the near-critical rank law.

## 6. Practical consequences

The next data collapse should use the full triplet `(P0,P1,P2)` and the conditional K cumulants, not another single exponent fit. The hard controls are:

1. matching symmetry (3.2);
2. Fisher angle `rho_KX^2 ~ L^-1/2`;
3. parity ladder (4.3);
4. metric-free odd ratio (4.6);
5. conditional-density CLT (5.2).

These tests use existing rank-sector histograms whenever available. They do not require a new microscopic observable and do not alter the original-U contract.
