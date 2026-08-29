# Issues #215/#275: a computable lattice proxy for `g[A_top,Q4 epsilon]`

Status: exact finite-lattice response, conditional CFT bridge, and a separate
mechanism conjecture. The three levels are not interchangeable.

## 1. The finite lattice coupling is a connected response

For one configuration let

\[
 q(\omega)=r(\omega)-1=A_{\rm top}(\omega).
\]

The previous insertion oracle supplies the matching-odd spin-four source

\[
 J_{D,4}(\omega)=\sum_v\chi_4(\ell_v)
 (I_{12}-I_{01})_v,
\]

with a separate version in which `chi4(ell)` is replaced by the local
landing-sector H4 mark. Every nonzero term has a canonical `ell`; simultaneous
`0->2` births cancel from `D`.

Introduce an actual finite-lattice source at registered phase `phi`:

\[
 \mathbb P_{p,h}(\omega)
 \propto \mathbb P_p(\omega)
 \exp\!\left[h\,\operatorname{Re}
 (e^{-4i\phi}J_{D,4}(\omega))\right].
\]

Differentiating the normalized finite sum gives exactly

\[
 \boxed{
 \left.\partial_h\langle A_{\rm top}\rangle_{p,h}\right|_{h=0}
 =\operatorname{Cov}_p\!\left(
 A_{\rm top},\operatorname{Re}(e^{-4i\phi}J_{D,4})
 \right).}
\]

This is the desired lattice coupling proxy: it is a genuine geometric/source
derivative of the exact global projector, not an exponent fit.

Define the birth-mass-normalized version

\[
 \gamma_{D,4}(N,\tau)=
 \frac{\operatorname{Cov}(A_{\rm top},J_{D,4})}{B_N},
 \qquad
 B_N=\left\langle\sum_v(I_{01}+I_{12})_v\right\rangle=M_N'(p).
\]

This normalization removes the total thermal/rank-pivotal density and retains
the complex spin-four phase and sign.

## 2. Why symmetry allows this response

Under primal/matching complement plus insertion reversal,

```text
A_top -> -A_top,
J_D4  -> -J_D4,
J_S4  -> +J_S4.
```

Therefore

```text
Cov(A_top,J_D4)  is complement-even and allowed;
Cov(A_top,J_S4)  is complement-odd and changes sign.
```

The latter is a built-in sign/null control when primal and matching backends
are combined. The `D` response is not killed by either Alexander parity or
spin: `chi4` is projective, C4-invariant, and changes sign under the Gaussian
`1+i` rotation.

This answers the selection question but not field identity. Symmetry says the
matrix element may be nonzero; it does not say that the only contributing
continuum state is `Q4 epsilon`.

## 3. Tiny exact coupling is nonzero

The oracle exhausts every configuration on axis `L=2`, Gaussian `(2,1)`, and
axis `L=4`. For axis `L=4`, translation invariance reduces the site sum to the
origin-root covariance times `N=16`; this is exact for a translation-invariant
global `q`.

At `p=1/2`:

| geometry | `Re Cov(A,J_D4)` | `Im Cov(A,J_D4)` | landing `Cov(A,J_D)` | `Re gamma` | `Im gamma` |
|---|---:|---:|---:|---:|---:|
| axis `L=2` | `1` | `0` | `0` | `1/3` | `0` |
| Gaussian `(2,1)` | `-49/128` | `21/16` | `0` | `-49/400` | `21/50` |
| axis `L=4` | `5013127/4194304` | `0` | `6977235/33554432` | `5013127/17240064` | `0` |

Every `D` covariance is identical on the matching backend, while every `S`
covariance reverses sign. All complement residuals vanish exactly.

What this proves:

- no finite-lattice symmetry forbids the line-resolved coupling;
- the registered axis phase has a positive tiny coupling;
- the Gaussian control carries a nonzero complex rotation character;
- the local landing version is already nonzero on axis `L=4`.

What it does not prove:

- survival or sign of the coupling as `N -> infinity`;
- dominance by an `x=21/4` field;
- the Q4/Jordan modular fingerprint.

Thus tiny nonzero closes the selection obstruction but supplies only a
finite-size sign/phase, not `g[A_top,Q4 epsilon]` itself.

## 4. What the archive can and cannot supply

Main `bedc94b` reconstructs the two essential-birth densities from the
committed axis `L=8` rank archive. At its balance root it already supplies

```text
B_N=M' = 8.3339594126152504344...,
D_N=f12-f01 = +0.02033467424339790896....
```

But the archive lacks

```text
ell/iota,
line or landing H4 at the births,
same-sample A_top * J_D4 cross moments.
```

The old marked-pivotal batches likewise retain H4 counts but not the global
rank coordinate or its product with the mark. Consequently neither archive
can reconstruct the connected coupling. This is an information boundary, not
a request for more retrospective fitting.

## 5. Scaling if the source flows to `x=21/4`

There are two useful normalizations and they must not be mixed.

### Canonically normalized continuum density

For a local field of dimension `x=21/4` on a fixed-shape torus with
`N=L^2`,

\[
 c_v\sim L^{-x}=N^{-21/8},
 \qquad
 \sum_vc_v\sim L^{2-x}=N^{-13/8}.
\]

This is the ordinary integrated irrelevant-field rule.

### The measured rank-birth counting measure

The unmarked birth mass is a thermal susceptibility. Assuming the percolation
value `nu=4/3`,

\[
 B_N=M'(p_c)\sim L^{1/\nu}=N^{3/8}.
\]

If the **relative** coupling isolates the Q4 correction, then

\[
 \boxed{\gamma_{D,4}=C_N/B_N\sim N^{-13/8}}.
\]

Therefore the raw observables in the proposed stream should scale as

\[
 C_N\sim N^{3/8-13/8}=N^{-5/4},
 \qquad
 C_N/N\sim N^{-9/4}.
\]

For the inherited top Jordan partner, multiply these laws by an affine
`A+B log N` factor. The modulus/phase dependence must also agree with the
existing Q4/Jordan oracle; a power alone is insufficient.

## 6. Minimal next sufficient statistics

An unbiased site-sum estimator does not require evaluating every root. For
each Bernoulli configuration choose one uniform root `V`, compute its two
counterfactual ranks, and multiply the root insertion by `N`:

\[
 \mathbb E_V[Nj_V\mid\omega]=\sum_v j_v(\omega_{-v}).
\]

Within every aligned batch retain

```text
q=r-1 and q^2,
I01, I12, direct-0->2,
ell, iota, Re/Im chi4(ell),
landing H4,
N*S_H4, N*D_H4,
q*(N*S_H4), q*(N*D_H4),
unmarked N*S birth mass.
```

For orientation pairs reuse the same configuration counter, random root and
delete-one batch. These raw sums reconstruct means, connected covariances,
the full common-field covariance, the even sign control, and `gamma_D4`
without saving configurations.

The threshold/permutation stream should additionally retain
`(K1,K2,ell,iota)` for birth-density tomography, but endpoint histograms alone
are not a substitute for the same-sample `qJ` cross moment.

## Claim layers

1. **Exact:** the exponential-tilt derivative equals the connected covariance;
   symmetry permits the `D` response; tiny couplings are nonzero with exact
   primal/matching sign controls.
2. **Conditional CFT bridge:** if `J_D4` flows predominantly to the thermal
   `Q4 epsilon` family, `gamma_D4` must follow `N^-13/8` and the independent
   modulus/Jordan fingerprint.
3. **Conjecture:** the line-resolved rank-one lifetime source has nonzero
   continuum overlap with the `Q4 epsilon`/energy-hull Jordan module and is the
   local mechanism behind the global matching-odd H4 response.

Reproduce with:

```bash
python scripts/rank_birth_atop_coupling_proxy.py \
  --output results/rank-birth-atop-coupling/latest.json
python scripts/rank_birth_atop_coupling_proxy.py --format markdown \
  --output results/rank-birth-atop-coupling/latest.md
python -m unittest discover -s tests -p 'test_rank_birth_atop_coupling_proxy.py'
```

