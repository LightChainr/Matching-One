# Orientation-amplitude smoke test from exact tiny tori

Status: exploratory exact finite-size evidence. The reference `p` below is used only to compare amplitudes; its last digits are irrelevant at these tiny sizes.

Take

\[
p_{ref}=0.592746050790.
\]

For each exact small matching polynomial `M_L(p)` from `exact_small_matching_polynomials.md`, evaluate

\[
A_M(L)=L_{phys}^{13/4} M_L(p_{ref})
\]

and

\[
B_M(L)=M'_L(p_{ref})/L_{phys}^{3/4}.
\]

The powers are motivated by the Mertens-Ziff scaling analysis of the ordinary matching root:

- `M_L(p_c)` is approximately `L^-13/4` if the root bias is `L^-4`;
- `M'_L(p_c)` scales as `L^(1/nu)=L^(3/4)`.

## Exact-polynomial evaluation

| geometry | parameter | `L_phys` | `M(p_ref)` | `L_phys^(13/4) M` | `M'(p_ref)/L_phys^(3/4)` |
|---|---:|---:|---:|---:|---:|
| axis | 2 | 2 | `+0.1585008563` | `+1.5079` | `1.8289` |
| axis | 3 | 3 | `+0.02496987709` | `+0.8873` | `1.7604` |
| axis | 4 | 4 | `+0.01033323611` | `+0.9353` | `1.7627` |
| diamond | 2 | `2 sqrt(2)` | `-0.04531724101` | `-1.3298` | `1.7531` |
| diamond | 3 | `3 sqrt(2)` | `-0.007869374679` | `-0.8625` | `1.7667` |

## What is striking

For all but the tiniest axis `L=2` point, the derivative metric factor

\[
M'_L/L^{3/4}
\]

is already nearly geometry-independent: roughly `1.75–1.77`.

By contrast, the scaled critical residual

\[
L^{13/4}M_L(p_{ref})
\]

has a clear orientation sign flip:

```text
axis:     +0.887, +0.935
diamond:  -1.330, -0.862
```

The closest physical-size comparison in the table, axis `a=4` versus diamond `d=3`, has amplitudes `+0.935` and `-0.862`, already equal in magnitude to about 8% despite both systems being very small.

This strongly suggests that the opposite-side root approach is primarily produced by a sign-changing amplitude in `M_L(p_c)`, not by a radically different thermal metric factor in `M'_L`.

## Root-bias interpretation

Linearizing the matching root,

\[
0=M_L(p_L^*)\approx M_L(p_c)+M'_L(p_c)(p_L^*-p_c),
\]

gives

\[
p_L^*-p_c\approx-\frac{M_L(p_c)}{M'_L(p_c)}.
\]

If

\[
M_L(p_c)\sim \pm A L^{-13/4},
\qquad
M'_L(p_c)\sim B L^{3/4},
\]

then

\[
p_L^*-p_c\sim \mp(A/B)L^{-4}.
\]

Thus an orientation-odd leading amplitude in `M_L` produces exactly the observed orientation-odd `L^-4` root bias.

## Current conjecture sharpened

The data motivate the following stronger statement to test:

> After physical-length normalization, the leading `L^-13/4` amplitude of the square-site torus matching function is predominantly spin-4 / orientation-odd, while the leading thermal derivative amplitude is orientation-even.

This is more specific than saying that "finite-size corrections depend on orientation."

## Falsification path

The next Pell pairs are decisive. For `(a,d)=(7,5),(17,12),(41,29),...`, measure directly at a common fixed `p_ref` near the threshold:

\[
A_A=L_A^{13/4}M_A(p_{ref}),
\qquad
A_D=L_D^{13/4}M_D(p_{ref}),
\]

and the derivative-normalized quantities.

The spin-4 version predicts:

1. `A_A` and `A_D` stabilize to opposite signs;
2. their magnitudes approach each other;
3. `B_A` and `B_D` approach the same positive limit;
4. the orientation gap of roots then follows an `L^-4` law without needing an assumed `p_c`.

A common shift in `p_ref` affects both geometries through the nearly equal derivative term. Joint fitting can separate that common thermal displacement from the orientation-odd residual, so the experiment need not depend sensitively on disputed final digits of `p_c`.

## Important caution

The exponent `13/4` is being used here because it is the Mertens-Ziff hypothesis corresponding to the observed `L^-4` root convergence; their direct finite-size slopes at modest `L` were not exactly asymptotic. This table should therefore be rerun for alternative correction/log models once larger orientation-paired data exist.
