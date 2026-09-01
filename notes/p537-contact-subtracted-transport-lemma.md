# Contact-subtracted transport and the quotient-reaching remainder

Status: mechanism synthesis after the exact N25 contact tensor and the held-out
N65 full-T quotient.  The two-size powers are numerical results; the collar
cancellation and pivotal bounds below are explicit theorem targets.

## 1. The subtraction exposes an L^-5 law

Let `T_contact` be the sum of the six contact-stage cells and define

\[
T_N^{\rm rem}=T_N^{\rm full}-T_N^{\rm contact},\qquad
J_N^{\rm rem}={A_N\over M_t}T_N^{\rm rem}.
\]

The exact N25 and held-out N65 values are

| N | `T_full` | `T_contact` | `T_rem` | `J_rem` |
|---:|---:|---:|---:|---:|
| 25 | `-8.3989031874e-5` | `-4.9488399165e-6` | `-7.9040191958e-5` | `-0.0051942129774` |
| 65 | `-7.4331423098e-6` | `-1.8958350633e-7` | `-7.2435588035e-6` | `-0.0015811275102` |

Consequently

\[
{T_{65}^{\rm rem}\over T_{25}^{\rm rem}}=0.09164399306,
\qquad q_T=2.50111483
\]

in the area variable `N`.  Equivalently,

```text
N^(5/2) T_rem = -0.247000600  at N25
                -0.246737626  at N65.
```

The centers differ by only `0.1065%`.  In original-U units the same statement
is `J_rem~N^-5/4`; its two-point power is `1.24478562`.  A paired delete-one
calculation gives

```text
SE(T_rem,65)       = 8.5562594e-7
corr(T_full,T_contact) = -0.05562
q_T                = 2.50111483 +/- 0.12362208
q_J                = 1.24478562 +/- 0.12362238
```

The `q_T=5/2` zero-parameter prediction from exact N25 is
`T_rem,65=-7.2512790e-6`; the observed difference is only `0.00902` standard
errors.  Thus the rational power is an unusually sharp two-size fingerprint,
although it is not by itself a proof of the asymptotic law.

The contact share falls from `5.892%` to `2.551%`.  The noncommuting local
contact operator is therefore a real, faster-decaying correction rather than
the leading asymptotic obstruction.

## 2. Candidate quotient-reaching bound

After subtracting the explicit contact OPE, the minimal useful bound is

\[
|J_L^{\rm rem}(p)|\le C
\sum_{r\le cL\atop r\ {\rm dyadic}}
{r^2\over L^2}\,
\pi_4(p;1,r)^2\pi_4(p;r,cL).                    \tag{1}
\]

This is not the retired blanket claim that every three-packet term vanishes.
It has three narrower ingredients.

1. **Contractible-collar quotient identity.**  If the two endpoint packets
   close inside their radius-`r` collars and the thermal/readout exploration
   does not cross the outer annulus, the axis/tilted fixed-M Schur difference
   equals the contact OPE already subtracted from `T_rem`.
2. **Pivotal domination.**  On the surviving event, the sum over the thermal
   site divided by `M_t` is a uniformly bounded signed pivotal kernel.  Arm
   separation then gives the three factors in (1); no extra `r^2` is charged
   for the thermal site.
3. **Near-critical uniformity.**  Quasi-multiplicativity, separation constants,
   the pivotal bound and `M_t~L^(3/4)` hold uniformly between `p_c` and the
   pooled root.

If `pi_4(a,b)<=C(a/b)^alpha`, equation (1) gives

\[
J_L^{\rm rem}=O(L^{-2\alpha}).                 \tag{2}
\]

With the square-percolation value `alpha=5/4`, equations (2) and
`M_t/A_N~L^-5/2` yield

\[
J_L^{\rm rem}=O(L^{-5/2})=O(N^{-5/4}),\qquad
T_L^{\rm rem}=O(L^{-5})=O(N^{-5/2}),             \tag{3}
\]

exactly the two-size fingerprint.  More importantly, the repository's strict
input `alpha>1` already gives `T_rem=O(N^(-9/4-delta))`, which is more than
enough for `T_N=o(N^-5/4)` once the three ingredients above are proved.

## 3. Root transport becomes a small theorem, not a new experiment

The matching-odd center law and slope lower bound

\[
M_L(p_c)=O(L^{-13/4}),\qquad \inf_I M_t\ge cL^{3/4}
\]

give `|p_L-p_c|=O(L^-4)` by the mean-value theorem.  This lies far inside the
near-critical window `L^-3/4`.  A uniform logarithmic derivative bound of
order `L^(3/4)` changes (1) only by `1+O(L^-13/4)` across that interval.
Numerically, reweighting N65 from the prescribed `p_ref` to the pooled root
changes full `J` by only `0.0155%`; this supports the scale separation but is
not the uniform theorem.

The next proof target is therefore precise: prove the contractible-collar
identity and bounded normalized pivotal density.  A new contact label, radius
scan or another local carrier sample does not address this gap.
