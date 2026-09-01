# P537 contact gauge-horizontal sidecar

This sidecar uses only the frozen exact N25 tables and the frozen N65 20M
contact-stage shards.  It does not change the canonical carrier, draw a new
configuration, add a descriptor, or overwrite the preregistered score.

For every source displacement `d`, it evaluates the finite identity

```text
T_C,d = T_C,d^hor + beta_d C_C,d,
C_C,d = p(1-p) sum_z <I_C Htilde>,
a_perp,d = a_d - <a_d> - beta_d S.
```

Thus `T_C^hor=sum_d <I_C Htilde u_z a_perp,d>` is invariant under the common
thermal-coordinate change `a -> a+cK+d`, whereas the original allocation obeys
`T_C -> T_C+c C_C`.  The N65 result retains all 63 displacement-specific
`beta_d` values and uses the joint 100-batch delete-one covariance of the new
source block and the independent P45 baseline.

## Aggregate results

```text
                         N25 exact midpoint       N65 estimate (joint SE)
T_full                  -8.39890318742124e-5     -7.43314230984601e-6 (8.49960e-7)
T_contact canonical     -4.94883991645077e-6     -1.89583506332088e-7 (6.18083e-8)
C_contact               -1.63571607585739e-4     -1.75852780207318e-5 (2.01852e-6)
sum beta_d C_contact,d  +1.89490872328009e-8     +2.75493310160022e-10 (4.58429e-11)
T_contact horizontal    -4.96778900368357e-6     -1.89858999642248e-7 (6.18363e-8)
T_remainder horizontal  -7.90212428705288e-5     -7.24328331020377e-6 (8.55635e-7)
```

The gauge-horizontal correction changes the contact allocation by `0.3829%`
at N25 and `0.1453%` at N65.  It is statistically resolved at N65 but tiny
relative to the full thermal response.  The contact identity closes to floating
zero in every N65 central/delete-one evaluation; the N25 midpoint residual is
`-2.29e-20`.

## Scale decision

The previously observed contact-subtracted power is preserved:

```text
representative              q_N from N25 -> N65       N^(5/2) T_rem at N25,N65
canonical                   2.50111483375              -0.247000600, -0.246737626
gauge-horizontal            2.50090370593              -0.246941384, -0.246728242
```

For the gauge-horizontal representative, the two scaled centers differ by only
`-0.08631%`.  Propagating the complete N65 joint jackknife gives
`q_N=2.50090 +/- 0.12363` (delta-method SE; 95% interval
`[2.25859,2.74321]`).  Therefore the `T_rem ~ N^-5/2` / `J_rem ~ N^-5/4`
fingerprint is not an artifact of the non-horizontal contact allocation.

This result does not make the contact-stage split unique.  It only removes the
specific common-thermal-coordinate ambiguity identified by `C_C` and shows
that doing so leaves the two-size remainder fingerprint essentially unchanged.
