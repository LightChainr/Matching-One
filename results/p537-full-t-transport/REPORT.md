# P537 full-T transport quotient

Status: `COMPLETE_SIGNED_RESPONSE_CONTRACTS`

The already-produced N65 global sufficient statistics give the complete
canonical-pair thermal response, after restoring the one C4-equivalent NN
source column omitted by the local-carrier producer:

```text
pooled root p        0.5927311266364432
T_t                 -7.43314231e-6  +/- 8.49959721e-7
J_65=A_65*T_t/M_t   -0.00162250989 +/- 0.00018553008
95% interval        [-0.00198614885, -0.00125887093]
```

The exact N25 full response is `J_25=-0.0055194314248394015`.  Therefore

```text
J_65/J_25 = 0.29396323 +/- 0.03361398
95% interval = [0.22807982, 0.35984664]
effective abs(J) power = 1.28130 +/- 0.11967
```

This is a substantial contraction of the complete normalized response, not
only of the contact-stage carrier.  It is consistent with `J_N -> 0` and hence
with `T_N=o(M_t/A_N)` over this first exact-to-MC size direction.  Two sizes do
not establish the asymptotic little-o statement.

The contact-completion carrier contributes `5.892%` of full `T_25` but only
`2.551%` of full `T_65`.  It is therefore a genuine transmitted microscopic
operator but is not saturating the complete response.  Its three leading
cells scale as the candidate `T_cell~N^-3`, which becomes
`J_cell~N^(5/4)N^-3=N^-7/4` after original-U normalization.  The remaining
asymptotic obstruction is the nonlocal remainder and its near-critical
uniform transport, not another local-contact descriptor.

## The pooled-root displacement does not generate the finite signal

Reweighting the same N65 sufficient statistics to the prescribed square-site
reference `p_ref=0.592746050790` gives

```text
J_65(p_ref)                 -0.00162276136 +/- 0.00018552167
J_65(root)-J_65(p_ref)      +2.51467029e-7 +/- 9.06647139e-8
absolute transport fraction 1.54986e-4
```

Thus the observed displacement from this reference to the pooled root changes
the response by only about `0.0155%` of its magnitude.  This removes finite
root motion as a numerical explanation of the N65 signal.  The reference is
not a rigorous enclosure of the mathematical `p_c`; a uniform near-critical
transport theorem remains the exact-`p_c` gap.

## C4 reconstruction and dependencies

The 63 retained directions alone give `J=-0.00150763301`.  The omitted
`y=z=+e1` NN column, reconstructed as the mean of `-e1,+e2,-e2`, contributes
`-0.000114876876`, about `7.08%` of the complete value.  It cannot be silently
dropped.

Uncertainty combines two independent delete-one groups: the new 20M source
block and the P45 100M baseline, each with 100 batches.  Baseline omissions
resolve the pooled root before reweighting the source table.  The root-minus-
reference contrast is paired within each group before their variance
contributions are added.

No configurations were generated or replayed.  The scorer, exact formulas,
input hashes and complete numerical output are in
[`RESULT.json`](RESULT.json).
