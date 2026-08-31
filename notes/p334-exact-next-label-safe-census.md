# P334: exact joint-safe probabilities separate role-within from role-between information

The complete vacant-label census is now available for every one of the
original40,000 paired checkpoints. It supplies exact `pi_joint_safe=h_joint/d`
without estimating a small per-prefix denominator from eight sampled quartets.
This enables the shared covariance analysis to distinguish genuine
within-safe/within-changing variation from variation between those two roles.
It does not itself score that decomposition.

## Exact definition and finite calculation

For each prefix, a label is safe for one orientation exactly when its next
rank equals that orientation's original checkpoint rank. A jointly safe label
preserves both ranks. R2 is always rank-preserving. Thus this is a probability
of the intersection under the old shared-label coupling, not the product of
the two marginal safety probabilities.

The R0 calculation uses the within-component contact-address span already
defined for the sampled-label readout. In R1, append the existing global
essential line to these winding generators. The span then gives the complete
after-insertion ambient image. Spanning-tree gauge changes add multiples of
the old line, so they cannot change whether a transverse new direction is
present. Untouched essential components are included by retaining the old
global line. No new R1 winding label is introduced. R2 needs no contact query.

Every remaining label is visited once, using the unchanged checkpoint UF;
there are no new site sets, random labels, suffixes or DP states.

| N | Prefixes | Vacant labels per prefix | Total positions | All labels jointly safe | No label jointly safe | Runtime |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 325 | 20000 | 132 | 2640000 | 2186 prefixes | 0 | 0.538s |
| 425 | 20000 | 173 | 3460000 | 2121 prefixes | 0 | 0.742s |

The sums of jointly safe counts are2356656 and3148562. A pooled safety rate
must not replace the saved individual-prefix rational probabilities in a
nonlinear role normalization. Prefixes with pi_joint_safe=1 have no changing
role; any weighted changing-role contribution there is zero, not a0/0 ratio.

## Why this is sufficient for the queued variance readout

For role c of exact probability pi_c, two independent same-role next-label
draws occur with probability pi_c². Their masked half-difference product has
conditional expectation `pi_c² B_c`, where B_c is the conditional response
covariance inside role c. Dividing that masked contribution by the exact
per-prefix pi_c gives `pi_c B_c`, the required within-role contribution.
The total response covariance minus the sum of those within-role terms is
the between-role part. A mixed-pair mask alone contains both contributions
and cannot be equated to between-role information. The coordinator implements
the shared-batch subtraction; this artifact only supplies the exact weights.

## Immutable delivery

Data lock: `e9dc7a1078b2c64b319f4a36ffc1c844e8426aa0`.
Code: `0e4db1b8ccae26f2953522ff7428162b12e9e8fa`.
Original prefixes: `9c495ab13e65f2bc93dc0849ee3b73f88724c4b1`.
Files: `results/p334-exact-next-label-safe-census/N{N}/N{N}.csv.gz`, one perN.
Columns are

```
N,batch,counter,k0,d,first_oldrank,second_oldrank,
first_safe_count,second_safe_count,joint_safe_count
```

Metadata and SHA manifests retain provenance. The forty original batches
remain the shared inference units. This deterministic finite census adds
no independent data population or intrinsic/universal notion of shared space.
