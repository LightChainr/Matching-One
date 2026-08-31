# One forced source, one unchanged global observer, one independent block

## Scientific decision

The finite checkerboard decimation proved at `20743651` forces
`Ctot_parent=Ctot_child+F4`, where F4 counts fully occupied square faces.
Its exact source transport is

`V_parent,end^cluster=2^(13/8)(V_child^cluster+V_child^F4)`.

The missing quantity is **V_child^F4 on the original global U**, not another
local response. This fixed experiment tests the narrow, falsifiable
zero-projection model: the forced extra source could be omitted at this
observer because its response vanishes at all four declared sizes.
We do not claim that four sizes exhaust all mechanisms or identify a field.
No nonzero sign follows from the topology theorem.

The [thermal-quotient derivation](../../notes/plaquette-source-thermal-quotient.md)
also removes a specific alternative explanation in advance. The single-site
part of F4 generates `4p^4(1-p)*partial_p`, a common thermal clock annihilated
by this U. Any resolved nonzero response must therefore involve the
centered two-, three-, or four-site part; it cannot be explained by that
common one-site density reparameterization alone.

This source was selected by a configuration identity, without fitting a
residual. It is an equilibrium plaquette fugacity, not a retimed or rescaled
version of the stopped P154 lag1 policy. The P154 and P334 decisions remain
unchanged.

## Exact paired estimator

For either translation-invariant homology observable O=q or E, write
F4=sum_f I_f. There are N translated elementary faces and P(I_f=1)=p^4.
Therefore

`J_O=Cov(O,F4)=N p^4 [E(O | I_0=1)-E(O)]`.

Use one uniform site permutation for the ordinary process and its
restriction to the N-4 non-face sites for the forced-face process. Each
orientation has its own face `{0,a,b,a+b}` in the declared cyclic labels.
Both orientations use the same replica permutation. The forced four sites
are occupied before the sweep; forced birth ranks count only remaining
sites and may be zero. Its canonical Bernstein polynomial has degree N-4,
whereas the ordinary polynomial has degree N.

For first/second rank births F1,F2, `q=-1+F1+F2`, `E=1-F1+F2`.
Reconstruct their value and thermal derivatives without a p-grid fit.
Differentiate the prefactor too:

`J_O,p=N[4p^3 delta_O+p^4 delta_O,p]`.

At the fresh pooled ordinary root Q=mean(q)=0, let
`Y=P4(E)`, `D=Q_p`, `A=N^(13/8)/2`. The response includes the moving root
and the changing thermal denominator:

```
p0dot = -mean(J_q)/D
Ddot  = mean(J_q,p) + Q_pp*p0dot
V_F4  = A/D * [P4(J_E,p) + Y_pp*p0dot - (Y_p/D)*Ddot].
```

No archived baseline or source statistic enters the score. Delete each
whole paired batch and recompute the root, jets, denominator, and response.
Save the full within-N covariance and original omission vectors. Different
N have independent fresh seed domains; their estimates are not paired.

## Fixed budget and stops

The machine-readable [contract](CONTRACT.json) fixes N65/85/130/170,
20 million replicas each, 100 equal batches per N, seeds, geometry order,
projector, root bracket, and familywise intervals. Freeze the complete
producer/scorer/contract in Git before any production counter is consumed.
Compilation and a tiny disjoint-seed structural smoke do not score this
experiment and are not scientific evidence.

Primary family: four V_F4 coordinates, two-sided normal intervals with
Bonferroni family level95%. Any interval excluding0 rejects the declared
zero-projection model. There is no forced winner if it survives.

The separate practical band is +/-0.50 in the declared bulk-fugacity U
units. It is a chosen finite-size resolution, not a universal or
data-estimated natural scale. If all four simultaneous intervals are
inside the band, stop prioritizing F4 as a main coupling at that resolution.
An interval wholly outside the band resolves a material coupling for that
size. Both decisions may coexist: a nonzero but small response is possible.

Publish an inconclusive outcome as such. No top-ups, extra sizes, sign
flips, alternate source normalizations, or new descriptors in this block.
An excluded model is not reinstated by a later exploratory coordinate.

## Lifecycle

observer: original root/slope-normalized global U; sector: ordinary
uncharged H4 orientation difference; source: equilibrium full-face F4;
geometry: four declared Gaussian quotient pairs; dependency groups:
`fresh_F4_N65`, `fresh_F4_N85`, `fresh_F4_N130`, `fresh_F4_N170`.
The theory motivating the source and old endpoint predictions are not
additional independent data in this experiment.
