# F3 activation-flux tomography

## Result

The existing paired N65 archive now resolves the directional projective character response into rank-one birth and exit of the first-born line at the rank-two completion boundary. This is a new view of the same 20 aligned batches, not a new evidence vote.

The birth and first-line-weighted exit entries below are signed character coordinates. Both net/gross columns use `|birth-exit|/(|birth|+|exit|)`: zero means cancellation and one means reinforcement/no cancellation.

| point | p | plateau H | birth H | rank-one exit H at tau2 | dH/dp | H net/gross | charged |d(v,w)/dp| | plug-in Wald chi2(2) | asymptotic p(2) | charged net/gross |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `half` | 0.500000000 | -0.00341247 | -0.0077467 | -0.0273776 | 0.0196309 | 0.5589 | 0.0169735 | 5.5959 | 0.06093 | 0.3214 |
| `h_cross` | 0.573633326 | 4.16334e-18 | 0.033879 | -0.0263673 | 0.0602462 | 1 | 0.0385765 | 3.7688 | 0.1519 | 0.6039 |
| `p_ref` | 0.592746051 | 0.00110298 | 0.0368225 | -0.0169815 | 0.053804 | 1 | 0.0418601 | 5.2515 | 0.07239 | 0.6395 |
| `upper_probe` | 0.650000000 | 0.00285719 | 0.0231637 | 0.0189922 | 0.00417144 | 0.09895 | 0.0325993 | 10.834 | 0.004441 | 0.6612 |

At the frozen exploratory H crossing, the H slope is `0.0602462`: first-birth `0.033879` minus first-line exit at rank-two completion `-0.0263673`. Their signs are opposite, so the two topological boundaries reinforce rather than cancel (`H net/gross=1`). This reinforcement appears in every fixed-p delete-one batch reconstruction.

At `p_ref`, H remains boundary-reinforced (`dH/dp=0.053804`, net/gross `1`). The point-estimate norm of the first-line-weighted exit projection onto the fixed `(v,w)` plane is larger than the birth projection, but their relative phase is not yet stable: its delete-one cosine spans `[0.083,0.547]` at the crossing and `[-0.054,0.487]` at `p_ref`. The crossing delete-one values hold the full-data root fixed and do not propagate root uncertainty. The safe claim is two-boundary reinforcement in H, not one resolved charged ray or verified C3 transport.

By `p=0.65`, H itself nearly cancels (net/gross `0.09895`) while the fixed charged-plane projection is nonzero under a plug-in Wald diagnostic (`chi2_2=10.834`, asymptotic `p=0.004441`). Because all p points share the same archive, this is a mechanism selector for a fresh child—not an independent significance claim or verified C3 transport.

## Exactness and dependence

For every orientation, batch, line, character and reported p, the Bernstein partition and rank/line transition identities close, including character derivative equals signed birth coordinate minus the first-line-weighted exit coordinate at rank-two completion. The maximum residual is `4.44e-15`; the maximum binomial normalization residual is `3.11e-15`.

`DIRECT_RANK2` atoms are reported but contribute exactly zero to the zero-sum projective characters: first `381`, second `404`. Full covariance-of-the-mean matrices are retained for H/A/D and u/v/w, and per-batch boundary coefficients retain arbitrary cross-p covariance for the paired contrast.

## What this changes

The next fresh N130 archive should freeze the complete A4 triplet and these two boundary fluxes together, and add a section-audited completion record. The present archive knows when rank two completes but not which marked complement representative completes it. That representative is not intrinsic: a fresh child must retain raw winding, section/basis, transporter and ambiguity metadata before testing a gauge-covariant stabilizer-C3 completion character. A larger deterministic basis-shear run cannot answer this question; a later physical defect/source insertion can test whether the charged plane is more than a projective re-expression of the same archive.

Claim boundary: This analysis exactly re-expresses one existing paired N65 archive. It can locate directional character response at the rank-one birth and rank-two completion boundaries, but it is not a new evidence vote, a continuum-field identity, an intrinsic completion-line state, or an independent physical twist-source experiment.
