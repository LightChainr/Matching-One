# The forced plaquette source has not yet resolved a global-U transmission

## Fixed-block result

The independent F4 experiment completed its entire frozen budget and gives
**INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP**. None of the four
familywise95% intervals excludes zero, and none fits wholly inside the
declared +/-0.5 practical band. Thus we neither establish transmission nor
establish a negligible response. The result does not narrow the physical
mechanism space by rejecting the zero-projection model.

| N | V_F4 | paired jackknife SE | simultaneous95% interval |
|---:|---:|---:|---:|
|65|0.0648901|0.2430704|[-0.5422282,0.6720084]|
|85|0.8085407|0.3815157|[-0.1443732,1.7614546]|
|130|0.0471853|1.3682577|[-3.3703195,3.4646901]|
|170|-0.7352727|2.2324582|[-6.3112958,4.8407505]|

These are derivatives with respect to the **bulk** source parameter in
exp(tF4), with F4 the number of fully occupied elementary faces. They are
not divided by N, not source-standardized, and not measurements of an
operator named spin4. The N85 point estimate does not override the fixed
four-coordinate family rule. No common sign or cross-size law was fit.

## What remains exact, and what this data did not supply

The finite checkerboard source identity remains

`Ctot_parent=Ctot_child+F4`,

so its endpoint U response is

`V_parent,end^cluster=2^(13/8)(V_child^cluster+V_child^F4)`.

At the configuration-law level the extra term is forced. The present
observer-level experiment asks whether that term can be omitted when
computing the original global U response. Its answer is unresolved at the
declared budget. We do not replace that answer with the existence of a
local correlation or the exact source identity itself.

The [thermal-quotient proof](plaquette-source-thermal-quotient.md) removes
the single-site density clock exactly. If a future independently justified
experiment resolved V_F4, it would necessarily involve its centered
multisite part. This statement alone supplies neither a sign nor a field
identity. The failed P154 lag1 and P334 contact-transfer decisions are not
rescinded or reinterpreted by this different source.

## One prescribed source; no archive training

The [protocol](../experiments/p337-f4-transmission-20260831/PROTOCOL.md),
producer, scorer, complete budget and stop rules were pushed at
`0f7a083770d31095e7b4d688d544637d8fc09658` before production. All80M
fresh permutation counters were generated afterward:20M at each N,
100 equal batches per N. The four sizes have distinct seed domains.
Each permutation is shared by both orientations and by the ordinary and
forced-face processes; these four sweeps are not four independent samples.

The exact identity

`Cov(O,F4)=N p^4 [E(O | one face full)-E(O)]`

allows a one-face paired estimator. The forced law has N-4 free Bernoulli
sites. Its correct canonical degree N-4, the derivative of p^4, the pooled
root motion and thermal-denominator response are all included. Every
baseline and every omitted-batch root comes from the new ordinary stream;
no old anchor, fitted residual descriptor, or selected source coordinate
enters. Full omission vectors and covariance are retained in the
[machine-readable result](../results/p337-f4-transmission-20260831/scored/score.json).

Raw batch statistics and receipts were committed as `f6006b61` before the
single scoring invocation. Score/report commit: `25ca3635`. No scientific
replay or repeated test suite was run. The only preproduction runtime smoke
was40 N170 permutations in a disjoint seed domain, with no scientific fit.

## Runtime and decision

NePnUn/N65,551oUR/N85,TVVfoB/N130,TgFr7R/N170 each used14 threads on
measured14.5CPU/25GiB containers. Production elapsed13.982/17.922/28.304/
34.510seconds respectively; all exited0. Producer peak RSS was4.2–5.1MB.
The runtime is not a reason to extend a frozen block after its result.

The +/-0.5 band was a declared finite-scale resolution, not a universal
natural scale. The fixed decision ends this block: no top-up, source
substitution, extra size, alternate sign convention, or descendant
descriptor search. The statistical outcome is unresolved, not a theorem
of zero. Further research should not count these80M replicas as evidence
for a measured microscopic-to-global transmission.

## Science card

- Changed mechanism space: no empirical zero-projection rejection in this
  block; the exact decimation source dictionary and density-clock quotient
  are separate mathematical results.
- Not proved: negligible coupling, continuum field identity, asymptotic
  exponent, a source winner, or rescue of the previously stopped sources.
- Observer/sector/source/geometry: original pooled-root global U; ordinary
  uncharged orientation contrast; equilibrium F4; Gaussian pairs at the
  four declared areas.
- Dependency: four independent fresh N groups; all four modes within each
  replica and all statistics within each batch remain paired. No old block
  is pooled and no diagnostic receives another independent evidence vote.
- Upgrade condition: an independently motivated, prospectively specified
  same-global-observer result with adequate resolution. This card does not
  authorize enlarging or retuning the completed experiment.
