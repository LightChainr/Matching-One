# Conditional C3 phase contract for #275

This is a bounded exact design review, not a new measurement or field identification. It follows #275 comment 5490477497 (2026-09-01): a real C3 readout has conjugate nontrivial Fourier modes, so one angle cannot identify H4 versus H8 with a free complex amplitude.

## Fixed conditional observation model

For the SAME phase-calibrated observable and physical normalizer, suppose

```
y_j(theta) = Re[A_s exp(i*s*(theta+2*pi*j/3))], j=0,1,2.
z(theta) = (2/3) sum_j y_j(theta) exp(-2*pi*i*j/3).
```

Then H4 gives `z=A4 exp(4i theta)`, while H8 gives `z=conj(A8) exp(-8i theta)`. This transformation law must be justified by the observable dictionary. It is not implied merely by naming a representation, nor automatically true for the original q/E, a twist-normalized trace, or an auxiliary contact source.

Let delta be a second physical rotation after absorbing the first phase into A.

## Exact decisions before spending compute

1. **One shared complex amplitude; known relative normalization.** Candidate vectors are `A*(1,exp(4i delta))` and `B*(1,exp(-8i delta))`. They are identical model spaces iff `12 delta=0 mod 2pi`. With two balanced isotropic blocks their squared principal separation is `sin^2(6 delta)`. Thus 15 or 45 degrees are optimal for this particular amplitude/noise contract, whereas 30,60,90 degrees remain aliases.
2. **Unknown nonzero SIGNED real gain at each size/geometry.** The relative phase is defined only modulo pi. The two candidate sets alias iff `sin(12 delta)=0`: 15 and 45 degrees now also fail. For example at 15 degrees H4 gain +1 and H8 gain -1 yield exactly the same two readouts. The maximal modulo-pi phase separation occurs at 7.5 or 22.5 degrees. Those are conditional angle targets, NOT a claim that an existing same-area Gaussian torus realizes them.
3. **Unknown POSITIVE gains.** Keep the phase equality AND positive-real-part condition. A phase difference of pi cannot be absorbed; aliasing returns to `12 delta=0 mod 2pi`.
4. **Independent unknown COMPLEX gains.** Every rotation aliases; each candidate fits arbitrary pairs of complex readouts. More angles cannot replace a physical amplitude/phase transport relation.

Two useful signed-gain zero relations, without dividing by a noisy signal, are

```
Im[z2*conj(z1)*exp(-4i delta)]=0,
Im[z2*conj(z1)*exp(+8i delta)]=0.
```

Both signals must be nonzero. Under noise, fit the original real vector with its full covariance and all allowed nuisance gains; these polynomial zeros are not a finite-sample confidence procedure.

## The exact norm-2 warning

For multiplication by `1+i`, delta=45 degrees and `exp(i12delta)=-1`. In an entire repeated norm-2 sequence, the H4 phase `(-1)^n` can be absorbed into independently signed per-size gains, while the aliased H8 phase is 1. More points on that same phase sequence do not fix this nuisance-induced nonidentifiability.

For `2+i`, `exp(i12delta)=11753/15625-i*10296/15625`, so the signed-gain phase separation is strictly nonzero. This is arithmetic only: multiplying a torus changes size and does not establish that a particular old C3 dataset obeys the required observer/normalizer/phase relation. It does not challenge frozen experiments with different observables or amplitude assumptions.

## Covariance-weighted comparison

Writing `X_s=2^-1/2 [I; R(s delta)]`, with effective spins 4 and -8, the cross Gram obeys `(X4^T X8)^T(X4^T X8)=cos^2(6delta) I`. For a positive-definite real covariance C, profile the alternate amplitude using

```
G = X4^T [C^-1 - C^-1 X8 (X8^T C^-1 X8)^-1 X8^T C^-1] X4.
```

The separation at amplitude a is `a^T G a`. Whitening preserves rank but changes efficiency. No empirical covariance or archive score is computed here. Singular covariance requires its actual statistical support, not an arbitrary ridge.

## Reproduction and actual checks

Requires Python and SymPy. Tested in CPython 3.13.5, a 4-CPU/4-GiB container. Run to a NEW file:

```
python experiments/p275-c3-phase-contract-20260901/phase_design.py --out /tmp/p275-phase-new.json
```

The script executes 66 exact checks: direct real-C3 Fourier evaluation, principal-angle Gram identities, joint ranks, the 15-degree signed-gain counterexample, Gaussian integer phase arithmetic, and non-isotropic positive-definite covariance controls. Local main runtime was about 24 seconds, peak RSS about 136 MiB. No MC, GPU, cloud task, new lattice enumeration, full-repository CI or physical field identification was performed. Source/input-only local replay reproduces the JSON.

The accompanying conversation bundle has 17 focused tests, including a separate standard-library Fraction/Gaussian matrix calculation and positive real-probability controls. Those extra test files and the separate N25 gauge review are NOT part of this two-file Git contribution.

This branch adds only this note and the deterministic script. It does not merge, alter the ledger, change priorities, or authorize a new angle/size experiment. The next necessary input is the admissible amplitude class for the SAME original observable; that input determines which existing second geometry can actually discriminate the two hypotheses.
