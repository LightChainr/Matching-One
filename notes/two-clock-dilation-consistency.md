# One dilation generator, two clocks: the mixed-context theorem

Issue #255 proposes that Gaussian-cover scale and annulus log-radius are two
clocks for one finite-dimensional dilation generator.  The sharp first result
is that this claim cannot be identified by comparing two separately fitted
spectra.  It *can* be identified by a small mixed-context experiment.

Remove the declared spin rotation from the Gaussian transfer and write

```text
A = exp(-t_A G),       t_A = log|m|,
U = exp(-t_U G),       t_U = log(R_out/R_in).
```

Then two exact consequences follow:

```text
[A,U] = 0,
A^p = U^q  whenever p t_A = q t_U.
```

The second identity avoids a matrix logarithm and all of its branch choices.
For example a norm-4 Gaussian clock has `t_A=log 2`; an annulus ratio two has
the same time and should give `A=U`.  A norm-2 Gaussian clock has
`t_A=(1/2)log 2`; two steps must equal the ratio-two annulus transfer.

## Why matching spectra are insufficient

Two latent realizations recovered independently are each defined only up to a
similarity transform.  Equal characteristic polynomials, equal discriminants,
and even the same rank-2 Jordan class do not align those two gauges.  The exact
oracle in `scripts/score_two_clock_generator.py` includes two Jordan matrices
with the right repeated eigenvalues but different nilpotent directions; their
commutator and the required power identity are nonzero.

This also sharpens the phrase “the same nilpotent direction.”  If

```text
G = g I + K,  K^2=0,
```

then in one shared basis

```text
A/e^(-g t_A) - I = -t_A K,
U/e^(-g t_U) - I = -t_U K.
```

The proportionality is meaningful only after the two experiments share source
and readout anchors.  A scalar logarithm or two separately reconstructed
Jordan blocks cannot establish it.

## The mixed-context theorem

Let the frozen sources `b_r` span the reachable state and frozen readouts
`c_l` span its dual observable space.  Measure both path orders:

```text
c_l^T A U b_r,         c_l^T U A b_r.
```

If every rectangle difference vanishes, reachability and observability imply
`AU-UA=0`.  Thus a two-source by two-readout block is already sufficient for a
rank-2 candidate.  This is the smallest useful experiment; another scalar
endpoint sequence is not.

For a non-scalar two-dimensional `A`, every commuting matrix lies in
`span{I,A}`.  Once the commutator vanishes, the commensurate-time power identity
fixes the clock law without naming either state.  If the rectangle is nonzero
in a morphism-sensitive row but zero in the unmarked endpoint sector, the
natural conclusion is deck/Smith/category memory rather than a third scalar
bulk correction.

## What the current archive says

The archive contains useful Gaussian low-rank summaries and correlated annulus
shell contrasts, but it does not contain a common two-dimensional source/readout
basis or both mixed path orders.  Therefore the one-generator hypothesis is
currently **unscored**, not failed.  This is a positive design result: the next
measurement should be the mixed rectangle, with norm-2 squared versus a
ratio-two annulus as the exact clock match.

The risky prediction is that the unmarked endpoint rectangle closes in rank
two, while a symmetry-allowed morphism-sensitive rectangle remains nonzero
until one deck/Smith state is added.  That outcome would separate a local bulk
generator from finite-cover memory with one compact experiment.
