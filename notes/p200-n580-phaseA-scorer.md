# Frozen N580 Phase A joint scorer

This scorer was frozen before the N580 Huawei target completed.  It accepts the
schema-compatible histogram, moments, and metadata files from the Phase A
runner and does not contain target-dependent mode selection.

It reuses the P50 intrinsic-center path.  In each full and delete-one batch
sample, it solves the root of the orientation-mean matching function and forms

```text
I_S  = N P4[S]
I_Du = N P4[D']/Mbar'
T_D  = N^(13/8) P4[D]
T_Su = N^(13/8) P4[S']/Mbar'.
```

At N580, the old binary64 recurrence evaluated at the historical `p=0.9`
bracket endpoint underflows.  The scorer therefore uses the repository's
existing arbitrary-precision binomial recurrence for the identical intrinsic
root and projector definitions.  This is a numerical stabilization, not a
change of coordinate or model.

The divisions by the intrinsic matching slope are the frozen thermal-width
coordinate used by P180.  No new width is fitted from N580.  The moments file
is consumed as a schema, batch, representation, counter-total, and exact
rank-gap identity gate; the four state coordinates themselves are reconstructed
from the full threshold histograms exactly as in P50/P180.

The 100 delete-one pseudovectors give the full target `4x4` covariance.  For
each frozen model, the scorer adds the committed N145/N290 prediction covariance
and evaluates one correlated quadratic form for

```text
x_N580 - x_model.
```

The order is fixed: ordinary q2 first, rank-2 Jordan second.  Marginal z values
are printed only as correlated diagnostics.  They are not four extra votes.

The test suite locks two independent paths: the committed N145 P50 histogram
reproduces the exact four-state vector archived by P180, and a deterministic
synthetic N580 hist/moments/metadata triple exercises the complete scorer before
any production target exists.
