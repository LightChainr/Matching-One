# Exact fixed-K motif covariance oracle

For a uniformly chosen `K`-subset `X` of `N` labelled vertices and two motif
embeddings `A` and `B`, the joint indicator moment is

`E[1_{A subset X} 1_{B subset X} | K] = binom(K, |A union B|) / binom(N, |A union B|)`.

Therefore an overlap-union histogram determines the exact rational covariance
of any two motif counts. Bilinearity then gives the mean and full covariance
matrix of arbitrary signed motif-count combinations. The implementation keeps
embedding multiplicity, rejects collapsed/out-of-range embeddings, and never
uses floating point.

The checked fixture exhausts every `K=0,...,6` subset of a labelled six-cycle.
It compares two equal-multiplicity edge families and two equal-multiplicity
triangle families. Both signed differences have zero fixed-K mean, while their
nonzero covariance matrix agrees exactly with direct subset enumeration.

This is an algebraic and synthetic regression oracle. It does not measure
covariance with the production orientation target, fit a control coefficient,
estimate achieved variance reduction or wall time, or satisfy the Issue #40
`>=2x` production promotion gate. Issue #40 remains open.
