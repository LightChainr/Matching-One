# P321 deep equal-area independent pool

## Result

All six planned independent campaigns are complete.  The pool retains the
full five-root covariance inside each run and combines disjoint RNG domains by
precision.  Effective samples per geometry are 10M at N=144 and 6M each at
N=576 and N=1296.

| N | square | rho 16/9 | rho 9/4 | rho 4 | rho 9 |
|---:|---:|---:|---:|---:|---:|
| 144 | .5927244101 | .5927066363 | .5926566762 | .5924926325 | .5914144868 |
| 576 | .5927612278 | .5927599405 | .5927399685 | .5927168333 | .5926205708 |
| 1296 | .5927593254 | .5927357442 | .5927312285 | .5927486666 | .5927279008 |

The fixed model

```text
p(N,rho) = pc + C_N(rho) N^-2 + D_N(rho) N^-3
```

gives `pc=.592742686993 +/- 7.12e-6` and `chi2=4.8973/4`
(`p=.2980`).  The preregistered conditional thermal-Q4/E4 relation gives
`chi2=3.6629/3` (`p=.3002`).  Adding the old 100k campaigns changes these to
`p=.2890` and `p=.3346`; the scientific decision is unchanged.  The E4 shape
is compatible, but this is not strong confirmation and it does not identify
the thermal one-point function.

## Why a naive D/C plot is invalid

For width normalization,

```text
C_width = C_N/rho^2
D_width = D_N/rho^3
D_width/C_width = D_N/(rho*C_N).
```

The point ratios are approximately `-149, -92, -5.6, -12.7, -6.0`, but the
first four `C_N` denominators have absolute z scores only
`1.58, .73, .29, 1.03`.  A symmetric delta-method error bar would therefore
turn a denominator singularity into an apparently measured large ratio.

The proof-carrying Fieller scorer instead finds:

- rho `1`, `16/9`, `9/4` and `4`: unbounded 99% ratio sets;
- rho `9`: a bounded 99% set `[-10.86,129.03]`, still containing zero;
- no geometry has a nonzero D/C ratio at 99% confidence.

Thus the five-point subleading curve is not identified.  This is not a
rejection of identity dressing, Jordan mixing or the E4 leading curve.  It
confirms the theoretical obstruction at `c2a5e2d`: a supplied
homology-resolved thermal one-point `F_t(tau)` can be scored, but the present
data must not be used to invent an arbitrary post-reveal E4-only correction.

## Integrity

- Six campaign counter intervals are pairwise disjoint.
- Repeated square histograms and moments are byte-identical inside every
  common-random campaign.
- The HZ N576-extra and Zy N1296-extra remote/local hashes agree file by file.
- Raw campaign hashes are retained without committing the approximately
  100 MB exploratory histograms.
- The portable result includes the complete 11 by 11 covariance of
  `[pc,C_N(5),D_N(5)]` and an old-100k sensitivity calculation.

## Scientific card

- **Mechanism space changed:** conditional E4 survives deeper data; an
  empirical five-geometry D/C selector does not yet exist.
- **Not proved:** E4 field identity, ordinary identity dressing, a Jordan
  extension, or a universal subleading curve.
- **Observer/sector/source/geometry:** individual `P2-P0` roots on five
  equal-area rectangle moduli, N=144/576/1296, full within-N covariance.
- **Dependency group:** six disjoint Huawei/local campaigns pooled by root
  precision; old 100k data are sensitivity-only.
- **Next lift:** compute the crossed/trivial thermal module amplitudes or
  `F_t(tau)`; if a purely empirical lift is still wanted, add an independently
  frozen fourth scale rather than merely plotting unstable ratios.
