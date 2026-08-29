# Deployable N325 fixed-p chiral response

This implements the charged norm-five response proposed in
`norm5-chiral-hecke-phase.md` without turning it into a general cover package.

## Exact common field

The parent is \(z=8+i\), \(N=65\).  Its two same-parent children are

\[
z(2+i)=15+10i,\qquad z(2-i)=17-6i,
\]

both of order 325.  The first child is nonprimitive, so a cyclic vertex index
is not a valid shared field convention.  Instead, every child site is labeled
exactly by

```text
(deterministic parent representative j, fiber k in Z/5).
```

The fiber representative is \(u_k=(k,0)\) modulo the corresponding
multiplier matrix, and the child coordinate is \(r_j+M_z u_k\).  The mapping
gate verifies 325 unique labels, exact projection back to the parent, and the
deck action \(k\mapsto k+1\) under translation by the parent period \((8,1)\).
The same Bernoulli bit indexed by \((j,k)\) drives both children.

There is an important geometry distinction.  The two same-parent children are
not reflections of one another.  The true reflected pair is

\[
(8+i)(2+i)\longleftrightarrow(8-i)(2-i).
\]

The engine accumulates that true reflection as a conjugacy null, which is zero
configuration by configuration.  It does not mislabel
\(R_{(8+i)(2+i)}-\bar R_{(8+i)(2-i)}\) as a null.

## Cheap unbiased marked row

For each replica, one root in the common 325-site field is selected by the
counter-stable schedule `replica mod 325`.  The existing primal/matching
landing-marked pivotal \(H_4\) reference is evaluated at that root.  Multiplying
by 325 gives an unbiased estimator of the full opposite-character marked sum:

\[
\widehat O_{\bar\chi}=
\frac{325}{2}\zeta_5^{-k_{root}}
[H^{piv}_{4,primal}-H^{piv}_{4,matching}].
\]

This preserves the fixed-root cost of #215/#225 instead of toggling all 325
sites in every configuration.  The stream saves the real four-vector

```text
[Re R_(2+i), Im R_(2+i), Re R_(2-i), Im R_(2-i)]
```

with its complete batch covariance, the true reflection null, and delete-one
covariance for the complex same-parent handed ratio.

## 20k smoke

The deterministic 20,000-replica, 20-batch smoke completed locally in about
19 seconds with eight workers.  The response vector was

```text
[-6.085 +/- 5.110, -5.536 +/- 3.379,
  0.929 +/- 2.371,  2.950 +/- 3.499].
```

The ratio phase was \(149.8^\circ\) with a delta-method standard error of
\(56.4^\circ\).  H4 is the nearest frozen phase in this smoke, but the run does
not discriminate: its phase uncertainty is comparable to the 65-degree
minimum target separation.  Straight \(1/\sqrt n\) projection gives about
18 degrees at 200k, so the frozen tenfold run has a realistic chance to answer
the three-way question without a huge production campaign.

Reproduce the smoke:

```bash
python3 scripts/norm5_chiral_fixedp_mc.py \
  --samples 20000 --batches 20 --workers 8 \
  --p 0.592746050790 --seed 2265325020000 --radius 1 \
  --output results/local-20260829/P226-norm5-chiral-fixedp-smoke/chiral_response.json
```

Run the frozen Huawei job:

```bash
python3 scripts/norm5_chiral_fixedp_mc.py \
  --samples 200000 --batches 100 --workers 16 \
  --p 0.592746050790 --seed 2265325000829 --radius 1 \
  --production-manifest experiments/p226_norm5_chiral_fixedp_production_20260829.json \
  --output /workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.json \
  --batches-output /workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.batches.csv
```

