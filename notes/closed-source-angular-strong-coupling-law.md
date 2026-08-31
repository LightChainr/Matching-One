# The N25 angular strong-coupling law approaches zero from below

The saved complete N25 `(K,g,q)` histograms determine an exact result
stronger than the previous `U=O(exp(-3t))` bound. With `lambda=exp(-t)`,
`A=25^(13/8)/2` and the unchanged projector `DeltaCos4=1152/625`,

```text
U(lambda)/A = -(625/1152) lambda^11
              +(390625/1152) lambda^13 + O(lambda^15).
```

In particular the leading coefficient is **negative**. The law eventually
approaches zero from below, and its derivative is eventually positive:

```text
U_t/A = +(6875/1152) lambda^11 + O(lambda^13).
```

This is deterministic symbolic analysis after the four-point readout,
not a prospective prediction or independent evidence. No coupling point,
random sample, configuration enumeration or fit was added. Both saved
populations remain the single dependency group
`p337-N25-axis-tilted-exhaustive-configuration-populations`.

## Exact computation and the moving root

The script [p337_closed_source_asymptotic.py](../scripts/p337_closed_source_asymptotic.py)
reads only the two locked integer histograms, not their finite-coupling
scores. It uses Python Fraction arithmetic through degree32, the declared
upper limit. For each geometry alpha it constructs

```text
Z_alpha(h,lambda)=sum count(K,g,q) h^K lambda^g,
Hq_alpha=sum q count(K,g,q) h^K lambda^g,
R1_alpha=sum_(q=0) count(K,g,q) h^K lambda^g.
```

The pooled root is obtained from the exact numerator

```text
Hq_axis Z_tilted+Hq_tilted Z_axis=0.
```

At lambda=0 this is `2(h^50-1)`, with h derivative100 at h=1. Formal
Newton iteration therefore selects the unique branch h0(0)=1. Its first
nonconstant coefficient is exactly1 at degree2:

```text
h0(lambda)=1+lambda^2-lambda^6+O(lambda^8).
```

The first shift follows directly from the common g=2 shell of25 isolated
occupied sites: `h^25-1=25 h lambda^2+...`. No fitted centering is used.
The pooled numerator vanishes identically through degree32 after the
computed substitution. Each geometry is normalized separately, and the
script differentiates in h **before** substituting h0(lambda). Finally

```text
P1_alpha=R1_alpha/Z_alpha,
Y_h=(P1_tilted,h-P1_axis,h)/DeltaCos4,
D=(Q_axis,h+Q_tilted,h)/2,
U/A=Y_h/D.
```

Thus h0 motion, partition derivatives, and the original D are all retained.
The output includes their full rational coefficients through32, together
with fixed-h1 and root-motion pieces for comparison. The same algebraic
implicit-function theorem guarantees a convergent local analytic germ;
the eventual signs are not a numerical inference from truncated values
at a finite coupling.

## The two precise cancellations before the leading angular response

The raw rank-one shells differ sharply. Axis first contributes at g=9,
with40 configurations and polynomial

```text
f9(h)=10(h^5+h^10+h^15+h^20).
```

These are the minimal straight winding-strip shells: four widths1..4,
five translations and two directions supply all40 saved configurations.
Tilted first contributes at g=13, with12850 configurations, so it cannot cancel or
alter the degree11 coefficient. It enters the next displayed order.

**First cancellation: normalized reciprocal shell.** The g=9 term in
axis P1 is `f9(h)/(1+h^25)`. Since `h^25 f9(1/h)=f9(h)`, this normalized
term is invariant under h->1/h and has zero h derivative at1. Hence
raw rank-one probability already starts at `20 lambda^9`, but its thermal
slope has no degree9 contribution. There is no degree10 slope either.
Taking only the raw derivative f9'(1)=500 would miss this cancellation.

**Second cancellation: root motion against the degree11 fixed-h slope.**
The raw g=11 axis shell has `f11(1)=1350, f11'(1)=16150`. The g=2
partition correction is `z2(h)=25h`. These and the fixed g=9 shell give
the following complete coefficient of `lambda^11` in axis P1_h:

| Contribution | Exact coefficient |
|---|---:|
| g=11 shell: `[f11/(1+h^25)]'` at1 | `-725/2` |
| g=2 partition: `[-f9 z2/(1+h^25)^2]'` at1 | `2875` |
| root shift: `h0_2 [f9/(1+h^25)]''` at1 | `-2500` |
| Sum | `25/2` |

The first two sum to `5025/2` at fixed h=1. Root motion cancels exactly
`200/201` of that coefficient, leaving `1/201`. This is not a choice of
thermal chart or an adjustable subtraction: it is the actual pooled-root
displacement required by Q=0. The cancellation is large but not complete.

Consequently

```text
P1_axis,h(h0)=+(25/2)lambda^11-150lambda^13+...,
P1_tilted,h(h0)=+(15375/2)lambda^13+...,
D(h0)=25/2-25lambda^2+... .
```

At the leading order, the `25/2` in the axis slope cancels the positive
constant `25/2` in D. The projector changes the sign because E=1-P1:

```text
[lambda^11](U/A)=-1/DeltaCos4=-625/1152.
```

The next coefficient also includes the lambda^2 correction to D; replacing
D by its limiting constant would change it. The answer therefore carries
both partition/root geometry and normalization, rather than simply the
smallest winding defect cost. The exponent11 is a property of these two
finite quotients and this original observer, not a universal exponent.

## A zero crossing, negative valley and eventual recovery are necessary

The already published [fixed four-coupling readout](closed-source-global-u-turnover-result.md)
has U(log16)>0. The exact leading term above is negative for all sufficiently
large t. Continuity forces at least one zero crossing at t>log16.
After a crossing into the eventual negative region, U must attain a finite
negative minimum before approaching0; analyticity makes U_t=0 there.
The derivative's positive leading coefficient proves eventual increase
toward0 from the negative side. Neither the zero nor the minimum is located
here, and neither is asserted to be unique.

This supplies an additional mechanism beyond the positive-coupling peak:
strong coupling ultimately leaves a small axis-dominated rank-one slope,
whose original angular sign is negative, after almost complete cancellation
by the moving root. The relatively large coefficient at degree13 is not
permission to estimate a zero by balancing two truncated terms; no
coupling scan, zero search or extension of the frozen four-point decision
has been performed.

The statement is **fixed N25 followed by t->infinity**. There is no exchange
with a thermodynamic limit, new continuum-field identification, or change
to the stopped F4-only source line. Exact integer input hashes, the code
commit, rational series and the decomposition above are retained in
[series.json](../results/p337-closed-source-asymptotic/series.json) and
[coefficients.csv](../results/p337-closed-source-asymptotic/coefficients.csv).

Run with the existing research environment, without installing packages:

```bash
/Users/lc/python-envs/research-py311/bin/python \
  scripts/p337_closed_source_asymptotic.py --order 32 \
  --output-dir /tmp/p337-closed-source-series-new
```

The output path must not already exist. This is a bounded exact-series
readout, with no external services and no general-purpose modeling layer.
