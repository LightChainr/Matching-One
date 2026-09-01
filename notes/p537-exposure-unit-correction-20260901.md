# P537 exposure-unit correction and mechanism-language boundary

**Status.** This is a fixed-blob, existing-data correction sidecar.  It does
not replace the frozen N25 or N65 results, create an independent evidence
block, or authorize a new simulation.  Priority is attention allocation, not
permission or a task lock.

## Confirmed source-coordinate mismatch

The fixed N25 scorer at `df4a64f68232eec5aa5b8c8a5d920062aaa7808e`
accumulates selected positive exposure as

```python
positive[i] += weight * wi * count * 2
```

while its signed source average has a separate `1/N` in
`sum_a16/(N*16)`.  Thus the reported positive exposure coordinate is `P`.

The fixed N65 scorer at
`bab37f21b061afa9e03286b6e9bf4560f879ecd5` accumulates

```python
exposure[st,mask-1] += 2*f*w*v[0]/N
```

and therefore stores the same positive mass as `P/N`.  The historical
`scale-fingerprint.json` directly formed `density25=n25/e25` and
`density65=n65/e65`.  Its cross-size exposure and conditional-density
decomposition mixed `P` with `P/N`.

`scripts/audit_p537_exposure_units.py` reads these immutable Git blobs,
verifies their Git blob SHA-1 and the required code markers, and then consumes
the historical scale result without changing it.

## Corrected finite-pair coordinates

For the `01`/double-contact cell, the signed masses remain

```text
K25 = -2.9372878646696404e-6
K65 = -9.227657800496906e-8
```

Using one common unweighted positive-mass coordinate gives

| quantity | N25 | N65 | two-point decay |
|:--|--:|--:|--:|
| exposure `P` | 0.03269409005269323 | 0.017724830091433354 | 0.6407331339 |
| conditional signed density `K/P` | -8.984155423612032e-5 | -5.206062767821259e-6 | 2.9808368721 |
| signed mass `K` | -2.9372878646696404e-6 | -9.227657800496906e-8 | 3.6215700060 |

The closure is

```text
0.6407331339 + 2.9808368721 = 3.6215700060.
```

Equivalently, using `P/N` at both sizes gives exposure decay
`1.6407331339` and density decay `1.9808368721`; the same signed-mass
decay is recovered.  The coordinate must be declared, but neither convention
is promoted to a universal exponent from two finite, non-dilation-equivalent
geometries.

## Frozen results that do not change

The correction changes only the positive-exposure coordinate and the density
formed from it.  It leaves intact:

- the N25 and N65 signed stage-by-contact matrices;
- the N65 determinant `-8.688216055121765e-14` and frozen decision
  `CONTACT_FUSION_COMPLETION_TRANSMITS`;
- the complete branch-only result
  `J65=-0.0016225098893862522 +/- 0.00018553008242164315` at
  `f9ba1ff690b07beefcc71e669f1f29581d4e264e`;
- the historical post-reveal shape score as a numerical diagnostic.

No raw TSV, baseline/source delete-one factors, or new configurations are
read by this sidecar.

## Interpretation corrections

### The additional `5/8`

The old `[3,29/8,3,3]` shape gives `Q=0.6364356783/4` and nominal
`p=.9589300602`, but it was constructed after the N65 reveal.  Once the
exposure units are aligned, the previous claim that the extra `5/8` is a
suppression in conditional signed strength is withdrawn.  The same two points
also do not identify it as an exposure-frequency or continuum-field exponent.

### The determinant and commutator name

A nonzero determinant establishes that the signed `stage x contact` table is
not a rank-one product.  Writing a raising matrix and solving a diagonal
operator from that same observed table is an algebraic commutator encoding,
not an independent identification of noncommuting microscopic operations.
Any physical commutator claim must define `F`, `B`, `FB`, and `BF` before
using the table and predict a quantity not used to solve those operations.

### The unselected response

Summing the selected N65 cells and converting with the fixed full-T factor
gives a point share of `2.550516%` of full `J65`.  The remaining point value
is `-0.001581127510162752`, but its only licensed name is
`complement_of_selected`.  It is not yet a spatially nonlocal mechanism,
unexplained causal fraction, or gauge-invariant operator loading.  Exact share
uncertainty still requires the shared full/selected delete-one factors.

## Reproduction

```bash
python3 scripts/audit_p537_exposure_units.py \
  --manifest analysis/p537_exposure_unit_audit_manifest.yaml \
  --output-json results/p537-exposure-unit-audit/latest.json \
  --output-md results/p537-exposure-unit-audit/latest.md
```

The generated result schema is `matching-one.p537-exposure-unit-audit.v1`.
The audit uses fixed Git blobs only, produces no new Monte Carlo evidence, and
does not overwrite any historical artifact.
