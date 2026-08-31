# The conditional clock shift has a direct/collective source allocation

This joins `3d760b86` C/L with `32270fa2` microscopic source means, both from
the same completed `0d1e586d` archive. Only their stored 20 batch vectors were
read: no raw reprocessing, network solve, or new simulation.

For each source `s`, use **the original orientation-specific R1 risk**:
`m_i,s=mean(Y_i,s)/r_i`, never a source-specific occurrence frequency. Apply
the same symmetric C/L identity. Then, within every original-batch LOO,

`sum_s C_s=C`, `sum_s L_s=L`, and `sum_s D_s=D`.

The source labels are original-checkpoint `H2` direct triggers, collective
completion (including triggers created by earlier safe insertions), and
unclassified original `Y` on whole-pair fallbacks. All fallback observations
remain in their original denominator and pairing.

| N / endpoint | L_direct ± SE | L_collective ± SE | L_unclassified ± SE |
|---|---:|---:|---:|
| 325 / canonical | +0.00114380 ± 0.00050986 | −0.00020310 ± 0.00018303 | +0.00002836 ± 0.00006657 |
| 425 / canonical | −0.00054200 ± 0.00055090 | −0.00004845 ± 0.00021125 | +0.00003690 ± 0.00015761 |
| 325 / integral | +0.00152158 ± 0.00074444 | −0.00131605 ± 0.00070585 | −0.00003820 ± 0.00012657 |
| 425 / integral | −0.00066690 ± 0.00085132 | +0.00051672 ± 0.00074444 | +0.00004470 ± 0.00019993 |

At N325 the canonical `L=+0.00096906 ± 0.00039612` points mainly toward the
original-direct source, with a smaller opposing collective contribution.
For the integral, larger opposing source point estimates leave a much smaller
net loading, `+0.00016733 ± 0.00008074`. Their errors are strongly dependent;
one cannot add source SEs in quadrature or call these source signs independently
established. N425 does not resolve the individual source signs either.

This is consistent with competition in how a completion is assigned to an
original direct gate versus collective paths; it is not a count of distinct
fields. The exact prefix-level source-integral identities derived in the
parallel theory note can interpret this allocation without another data run.

The result retains every source's weighted `Y`, conditional `m`, C/L/D values
and common LOO vectors. A single 45-coordinate covariance joins the earlier
27-coordinate C/L-plus-four-state block with the three sources' 18 C/L/D
readouts. Its rank is at most 19; it is **never inverted**, and no additional
omnibus or independent-source significance test is performed.

Scientific card: this links the mean conditional-clock loading to named
microscopic sources while preserving the covariance responsible for their
small total. It does not alter the main 99.85% four-state variance finding,
which concerns the gated R1 integral only, not full `A_top` or global variance.
The next full-observable step needs the R0/R1/R2 joint terms, not extrapolation
of this restricted layer.

Outputs: `results/p334-source-loading-crosswalk/{score.json,REPORT.md}`.
Reproduce: `python3 scripts/p334_source_loading_crosswalk.py`.
