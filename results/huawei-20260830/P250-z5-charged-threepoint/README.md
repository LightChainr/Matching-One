# P250 N325 Z5 charged three-point production

## Acquisition and integrity

- Scientific freeze: `95d170e`; execution-only partition: `5909987`;
  runner commit recorded by every shard: `846e000`.
- Common contract: seed `25011312220260831`, p `0.59274605079`, radius 1,
  exact transported Z5 generator, one field shared by both hands.
- Replica coverage: Zy `[0,340000)`, XP `[340000,670000)`, HZ
  `[670000,1000000)`; 34+33+33 retained 10k batches.
- All three remote stderr files are empty.  Remote/local SHA256 checks pass
  for every response, batch table, provenance file, stdout/stderr and exact
  gate.  The three exact gates have the same SHA256 `54b9111859...`.
- The merge oracle verified exact disjoint coverage and recomputed the full
  primary 8x8 and joint 24x24 covariance from all 100 batch rows.

## Frozen reveal

The zero-parameter cross-product closure was scored first:

`C113_plus*C122_minus-C113_minus*C122_plus =
(1.287 + 0.808 i) 1e-13`, chi-square `0.840/2`, p `0.657`.

The conditional single-field cubic templates then gave:

| frozen model | chi-square/df | p |
|---|---:|---:|
| H12 | 1.744/4 | 0.783 |
| H8 | 3.103/4 | 0.541 |
| H4 | 5.214/4 | 0.266 |

The nonneutral joint-zero control gives `10.595/8`, p `0.226`; maximum
configurationwise DFT conjugacy residual is `2.48e-17`.

All three H4/H8/H12 targets survive.  H12 ranks first, but this is not a
field identification.  A post-reveal, explicitly non-primary signal check
puts the complete eight-real primary vector at chi-square `8.564/8`, p
`0.380` against zero.  Thus the closure survival is presently compatible
with a nearly null charged cubic vector; it does not yet demonstrate a
nonzero common fusion eigenphase.  The one-million-replica result identifies
the current compact three-anchor local-H4 cubic as signal-limited, not a need
to choose a fourth phase model.

## Five-line scientific card

1. **Question:** do the two primitive neutral Z5 cubic channels share one
   handed eigenphase, and does that phase select H4, H8 or H12?
2. **Observable:** common-field, three-anchor, all-five-fiber local-H4 DFT;
   C113/C122 in both norm-five hands with conjugate and nonneutral controls.
3. **Result:** closure survives (`p=0.657`); H12/H8/H4 all survive joint GLS
   (`p=0.783/0.541/0.266`); exact/control gates pass.
4. **Interpretation:** the eight-real cubic vector itself is not detected
   (`p=0.380` versus zero), so model ranking is not field identification.
5. **Next selector:** do not simply add replicas to this compact cubic row;
   freeze a charged two-point-normalized/separation-aware observable that can
   establish nonzero cubic support before reusing the phase templates.

The 20k smoke is not combined with production.  No phase, anchor, p, radius,
or score was changed after seeing either smoke or production data.
