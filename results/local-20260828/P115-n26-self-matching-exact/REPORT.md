# Issue #115: N=26 exact C4 self-matching target

## Protocol order

The two N=26 targets were generated without target coefficients and pushed in
commit `1add3384598da2ce749a429af3b55c7a4ebc20a0` before exhaustive
enumeration.  The canonical channel is `either`.

- Geometry-only shortest nontrivial support: `s=5`, with an explicit length-5
  winding certificate.  Frozen target: `M(p)=2 I_p(5,5)-1`; vector SHA256
  `df1837531d6a6b8296baf1422b23c19b5c5b1dbe70eea4bc180c816f5a0c0cad`.
- Thirteen antipodal two-site orbits.  Frozen majority target:
  `M(p)=2 I_p(7,7)-1`; vector SHA256
  `6fe1e48cdcfa78ef9c7df5a1d6abc00a7eb780bbfd26cf642a0dd08211bcf3cb`.

The C++ kernel was committed and pushed as `68f6d81daea30475b58bd109625e044260d7f547`
before the authoritative runs.  Its N=10 coefficients agree exactly with the
independent Python oracle for all five wrapping channels.

## Frozen scores

Both fixed laws fail coefficient-by-coefficient.  No generalized-Beta fit was
performed.

| Frozen law | First difference | Number of different coefficients | Result |
|---|---:|---:|---|
| geometry `Beta(5,5)` | `k=5`: observed `-65624`, target `-65528`, difference `-96` | 16 | fail |
| antipodal majority `Beta(7,7)` | `k=5`: observed `-65624`, target `-65780`, difference `+156` | 16 | fail |

The exact N=26 integer Bernstein vector, for `k=0,...,26`, is

```text
[-1,-26,-325,-2600,-14950,-65624,-227292,-631540,-1414335,
 -2536768,-3567343,-3729856,-2444442,0,2444442,3729856,
 3567343,2536768,1414335,631540,227292,65624,14950,2600,325,26,1]
```

Its canonical compact SHA256 is
`32930e6e63c2d0ee22070390068d1ca86ec4a2465da43154e05febe4296dfdc8`.

## Exact structure after scoring

The polynomial is anti-palindromic and identical across `cross`, `both`,
`either`, `direction_0`, and `direction_1`, even though the five individual
wrapping-count rows differ.  In the power basis,

```text
M(p) = -1 + 156 p^5 - 338 p^6 + 260 p^7 - 260 p^8 - 338 p^9
       + 1144 p^10 + 3536 p^11 - 13702 p^12 + 15628 p^13
       - 3016 p^14 - 10088 p^15 + 11492 p^16 - 5798 p^17
       + 1482 p^18 - 156 p^19.
```

SymPy 1.14.0 gives

```text
M(p) = -(2p-1) Q18(p),
```

where the explicit degree-18 factor is retained in `score.json`.  For
`F=(1+M)/2`, the first nonzero term is `78 p^5`; hence `F` has the exact
factor `p^5`, and central antisymmetry gives the corresponding `(1-p)^5`
tail factor for `1-F`.

At occupation five, the N-site Bernstein numerator of `F` is `78`.  The
geometry `Beta(5,5)` law predicts `126`, while antipodal-majority `Beta(7,7)`
predicts `0`.  The raw graph has 156 `either`-wrapping five-site masks, split
as 78 in each quotient direction.  This is the earliest exact obstruction to
both proposed mechanisms; it is descriptive structure, not a refitted target.

The two discrepancy polynomials both contain the forced boundary/symmetry
factor `p^5 (1-p)^5 (2p-1)`; their remaining degree-8 factors are recorded
exactly in `score.json`.

## Independent reproduction

The authoritative 10-thread run exhausted all `67,108,864` masks in 2.55
kernel seconds (3.01 wall seconds).  A separate single-thread run exhausted
the same mask space in 15.65 kernel seconds (15.89 wall seconds).  After
removing only `threads` and `elapsed_seconds`, every JSON field and every
integer coefficient is identical.

See `commands.txt`, `environment.txt`, `metadata.json`, and
`checksums.sha256` for the exact provenance.

## Validation

- N=10 Python/C++ all-channel oracle and thread-invariance tests: 2/2.
- Geometry freeze, N=10 Beta(3,3), fixed-score, and reproduction tests: 8/8.
- Repository-wide local suite: 165/166.  The sole failure is the unchanged
  `test_exact_matching_zero_map` high-precision root audit (`ROOT_FAILURE`);
  neither that script nor its test differs from `origin/main` in this branch.
