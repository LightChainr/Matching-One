# Claim-wording erratum for the r=1 bounded-summary package

Date: 2026-09-02
Scope: additive correction on `research/p429-r1-noncompression-20260902`.
Does not reopen the search, does not change any rational, incidence list,
or census count.

## Checklist

| Item | Was | Now |
|---|---|---|
| Manuscript verdict token | `NO_COMPRESSION_WITNESS_FOUND` | `BOUNDED_SUMMARY_INSUFFICIENT` |
| Frozen search-protocol token | `NO_COMPRESSION_WITNESS_FOUND` | **unchanged** in `notes/p429-r1-search-contract-20260902.md` (locked before search) |
| Size claim | “n=7 remains smallest in n” / “minimal in n” | “smallest witness found in the declared enumerated families” |
| Global minimality among all plane TT networks | implied by the size claim | **not claimed** |
| `verify_witness.py` | “independent” check, easy to read as a second implementation | search-independent hard-coded witness verifier, using the same stdlib primitives as the search library; **not** a fully independent implementation |
| Cut-network minimality | “not proved” | `UNRESOLVED` (same meaning, explicit token) |
| Embedding | `TWO_PORT_EMBEDDING_SUFFICES`, citing parallel-gadget §6 for all plane TT gadgets | `GENERAL_REALIZATION_LEMMA`; §6 is not a surjectivity proof |
| Named HNF occupation of the n=7 pair | not constructed | still not constructed |

## Why the verdict token changed

The frozen protocol allowed only two search tokens:
`NO_COMPRESSION_WITNESS_FOUND` and `BOUNDED_SEARCH_CLOSES_UNDER_SUMMARY`.
The first records that a declared-class split was found. Manuscript wording
uses `BOUNDED_SUMMARY_INSUFFICIENT` for the same mathematical outcome: the
predeclared tuple `(S(z), n, H2, b2, r=1 neighbourhood)` is not sufficient
for the frozen depth-2 language on the tested class. No new split is
involved.

Machine JSON keeps `search_protocol_token` for provenance and sets
`verdict` to the manuscript token.

## Files touched by this erratum

- `notes/p1-plane-tt-realization-lemma-20260902.md` (new)
- `notes/p1-n7-torus-embedding-20260902.md`
- `notes/p429-r1-claim-boundary-20260902.md`
- `notes/p429-r1-noncompression-certificate-20260902.md`
- `notes/p1-theorem-scope-skeleton.md`
- `notes/p429-r1-search-contract-20260902.md` (provenance header only)
- `results/p429-r1-bounded-summary/*`
- `research/summary_search/verify_witness.py` (docstring)
- `tests/test_p429_r1_noncompression_witness.py` (docstring)
- generator scripts that emit the verdict string (no hunt re-run)

## Still forbidden, unchanged

- r=2 non-compression claim
- Euclidean / scalar-encoding / continuum / CFT statements
- all-graphs minimality of the cut network
- “every HNF torus realises this pair”
- new graph search, new descriptor, Monte Carlo
