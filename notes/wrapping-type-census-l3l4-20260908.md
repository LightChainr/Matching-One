# #640 wrapping-type census, L=3 and L=4 (seconds; L=5 not in this note)

Exact 4×4 tables live in `results/wrapping-type-census/`. Tripwire: collapsed `1{primal either}-1{matching either}` matches committed Bernstein integers on axis L=3, diamond L=3, axis L=4.

Classifier: `torus_homology.wrapping_channels`. Label `both` is `dir0 ∧ dir1` and therefore **lumps rank-2 cross with rank-1 spirals**. The `cross` bit exists on the same object and is not a column of these tables.

## What the three questions did

1. **"Wraps" is several indicators.** On axis L=3, `dir0` and `dir1` appear as separate blocks (k=3: 3+3; k=4: 18+18; k=5: 36+36). A single wrap bit would hide that. `dir0×dir0 = dir1×dir1` at every k on all three geometries.
2. **Primal and matching wrapping are not exclusive.** The joint table is not supported only on complementary cells. At k=0 the unique config is primal `neither` vs matching `both`. At k=9 it is the reverse.
3. **The L=3 sign change k=5→6 is two-term cancellation, not a single term crossing zero.**
   - k=5: `neither×both = 45` (D=−1) and `both×neither = 9` (D=+1), net `−36`. Same-direction diagonal (`dir0×dir0`, `dir1×dir1`) is D=0.
   - k=6: `both×neither = 36` and a 6-count `both×both`; no `neither×both`. Net `+36`.

## Five-cell support (all three geometries)

The 4×4 is empty except

```text
neither × both     D = −1
both    × neither  D = +1
dir0    × dir0     D =  0
dir1    × dir1     D =  0
both    × both     D =  0
```

No `neither×neither`, no `dir0×dir1`, no `dir×both`. That is the Mertens–Ziff torus pairing (PRE 94, 062152 (2016), §II, eqs. (9)–(11), (19)): no black wrap ⇒ white cross-wraps; single-direction wraps pair; cross-wrapping is exclusive; spiraling counts match. The D-carrying cells are only the exclusive both-axes mismatch, so

```text
a_k = #(both × neither) − #(neither × both)
```

bit-for-bit. Axis L=3 totals: n×b = 259, b×n = 91, dir0 = dir1 = 78, b×b = 6. Axis L=4: 36559 / 9045 / 9406 / 9406 / 1120. Diamond L=3: 153999 / 33739 / 34065 / 34065 / 6276.

`both×both` is spiral×spiral (D=0), not two independent H and V components. Rank-2 sits in `both×neither`. Complement symmetry `table_k[a][b] = table_{N-k}[b][a]` **fails** as soon as single-direction wrap appears (NN vs NN+NNN). That is the generating fact behind `M(p)+M(1-p) ≠ 0`.

This is compatible with PR #646's degenerate both-same picture (x/y D-neutral, both-two empty at these sizes) but does not prove it: the wrapping-type `both` column does not split both-same vs both-two. Axis L=5 (#651) is the first size that can put mass on a sixth cell.

Diamond L=4 (`2^32`) and axis L=5 (`2^25`) are not in this note.

Does not settle #635. Integers only in the JSON tables; floats are absent.
