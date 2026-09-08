# #640 wrapping-type census, L=3 and L=4 (seconds; L=5 not in this note)

Exact 4×4 tables live in `results/wrapping-type-census/`. Tripwire: collapsed `1{primal either}-1{matching either}` matches committed Bernstein integers on axis L=3, diamond L=3, axis L=4.

## What the three questions did

1. **"Wraps" is several indicators.** On axis L=3, `dir0` and `dir1` appear as separate blocks (k=3: 3+3; k=4: 18+18; k=5: 36+36). A single wrap bit would hide that.
2. **Primal and matching wrapping are not exclusive.** The joint table is not supported only on complementary cells. At k=0 the unique config is primal `neither` vs matching `both`. At k=9 it is the reverse.
3. **The L=3 sign change k=5→6 is two-term cancellation, not a single term crossing zero.**
   - k=5: `neither×both = 45` (D=−1) and `both×neither = 9` (D=+1), net `−36`. Same-direction diagonal (`dir0×dir0`, `dir1×dir1`) is D=0.
   - k=6: `both×neither = 36` and a 6-count `both×both`; no `neither×both`. Net `+36`.

Diamond L=4 (`2^32`) and axis L=5 (`2^25`) are not in this note.

Does not settle #635. Integers only in the JSON tables; floats are absent.
