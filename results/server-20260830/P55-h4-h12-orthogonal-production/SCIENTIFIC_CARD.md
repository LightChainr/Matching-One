# Scientific card

**Question:** Does an H12 completion survive the two exact opposite-alias
Gaussian rows?

**Acquisition:** 600M paired replicas per size, 300 batches, one frozen source
commit and exact contiguous counter coverage across three hosts.

**Result:** H4-only p `0.329`; zero p `0.425`; `A12=0.147 +/- 0.191`,
z `0.771`, p `0.441`.

**Meaning:** H12 is unresolved.  H4 remains compatible, but this target does
not establish it because zero remains compatible as well.

**Stop:** do not run the held-out alias; first strengthen the independently
motivated matching-odd source coordinate.
