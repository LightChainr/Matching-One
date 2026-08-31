## Independent F4 → original global U block

Frozen source/scorer/protocol: [0f7a0837](https://github.com/LightChainr/Matching-One/commit/0f7a083770d31095e7b4d688d544637d8fc09658).
[Protocol](https://github.com/LightChainr/Matching-One/blob/0f7a0837/experiments/p337-f4-transmission-20260831/PROTOCOL.md)
and [exact thermal quotient](https://github.com/LightChainr/Matching-One/blob/0f7a0837/notes/plaquette-source-thermal-quotient.md).

This is the one source forced by `Ctot_parent=Ctot_child+F4`, now tested
on the unchanged root/slope-normalized global U. Its one-site density-clock
part cancels exactly; a nonzero result must involve its centered multisite
part. No sign or continuum field identity is assumed.

Fresh N65/85/130/170:20M permutations each,100 paired batches per N,
independent seeds across N. Ordinary/one-full-face ensembles and both
orientations share each permutation. The conditional Bernstein degree is
N−4, not N. Fresh ordinary streams supply every root/jet; no old anchor is
pooled. Four Bonferroni95% coordinates test zero projection; a separate
fixed +/-0.5 band sets the finite-resolution stop. No top-ups or source scans.

**Completed: unresolved, fixed block stopped.** All80M permutations are
archived at`f6006b61`; [frozen score/report](https://github.com/LightChainr/Matching-One/blob/25ca3635/results/p337-f4-transmission-20260831/scored/REPORT.md)
is committed at`25ca3635` with full paired omissions/covariance.

|N|V_F4 ± SE|four-coordinate simultaneous95% interval|
|---:|---:|---:|
|65|0.06489 ±0.24307|[-0.54223,0.67201]|
|85|0.80854 ±0.38152|[-0.14437,1.76145]|
|130|0.04719 ±1.36826|[-3.37032,3.46469]|
|170|-0.73527 ±2.23246|[-6.31130,4.84075]|

Zero projection is NOT_EXCLUDED; the +/-0.5 practical band is also
unresolved. No measured global transmission, negligible-response result,
or sign is claimed. **No top-up or source/size substitution.**

Execution lanes:NePnUn/551oUR/TVVfoB/TgFr7R respectively,14 workers each;
all exited0 and outputs are local, with14–35s producer times. The stopped
P154 lag1 and P334 forecast decisions are unchanged. No cross-team DMs.
