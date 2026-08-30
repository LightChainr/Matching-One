# P250 radius-six raw lock

The authorized HZsCM6 production completed before any scientific score was
opened.  The copied response, 400-row batch archive, exact gate, timestamps,
and logs match the remote files byte for byte under SHA256.

```text
host: Huawei-CodeBuddy-HZsCM6
runner/scorer: 44a3a36
authorization: 18d188d
samples: 1,200,000
batches: 400 x 3,000
seed: 25060610120261250
counters: [0, 1,200,000)
started: 2026-08-30T09:53:16Z
completed: 2026-08-30T10:06:24Z
elapsed: 788 s
exit: 0
stderr: empty
```

The response reports the frozen schema, runner commit, 13-point geometry,
28 complex coordinates per batch, and the exact authorized run domain.  The
CSV contains exactly 400 consecutive batches; the first begins at counter 0,
the last begins at 1,197,000, and the sample sum is 1,200,000.

The archive also preserves an earlier zero-second protocol rejection.  That
invocation used the runner defaults rather than the authorized 1.2M argument,
so `validate_manifest` rejected it with `run differs from manifest for
samples`.  It exited 1 before generating any scientific samples.  It is kept
as evidence that the execution gate worked, not pooled with production.

At this commit `score_status=NOT_RUN_AT_RAW_LOCK`.  The new batch hash is
`8734a15e...a42cf`; the response hash is `b061d6dc...f4e33`.  The next and only
allowed scientific action is the scorer frozen at `44a3a36`, with alpha,
pivots, rank ladder, and bridge lock unchanged.
