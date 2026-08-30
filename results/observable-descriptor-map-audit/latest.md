# Observable descriptor map audit

This is a complete audit of the finite map registry, not a scan of every repository artifact.

| item | exact count |
|---|---:|
| valid descriptors | 200 |
| ordered descriptor pairs | 40000 |
| registered exact maps | 952 |
| fail-closed pairs | 39048 |
| inverse checks | 952 |
| composable paths | 5720 |
| inverse/composition failures | 0 |

Connected-component size histogram: `{'1': 24, '2': 8, '4': 24, '8': 8}`.

Blocked reasons: `{'cannot map a scalar value to an orientation contrast': 12800, 'no exact topology map': 26248}`.

Registered-edge SHA-256: `dfc984a89796bb1e6ff5bf0344554c8e7ad37451c7ac90983a8362f05bf438fc`.

## Interpretation boundary

This exhausts the current finite descriptor registry only. It does not scan all scorers or artifacts, register a new topology identity, or prove repository-wide descriptor adoption.
