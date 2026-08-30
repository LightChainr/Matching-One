# P333/P321/P370 generic-Q detach intertwiner gate

The frozen generator is the standard FK detach: split a non-singleton site with coefficient 1; detaching an existing singleton multiplies the state by `Q`. The exact calculation uses `Q=1+epsilon` and `X=X0+epsilon V` modulo `epsilon^2`.

| width | dim V_w | affine Hom jet | + endpoint | + radical Gram | + source | full Gram + source | decision |
|---:|---:|---:|---:|---:|---:|---:|---|
| 2 | 2 | 2 | 0 | 0 | 0 | 0 | `reopened_unique` |
| 3 | 5 | 2 | 0 | empty | empty | empty | `remains_empty` |
| 4 | 14 | 2 | 0 | empty | empty | empty | `remains_empty` |

## Exact decision

The standard generic-Q detach first makes the algebra sharper: its affine-Hom jet has two moduli at every tested width, and endpoint normalization uniquely selects X0=T,V=0. That selected line survives at translation-degenerate width two but fails the radical Gram equation at widths three and four, so the full physical intersection does not reopen. This is also a zeroth-order obstruction: every full-Q jet projects to the already inconsistent join/endpoint/radical-Gram/source system, hence no detach velocity V can repair it. A physical confluence requires a larger marked or direct-sum module, not merely the missing scalar loop weight.

- Width 2: `reopened_unique`; particular radical-velocity rank 0; all tangent radical velocities zero=True; inherited zeroth-order no-go=False.
- Width 3: `remains_empty`; no final velocity exists; inherited zeroth-order no-go=True.
- Width 4: `remains_empty`; no final velocity exists; inherited zeroth-order no-go=True.

The stronger unprojected first-jet Gram equation is reported only as a secondary check. It cannot rescue a failure of the frozen radical condition.

## Boundary

- This is exact rational first-jet algebra at Q=1 in widths 2,3,4 only.
- The detach operator is the minimal FK connectivity detach; no occupied-edge transfer matrix, marked cluster fugacity, or higher Q jet is included.
- A nonempty microscopic intersection would not by itself identify an LCFT field or the formal K of PR #393.
- An empty intersection rules out this declared minimal lift only; a larger marked or direct-sum module can carry additional extension data.
