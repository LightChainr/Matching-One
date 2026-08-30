# P155 general-period thermal-null pilot

The exact microscopic counterterm does **not** expose a conditioned second
response direction at the first injective `R=8` Gaussian checkerboard sizes.
The frozen 20k pilot therefore stops without the allowed 100k expansion.

Rows below are `(global cross, O_alpha*)`; columns are `(t,lambda)`, with
`O_alpha*=O_local_H4+(3/64)epsilon_cell`.

```text
N260, R8 = [[6.17370, 4.07470],
            [0.07729, 0.00412]]
lambda z = 0.218, determinant z = -1.984
singular values = (7.39745, 0.03914), condition = 189.0

N340, R8 = [[6.64540, 4.41160],
            [0.03449,-0.00204]]
lambda z = -0.088, determinant z = -0.868
singular values = (7.97649, 0.02078), condition = 383.9
```

Neither size passes the frozen `|lambda z|>=3`, `|det z|>=3`,
`condition<=50` gate.  Even a `sqrt(5)` 100k variance projection cannot bring
the lambda gate close to three, and both central condition numbers already
fail.  No expansion was run.

The radial diagnostics sharpen the negative result.  At R4 the lambda entry is
individually nonzero at both sizes (`z=4.20,3.47`), but the full rows are nearly
thermal-parallel: condition numbers are `4015` and `302`, with determinant z
`-0.085` and `-1.27`.  At R2 neither lambda response is resolved and condition
numbers remain above `184`.  Thus the N10 exact staggered matrix element is not
simply missing because every finite-R local response vanished; rather, the
surviving intermediate-radius response lies on the already-resolved thermal
line and then becomes compatible with zero by injective R8.

This is evidence that the frozen microscopic zero is UV/landing-scale
sensitive for this readout.  It is not an exclusion of a distinct odd RG
direction: a different local operator could still couple to one.  No
generalized eigensystem, exponent or CFT identity is reported.

Production ran on Huawei HZsCM6 (`aarch64`, 16 cores, GCC 10.3.1), commit
`0b25b15`, seed `15583020260830`, counters
`[15500000000,15500020000)`.  Kernel sampling time was 0.202 seconds.

