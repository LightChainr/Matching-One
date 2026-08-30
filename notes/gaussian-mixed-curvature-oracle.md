# Exact norm-2/norm-5 mixed-curvature oracle

This certificate resolves the bounded formal-algebra part of parent issue
#138.  Define the commuting mixed difference

`H25[F](N) = F(10N)-F(5N)-F(2N)+F(N)`.

For an affine logarithmic law `F(N)=A+B log N`, every term cancels exactly.
For `F(N)=A+B log N+C(log N)^2`, formal expansion in independent symbols
`log N`, `log 2`, and `log 5` leaves only

`H25[F](N) = 2 C log(2) log(5)`.

This is a mixed-scale curvature: neither a constant nor a rank-2 affine log
cocycle can produce it, whereas a quadratic logarithm does.

For comparison, an ordinary power `Y(N)=N^(-beta)` has normalized response

`H25[Y](N)/Y(N) = (1-2^(-beta))(1-5^(-beta))`.

The checked contract freezes the exact rational values `2/5`, `18/25`, and
`217/250` for integer `beta=1,2,3`.

## Boundary

This is a formal transform oracle only.  It reads no observations and does
not identify a Jordan block, logarithmic rank, exponent, or physical channel.
Those interpretation and data-selection parts of parent issue #138 remain
open.
