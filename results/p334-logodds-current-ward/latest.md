# Log-odds projective-current Ward identity

The exact finite-volume identity is

`dA_O/deta = p(1-p) (J_O,birth - J_O,exit)`.

It holds as a polynomial for every recorded line orbit and through the complete
eta jet of order 6 at the frozen p_ref.

## gaussian-3-plus-2i (N=13)

- all exact gates: `True`
- empty/full rank-one amplitude vanishes: `True`
- integrated net-current sum rule: `True`
- exact unique orbit/total stationary points: `True`
- coordinate-free shares at p_ref: axis_orbit=0.755739917417081006, diagonal_orbit=0.244260082582918994

## gaussian-4-plus-1i (N=17)

- all exact gates: `True`
- empty/full rank-one amplitude vanishes: `True`
- integrated net-current sum rule: `True`
- exact unique orbit/total stationary points: `True`
- coordinate-free shares at p_ref: axis_orbit=0.764844997919214790, diagonal_orbit=0.235155002080785210

## Consequence

A zero of an orbit-resolved birth-minus-exit current is a stationary point of that orbit's finite H4 character amplitude. A zero of the total current is the corresponding total-amplitude stationary point and a pole of signed orbit shares. These locations survive every regular scalar reparameterization of the thermal coordinate.

## Boundary

- This is an exact finite-volume Ward/continuity identity, not a continuum Ward identity.
- It upgrades orbit shares and net-current zeros to coordinate-free finite observables; it does not assert their asymptotic limits.
- The N13 and N17 coefficient tables are reused without Monte Carlo or path enumeration.
