# Finite-abelian twist tomography of the intrinsic rank source

Status: `exact_finite_abelian_twist_transform`.

## Exact transform

|{alpha in Hom(Z^2,A): alpha|Lambda=0}|=n^(2-rank Lambda).

`S_n=sum_alpha T_alpha=n^2 P0+n P1+P2.`

`S_n/n=Z_top(s=-log n)=n P0+P1+n^(-1) P2.`

Thus order-2 and order-3 twist averages, together with `S_1=1`, invert to:

- `P0 = (S_3-2 S_2+1)/2`
- `P1 = S_2-1-3 P0`
- `P2 = 1-P0-P1`

## Projective refinement

For nonzero alpha in F_q^2, T_alpha=P0+L_ker(alpha), where L_line is the probability of a rank-one primitive winding line reducing to that projective class modulo q.

This converts the finite-field construction from proposed saturation tomography into modular winding-line tomography: saturation is already exact by #269.

## Executable gates

- finite-abelian enumeration rows: 30
- prime projective audits: 4
- prime line-tomography audits: 3
- all gates pass: `True`

## Consequence

- Twist averages at group orders 2 and 3 reconstruct the entire unmarked rank-source functional, with normalization S_1=1.
- Individual prime twists are a modular projective-line tomography of the #334 first/plateau winding direction.
- Integral saturation makes r_q=r for every prime: finite-field ranks cannot reveal an additional Smith/index state on these carriers.
- The aggregate transform depends only on |A|, while individual twist sectors retain line-incidence information.

## Boundary

- The result is an exact lattice/cohomology transform, not yet an identification with a local CFT field.
- A finite set of primes resolves winding lines only modulo those primes; it does not reconstruct an unbounded integral line without an external size bound.
- The projective refinement concerns the rank-one carrier line, not a new saturation index, because #269 proves saturation.
