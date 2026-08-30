# Orbit-resolved source/sink phase diagram

Status: `exact_orbit_source_sink_phase_diagram`.

The 76/24 split is not birth-dominated. At p_ref each orbit net is only a small residual of two large positive birth/exit currents. Birth and exit have similar axis/diagonal composition, and their small composition skew reverses between N13 and N17, making both orbit nets reverse together while their reinforcing signed share stays close. The signed-share slope itself reverses and roughly doubles, so 76/24 is not a geometry-independent constant.

## N13: gaussian-3-plus-2i

| exact curve | interior roots in (0,1) | locations |
|---|---:|---|
| axis_birth | 0 | none |
| axis_exit | 0 | none |
| diagonal_birth | 0 | none |
| diagonal_exit | 0 | none |
| axis_net | 1 | 0.600766413703138930427450372989 |
| diagonal_net | 1 | 0.559816490590112258706163663977 |
| birth_character_total | 0 | none |
| exit_character_total | 0 | none |
| total_net | 1 | 0.604189609835120171857440404708 |

At `p_ref`, axis signed share = `0.755739917417081006371404207413` with slope `-28.3067163546693481222514871725`.

Activity residual fractions `|birth-exit|/(birth+exit)`: axis `0.0250814086426450741951889288385`, diagonal `0.0772372831586056705241752447749`.

Axis composition: birth `0.913676469271463056475843148820`, exit `0.896080912852695726473173675327`.

Thus the signed-share zero/pole map is: axis zero at `0.600766413703138930427450372989`, diagonal zero at `0.559816490590112258706163663977`, and common pole at `0.604189609835120171857440404708`.

Nearest roots to `p_ref`:

- axis_net: `0.600766413703138930427450372989`; `p_ref-root=-0.00802036291313893042745037298934`; slope `-11.9376375733762853306214721146`.
- diagonal_net: `0.559816490590112258706163663977`; `p_ref-root=0.0329295601998877412938363360234`; slope `-0.935779767199983341635757875766`.
- total_net: `0.604189609835120171857440404708`; `p_ref-root=-0.0114435590451201718574404047083`; slope `-11.1183939571956699259684489046`.

Phase intervals:

- `(0, 0.559816490590112258706163663977)`: axis positive, diagonal positive, total positive; orbit contributions **cancel**.
- `(0.559816490590112258706163663977, 0.600766413703138930427450372989)`: axis positive, diagonal negative, total positive; orbit contributions **reinforce**.
- `(0.600766413703138930427450372989, 0.604189609835120171857440404708)`: axis negative, diagonal negative, total positive; orbit contributions **cancel**.
- `(0.604189609835120171857440404708, 1)`: axis negative, diagonal negative, total negative; orbit contributions **cancel**.

## N17: gaussian-4-plus-1i

| exact curve | interior roots in (0,1) | locations |
|---|---:|---|
| axis_birth | 0 | none |
| axis_exit | 0 | none |
| diagonal_birth | 0 | none |
| diagonal_exit | 0 | none |
| axis_net | 1 | 0.588553491659106509175755689426 |
| diagonal_net | 1 | 0.606810393541698534611658667147 |
| birth_character_total | 0 | none |
| exit_character_total | 0 | none |
| total_net | 1 | 0.586722083251881593899454930032 |

At `p_ref`, axis signed share = `0.764844997919214789577486300717` with slope `55.5339941988900181529328636280`.

Activity residual fractions `|birth-exit|/(birth+exit)`: axis `0.0145781334309210676639707490300`, diagonal `0.0365181005705846071197749362461`.

Axis composition: birth `0.885660450193297078642394926669`, exit `0.895611586626106671319832305845`.

Thus the signed-share zero/pole map is: axis zero at `0.588553491659106509175755689426`, diagonal zero at `0.606810393541698534611658667147`, and common pole at `0.586722083251881593899454930032`.

Nearest roots to `p_ref`:

- axis_net: `0.588553491659106509175755689426`; `p_ref-root=0.00419255913089349082424431057375`; slope `-14.7400148493281370545377173581`.
- diagonal_net: `0.606810393541698534611658667147`; `p_ref-root=-0.0140643427516985346116586671466`; slope `-1.37746164349983798864423886120`.
- total_net: `0.586722083251881593899454930032`; `p_ref-root=0.00602396753811840610054506996756`; slope `-13.4002409713082302649979900745`.

Phase intervals:

- `(0, 0.586722083251881593899454930032)`: axis positive, diagonal positive, total positive; orbit contributions **cancel**.
- `(0.586722083251881593899454930032, 0.588553491659106509175755689426)`: axis positive, diagonal positive, total negative; orbit contributions **cancel**.
- `(0.588553491659106509175755689426, 0.606810393541698534611658667147)`: axis negative, diagonal positive, total negative; orbit contributions **reinforce**.
- `(0.606810393541698534611658667147, 1)`: axis negative, diagonal negative, total negative; orbit contributions **cancel**.

## Next falsifiable prediction

On the next independently chosen two-orbit Gaussian quotient, the axis and diagonal net-flux zeros should remain a close ordered pair; between them the two chi4-weighted orbit contributions reinforce, while outside them they cancel. A quotient that lacks this paired-zero window falsifies the common-activity/counterflow classification.

## Boundary

- All root counts and isolating intervals come from exact rational polynomial arithmetic.
- Decimal root locations and slopes summarize exact rational brackets; they are not fitted values.
- The next-quotient statement is a mechanism prediction, not evidence already supplied by N13/N17.
- No new size enumeration, Monte Carlo sample, Huawei production, PR, or merge is used.
