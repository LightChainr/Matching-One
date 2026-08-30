# Scientific card: N680 lands closest to the frozen two-mode recurrence

- Heldout result: `A_H=-0.00216756 +/- 0.00055693` from 120M fresh
  samples/shape.  The exact H4 child direction is negative:
  `Delta K_A=-0.00345609 +/- 0.00088800` (`-3.892 sigma`).
- Frozen point forecasts: the residual is `-0.586 sigma` to the two-mode
  recurrence, `+1.563 sigma` to the free single transfer, `-2.048 sigma` to the
  frozen nominal single mode, and `+12.110 sigma` to scale-neutral.  The
  preregistered closest forecast is therefore the two-mode recurrence.
- Predictive scores: after source-fit uncertainty, residuals are `-0.295`,
  `+0.747`, `-2.024`, and `+7.159` predictive SE in the same frozen order.
  Thus N680 strongly closes scale-neutral, tensions the single frozen mode,
  and ranks two-mode ahead of free-single without decisively separating the
  first two.
- Power outcome: N680 variance was 1.99 times the N340-scaled projection, so
  the actual two-mode/free-single forecast separation was 2.15 measurement SE,
  not the designed 3.03.  This is a missed power target, not permission to add
  samples or refit the models after reveal.
- Orthogonal control: projective scalar is
  `-0.00031165 +/- 0.00078511` (`-0.397 sigma`), again null.
- Mechanism-space update: the fourth independent generation follows the
  overshoot/return shape predicted by an opposite-sign correction more closely
  than either one persistent transfer or no decay.  The correction identity
  should now be refit with four-generation covariance and one residual degree
  of freedom.
- Does not prove: a unique second eigenvalue or continuum correction exponent.
  The heldout ranking is real; two-mode versus free-single remains unresolved.
- Dependency: all four forecasts were frozen at `ba4ca6f` from `4024a7c`.
  N680 uses a new seed/counter block and full paired batch covariance.
