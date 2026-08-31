# Thermal dipole, center displacement and birth-clock spread

This analysis uses the original common-label fork block. It combines the final
thermal moment readout `9059776d` with the forthcoming full-population continuous
birth-clock moment readout. It does not sample new paths or choose new extrema.

Let the two continuous uniform-order-statistic birth times be tau1<=tau2,
`C=(tau1+tau2)/2`, and `W=tau2-tau1`. All unperturbed means are full-population
means in one physical orientation. H is the derivative under the common label
policy and one of the fixed plus/minus marks.

## An exact three-part dipole

The rank-one interval gives

```
I0 = integral H(E(p)) dp       = -H(W)
I1 = integral p H(E(p)) dp     = -H(CW)
dipole_at_pref = I1-pref I0
               = -muW H(C) -(muC-pref) H(W) -H(Cov(C,W)).
```

Thus the centered first thermal moment separates into a center-displacement
term, a mean-lifetime term and a connected spread term. Their errors must be
propagated jointly; subtracting separately reported errors is invalid.

The moment derivatives are

```
H(Cov(C,W)) = H(CW)-muC H(W)-muW H(C)
H(Var C)   = H(C^2)-2 muC H(C)
H(Var W)   = H(W^2)-2 muW H(W).
```

Compute these products separately for first and second orientation. Only then
form S=(first+second)/2 and D=(first-second)/delta_cos4. In every delete-one-batch
replicate, pool the retained raw means before recomputing the products. This
preserves all shared source and mean uncertainty.

## Marginal spread versus joint fluctuation

The apparently connected C/W quantity obeys a further exact identity:

```
Cov(C,W) = (Var(tau2)-Var(tau1))/2.
```

Its derivative measures unequal broadening of the two birth marginals. It
does not, by itself, identify information in their copula. The distinct joint
direction is

```
Cov(tau1,tau2) = Var C-Var W/4.
```

Both are recoverable from the same five supplied full-population moments.
If v1=Var(tau1) and v2=Var(tau2) remain positive in all original-batch deletions,
retain the dimensionless derivative

```
H(rho12) = H(cov12)/sqrt(v1 v2)
           -rho12/2 [H(v1)/v1+H(v2)/v2].
```

This is a normalized response coordinate, not a fitted model or inferred field
count. Its source parameter remains the original common-label mark strength.

## Continuous moments

Given integer ranks k1<=k2 and n1=N+1, n2=N+2, the moment reader integrates over
the common uniform order-statistic clock:

```
E(C|ranks)   = (k1+k2)/(2 n1)
E(W|ranks)   = (k2-k1)/n1
E(CW|ranks)  = [k2(k2+1)-k1(k1+1)]/(2 n1 n2)
E(C^2|ranks) = [(k1+k2)^2+3 k1+k2]/(4 n1 n2)
E(W^2|ranks) = (k2-k1)(k2-k1+1)/(n1 n2).
```

These are not the squares of the conditional mean clock times. The reader
supplies them once; this coordinator consumes only its final committed batch
vectors and the saved thermal moments.

Status: formulas fixed before receiving the final baseline-moment archive.
