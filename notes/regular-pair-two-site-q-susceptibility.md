# Q activation is already interacting on a fixed four-path occupation

The [two-site lattice witness](local-pair-two-insertion-geometry.md)
does more than expose the pole of the uncompleted pair tensor. Its
canonical regular completion has a finite, explicitly nonadditive
Q-activated colour susceptibility. This statement concerns a fixed
occupation's marked colour contraction, not an occupation-averaged
connected correlator or a continuum exponent.

## 1. The completed two-insertion kernel

Let `Kreg=K2bar+K0bar` with the singlet coefficient identically one.
Summing exact colour-equality patterns gives

```
<K2bar,K2bar> = Q(Q-3)(3Q^2-9Q+8) / [8(Q-1)(Q-2)],
<K0bar,K0bar> = (2Q^2-4Q+3) / [2Q(Q-1)],
<K2bar,K0bar> = (Q-3) / [4(Q-1)].
```

The three Q1 residues are respectively `1/2,1/2,-1/2`.
The cross term occurs twice, so the pole cancels exactly. The result is

```
Hreg(Q) = ||Kreg||_F^2
        = (Q-1)(3Q^3-12Q^2+20Q-24) / [8Q(Q-2)].
Hreg(1)=0,                 Hreg'(1)=13/8.
```

One can read the derivative directly from colour-equality patterns:
two-pair assignments contribute `3/8`, one-pair assignments `-27/4`,
and four-distinct assignments `8`, totaling `13/8`. These are rational
diagram-continuation contributions, not probabilities at Q1.

## 2. An exact nonlinear Q source, with no fitted coefficient

Retain the fixed occupied set of the 8x8 witness. Replace only its two
vacant marked vertices by `1+lambda_x*Kreg` and `1+lambda_y*Kreg`.
The unmarked colour weight is `Q^74`; 70 exterior components are
spectators, and the two marks share the same four component colours.
Normalize by the unmarked fixed-occupation weight. The exact ratio is

```
F_A(Q;lambda_x,lambda_y)
 = 1 + b(Q)(lambda_x+lambda_y) + c(Q)lambda_x*lambda_y,
b(Q)=(Q-1)/Q^3,            c(Q)=Hreg(Q)/Q^4.
```

Thus `F_A(1;lambda_x,lambda_y)=1` for every finite lambda, but

```
Xi_A(lambda_x,lambda_y)
 = d_logQ log F_A(Q;lambda_x,lambda_y)|Q1
 = lambda_x + lambda_y + (13/8)lambda_x*lambda_y.          (1)
```

In particular `d_lambda_x d_lambda_y Xi_A=13/8`. This is the mixed
two-mark coefficient of the **log of the conditional colour weight**.
No sampling, activity value, field fit or assumption about global U was
used. The number 13/8 here is a finite four-port contraction coefficient;
it is not the exponent in the project's N-dependent normalization.

The site-average convention used by the N25 first-response experiment
would set `lambda_x=lambda_y=epsilon/N`, giving a contribution
`(13/8)*epsilon^2/N^2` from this unordered pair. Differentiating twice
with respect to the common epsilon produces `13/(4N^2)`. Keeping the
independent-lambda convention in (1) avoids a factorial ambiguity.

## 3. Which local-mechanism simplification this excludes

Multiplying two independent single-mark closed weights would give

```
[1+b(Q)lambda_x]*[1+b(Q)lambda_y].
```

Its two-mark coefficient is `b(Q)^2=O((Q-1)^2)`, and its first Q
derivative at Q1 vanishes. The actual shared-component contraction
instead has `c(Q)=(13/8)(Q-1)+O((Q-1)^2)`.
Hence replacing the completed colour interaction by independent
first-insertion factors loses a nonzero term already at first order
in Q-1. Equivalently, a purely additive activated occupation source
captures the first lambda tangent but cannot represent the complete
two-mark Q-activated family, even on this explicit original-lattice
configuration. The missing coefficient is fixed by the colour sewing,
not by a new contact regression.

This exclusion is coefficientwise and uses shared exterior components.
It does not claim long-range propagation, a universal two-point
amplitude, a nonzero summed homogeneous mixed response, or identification
of a logarithmic CFT field. Those questions involve further contractions
and occupation sums. The [single original-U activation calculation](regular-pair-activation-original-u.md)
retains its specified first-source scope and does not approximate (1)
by the square of a site-average mark.
