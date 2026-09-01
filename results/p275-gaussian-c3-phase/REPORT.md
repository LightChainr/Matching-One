# Paired Gaussian real-C3 phase gate

Status: `H8_SELECTED_H4_STOP` at the frozen sample count; no top-up.

The production used two equal-area Gaussian quotients,

```text
g1 = 8+i,       g2 = 7+4i,
N1 = N2 = 65,   delta = atan2(5,12) = 22.619864948 degrees,
```

and applied the same counter-derived 130-bit field to both deterministic bond
orders for each of 2,000,000 paired replicas.  Only the 100 batch means of
`(Re z1, Im z1, Re z2, Im z2)` were retained.  The common intrinsic `tau=i`
Pinson--Arguin baseline was subtracted before scoring.

## Frozen result

```text
mean z1 = +0.00208049662712 - 0.00262405113148 i
mean z2 = -0.00110625337288 + 0.00219494722788 i

H4: chi2 = 73.6412 / 1 df, p = 9.36904e-18, signed gain = +0.76842
H8: chi2 =  1.1122 / 1 df, p = 0.291603,   signed gain = +0.73660
```

At the frozen `alpha=0.01`, H4 fails and H8 survives, so the preregistered
decision is

```text
H8_SELECTED_H4_STOP
```

No additional replicas are authorized or drawn.

The raw projective ratio is

```text
z2/z1 = -0.7188453 + 0.1483591 i,
|z2/z1| = 0.7339952,
arg(z2/z1) = 168.3387 degrees.
```

Modulo the sign allowed by the real-gain contract, its angular residual is
`-10.70 degrees` from the H8 phase and `+77.86 degrees` from the H4 phase.

## Scientific meaning

This is a positive sector split, not a reversal of the ordinary global-channel
H4 evidence.  The phase-calibrated primitive real-C3 observer selects its H8
alias under a physical equal-modulus rotation, while the quotient-prism/global
`A_top` channel can still select H4.  Consequently a real-C3 response must not
be used as another vote for the field already identified in the global
channel: the observer projection changes the surviving harmonic sector.

The result identifies the frozen finite signed-real transport law.  It does
not identify a continuum primary, prove an OPE coupling, or apply to E_top,
rho-child characters, or original square-site U without an explicit map.

Input/output hashes:

```text
edaa8647f6e3e642e008d389fab6d01c03fbd6042a1a3f6aaa2c76f4923f937a  paired-batches.csv
09245833bae30897ba488219bc4d7ec08c791d616df23ca7c67a2220252b2499  paired-batches.csv.metadata.json
ab0896a5ac60a4b4c4d2464bcd36f91a5c5d6e9d7dd2320b7636995904cf836b  RESULT.json
```
