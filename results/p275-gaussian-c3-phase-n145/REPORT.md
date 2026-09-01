# Held-out N145 three-way primitive-C3 phase gate

Status: `H0_EVEN_CHARACTER_SELECTED`; no top-up.

The held-out pair is

```text
g1 = 12+i,       g2 = 9+8i,
N1 = N2 = 145,   exp(i delta)=(4+3i)/5,
delta = 36.8698976458 degrees.
```

Five million fresh paired common-field replicas were divided into 100 batches
with seed `20260901277`.  The observer, `tau=i` baseline, normalizer and four
retained coordinates are identical to the N65 gate.  No N65 random stream is
reused.

## Frozen three-model decision

```text
mean z1 = +0.000850646627122 - 0.00147301676729 i
mean z2 = -0.000917453372878 + 0.00161124610575 i

H0: chi2 =  0.00154675 / 1 df, p = 0.968628,    gain = -1.09016
H4: chi2 = 14.56361649 / 1 df, p = 0.000135506, gain = +1.06343
H8: chi2 = 53.43005623 / 1 df, p = 2.67966e-13, gain = -1.05446
```

Exactly one model survives at the frozen `alpha=0.01`, hence

```text
H0_EVEN_CHARACTER_SELECTED
```

The convention-invariant projective statistic independently exposes the same
geometry:

```text
X = 0.99992606 +/- 0.00766449

H0 prediction = +1
H4 prediction = +0.42197248
H8 prediction = -0.64387845
```

No extra sample is authorized or drawn.

## Mechanism update

The N65 H4/H8 gate selected H8 only within its frozen two-model set.  At that
angle `-8 delta` is within one degree of `pi`, so signed-scalar H0 also passed
post reveal.  The new angle separates all three lines and shows that the
surviving object is embedding-even, not a local H8 harmonic.

The primitive complex response keeps one projective direction under physical
rotation and changes only by a signed real amplitude.  It is therefore a
nontrivial homology-sector character with spin-zero embedding transport.  This
also changes how the older norm-2 sign alternation should be read: a negative
cover ratio alone cannot identify H4 unless positivity of the radial gain is
independently justified.  The global `A_top`/quotient-prism H4 result is not
altered; the two observers live in different transport sectors.

The algebraic operation is in fact exact.  Multiplication by the Gaussian unit
`i` is an automorphism of every quotient `Z[i]/(g)` and maps period-basis
homology `(m,n)->(-n,m)`.  It exchanges the unoriented target lines `l0,l1`, so
`P_l0=P_l1` at every finite size and every bond probability.  Since the
`tau=i` baseline has the same equality,

```text
z = omega * [(P_l2-pi_l2) - (P_l0-pi_l0)] in omega R.
```

The signed-real H0 line is therefore a theorem for this Gaussian family, not a
new fitted field.  This closes further Gaussian-angle Monte Carlo with the same
character; local spin can only be tested after leaving this unit-rotation
symmetry class or changing the observer.

Input/output hashes:

```text
0e08402c3772baf4c277fbd22cfb973abe509c6408b2eb211ea2a73d5d9e904f  paired-batches.csv
a324ffdcd1f037f9a5f9890b406095a70ab41ff7ccf3c44bbde0d630e4095bf6  paired-batches.csv.metadata.json
88df1e3b604b19a5a42be5df81bceca1d40b58726677b9713b130479a5b7602b  RESULT.json
```
