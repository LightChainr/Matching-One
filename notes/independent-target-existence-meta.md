# Independent matching-odd target synthesis

## Question

Issue #57 by itself does not reject a zero child effect. That fact is useful for describing the power of the norm-5 block, but it should not be promoted into a project-wide statement that the matching-odd signal has not been established.

The clean existence question can be asked without fitting an exponent and without combining correlated derivative/root views: use the two independent prospective matching-odd target blocks already present in the evidence ledger.

## Selected blocks

Only

- `issue43_n185_n265_deltaM`;
- `issue57_norm5`.

Both are scored primary matching-odd blocks. Their `raw_data_group` values are distinct. The target samples were generated on independent streams, so the cross-block target covariance is zero.

No P49 derivative, P45 root, P50 full-curve, Krawtchouk or other deterministic view is included.

## Result

The zero-effect target-only goodness-of-fit is

```text
N185/N265: 29.4093843121 / 2
N325/N425:  1.7763512394 / 2
joint:     31.1857355515 / 4
p:          2.8055953e-6
```

The block-specific fixed H4 predictions give

```text
N185/N265: 3.0459757733 / 2
N325/N425: 0.4163037640 / 2
joint:      3.4622795373 / 4
p:          0.48363695
```

The independent predictive NLPDs add to

```text
H4  = -35.7946059274
zero= -21.9881871377
Delta NLPD(H4-zero) = -13.8064187897 nats.
```

The NLPD difference is the proper model-ranking summary across the two registered predictive distributions. The raw chi-square difference is not advertised as a calibrated likelihood-ratio statistic because the predictive covariance under H4 includes propagated source uncertainty while zero does not.

## Interpretation

The scientifically correct boundary is:

> The norm-5 children alone are compatible with zero, but the independent N185/N265 and N325/N425 matching-odd target blocks jointly reject the global zero-effect description strongly, while the corresponding fixed H4 predictions are jointly compatible.

This changes the execution priority. Signal existence is not the principal bottleneck. The high-information questions are now mechanism and representation:

1. whether a genuinely local pivotal/four-arm observable carries a related H4 response;
2. whether the finite-size thermal block is analytic or logarithmic/Jordan;
3. whether noncyclic Smith/quotient arithmetic contributes at the current sizes;
4. whether the same field has the predicted torus-modulus `E4/g2` fingerprint.

## Evidence boundary

This is a post-reveal synthesis of two evidence rows that already exist. It must not be inserted as a third additive primary vote in the prequential ledger. Its purpose is to prevent the low-power zero test of one block from being mistaken for the global existence question.
