# Three-generation same-lineage H4 recurrence

No new simulation is used. Exact alternating H4 geometry signs are divided out before analysis.

With `lambda0=2^(-13/8)=0.324209889`, the exact three-point identity gives `lambda1=0.212404 +/- 0.236438` (delta method). The point lies in `(0,1)`, but is only `0.898` SE above zero; the 95% delta interval is `[-0.251,0.676]`.

The amplitudes are `c0=-0.068889` and `c1=+0.052838`. They have opposite signs. The correction/leading magnitude ratio falls `0.767 -> 0.502 -> 0.329 -> 0.216` from N85 through predicted N680.

Two-mode N680 prediction: `A_H=-0.0018413 +/- 0.000955`; its exact child geometry makes the pair negative, `-0.0029359`.

| model | in-sample q/df | N680 A_H | prediction SE |
|---|---:|---:|---:|
| frozen single lambda0 | 11.747/2 | -0.0010272 | 8.57e-05 |
| free single lambda | 1.271/1 | -0.0030378 | 0.00102 |
| scale-neutral | 17.686/2 | -0.0089120 | 0.00076 |
| fixed-lambda0 plus correction | interpolation/0 | -0.0018413 | 0.000955 |

The useful new coordinate is a point-estimate fast, opposite-sign correction that naturally produces the N170 overshoot and N340 return. Its sign/rate is not yet resolved: about 80% of lambda1 variance comes from N340. N680, not the zero-df interpolation, is the clean discriminator. No exponent is fitted across unrelated geometries.
