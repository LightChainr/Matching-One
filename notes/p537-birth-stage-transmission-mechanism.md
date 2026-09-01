# Birth-stage gating as the missing microscopic-to-global transmission map

## Mechanism extracted from the exact collar

The radius-one collar result supplies a two-stage transmission law rather than
another undifferentiated source response.  In the frozen physical basis

\[
L=\begin{pmatrix}
L_{first,absent}&L_{first,present}\\
L_{second,absent}&L_{second,present}
\end{pmatrix},
\]

the exact N25 matrix is approximately

\[
10^7L=
\begin{pmatrix}
2.710045&1.645514\\
-1.697815&-0.038140
\end{pmatrix}.                                               \tag{1}
\]

The source-presence retention factors are

\[
r_1={L_{first,present}\over L_{first,absent}}=0.6071907\ldots,
\qquad
r_2={L_{second,present}\over L_{second,absent}}=0.0224641\ldots. \tag{2}
\]

Their ratio is `27.0294`, while the diagonal-rescaling-invariant cross-ratio

\[
\chi={L_{first,absent}L_{second,present}\over
           L_{first,present}L_{second,absent}}
     ={r_2\over r_1}=0.0369967\ldots                         \tag{3}
\]

is far from the rank-one value one.  Calibrating the common-amplitude model
on the first birth leaves the exact second-stage interaction

\[
I_2=L_{second,present}-r_1L_{second,absent}
   ={\det L\over L_{first,absent}}
   =+9.9275800\ldots\times10^{-8}.                           \tag{4}
\]

Equations (2)--(4) identify the minimal second channel: `axial2` source
presence largely preserves first-birth transmission but almost extinguishes
the second-birth response through a signed cancellation.

## Relation to the earlier temporal-transmission null

The prospective #154 lag-one experiment found no large aggregate completion
response and correctly stopped that source as the main H4 explanation.  The
collar result suggests a sharper reading of such a null: a weak completion
mean need not imply that the microscopic source fails to reach the completion
sector.  A real rank-two coupling can be present while its source and
root-Schur pieces cancel specifically at the second birth.

This does not retroactively identify the #154 lag-one source with the N25
`axial2` source; the random streams, scales, and source definitions differ.
It supplies a concrete transmission mechanism to test:

```text
microscopic pair source
    -> radius-one four-arm collar
    -> outer double join (J_B=J_W=1)
    -> first/second birth gate
    -> signed P4-Schur readout.
```

The critical observable is no longer a raw completion amplitude.  It is the
within-stream cross-ratio `chi`, which cancels independent row and column
normalizations and asks whether completion is attenuated *relative to entry*.

## One-shot prospective discriminator

For any future production block that already records the two rank births,
add only the radius-one corner word and the two outer join bits.  Freeze the
same four cells as in (1) and report

\[
                         \chi_N={L_{1,A}L_{2,P}\over L_{1,P}L_{2,A}}. \tag{5}
\]

The competing mechanisms make sharp predictions:

- common pure-thermal amplitude: `chi_N -> 1`;
- collar birth-stage gate: `chi_N` remains separated below one, with the
  second-birth present cell strongly attenuated;
- finite-N contact artifact: `chi_N` moves rapidly toward one when the collar
  is embedded in larger tori or when the source orbit is moved outward.

No descriptor search is needed.  One held-out size and the predeclared four
cells decide which route survives.

## Boundary

The values above are exact signed P4-Schur contributions for the N25
radius-one collar, not transition probabilities.  The proposed cross-size
behavior is a mechanism hypothesis, not a result.  The exact result is the
finite rank-two matrix and its stage-selective cross-ratio.
