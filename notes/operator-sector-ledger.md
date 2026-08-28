# Operator / sector evidence ledger

Purpose: keep exact identities, finite-size observations, scaling deductions and LCFT/operator identifications at separate evidence levels.  A failed downstream interpretation must not erase upstream numerical facts.

Last updated: 2026-08-28.

## Evidence labels

- **E0 exact** — finite combinatorial/algebraic identity or exact arithmetic check.
- **E1 observed** — direct Monte Carlo/transfer result with declared uncertainty.
- **E2 predictive finite-size law** — survives independent seed and at least one frozen/held-out prediction.
- **E3 mechanism deduction** — follows from E0–E2 plus explicit scaling assumptions.
- **E4 operator identification** — assignment to a specific CFT/LCFT field; requires independent predictions beyond exponent matching.

## A. Finite matching involution

### A1. Periodic matching identity

Status: **E0 exact**.

Black primal / white matching topology gives the finite periodic matching relation used throughout the project. Tiny exhaustive tests agree to arbitrary-precision roundoff. The five wrapping-difference channels in the current complementary construction are configuration-identical, not independent estimators.

Consequences:

- wrapping-only GLS is structurally rank one;
- `S=(R_G+R_hat)/2` and `D=(R_G-R_hat)/2` are the natural matching-parity combinations.

Primary falsifier: an exact counterexample on a correctly implemented periodic quotient. None known.

## B. Matching-odd orientation sector

### B1. Nonzero same-N angular effect

Status: **E2 predictive finite-size evidence**.

Evidence:

- independent 100M P31 seed reproduces signs at N=65,85,130,145,170;
- frozen H4 `N^-13/8` model predicts held-out N=145,170 better than zero and tested alternatives;
- prospective Gaussian-doubling #38 passes the parameter-free `-2^-13/8` ratio on two fresh lineages.

Robust statement:

> there is a matching-odd orientation-dependent correction in the odd square-lattice harmonic class whose finite-size behavior over the tested range is accurately described by `DeltaM ~ DeltaCos4 N^-13/8`.

Not yet established:

- unique asymptotic exponent 13/8;
- H4 rather than H12/H20;
- a specific LCFT field.

Next falsifiers: #43, #49, #50, #55.

### B2. H4 angular identification

Status: **E2/E3 boundary, not unique**.

Support:

- H4 explains P31/P32;
- H8 is unresolved and does not improve held-out prediction;
- pi/4 doubling selects odd `m` harmonics.

Open alternative:

- H12 is also odd under pi/4 rotation and is not excluded by doubling alone.

Decisive tests:

- #55 N=305/325 opposite H12/H4 alias ratios;
- #36 stage-2 angular regression;
- N=1105 multi-angle commuting-square test only after power/provenance gates.

## C. Matching-even orientation sector

### C1. Leading `N^-1` / `L^-2` angular correction

Status: **E1 with strong finite-size consistency; prospective confirmation still desirable**.

P31 100M `S=(R_G+R_hat)/2` data are compatible with

`P4[S] ~ A_I N^-1`,

with an inverse-variance common amplitude near `A_I=0.0106 +/- 0.00094` and acceptable five-size scatter.

Candidate interpretation: square-lattice identity-family spin-4 field with `x=4` and matching-even parity.

This candidate is motivated by standard square-lattice CFT irrelevant-operator classifications but is not yet identified in percolation.

Falsifiers/controls:

- #43 simultaneously scores the frozen S and M sectors on new N=185/265;
- #42 square-bond self-dual control;
- #44 C4 self-matching site triangulation control.

## D. Thermal-family spin-4 candidate

### D1. Algebraic existence

Status: **E0 algebraic within the ordinary Virasoro quotient**.

For `c=0,h=5/8`, after the level-2 null relation

`(L_-2 - 2/3 L_-1^2)|h>=0`,

a non-null level-4 quasiprimary survives:

`Q4=(40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4)|h>`.

Checks committed in `scripts/virasoro_level4_candidate.py`:

- `L1 Q4=0`;
- not in the level-2 null descendant span;
- not a simple `L_-1` image;
- nonzero formal Shapovalov norm.

Bulk chiral pairing `(h+4,h)` or conjugate gives

`x=21/4`, spin `+/-4`.

### D2. Identification with the measured matching-odd field

Status: **E4 conjecture, now high-value but not proved**.

Why attractive:

- `L^(2-x)=L^-13/4=N^-13/8` matches the observed central residual law;
- division by the thermal slope `M'~L^3/4` gives a root bias `L^-4`;
- field is spin 4 and therefore changes sign under 45-degree rotation;
- matching/thermal parity pattern is compatible with P33 central-reflection data.

What is still missing:

1. proof/derivation of matching-involution parity of the physical LCFT scaling field;
2. exclusion of lower/equal-dimension competing matching-odd singlet fields;
3. treatment of the logarithmic partner/Jordan module;
4. prospective derivative-spectrum and control-model predictions.

Primary issues: #37, #44, #48, #49.

## E. Matching-parity derivative spectrum

Status: **E3 frozen prediction; prospective test pending**.

For an eigenfield with matching parity `eta`, the linearized involution predicts at the intrinsic center

`S^(n) != 0 iff (-1)^n=eta`,

`D^(n) != 0 iff (-1)^n=-eta`.

After the same-N spin-4 projector, the two-field model freezes:

```text
matching-even x=4:
  P4[S]   ~ N^-1
  P4[D']  ~ N^-5/8

matching-odd x=21/4:
  P4[D]   ~ N^-13/8
  P4[S']  ~ N^-5/4
```

Retrospective P33 analysis is development only. Prospective score is reserved for new statistics/geometries under #48/#43.

## F. Root mechanism

### F1. Local residual-to-root closure

Status: **E2 observed mechanism at finite size**.

P35 finds

`C=-DeltaRoot * mean(M') / DeltaM`

within roughly `3e-4` of one across all five P33 sizes. Thus direct and linearized roots agree and the local root lies in the linear regime.

### F2. Asymptotic `N^-2` angular-normalized root law

Status: **not yet E2**.

Correct invariant under the H4 model is

`A_p=-N^2 DeltaRoot / DeltaCos4`,

not bare `N^2 DeltaRoot`.

Frozen prediction from independent source amplitudes:

`A_p=A_M/B ~= 0.4510`.

Existing 10M radial data drift; cross-size covariance audit shows the drift is not an artifact of diagonal covariance. High-stat #45/#49 remains decisive.

## G. Historical post-L^-7 annihilator

Status: **open**.

Mertens-Ziff small-L data suggested an accelerated root near `L^-7`. If the leading central amplitude is

`L^-13/4 (1+c L^-q+...)`,

the accelerated root scales as `L^-(4+q)`; genuine `L^-7` means `q=3`.

Thermal-family quasiprimary counting gives the next ordinary spin-4 pair at total descendant level 10, i.e. relative `q=6`, not 3. Therefore the historical 7 cannot simply be labelled the next ordinary thermal spin-4 descendant.

Alternatives: another conformal family, nonlinear/composite correction, logarithmic mixing, second-order correction sector, or preasymptotic effective power.

The corrected critical-Potts branch supplies a specific standard-spectrum
candidate: the diagonal Potts singlet `V_<1,4>` has `x=33/4` and spin zero.
Under the conditional interchiral/OPE parity rule `eta_s=(-1)^(s-1)`, it is
matching odd.  Its contribution `L^(2-x)=L^-25/4` is exactly three powers below
the leading `L^-13/4` term and therefore produces an `L^-7` accelerated root.
The spectrum membership and arithmetic are E0; the matching parity is E3 under
the stated automorphism assumption; identifying the historical signal with a
nonzero `V_<1,4>` coupling remains E4.

Decisive test: #47 modern threshold-rank replay with held-out exponent challenge.

## H. Logarithmic partner

Status: **theoretically allowed; not numerically required at current precision**.

Because percolation is a `c=0` LCFT and the thermal Kac field participates in logarithmic structure, terms such as

`N^-13/8 (A+B log N)`

must remain admissible. P32's log extension did not improve held-out prediction, so current data do not require nonzero B.

Gaussian doubling is especially diagnostic because a log partner creates a coherent nonzero pure-power doubling residual proportional to `B log 2`. #50 preserves this as a dedicated Jordan diagnostic.

## I. Universal matching crossover / kappa3

Status: **universal-object program; exact value unresolved**.

`kappa3=Mcal'''(0)/Mcal'(0)^3` removes the nonuniversal thermal metric factor. Exact-threshold square-bond/triangular controls are compatible with values near `-5/3` but do not establish the rational value.

The analytic target is now formulated as a torus homology-projector / integrated-thermal-operator cumulant ratio in the `Q->1` Potts/FK theory (#54).

Do not perform rational reconstruction until:

- same-modulus universality is demonstrated across microscopic realizations;
- derivative covariance/systematic finite-size corrections are controlled;
- ideally an independent continuum/TCSA/TL calculation exists.

## J. Current strongest research chain

The strongest currently defensible chain is:

1. exact finite matching involution/topology;
2. independently reproduced same-N angular matching-odd signal;
3. held-out H4 `N^-13/8` predictive success;
4. prospective no-fit Gaussian-doubling success;
5. thermal-even central character in the full-curve pilot;
6. exact existence of an `x=21/4`, spin-4 thermal-family quasiprimary candidate;
7. local residual-to-root closure.

The weakest links still separating this chain from an operator identification are:

- H4 versus H12;
- physical matching parity of the LCFT field;
- clean full-curve prospective replay/root scaling;
- self-dual/self-matching controls;
- logarithmic-module interpretation.

All new computation should target one of these weak links rather than merely increasing significance of an already-established center signal.
