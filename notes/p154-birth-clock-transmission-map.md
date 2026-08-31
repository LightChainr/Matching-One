# P154: the missing high-gain transmission is an orientation clock

**Outcome:** a one-step source has four birth-conditioned clock enrichments.
Their common part is invisible to original U; three contrasts remain. The
previous two-rigid-birth model kept only the geometry-even relative clock,
which is weakly transmitted. The geometry-contrast, common-birth clock has
about **915 times its nominal U gain at N340**. This supplies a concrete
source → births → U map without fitting another archive descriptor.

It also gives two different, amplitude-free channel fingerprints: a locally
flat orientation/common-birth clock makes entry and completion move together;
an orientation/relative-birth clock makes them oppose. These are **conditional
mechanism predictions**, not findings that the actual source follows either
clock. The current numerical rectangles are not changed or relabelled as
physical theories. No new prospective observations were inspected.

## 1. The microscopic coefficient is a birth-selected mark

Use the fixed lag=1 source from `4daae57e`, in bulk cluster-count units.
For birth j=1,2 and geometry g, let b_jg(k)=Pr(K_j=k), F_jg its canonical
CDF, and z_(k−1)=s−E[s|k−1,rank,g]. Rank centering gives

```text
J_jg(p) = sum_k Bin(N,k;p) b_jg(k) m_jg(k),
m_jg(k) = E[z_(k−1) | K_j=k],
F'_jg(p) = N sum_k Bin(N−1,k−1;p) b_jg(k).
```

The normalized birth posterior is
`nu_jg,p(k)=N Bin(N−1,k−1;p)b_jg(k)/F'_jg(p)`. Hence exactly

```text
alpha_jg(p) := N J_jg/F'_jg
             = E_nu[(Np/K_j) m_jg(K_j)].
```

This is a source-enrichment coefficient among impending births, rather than
the magnitude of the bulk source. A 1/N clock scale requires bounded marks
and concentration of K_j/N away from zero; centering alone does not imply it.
Both direct 0→2 events enter both birth posteriors. The full proof, derivative
formula, and non-rigid interpretation of the pooled relative mark Xi are in
[the birth-posterior note](p154-lag-one-birth-posterior-clock.md).

## 2. Four clocks modulo one invisible clock

Let eta_g=(+1,−1) for the two geometries and tau_j=(−1,+1) for first and
completion birth. At each p make the exact four-coordinate decomposition

```text
alpha_jg = a00 + eta_g a10 + tau_j a01 + eta_g tau_j a11.
```

Each coefficient is one quarter of the corresponding signed sum of the four
alpha values. These are geometry/birth **contrasts**, not continuum parity
assignments, independent fields, or three mutually exclusive theories.
With q=F1+F2−1 and E=1−F1+F2, their actual source derivatives satisfy

```text
Jq_g = [(a00+eta_g a10) q'_g + (a01+eta_g a11) E'_g]/N,
JE_g = [(a00+eta_g a10) E'_g + (a01+eta_g a11) q'_g]/N.
```

Take all unmarked jets at the original pooled root. Write

```text
A=N^(13/8)/2, D=mean(q'), B=P4(E'), H=P4(E''), T=mean(q''),
P4(f)=(f_first−f_second)/delta_cos4.
```

For a mode with baseline pair `(f_q,f_E)`, define

```text
C[f] = A/(N D) * {P4(f_E') − mean(f_q) H/D
                  − (B/D)[mean(f_q')−mean(f_q) T/D]},
L[f] = A/(N D) * {P4(f_E) − (B/D)mean(f_q)}.
```

The complete original derivative is the exact first-jet map

```text
v = sum_(m=10,01,11) [C_m a_m + L_m a'_m],

m10: (f_q,f_E)=(eta q', eta E'),
m01: (f_q,f_E)=(E',q'),
m11: (f_q,f_E)=(eta E',eta q').
```

For m00 the pair is `(q',E')`, and both coefficients vanish identically.
Thus even a p-dependent common clock disappears. This common-temperature
null was already stated in Issue154; the new result here is the explicit
remaining mode map and its markedly unequal gains. Discarding a'_m because
one has estimated a pointwise clock ratio is not justified.

One compact derivation uses x=mean(q_t(p)) as the shared thermal coordinate
and R_t(x)=P4(E_t(p_t(x))). Then U=A R_t'(0) and, with D(p)=mean(q'(p)),

```text
v = (A/D) d_p {P4(JE) − [P4(E')/D] mean(Jq)} at the root.
```

Inserting each mode yields C and L. All entry/completion pieces use the same
full-source root movement and slope movement.

## 3. Nominal gains: orientation coherence is the efficient route

The following are a first point evaluation using **only already published
unmarked source-root jets** in
`7da1eeb0:results/norm4-source-endpoint-1m/latest.json`, `by_N[N].source`.
We did not read or fit source responses to choose a direction, reconstruct
paths, extract another old feature, or inspect new production. N85/N340 are
the upcoming experiment endpoints; N260 connects the previous planning note.

| N | C10, orientation/common-birth | C01, common relative-birth | C11, orientation/relative-birth | L11, thermal variation of that relative clock |
|---|---:|---:|---:|---:|
|85|60.73277|0.25324|5.17012|10.07563|
|260|219.78424|0.14951|12.19155|23.70077|
|340|246.38816|0.26940|12.32343|23.96402|

All gains multiply dimensionless alpha-mode values or their derivatives
with respect to p. C and L already include 1/N. At N340 a locally flat
pure a10 of **0.0020293** would give v=0.5 at these baseline centers. By
comparison pure a01 would require 1.85596. Neither amplitude has been
observed or forecast here; these are inverse sensitivities.

Why the difference? `P4(eta E'')=2 mean(E'')/delta_cos4`: an orientation
clock exposes the large ordinary E curvature to the angular readout. The
geometry-even relative clock instead depends on already small angular jets.
The exact map specifies how a microscopic source could exploit that
curvature; it does not establish that this source actually does.

The prior rigid two-birth hypothesis corresponds to a01 only, with
`delta=tau2−tau1=−2a01/N` and `C01=−2K_rel/N`. Its small nominal v therefore
does not bound the orientation-clock mechanisms. No source rescaling or
new arbitrary gain was introduced to obtain this distinction.

## 4. Two channel fingerprints for the existing readout

Consider two **first-jet restrictions at the prescribed root**, not a claim
that a nonzero clock can be constant throughout the complete finite-p range:

- Coherent orientation clock: a01=a11=a01'=a11'=a10'=0; a10 is free.
- Cancelling orientation/relative clock: a10=a01=a10'=a01'=a11'=0; a11 is free.

A common a00 and its derivative may be present in either case; they cancel
from each normalized readout. Under either restriction the vector
`(v_entry,v_completion)` lies on the corresponding unmarked gain line.
No amplitude estimate is needed to predict its line direction.

| N | Coherent line completion/entry | Cancelling line completion/entry |
|---|---:|---:|
|85|+1.186099|−1.186100|
|260|+1.117460|−1.117456|
|340|+1.105303|−1.105299|

For example, N340's coherent gain vector is `(117.03215,129.35600)`;
the cancelling vector is `(−117.03236,129.35580)`. Thus large, nearly
equal-and-opposite readout responses have a concrete temporal mechanism
that carries much less net U, while same-sign responses add coherently.
Negative common amplitudes reverse both signs but preserve these ratios.

For formal use, impose the **line restriction**, not a noisy response ratio:
`C_entry*v_completion−C_completion*v_entry=0`, using the appropriate
mode's unmarked coefficients and their joint uncertainty. The table is a
nominal calibration, not an uncertainty-propagated prospective band.
Both lines contain zero; neither is an exhaustive theory. Non-flat clocks,
mode mixtures, and direct-event currents can break these simple fingerprints.

This note supplies theory before seeing new output. It does not add a fourth
rescue template, change the other team's frozen scoring, or pool discovery
with fresh evidence. If the original numerical restrictions fail, they stay
failed. These independently stated jet restrictions can be considered on
their own terms with the already collected channels and shared covariance.

## 5. A separate event-current ambiguity survives

The [exact current note](p154-lag1-current-commutator.md) proves that the
three event currents enter a rank-gradient pairing with
`w02=w01+w12`. The signed direction `(-1,+1,-1)` leaves all endpoint rank
responses unchanged. A purely geometry-antisymmetric direct0→2 current
can generate exactly cancelling entry/completion with **zero net v**.
Consequently the cancelling-clock fingerprint is not unique evidence for
that clock. The existing 01/02/12 records, after the common root/slope
feedback is removed, already distinguish event support without collecting
another observable.

## 6. Finite-step correction to a constant-mark approximation

For an auxiliary constant birth mark m, put
`S_j(p)=sum_k Bin(N,k;p)b_j(k)` and `T_p=(1−p)d_p`. Exactly

```text
[N+T_p] S_j=F_j',  S_j(0)=0,
S_j=F_j'/N − T_p F_j'/N² + T_p² S_j/N².
```

The differential equation needs the boundary condition: `(1−p)^N` is a
homogeneous solution. For the real source always retain the exact remainder
`R_j=sum_k Bin(N,k)b_j(k)[m_j(k)−m]`, so `J_j=m S_j+R_j`.
The leading common `m F_j'/N` is screened by U; the next term is a
thermal-shape correction. In the shared x coordinate, for `c(p)=1−p`,

```text
L_U[c F1'',c F2'']
   = A * [(2c T+c' D) R''(0) + c D² R'''(0)].
```

Here L_U is the original source-to-v functional **without** an added 1/N,
and R is the unmarked projected E curve in x. Therefore the first common
constant-mark contribution is minus m/N² times this expression, plus the
explicit higher-step and nonconstant-mark remainders. It is not exactly a
shared canonical clock. A nonzero constant mark over the entire finite
birth support is incompatible with rank centering at its last possible
birth; only a controlled critical-window approximation is plausible.
We have not assumed the remainders small, assigned them an exponent, or
fitted m from the archive.

The [quotient companion](p154-original-u-clock-quotient.md) also gives the
equivalent exact auxiliary remainder with denominator N(N−1), and valid
ordered-birth examples where a real nonclock deformation changes intrinsic
curvature while original U stays zero. These specify what a null U result
cannot exclude; they do not request additional production readouts.

## Scientific card and reproducibility

- Changed mechanism space: the weak geometry-even two-clock channel no
  longer stands in for all temporal transmission. A high-gain orientation
  clock and a cancelling orientation/relative clock have distinct first-jet
  predictions for the already fixed original-U readout.
- Does not prove: actual source amplitudes, new field identity, physical
  energy coupling, or uniqueness of a cancelling mechanism.
- Observer/sector/source/geometry: original pooled-root slope-normalized
  P4(E') U; lag1 rank-centered bulk CB+CW; original N85/260/340 pairs.
- Dependence: deterministic theory; numeric gains reuse unmarked jets from
  the existing norm4 archive, not another independent evidence block.
- Next observation that changes weight: the fresh entry/completion/net
  result, interpreted against specified jet restrictions and the already
  retained 01/02/12 event support. No added descriptor or lag scan.

Run `scripts/p154_birth_clock_transmission_map.py --output NEW_JSON_PATH`
with the managed research Python. The saved
[`gains.json`](../results/p154-birth-clock-transmission-map/gains.json)
contains the exact source commit/hash, baseline jets, all value/derivative
gain vectors, and explicit zero-new-sample/no-source-fit flags. Gain
uncertainty was not newly propagated. This local calculation took well under
a second and required no cloud resource or package installation.
