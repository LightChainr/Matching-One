# Integrable Potts field theory route to the universal charge scaling function

2026-09-14.  Literature-grounded continuum programme for #767.  This note does NOT claim that the required percolation magnetic-sector NLIE has already been written in the needed form.  It narrows the problem to one sharply identifiable finite-volume sector.

## 1. The lattice transfer has already identified the continuum object

The fixed-width safe transfer gives two microscopic excitation energies

\[
I^0_{4,w}(p),\qquad I^0_{8,w}(1-p).                           \tag{1.1}
\]

At criticality both converge to the SAME magnetic primary gap

\[
wI^0\to2\pi x_m,
\qquad x_m=5/48.                                              \tag{1.2}
\]

Near criticality, universality plus complementary thermal orientation suggests one universal magnetic finite-size energy scaling function

\[
\mathcal E_m(X),                                               \tag{1.3}
\]

with primal and matching sectors sampling opposite thermal arguments.  Therefore the leading charge scaling function is not a difference of two unrelated continuum sectors but

\[
\boxed{
\mathcal F(X)=\mathcal E_m(X)-\mathcal E_m(-X).}             \tag{1.4}
\]

Equation (1.4) automatically implies

\[
\mathcal F(-X)=-\mathcal F(X),                               \tag{1.5}
\]

which is precisely the odd-derivative/even-cancellation pattern seen in the safe transfer.

So the continuum task is:

> compute the finite-volume energy of the **magnetic Potts sector** under the thermal perturbation, continuously at `q=1`.

## 2. Why integrable Potts field theory is relevant

The scaling `q`-state Potts field theory for `q<=4` under its thermal perturbation is integrable.  Chim--Zamolodchikov supply the kink S-matrix; Dorey--Pocklington--Tateo develop finite-size TBA/NLIE descriptions and explicitly note that the continuous-parameter NLIE is particularly suitable for the `q->1` percolation limit.

Delfino--Cardy take the same exact Potts scattering theory to `q->1` and compute order/disorder/energy form factors and percolation amplitude ratios.  Thus the bulk massive field theory is already known to survive the percolation continuation in a useful sense.

What is NOT already supplied by these facts is the particular finite-volume magnetic/topological energy (1.3).

## 3. A useful template from the three-state Potts excited TBA

Lencses--Takacs (2014) construct excited-state TBA equations for the scaling three-state Potts model.  In the ferromagnetic phase they generate additional finite-volume sectors by inserting a symmetry twist; the TBA pseudoenergy acquires a twist parameter, and the UV limits are matched to CFT sectors.  They also use Kramers--Wannier duality to relate the two thermal signs.

This gives a concrete template:

1. finite-volume Potts sectors can be selected by TBA/NLIE twists/sources;
2. the same thermal field theory can be followed through both signs of the thermal coupling;
3. the UV conformal weight of the selected sector is a decisive identification check.

The three-state twist itself cannot simply be copied to noninteger `q=1`; its role is methodological.

## 4. The percolation magnetic UV fingerprint is c_eff=-5/4

For a CFT sector whose lowest state has total scaling dimension `x`, the cylinder energy can be written

\[
E(R)
=-\frac{\pi c_{eff}}{6R}+\cdots,
\qquad
c_{eff}=c-12x.                                                \tag{4.1}
\]

Percolation has

\[
c=0,
\qquad
x_m=5/48.                                                     \tag{4.2}
\]

Therefore the required magnetic-sector finite-volume equation MUST have ultraviolet effective central charge

\[
\boxed{
c_{eff}^{(m)}
=0-12\frac5{48}
=-\frac54.}                                                   \tag{4.3}
\]

This is a hard sector-dictionary criterion.  A candidate twist/source with a different UV value is not the homology-safe magnetic sector, regardless of how attractive its infrared interpretation looks.

## 5. How to search for the q->1 magnetic twist

There are two plausible routes.

### Route A: analytic continuation of a generic-q twist

Start from a Potts loop/RSOS/NLIE representation in which noncontractible-cluster or symmetry sectors are parameterized continuously.  Continue the sector parameter together with `q->1` and impose (4.3).

The periodic TL interpretation on the lattice is useful here: the safe/open sector is the sector carrying the magnetic excitation.  Its loop/noncontractible-weight parameter should have a direct continuum twist analogue.

### Route B: excited-state analytic continuation

Start from the continuous-q ground-state NLIE and generate an excited solution by contour/source continuation, as in standard excited-state TBA constructions.  Choose the source pattern so the UV dilogarithm result is (4.3).

This may be technically cleaner if no explicit noninteger Potts symmetry twist is available.

In either route, (4.3) removes much of the ambiguity.

## 6. The desired output

Let

\[
r=m_{phys}R                                                    \tag{6.1}
\]

be the usual dimensionless massive finite-volume variable.  Solve for the magnetic excitation energy above the vacuum,

\[
\mathcal E_m^{(+)}(r),\qquad
\mathcal E_m^{(-)}(r),                                       \tag{6.2}
\]

on the two signs of the thermal perturbation, with a common mass normalization.

Kramers--Wannier / Potts duality should identify these as two branches of one function with reversed thermal coordinate.  The charge function is

\[
\boxed{
\mathcal F(r)
=\mathcal E_m^{(+)}(r)-\mathcal E_m^{(-)}(r).}                \tag{6.3}
\]

After matching the lattice thermal metric,

\[
X\propto(p-p_c)w^{3/4},                                      \tag{6.4}
\]

(6.3) should be the continuum limit of

\[
w\Theta_w(p).                                                \tag{6.5}
\]

## 7. Immediate lattice checks on a candidate NLIE solution

A successful continuum construction should satisfy ALL of the following.

### UV magnetic gap

\[
\mathcal E_m(0)=2\pi\frac5{48}.                              \tag{7.1}
\]

Equivalently its sector effective central charge is `-5/4`.

### Dual oddness

\[
\mathcal F(-X)=-\mathcal F(X).                               \tag{7.2}
\]

### First derivative

After fixing the thermal metric,

\[
\mathcal F'(0)                                                \tag{7.3}
\]

must reproduce the safe pivotal / Perron slope amplitude.  In the raw square-site logit coordinate the current widths give a sequence trending near `0.8`.

### Third derivative

The transfer gives a stable scaled third derivative near

\[
\mathcal F'''(0)\approx0.064                                 \tag{7.4}
\]

in the same uncalibrated coordinate.  Ratios after metric fixing provide a sharper comparison than either derivative separately.

### No leading second derivative

The universal contribution must obey

\[
\mathcal F''(0)=0.                                           \tag{7.5}
\]

The lattice `Theta_hh` is indeed parametrically smaller than the individual sector second derivatives.

### Infrared consistency

For large positive/negative thermal mass, the magnetic sector energy should match the appropriate kink/particle/disorder excitation of the two Potts phases.  This is where the known exact S-matrix and order/disorder form factors provide independent checks.

## 8. Why this is better than fitting an arbitrary crossover function

The charge-neutral description in #767 currently allows an unknown odd function `F(X)`.  The integrable-field-theory route would replace that by a specific finite-volume spectral problem with

- known bulk S-matrix;
- known thermal mass-coupling relation in Potts field theory;
- a fixed UV magnetic conformal weight;
- exact duality between thermal signs;
- numerical TBA/NLIE methods.

That is a much narrower and more falsifiable problem than inventing a phenomenological interpolation.

## 9. Main unresolved technical point

The current literature checked in this pass does NOT hand us a ready-made generic-`q`, magnetic-sector excited NLIE whose `q->1` limit visibly has `c_eff=-5/4`.

- Dorey--Pocklington--Tateo provide the continuous-q thermal finite-size framework and emphasize its usefulness at percolation.
- Lencses--Takacs show explicitly how twist/excited sectors work in the three-state model.
- Delfino--Cardy show that the `q->1` massive Potts form-factor continuation is meaningful for percolation observables.

The missing bridge is the magnetic-sector finite-volume continuation at noninteger q.

This is now a precise literature/derivation target rather than a vague “maybe TBA helps” suggestion.

## 10. Claim boundary

The integrability of the thermally perturbed Potts scaling theory and the cited finite-size/form-factor results are literature facts.  The identification (1.4) is a universality/duality synthesis strongly supported by the lattice transfer parity data.  The existence and exact form of a `q=1` magnetic-sector NLIE with UV `c_eff=-5/4` remain to be established.
