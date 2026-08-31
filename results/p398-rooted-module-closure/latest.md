# P333/P398: minimal closure by coupling rooted and old landing marks

The seven new rooted derivative marks alone fail the complete inherited P333 intersection. The sole remaining obstruction is one rational C4 charge-one doublet. Coupling the **already-defined two-column landing doublet** to the seven rooted marks gives a unique full affine/endpoint/radical-Gram/source solution: `X0=diag(T,R), V=0`.

| candidate | dim W | rank G0 | dim radical | radical C4 dimensions (trivial, charge1, sign) | Gram-skew rank | full inherited intersection |
|---|---:|---:|---:|---|---:|---|
| rooted7 | 21 | 15 | 6 | 3, 2, 1 | 2 | empty |
| rooted7_plus_existing_landing_charge1 | 23 | 19 | 4 | 3, 0, 1 | 0 | unique |

## Why this is the full intersection, not just a canonical-point check

1. Write X0=[[A,0],[B,C]] and V=[[a,0],[b,c]]. Mark transport fixes C=R,c=0.
2. The order-zero ordinary affine block has A=alpha*T; fixing the ordinary source sets alpha=1.
3. Each row of B is a common left invariant; fixing the ordinary source forces B=0.
4. P_(i+1)T=TP_i cancels the inhomogeneous ordinary jet term, so a=beta*T; V(source)=0 sets beta=0.
5. R E_i=E_(i+1)T cancels the lower jet term; each row of b is a common left invariant and V(source)=0 forces b=0.
6. Thus the source-transport-normalized affine first jet is uniquely X0=diag(T,R), V=0. Endpoint and radical are checked at that point; no surviving modulus can alter its Gram residual.

The exact ordinary shifted-Hom rank is 195/196 and the common-left-invariant rank is 13/14. Hence the source-normalized affine point is unique before imposing Gram. For rooted7 the Gram restriction has coefficient rank 0 and augmented rank 1; its explicit radical-vector witness is in `latest.json`. For the coupled nine-mark lift every inherited gate is zero exactly, and the surviving four-dimensional radical still has nondegenerate leading Gram form: success is not a vacuous radical-exhaustion trick.

## Minimum coupled mark

The rooted7 residual is a nondegenerate two-dimensional rational C4 charge-one block. A one-dimensional rational C4 mark is trivial or sign and cannot pair with it; the existing two-column landing response removes precisely this residual block.

The two new-to-this-module columns are `G1[:,landing0]-G1[:,landing2]` and `G1[:,landing1]-G1[:,landing3]`, already tested in the old landing branch. Their isolated width-four failure and the rooted7 failure do not survive their coupling. The seven rooted marks are not a standalone replacement for the old response family.

Here coupling means sharing the ordinary Gram/source block and configuration-dependent emission. It is not a fitted off-diagonal exchange between the two mark families or evidence for a physical interaction.

## Semantics and boundary

This uses the prior P333 retained-response accumulator convention: joins act identically on emitted marks, while detach emits the configuration-dependent row `C0^T P_i`; it does **not** pretend the seven reference vectors are an invariant active connectivity submodule. Under active-reference semantics their join closure is 13 and join+detach closure is all 14, with a one-coordinate witness already at `J0 AP_sum`. That separate saturation is recorded, not used as the full-intersection decision.

The inherited Gram gate is the first-jet form restricted to `ker G0`. Unprojected G0/G1 self-adjoint residuals are explicitly nonzero in the artifact; the result does not establish a stronger all-Gram/all-Q closure, a physical transfer module, or a continuum Jordan/field identity.

## Reproduce

```bash
python3 scripts/p398_rooted_module_closure.py
python3 -m unittest discover -s tests -p test_p398_rooted_module_closure.py
```
