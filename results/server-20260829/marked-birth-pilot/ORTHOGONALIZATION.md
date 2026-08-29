# Held-out batch-Gram orthogonalization

N65 trains `alpha=0.00626022695459707404599426+(-0.0209896457381260144007131)i` using the centered complex batch Gram normal equation.
N130 is held out; P50 N145 is external direction only.

The projection does not remove thermal growth:

- raw connected q2 transfer: `{'re': '1.26439187314319679165298', 'im': '-0.00864665042631079889628515', 'abs': '1.26442143822151191369753'}`;
- orthogonal connected q2 transfer: `{'re': '1.26446925951222720883551', 'im': '0.0012460117256518300614657', 'abs': '1.2644698734238869407709'}`.

The reason is stronger than a weak regression. Gate algebra gives exactly, orientation by orientation,

```text
Cov(q,J_D)=J_S/2+(p-1/2-<q>)J_D,
Cov(q,J_S)=J_D/2+(p-1/2-<q>)J_S.
```

Thus the connected rank response contains no independent matrix element beyond source means.
A true field Gram additionally needs two-root `J_D*conj(J_S)` and `|J_S|^2`, and a non-tautological coupling needs an independent global observer.

## Scientific card

1. MECHANISM SPACE: the apparent thermal growth is an exact rank-gate contact channel, not free evidence for a Q4 matrix element.
2. NOT PROVED: the noisy mean-J_D q2 exponent remains a candidate; batch Gram is not the missing field-level Gram.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: A_top | Alexander odd | J_D/J_S gate doublet | q2 Gaussian parent and held-out child.
4. DEPENDENCY GROUP: orthogonal and raw scores are deterministic transforms of the same 634040d pilot.
5. UPWEIGHT OBSERVATION: mean-J_D radial transfer or an independent-observer coupling after two-root field orthogonalization.
