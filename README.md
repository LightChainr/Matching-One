# Matching One

Exploratory project on the **square-lattice site percolation threshold**

$$
p_c^{\mathrm{site}}(\mathbb{Z}^2) = 0.59274605079210(2)
$$

and its exact partner on the matching lattice

$$
p_c^{\mathrm{site}}(\mathbb{Z}^2) + p_c^{\mathrm{site}}(\mathrm{NN+NNN}) = 1.
$$

The second number is $1-p_c \approx 0.40725394920790$. The two constants are the same problem.

There is **no known closed form**. Bond percolation on the same lattice is exactly $1/2$ by duality. Site percolation is not self-dual; that missing closed form is the object of this repo.

## What this is not

Not a claim that $p_c$ equals a simple combination of $\pi$, $e$, or $\varphi$. Those guesses fail at $10^{-3}$ or worse against a 14-digit value. The working stance is:

1. Record the constant and the exact matching identity.
2. Systematically exclude low-complexity closed forms (PSLQ / integer relations).
3. Prefer lattice-native questions (matching, critical polynomials, site-to-bond maps) over numerology.

## Working value

| Quantity | Value | Source |
|---|---|---|
| $p_c$ square site | `0.59274605079210(2)` | Jacobsen, graph polynomial / transfer matrix (2015) |
| $1-p_c$ matching site | `0.40725394920790(2)` | exact matching relation |
| $p_c$ square bond | `1/2` | Kesten |
| $p_c$ triangular site | `1/2` | self-matching |

Use Jacobsen's figure unless a later independent computation moves the last digits.

## Repo layout

```
constants/     accepted numerical values and citations
notes/         working notes, excluded expressions
scripts/       small reproducible checks
```

## Status

Initialized 2026-08-28. First issues track exclusion work and matching-lattice notes.
