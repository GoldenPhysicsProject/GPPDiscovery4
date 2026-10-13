# symbolic/ — symbolic-check tooling

**Sources:** Codex's 2026-09-30 check trio (via GPP-bridge
`research/codex/scripts/`) and RH sprint wave 70 (Turbo's sprint agents).

**Epistemic status:** computer-algebra verification (PROVED at the symbolic
level, as marked in each file) backed by independent high-precision numeric
evaluation. These are the program's "trust, but verify symbolically" layer
for closed forms that numerics alone cannot certify.

| File | What it checks |
|---|---|
| `check_euler_source_2026_09_30.py` | Euler-source identities for the ground-transform checks; exports `prime_powers` used by the other two. |
| `check_ground_transform_2026_09_30.py` | Ground-transform φ-score evaluations (imports the Euler source module — keep the trio together). |
| `check_neutral_threshold_2026_09_30.py` | Neutral-threshold φ-score evaluations, the δ = 1/2 half-density threshold family. |
| `a6_antiderivative_check_wave70.py` | Sympy verification of the A6 antiderivatives Ψ_r, Ψ_f and the closed forms A_rise = 2 + (2√2−4)/log 2, A_fall = −2 + (4√2−4)/log 2, against 30-digit mpmath quadrature. |

Requires numpy/scipy (trio) and sympy/mpmath (wave-70 check).
