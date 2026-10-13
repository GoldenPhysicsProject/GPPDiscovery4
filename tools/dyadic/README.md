# dyadic/ — dyadic prime-sum harness

**Source:** RH research sprint, waves 21 and 26 (Turbo's sprint agents);
banked in the sprint ledger and PROGRAM.md.

**Epistemic status:** the tiling identities are EXACT; the wave-26 remainder
bound is PROVED by interval certification (exact rational activity decisions,
interval enclosures at kinks) — C₁ ≤ 1.20069 for the compensated dyadic
constant. The wave-21 script is the numerical companion that quantifies why
naive dyadic induction cannot close.

| File | What it does |
|---|---|
| `dyadic_tiling_wave21.py` | Verifies the exact dyadic tiling identities (w_x + w_{2x} = 1) and the deficit identity U(2x) = −U(x) − D(x); the IH quantity enters with a minus sign, which is the induction barrier. |
| `dyadic_tiling_wave21_results.json` | Recorded output of the wave-21 run. |
| `dyadic_cert_wave26.py` | Interval-certified bound on the windowed dyadic remainder B_w(x); sup over x ≥ 2 attained at kinks/endpoints, evaluated with exact Fraction arithmetic. |
| `dyadic_cert_wave26_results.json` | Recorded output of the wave-26 certification run. |

Related: `experiments/2026-10-12_dyadic_antilocking.py` reproduces the
wave-30/31 F₁/F₂ anti-locking signature (Corr ≈ −0.983) at small scale.
