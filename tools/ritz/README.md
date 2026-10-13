# ritz/ — fixed-window Ritz machinery

**Source:** Claude's Ritz scripts, ingested 2026-10-08 (Ritz tests in the RH
sprint's fixed-window program); staged verbatim from the Discovery4 staging
bundle (2026-10-09), with `sys.path` made self-contained.

**Epistemic status:** NUMERICAL. Prime vs zero matrices agree to ~1e-6 in
these tests; the F_θ discriminator follows the law λ_min − 1 ≈ −0.052θ².
Nothing here is a proof; this is the program's standard falsification rig —
every candidate positivity argument must survive it, and the shifted product
F_θ(s) = Ξ(s−θ)Ξ(s+θ) is the adversarial control.

| File | What it does |
|---|---|
| `theta_control_ritz_4_spk4.py` | Base harness: builds the pole-neutral Ritz form E, the prime form Ŵ, and the neutral basis; run first. θ-scan prints λ_min − 1 per θ. |
| `ritz_analysis.py` | Analysis layer: `setup`, `trigbasis`, loss profiles per prime jump. |
| `test23.py` | Tests 2–3: loss-of-positivity prime profiles and the degree-two residual prediction. |

**Recorded results:** the original 2026-10-08 runs are summarized in the
staging bundle's `results/ritz-numerical-results-2026-10-08.md`; a fresh
reproduction from this directory is recorded in `REPRODUCTION-2026-10-13.md`.

Requires numpy (Legendre quadrature via `numpy.polynomial.legendre`).
