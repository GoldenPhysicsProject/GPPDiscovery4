# Provenance — 2026-10-13 Weil controls batch

Turbo's numerical workbench output for 2026-10-13, deposited from
`~/workspace/ym-k-experiments/` into the Discovery4 workbench so the rigs and
their results live with the harness they build on. All results **NUMERICAL**.
Nothing in this folder is a claim about RH.

## Rig dependency

Every script here runs on the Ritz/Weil reduced-form harness in the parent
folder (`../ritz_analysis.py`, `../theta_control_ritz_4_spk4.py`), imported via
`sys.path` — the same harness whose reproduction against the 2026-10-08
baseline is exact (see `../REPRODUCTION-2026-10-13.md`). Controls throughout:
the shifted product F_θ(s) = ξ(s−θ)ξ(s+θ) (θ = 1/4 working value) and the
Davenport–Heilbronn function (FE, no Euler product, first off-line zero at
0.8085171825 + 85.6993484854i). Where the other seats produced Arb interval
certificates for the same quantities, that is noted per experiment.

## Contents, with verdicts

- `haar_plancherel_split.py` — Sector split of the finite-model energy form E.
  Verdict: E is a termwise Haar square (archimedean and prime jump-class sums
  each PSD-singular, min eig ~1e-15; ρ(y) > 0 on the whole panel). All
  θ-discrimination is Plancherel-side (reduced spectrum vs the W metric).
  "Haar builds; Plancherel judges."

- `ground_vector.py`, `ground_vector_bundle.py`, `ground_vector.md`,
  `ground_vector_L*_th*.csv` — Is there a single ground/null direction at
  θ = 0? Verdict: NO. The bottom is a near-degenerate bundle at λ = 1
  (7/12/12 eigenvalues within 1e-6 at L = 1/2/4; neutral dim ~29). Per-vector
  profiles are basis artifacts.

- `bundle_mechanism.py`, `bundle_mechanism.md`,
  `bundle_mechanism_results.txt` — What is the bundle? Verdict: Claude Code's
  prolate reading confirmed. Bundle vectors' Φ-weighted Fourier transforms
  have their first eight minima on the ζ ordinates (max deviation 0.13) —
  the Connes–Consani–Moscovici signature. Bundle dimension is flat in L
  (12 for L = 1.5–4); the decisive variables are Φ-width / bandwidth, not L.

- `bundle_counting.py`, `bundle_counting.md`, `bundle_counting_results.txt`,
  `bundle_counting_data.npz` — Claude Code's counting law for the bundle.
  Verdict: CONFIRMED with measured Ω (≈ kmax + 12 ± 2): over Ω ∈ [30, 61],
  rank(D − bundle) = 2N(Ω) + 2, i.e. bundle = D − 2N(Ω) − 2 (±1); crossings
  at γ₅…γ₉ step the rank by exactly +2 each. Reported to the seats
  (GPP-bridge PR #11, 2026-10-13).

- `dh_control.py`, `dh_results.md`, `dh_results.json`,
  `dh_coefficients.txt` — Davenport–Heilbronn control. FE verified first
  (worst residual 4.8e-18); identity ‖E_DH − Z_DH‖/‖Z_DH‖ = 4.7e-7 at L = 4.
  Verdict: DH saturates like ζ at working windows (bottom ≈ 0 for L ≥ 2), so
  window saturation is FE + low zeros, NOT an Euler-product certificate; the
  rig separates DH only at small L or with a modulated carrier aimed at
  t ≈ 85.7. The modulated separation was later CERTIFIED by Codex's exact
  rational witnesses + Claude Code's Arb certificates (−0.64137 at h=3,k=3;
  −1.11918 at h=7/2,k=4). `dh_coefficients.txt`: b(n) for −f′/f, n ≤ 500
  (b(2)=0.1969, b(3)=−0.3121, b(5)=0; mixed signs on [1,403] = no Euler
  product).

- `ftheta_sector_split.py`, `ftheta_sector_results.md` — Codex's assigned
  test (PR #11, 2026-10-13): where does the completed archimedean sewing
  break for F_θ? Verdict: nowhere locally. Both completed sectors remain
  termwise Haar squares (full-space sector bottoms ≥ −1.2e-15 at all θ);
  the defect λ_min − 1 = −0.0520θ² exists only in the complete sewn form,
  carried by one direction with energy 76.7% archimedean / 23.3% prime
  (constant in θ). DH shows the opposite anatomy (prime sector alone
  indefinite, min eig −0.436; completed form saturates). Reported to the
  seats (PR #11 comment, 2026-10-13).

- `zeta_ordinates_140.npy` — ζ zero ordinates used by the bundle scripts.

## Reporting trail

Bridge notes: `MuseResearch/turbo-note-2026-10-13-weil-window-scan.md`
(+ corrections), `turbo-note-2026-10-13-haar-positivity-final-form.md`
(+ codas), GPP-bridge `CONVERSATION.md` entries of 2026-10-13, and PR #11
comments of the same date. Two claims from the morning were corrected the
same day and the corrections stand in those records: the windowed ζ bottom
is positive-but-exponentially-small (prolate), not exactly 0; and the
"purely odd defect direction / n ≤ 5 dominance" observations were
rig-specific, not universal.
