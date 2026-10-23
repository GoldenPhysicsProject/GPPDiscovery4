# F_θ archimedean-sewing sector test — results

Numerics, 2026-10-13. Test: explicit θ branches in the paired
entropy/zero-count rig; isolate where completed Arch sewing breaks.
Script: `ftheta_sector_split.py` (rerunnable, ~80 s; copy-and-adapt of
`haar_plancherel_split.py` / `ground_vector.py`, originals untouched; the DH
part is copied verbatim from `dh_control.py` with one extension noted below).
**Label: NUMERICAL throughout. Nothing here is a claim about RH.**

## Setup (exact)

- Rig: `ritz_analysis.setup` + `neutral_basis` (GPPDiscovery4 `tools/ritz/`),
  kernel machinery from `theta_control_ritz_4_spk4.py` — identical to the
  ζ-baseline runs (`haar_plancherel_split.py`, `ground_vector.py`).
- Basis: `trigbasis(30, 0.75, (0,1,2))` — 243 functions {x^j, x^j cos kx,
  x^j sin kx, j=0,1,2, k = 0.75…30}. Window **L = 4.0**, grid h = 0.01
  (801 points), prime cutoff **ycut = 6.0** (jumps at y = log n, n ≤ 403).
- Form: E = Σ_jumps M_j, M_j = D_jᵀ diag(Φ(u)Φ(u−y_j)·h·wt_j) D_j, D_j the
  shift-difference of the basis. Reduced form: eigenvalues of BᵀEB on the
  pole-neutral space (B whitens the pole weight Wm, kills the pole-moment
  constraints Cm); λ_min − 1 is the defect. ζ: λ_min − 1 = 0.
- Sectors (from the labeled per-jump list):
  - **arch**: 40-point Gauss panel of the archimedean weight ρ(y)·2cosh(θy),
    ρ(y) = e^{−y/2}/(1−e^{−2y}).
  - **prime**: jumps at y = log n with weight Λ(n)/√n · 2cosh(θ log n)
    (F_θ's two shifted branches enter exactly through this cosh factor and
    through the completed kernel Φ_θ = (Φe^{−θ·}) ∗ (Φe^{+θ·})).
  - **pole**: the weight W (ζ: (C/2)cosh(u/2)Φ; F_θ: two-branch,
    coefficients ξ(2θ), ξ(1+2θ), pole functions cosh/sinh((1/2∓θ)u)) — it is
    the metric of the reduced form, not a jump sector.
- θ values: 0 (ζ), 0.05, 0.10, 0.15, 0.20, 0.25 (F_θ); headline θ = 1/4.
- Sanity gates, all passed: sector sums reproduce E exactly
  (‖E − (E_a+E_p)‖/‖E‖ ≈ 1–2e-16 at every θ); W-whitening residual
  ‖BᵀWB − I‖ ≈ 0.7–3.6e-6; DH Mellin gate max rel err 2.4e-15.

## Sector table (L = 4.0)

Full-space bottom eigenvalues (termwise Haar-square check; PSD ⟺ ≥ 0) and
neutral-space quantities. Sector bottoms are min eigenvalues of the
sector-restricted reduced forms BᵀE_aB, BᵀE_pB (Rayleigh quotients against
the W-metric; 1 = threshold). They do **not** sum to λ_min.

| quantity | ζ (θ=0) | F_{1/4} (θ=0.25) |
|---|---|---|
| min eig E_arch (full space) | −1.65e-15 | −7.5e-16 |
| min eig E_prime (full space) | −1.2e-16 | −1.6e-16 |
| min eig E (full space) | −2.1e-15 | −9.5e-16 |
| neutral dim | 30 | 42 |
| sector bottom, arch (vs W) | 0.2745 | 0.1338 |
| sector bottom, prime (vs W) | 5.04e-05 | 6.43e-03 |
| **λ_min of sewn form** | **0.99999987** | **0.99675** |
| **λ_min − 1** | **−1.3e-07 (≈ 0)** | **−3.25e-03** |
| bottom-dir energy share, arch | 29.9% | 76.7% |
| bottom-dir energy share, prime | 70.1% | 23.3% |
| arch-cheapest dir: its λ (sewn) | 1.000265 | 1.000145 |
| prime-cheapest dir: its λ (sewn) | 1.7443 | 1.4185 |
| sewing functional E_a[c_a]+cross+E_p[c_p] | 2.745 | 2.419 |

Notes on the table:

- Both sectors are PSD in the full space for F_θ exactly as for ζ
  (bottoms ≥ −1.2e-15, i.e. zero at float precision): each F_θ jump class is
  still a termwise Haar square. **Locally, nothing breaks anywhere.**
- One-sided sector bottoms against W are ≪ 1 for **both** ζ and F_θ: each
  sector alone has near-null neutral directions (its mass is masked by the
  other sector there — along the prime-cheapest direction the archimedean
  expectation is 1.41–1.74, and vice versa the arch-cheapest direction
  carries prime expectation 0.73–0.87). A one-sided sector bottom < 1 is
  therefore **not** a defect location; it does not discriminate ζ from F_θ.
- The pair of sector-cheapest directions also cannot assemble the defect:
  the sewing functional E_a[c_a] + (c_aᵀE_p c_a + c_pᵀE_a c_p) + E_p[c_p]
  = 2.42 > 2 at θ=1/4 (2.74 for ζ). The F_θ bottom direction is a genuinely
  different, third direction.
- The θ=0 bottom direction's shares (29.9/70.1) are **not** comparable to
  the θ>0 profile: at θ=0 the bottom is the known near-degenerate prolate
  bundle and any basis vector in it is a basis artifact. For θ ≥ 0.05 the
  bottom is isolated and its profile is θ-independent (shares stable to
  ±0.1%: arch 76.6–76.8%); the θ→0⁺ limit of the defect direction is
  archimedean-heavy, 76.7% / 23.3%.

## θ sweep (full sewn form)

| θ | λ_min − 1 | (λ_min − 1)/θ² | sector bottoms arch / prime | bottom-dir shares arch / prime |
|---|---|---|---|---|
| 0.00 | −1.34e-07 | — | 0.2745 / 5.04e-05 | 29.9% / 70.1% (bundle) |
| 0.05 | −1.294e-04 | −0.0518 | 0.1351 / 6.26e-03 | 76.7% / 23.3% |
| 0.10 | −5.191e-04 | −0.0519 | 0.1350 / 6.28e-03 | 76.6% / 23.4% |
| 0.15 | −1.169e-03 | −0.0520 | 0.1347 / 6.31e-03 | 76.6% / 23.4% |
| 0.20 | −2.081e-03 | −0.0520 | 0.1343 / 6.36e-03 | 76.7% / 23.3% |
| 0.25 | −3.254e-03 | −0.0521 | 0.1338 / 6.43e-03 | 76.7% / 23.3% |

Fits: λ_min − 1 ≈ c₂θ² with **c₂ = −0.0520** (one-parameter; the banked
pole-neutral control law). Two-parameter fit c₂θ² + c₄θ⁴: c₂ = −0.0519,
c₄ = −0.0022. The defect is quadratic to the precision of the rig, the
sector bottoms are θ-independent to ~1%, and the defect direction's sector
profile is constant across the sweep.

## DH control (same L, h, ycut; copied from `dh_control.py`)

DH is entire — no pole weight W, no neutral space — so only the full-space
PSD check and the arch-scale (A = arch-only Gram) generalized bottoms
transfer. Its prime weights b(n) (from −f′/f) have **mixed signs**
(range [−7.68, +7.92] on n ≤ 403), so its prime sector is not a Haar square.

| quantity | DH |
|---|---|
| min eig E_arch (full space) | −1.4e-15 (≈ 0, PSD) |
| min eig E_prime (full space) | **−0.436** (indefinite) |
| min eig E (full space) | −2.3e-15 (≈ 0) |
| μ(E) vs arch scale | −1.1e-10 (≈ 0, saturates) |
| μ(E_prime only) vs arch scale | **−1.0000** |

So DH's anatomy is the opposite of F_θ's: its prime sector **alone** is
genuinely indefinite (sector-local break, visible already in the full
space), while its completed form still saturates at this window
(μ(E) ≈ 0 — its first off-line zero is at height ≈ 85.7, unresolved here,
consistent with the earlier DH verdict). The sector test therefore *can*
see a sector-local break when one exists — and for F_θ there is none.

## Verdict (NUMERICAL)

In this rig, the F_θ negative mass lives in **no single sector and in no
cross-sector combination of the sectors' own cheapest directions — only in
the complete sewn form**. Termwise, both completed sectors of F_θ remain
Haar squares exactly as for ζ (full-space sector bottoms ≥ 0 at float
precision for every θ), and the two sector-cheapest neutral directions
cannot assemble a defect (sewing functional 2.42 > 2 at θ = 1/4; each such
direction has sewn λ ≥ 1.0001). The failure appears only in the sum
E_arch + E_prime measured against the two-branch pole weight W on the
pole-neutral space, as a clean quadratic defect λ_min − 1 = −0.0520 θ²,
carried by a single isolated direction whose energy is **76.7%
archimedean / 23.3% prime**, constant across θ ∈ [0.05, 0.25]. In sector
terms: F_θ's local pieces are intact everywhere (consistent with the
prime-local TFD/Haar-square picture); the completed sewing breaks on an
**archimedean-heavy direction of the completed form** — the shifted
two-branch Γ-factor/pole completion is where the −0.052θ² defect
concentrates. The DH control confirms the test's resolving power by
showing the contrary anatomy (a genuinely indefinite prime sector) for a
function whose failure is local rather than sewn.

## Conditioning / noise floor / what the rig could not resolve

- **Clean regime**: all ζ/F_θ numbers at L = 4.0, inside the rig's clean
  window (the small-L zone, L ≲ 1, is conditioning-flagged in this rig
  family and was not used). Dense-eigensolver backward-error scale for the
  reduced forms ≈ 2e-13; the W-whitening residual ‖BᵀWB − I‖ ≈ 1–4e-6 is
  the practical accuracy scale for λ — ζ's λ_min − 1 = −1.3e-7 sits
  *below* it (numerical zero, consistent with the banked 1.3e-7 tolerance),
  while the smallest F_θ defect used, −1.29e-04 at θ=0.05, is ~40× above
  it. DH's arch Gram is genuinely ill-conditioned (kept 43 of 243,
  cond 2.4e10, ridge 1e-11 as in `dh_control.py`); μ(E) = −1.1e-10 should
  be read as "zero at DH's noise floor", and μ(E_prime) = −1.0000 and
  min eig E_prime = −0.436 are far above any floor.
- **Sector bottoms are not cross-θ comparable in fine detail**: the
  neutral space changes with θ (dim 30 at θ=0 with 2 pole constraints,
  dim 42 at θ>0 with 4), so sector bottoms are quoted as qualitative
  anatomy, not as a θ-family.
- **Not resolved — a sector split of the coefficient c₂ itself.** The
  shares above decompose the defect *direction's energy*, not the
  eigenvalue defect: c₂ is a property of the sewn form, and splitting it
  sector-wise by perturbation theory would need a θ-independent neutral
  reference space (the θ=0 bottom is a near-degenerate bundle and the
  neutral spaces are θ-dependent), which this rig does not provide.
  What is established is where the defect direction's mass sits
  (arch-heavy, 3.3 : 1) and that no sector or sector-pair carries the
  defect alone.
- The "archimedean sector" is the 40-point Gauss-panel discretization of
  the continuous archimedean weight used throughout this rig family;
  panel-refinement was not varied in this run (the ζ baseline and all
  banked control laws use the same panel, so the comparison is matched).
