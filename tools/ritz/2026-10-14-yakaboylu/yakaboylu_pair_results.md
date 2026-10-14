# Yakaboylu reflection-pairing check — results

Turbo, 2026-10-14. Script: `yakaboylu_pair_check.py` (runs in ~13 s).
**Label: NUMERICAL throughout. Nothing here is a claim about RH.**

## Question

Yakaboylu (arXiv:2408.15135v17) reduces RH to positivity of a reflection
pairing over the zeros; Turbo's audit of v17 (bridge
`MuseResearch/turbo-note-2026-10-14-yakaboylu-v17-audit.md`) predicts
that for the shifted product F_θ, the Weil zero-side pairing splits per
mirror pair into a symmetric (positive) part and an antisymmetric
(defect) part. This check verifies that decomposition with the existing
rigs and quantifies the antisymmetric mass as a function of θ.

## Setup and conventions (the rigs' own)

- Basis `trigbasis(30, 0.75, (0,1,2))` (243 real functions of u), window
  L = 4.0, grid h = 0.01 on u ∈ [−L, L], ycut = 6.0 — identical to
  `weil_window_scan.py` / `ftheta_sector_split.py`.
- Scalar transform (this is `ritz_analysis.F_matrix`, per basis function),
  for a real test function f = Σ c_j g_j:
  ĝ_θ(z) = ∫ Φ_θ(u) f(u) e^{(z−1/2)u} du, evaluated as the rig's
  rectangle sum. Φ_θ is the F_θ archimedean profile (`make_labeled`);
  Φ_0 is ζ's profile.
- Zero pairs from `zeta_ordinates_140.npy` (γ_k, k = 1..140,
  γ ∈ [14.1347, 302.6967]):
  a_k = 1/2 + θ + iγ_k, b_k = 1/2 − θ + iγ_k = 1 − conj(a_k).
- Pair contribution, m/d split, masses:
  pair(a,b) = ĝ(a)ĝ(1−a) + ĝ(b)ĝ(1−b) = 2 Re[ĝ(a) conj(ĝ(b))],
  m = (ĝ(a)+ĝ(b))/2, d = (ĝ(a)−ĝ(b))/2,
  S = Σ_pairs 2|m|², A = Σ_pairs 2|d|², pair sum = S − A.

Test functions (coefficient vectors held fixed across the θ sweep;
Wm-normalized neutral Ritz directions of the θ = 0.25 rig, c′Wm c = 1):

- **Carrier**: bottom neutral direction at θ = 0.25, L = 4. Completed
  form value λ_min − 1 = **−3.253923e-03** (reproduces the banked
  −3.254e-03; this is the direction where the F_θ completed Weil form
  is negative).
- **Control**: top neutral direction at θ = 0.25. Completed form value
  λ_max − 1 = **+9.393806e-01** (strongly positive).

## Identity verification (θ = 0.25, first 50 pairs)

pair_direct = ĝ(a)ĝ(1−a) + ĝ(b)ĝ(1−b) with ĝ(1−a), ĝ(1−b) evaluated
directly (no conjugation shortcut), against both closed forms:

| function | max \|pair_direct − 2Re[ĝ(a)conj(ĝ(b))]\| | max \|pair_direct − 2(\|m\|²−\|d\|²)\| | max \|ĝ(1−a) − conj(ĝ(b))\| |
|---|---|---|---|
| carrier | 3.388e-21 (rel 3.3e-18) | 5.432e-20 (rel 5.3e-17) | 0 exactly |
| control | 4.337e-19 (rel 2.3e-18) | 5.551e-17 (rel 2.9e-16) | 0 exactly |

Relative errors are against max(|ĝ(a)|²+|ĝ(b)|²) over the 50 pairs.
**The identity pair = 2(|m|² − |d|²) holds to machine precision**, as it
must: it is the parallelogram identity plus the reality relations
ĝ(1−a) = conj(ĝ(b)), ĝ(1−b) = conj(ĝ(a)), which hold exactly in the
discrete sum (Φ_θ, f real). The pair identity itself is θ-independent
algebra; nothing about it is special to θ = 0.25.

## θ = 0 sanity

At θ = 0, a = b = 1/2 + iγ pointwise, and **max |d| = 0 exactly** over
all 140 pairs for both test functions. The antisymmetric channel is
identically closed on the line.

## S / A table (140 pairs)

Carrier (bottom neutral direction):

| θ | S | A | S−A | A/S | A/θ² |
|---|---|---|---|---|---|
| 0.00 | 3.87478956e-02 | 0 | 3.87478956e-02 | 0 | — |
| 0.01 | 1.66597938e-05 | 2.59384698e-06 | 1.40659468e-05 | 1.557e-01 | 2.593847e-02 |
| 0.02 | 1.66342397e-05 | 1.03759469e-05 | 6.25829282e-06 | 6.238e-01 | 2.593987e-02 |
| 0.05 | 1.64567933e-05 | 6.48741288e-05 | −4.84173355e-05 | 3.942e+00 | 2.594965e-02 |
| 0.10 | 1.58435452e-05 | 2.59846242e-04 | −2.44002697e-04 | 1.640e+01 | 2.598462e-02 |
| 0.25 | 1.24565258e-05 | 1.63941885e-03 | −1.62696232e-03 | 1.316e+02 | 2.623070e-02 |

Control (top neutral direction):

| θ | S | A | S−A | A/S | A/θ² |
|---|---|---|---|---|---|
| 0.00 | 2.97060740e+00 | 0 | 2.97060740e+00 | 0 | — |
| 0.01 | 4.66412141e-01 | 7.04234203e-07 | 4.66411437e-01 | 1.510e-06 | 7.042342e-03 |
| 0.02 | 4.66429960e-01 | 2.81703419e-06 | 4.66427143e-01 | 6.040e-06 | 7.042585e-03 |
| 0.05 | 4.66554708e-01 | 1.76107247e-05 | 4.66537098e-01 | 3.775e-05 | 7.044290e-03 |
| 0.10 | 4.67000522e-01 | 7.05038050e-05 | 4.66930018e-01 | 1.510e-04 | 7.050380e-03 |
| 0.25 | 4.70133618e-01 | 4.43322879e-04 | 4.69690295e-01 | 9.430e-04 | 7.093166e-03 |

θ² fit of A (through origin, over the five nonzero θ):
carrier α = 2.622412e-02, control α = 7.092022e-03 (both fits dominated
by the θ = 0.25 point). The cleaner small-θ coefficients are the
θ = 0.01 values of A/θ²: **α₀ = 2.5938e-02 (carrier)**,
**α₀ = 7.0423e-03 (control)**. A/θ² drifts up by only ~1.1% (carrier)
and ~0.7% (control) between θ = 0.01 and θ = 0.25, i.e.
A = α₀θ² + O(θ⁴). Successive log-log slopes dlogA/dlogθ:
carrier 2.000, 2.000, 2.002, 2.010; control 2.000, 2.000, 2.001, 2.007.
**The antisymmetric (defect) mass is O(θ²) for small θ**, for both a
defect-carrying and a positive test function — the θ² law is in the pair
algebra (d is odd in θ), not in the choice of f.

## Reading of the numbers

- S ≥ 0 by construction (sum of squares); all negativity of any pair sum
  comes from A. For the **carrier at θ = 0.25, A/S = 131.6**: the pair
  sum S − A = −1.627e-03 is −A × (1 − 0.0076), i.e. **99.2% of the
  negative mass is antisymmetric** (mirror-difference) mass. The
  carrier's S is itself tiny (1.25e-5 vs the control's 0.470) because
  the carrier is tuned so ĝ nearly vanishes at the pair points — both
  m and d are small there — but d decays only like θ while m is
  suppressed further, so A wins as θ grows.
- The carrier's fixed-f pair sum changes sign between θ = 0.02
  (S−A = +6.26e-06) and θ = 0.05 (S−A = −4.84e-05); interpolating
  A = α₀θ² against S ≈ 1.66e-5 puts the crossing at θ* ≈ 0.025. (This
  is a zero-side, fixed-f statement, not the completed form's window
  crossing.)
- For the **control**, A/S ≤ 9.4e-04 at every θ: the pair sum is
  S-dominated and positive throughout, as expected for a direction
  where the completed form is positive. The d-channel is present
  (A = α₀θ² ≠ 0) but negligible.
- Rig cross-check at θ = 0.25: the rig's own zero-side matrix
  (`zero_side_matrix`, same 140 zeros, both half-planes) gives
  c′Zc = −3.25392466e-03 (carrier) and +9.39380589e-01 (control),
  against 2(S−A) = −3.25392464e-03 and +9.39380589e-01 from the pair
  sums — relative differences 6.1e-09 and 6.2e-10 (absolute ~2e-11;
  floating-point reassociation between the matrix and per-pair routes).
  The factor 2 is the lower half-plane: Z sums over all zeros, the pair
  sum covers the upper half-plane. So the m/d split is a decomposition
  of the rig's actual zero side, not of a toy surrogate.

## Caveats

- **Zero side only.** S − A is the zero-side pair sum, truncated at
  140 pairs (γ ≤ 302.70). It is NOT the full Weil form: the
  archimedean/prime-side balance via the explicit formula (the rigs'
  E − W completed form) is a separate, existing computation and is not
  re-derived here. (For these neutral directions the two agree at
  θ = 0.25 through the rig's E ↔ Z identity, per the cross-check above
  and test 1 of `ritz_analysis.py`; that identity is the rig's
  business, quoted, not reproved.)
- **Truncation**: at θ = 0.25 the first 50 pairs already saturate S
  and A to all printed digits (S50 = S, A50 = A at 8 significant
  figures, both functions), so the 140-pair truncation error in S and
  A is below print precision there; the same transform decay applies
  at the other θ values but was not separately bounded.
- **θ = 0 row scale**: the θ = 0 row uses ζ's own profile Φ_0
  (max |Φ_0| ≈ 0.894 on the grid) while every θ > 0 row uses the F_θ
  product profile (max |Φ_θ| ≈ 0.320, nearly θ-independent). S(0) is
  therefore not on the F_θ family's scale and is not its θ → 0 limit;
  only A = 0 / d = 0 at θ = 0 is quoted from that row, and that is
  profile-independent. Within the F_θ family S is smooth
  (1.666e-05 → 1.246e-05 over θ ∈ [0.01, 0.25] for the carrier).
- Quadrature is the rig's rectangle rule at h = 0.01. The identity
  check is algebraic in the computed ĝ values, so quadrature error
  does not enter it; S and A inherit the rig's quadrature accuracy.
- F_θ has off-line zeros by construction; it is a control. This check
  concerns the pair algebra and its θ-scaling only.
