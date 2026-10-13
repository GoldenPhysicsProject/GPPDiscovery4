# Ground neutral vector of the Ritz/Weil reduced form — NUMERICAL

Script: `ground_vector.py` (+ supplement `ground_vector_bundle.py`).
Machinery: `ritz_analysis.setup` + `neutral_basis`, trig basis
`trigbasis(30, 0.75, (0,1,2))` (243-dim), h=0.01, ycut=6.0.
CSVs: `ground_vector_L{1,2,4}_th0.csv`, `ground_vector_L4_th025.csv` (u, v).

## Headline

There is **no isolated ground direction**. At θ=0 the bottom of the reduced
spectrum is a **near-degenerate bundle** at λ=1: 10–12 of the 30 neutral
directions have |λ−1| < 1e-6 at every L (6–11 within 1e-7). "Everything = 0"
in this rig is not a single null vector; it is a fat near-null subspace — the
signature of an exact operator identity E ≡ W on the pole-neutral space,
discretized. The single vector U[:,0] returned by the eigensolver is an
arbitrary orthonormal-basis choice inside that bundle, and its profile is
therefore **not a stable invariant** (details below).

## Step 1 — bottom vector at θ=0 (per-L, as requested)

| L | λ_min−1 | λ_2−1 | neutral dim | mass center |u|<0.5 | mass log2-band | mass 0.85–1.5 | mass tail >1.5 | ‖even‖/‖v‖ | ‖odd‖/‖v‖ |
|---|---------|-------|------|------|------|------|------|------|------|------|
| 1.0 | −2.02e-08 | −1.77e-08 | 28 | 0.000 | 0.002 | 0.998 | 0.000 | 0.707 | 0.707 |
| 2.0 | −2.37e-08 | −3.82e-10 | 30 | 0.000 | 0.000 | 0.860 | 0.140 | 0.942 | 0.335 |
| 4.0 | −1.33e-07 | −1.40e-09 | 30 | 0.000 | 0.000 | 0.462 | 0.538 | 0.174 | 0.985 |

Mass lives entirely at |u| ≳ 0.85 (nothing at center, nothing at ±log 2),
migrating outward with L. Even/odd content flips between L values.
Shape correlations on common support: L1–L2: −0.337, L2–L4: +0.495,
L1–L4: +0.373 — **no L-independent profile stabilizes**.

## Step 2 — closed-form candidate overlaps (L=4, θ=0)

| candidate | corr |
|---|---|
| Φ(u) | +0.019 |
| Φ(u)·cosh(u/2) | +0.020 |
| Φ′(u) | −0.000 |
| best-fit Gaussian | σ → 3.0 (boundary of search grid), |corr| = 0.054 |

No candidate matches. Consistent with the degeneracy diagnosis: U[:,0] is a
random direction in the bundle, so overlap ≈ 0 is expected, not informative
about the bundle itself.

## Step 3 — θ=0.25 contrast (L=4)

λ_min−1 = **−3.254e-03** (λ_2−1 = −3.252e-03: a near-degenerate split pair).
Bottom vector: corr with the θ=0 vector = **+0.127**; **purely odd**
(‖odd‖/‖v‖ = 1.0000, ‖even‖ = 0.0002); **100% of mass in the tails |u|>1.5**;
v(0) ≈ 0. E-expectation 0.996746 vs W-expectation 1.000000.

## Step 4 — jump fingerprint

Single ground vector (L=4, θ=0): E-exp = 0.99999987, Wm-exp = 0.99999983,
(E−W)-exp = +4.0e-08. Split: **arch 29.92%, prime 70.08%** of E.
Top prime jumps: n=3: 0.4577, n=4: 0.1222, n=2: 0.0786, n=5: 0.0422,
n=7: 1.2e-04, then a cliff (n=8: 8.7e-07, n=11: 2.3e-12).

Bundle-averaged (12 near-null dirs, normalized, L=4): **arch 46.5%,
prime 53.5%**; n=3: 0.2588, n=2: 0.2314, n=4: 0.0363, n=5: 0.0088,
n=7: 1.9e-05. The robust bundle-level facts: only prime powers n ∈ {2,3,4,5}
carry weight; n=3 and n=2 dominate; everything n ≥ 7 is ≤ 2e-5.

## Honest read

The null direction is **not a recognizable closed form — because it is not a
direction at all**. The exact identity behind "everything = 0" shows up in
this discretization as a ~12-dimensional near-null bundle at λ=1 (a discrete
shadow of E ≡ W as operators on the neutral space), so any single-vector
extraction is an artifact of the eigensolver's basis choice inside the
bundle: no stable profile, no parity, no closed-form overlap. What *is*
stable and formalization-shaped is the bundle itself plus its fingerprint:
the balance closes using only the smallest prime powers (2, 3, 4, 5), with
the archimedean class contributing roughly half at bundle level. A concrete
target would be: prove the reduced form is the identity on the neutral space
modulo the span of directions whose jump support is confined to n ≤ 5 plus
arch — i.e., characterize the cokernel/kernel of (E − W) restricted to the
neutral space, not a ground-state profile. Label: NUMERICAL throughout.
