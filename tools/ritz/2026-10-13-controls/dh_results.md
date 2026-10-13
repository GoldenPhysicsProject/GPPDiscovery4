# Davenport–Heilbronn control for the Ritz/Weil rig — results

**Turbo, 2026-10-13. Label: NUMERICAL throughout** (finite basis, finite grids, mpmath root-finding).
Script: `dh_control.py` (this directory; raw numbers in `dh_results.json`).
Question posed: the rig's reduced form saturates exactly for ζ (banked: λ_min − 1 = 0 within 1.3e-7 for L ≥ 1.5, pole normalization) and fails for the shifted product F_θ (−0.052θ² law). Does that saturation use only the **functional equation** (DH has one too) or the **Euler product / local unitarity** (DH lacks one and has off-line zeros)?

## 1. The DH function and its functional equation — VERIFIED

f(s) = Σ a(n) n^(−s), a 5-periodic: a = (0, 1, x, −x, −1) for residues (0,1,2,3,4).
With χ mod 5, χ(2)=i: Gauss sum τ(χ) = −1.1755705045849463 + 1.9021130325903070i (|τ|² = 5 exactly),
root number ε(χ) = τ(χ)/(i√5) = 0.8506508083520399 + 0.5257311121191336i.
The FE F(s) = F(1−s) for F(s) = 5^(s/2) π^(−(s+1)/2) Γ((s+1)/2) f(s), together with real coefficients,
forces β/α = ε(χ) for f = αL(s,χ) + βL(s,χ̄), i.e.

**x = +0.28407904384041227** (= i(1−ε)/(1+ε); a first derivation with ε̄ in place of ε gave −x and an FE failure at 1e-2 — caught by the verification gate, corrected, re-verified).

FE residuals |F(s) − F(1−s)| at 9 off-line points (dps 50): worst **4.8e-18** — far under the 1e-8 gate. PASS.
f is entire (Σa = 0), so: **no pole weight W, no pole constraints** in the rig adaptation (documented choice; the ζ rig's W/cons exist only because ξ has the s(s−1) zeros at 0,1).

## 2. Zeros actually found

Strip Re ∈ [0.2, 1.3], Im ∈ [0.1, 40]: argument-principle count = **15.000**; off-line half-rectangle counts Re∈[0.55,1.3] and Re∈[0.2,0.45] both **0.00**. findroot located exactly those 15 (max |f| residual 9.7e-15), **all on the critical line** (Re = 0.5000000000 to 10 digits):

5.0941598446, 8.9399144081, 12.1335454258, 14.4040031123, 17.1302394006, 19.3088001742, 22.1597077651, 23.3453701120, 26.0949673462, 27.9237988216, 30.1594179063, 31.9644999712, 33.6998622816, 35.8908547201, 37.4554623430 (plus 40.1625784769 just above the band).

**There is NO off-line pair inside the strip** — the expectation of one near Re ≈ 0.81 at these heights is not borne out. The off-line zeros are higher: argument count in Re∈[0.55,1.3], Im∈[40,200] = **4.00**, all located (each implies its FE mirror at 1−Re and conjugates):

| zero | |f| |
|---|---|
| **0.8085171825 + 85.6993484854i** | 7.6e-15 |
| 0.6508300806 + 114.1633427308i | 5.0e-16 |
| 0.5743560505 + 166.4793059132i | 3.7e-15 |
| 0.7242576946 + 176.7024612429i | 9.3e-15 |

The famous Re ≈ 0.8085 pair is real — at height **85.70**, not below 40.

## 3. Prime-side coefficients b(n)

−f′/f = Σ b(n) n^(−s), b = (a·log) ∗ a^(−1) (Dirichlet convolution, exact divisor sums; an early pair-enumeration bug that silently zeroed a^(−1) was caught by hand-checking b(4)).
Checks: b(p) = a(p)·log p at primes (b(2)=0.19691, b(3)=−0.31209, b(5)=0, b(11)=2.39790, b(19)=−2.94444).
On [1, 403]: b ranges over **[−7.68, +7.92]** — mixed signs, support everywhere (b(n)=0 only when 5|n). Against ζ's Λ(n) ≥ 0: the DH prime side is **not** a sum of positive jumps; no Euler product, as expected.

## 4. Kernel and the explicit-formula identity — VERIFIED

Φ_DH(u) = 2e^(3|u|/2) Σ_n a(n)(n/√5) exp(−πn²e^(2|u|)/5) (theta-series inverse Mellin).
Mellin spot checks ∫Φ_DH e^((s−1/2)u) du vs F(s): max rel err **2.4e-15**.
Archimedean density ρ_DH(y) = e^(−3y/2)/(1−e^(−2y)) = Σ over Γ((s+1)/2) poles s=−(2k+1) of e^((pole−1/2)y) — the same rule that reproduces ζ's ρ(y) from Γ(s/2).

**Identity check (L=4, trig basis 243, ycut 6): ‖E_DH − Z_DH‖/‖Z_DH‖ = 4.715e-07** — the DH explicit formula closes in the rig at the same level as ζ's 2026-10-08 test 1 (~4e-7). Adding the four high off-line zeros to Z changes this not at all (identical digits): their contribution through the kernel smoothing is below the floor.

## 5. The reduced-form scan (arch normalization — see note)

f entire ⇒ no pole scale exists. The reduced quantity used is μ(L) = min generalized eigenvalue of (E, A_arch), A_arch = archimedean-only Gram (the positive trivial scale), computed **identically for ζ** as the control column. Banked pole-normalized references: ζ λ_min−1 = 0 (±1.3e-7, L≥1.5); F_0.25 λ_min−1 = −3.254e-3.

| L | μ_DH(L) | μ_ζ(L) | kept dim DH / ζ (of 243) |
|---|---|---|---|
| 0.75 | **−9.6267e-01** | 1.0000129 | 34 / 26 |
| 1.00 | **−2.2769e-01** | 1.0000329 | 37 / 28 |
| 1.50 | **−1.0961e-04** | 1.0000334 | 41 / 29 |
| 2.00 | −3.19e-11 | 1.0000334 | 43 / 29 |
| 3.00 | −1.31e-10 | 1.0000334 | 43 / 29 |
| 4.00 | −1.12e-10 | 1.0000334 | 43 / 29 |

Conditioning, stated plainly: the arch Gram is low-rank (kept 43/243 for DH at L≥2, cond ≈ 2.4e10 after the 1e-11 cut). Empirical noise floor for μ: Z_40 is PSD by construction (on-line zeros only) yet its computed bottom is −3.9e-10, so |μ| ≲ 4e-10 is zero. μ_DH at L ≥ 2 (−3e-11…−1.3e-10) is **at that floor**; μ_DH at L ≤ 1.5 is 3–9 orders above it — real.

Reading of the columns: for ζ, P_ζ ⪰ 0 forces μ_ζ ≥ 1 structurally, and it sits at 1 + 3.3e-5 (bottom direction almost prime-free). For DH the prime part is indefinite, and the bottom plunges at small windows — but the plunge **decays with L and is gone by L = 2**.

## 6. Verdict

**At the rig's working windows (L ≥ 2), DH saturates like ζ: bottom = 0 within ~1e-10.** The mechanism is visible in the data: DH's low zeros (all of them, up to height 40) are on the line, so its zero-side form is a sum of squares and the explicit-formula cancellation closes (identity verified at 4.7e-7); its off-line zeros sit at height ≥ 85.70, above the Mellin content of L ≤ 4 test functions, and are invisible to the form (adding them changes nothing). DH separates from ζ only at small windows (μ = −0.96, −0.23, −1.1e-4 at L = 0.75, 1, 1.5, vs ζ's flat 1.00003) — the same small-window zone already flagged unquotable for ζ.

So: **the exact saturation is not, by itself, an Euler-product certificate.** It is what FE + on-line low zeros produce. The contrast with F_θ is one of *height*, not of kind: F_θ's off-line zeros are θ-shifts of ζ's low zeros, present at every height, so F_θ fails at every window (−0.052θ²); DH's first off-line zero is at 85.70, so DH passes every window the rig can currently build. A discriminator that sees DH's defect needs test functions with Mellin content at t ≈ 86 (much larger effective window / oscillatory basis), not a better bottom-eigenvalue at L ≤ 4.
