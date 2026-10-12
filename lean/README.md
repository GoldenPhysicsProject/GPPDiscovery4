# Turbo's Lean workbench

Kernel-checked Lean 4 (v4.33.1, Mathlib v4.33.1) formalizations produced by Turbo.
Every file here carries its verification status honestly: **LEAN-CHECKED** means the
module built with `lake build` exit 0, zero `sorry`/`admit`, zero custom axioms.

## Contents

### `ArithmeticFock/Basic.lean` — algebraic bosonic Fock core (LEAN-CHECKED)

Snapshot (2026-10-11) of the module live at
`RhWaveLean/ArithmeticFock/Basic.lean` in the RH workbench. 413 lines, 24 theorems.

- Occupation configurations `ι →₀ ℕ`, bosonic Fock space `Fock ι := Config ι →₀ ℂ`
  (algebraic, uncompleted — no Hilbert completion yet).
- Creation/annihilation operators in the normalized algebraic convention
  (`a†|n⟩ = |n+1⟩`, `a|n⟩ = nᵢ|n−eᵢ⟩`); all four CCR relations proved:
  `comm_create_annihilate`, `comm_create_annihilate_ne`,
  `comm_annihilate_annihilate`, `comm_create_create`.
- Number operators and eigenvalues (`numberOp_eigenvalue`, `weightedNumber_eigenvalue`).
- Prime-tower identification: `toNat`/`fromNat` round trips between prime-indexed
  configurations and positive naturals via factorization (`fromNat_toNat`,
  `toNat_fromNat`), `sum_mul_log_eq_log_toNat`, and the log-weighted number
  operator acting as multiplication by `log n` (`weightedNumber_log_eigenvalue`).

Verification: `lake build RhWaveLean.ArithmeticFock.Basic` → exit 0
(verified 2026-10-11, twice, independently of the authoring agent).
Axioms: standard only (`propext`, `Classical.choice`, `Quot.sound`).
Build notes: Mathlib 4.33.1 API adaptations were needed (`Finset.prod_map` vs
`prod_image`, `∏ x ∈` binder syntax, `Finsupp.linearCombination_single` arity).

To rebuild: this file is a snapshot; the buildable original lives in the
`prove-the-riemann-hypothesis/lean/` workbench (`lake build RhWaveLean.ArithmeticFock.Basic`).

## Related work (not in this repo)

- **GPPVerify stub-closing campaign** (branch `turbo/stub-closing` on
  `GoldenPhysicsProject/GPPVerify`, workbench-only until Daniel reviews):
  making the shadow–τ identification claims literally true (dimension-aware
  antiholomorphic correspondence, equator corollary) and closing `open_` stubs
  with real proofs. In progress as of 2026-10-12.
