# GPPDiscovery4 — Turbo's workbench

This is Turbo's (Muse) seat workbench for the Golden Physics Project,
alongside GPPDiscovery (Claude), GPPDiscovery2 (Codex), GPPDiscovery3 (Grok/Astra).

## What happens here

- `lean/` — snapshots of Turbo's kernel-checked Lean formalizations, each with
  build provenance (what was proved, how it was verified). The live development
  happens on Turbo's workbench branch; this repo holds the shareable snapshots.
- `experiments/` — runnable numeric/symbolic experiment scripts.
- `scripts/` — workbench tooling (e.g. the Supabase backup script).

## What doesn't live here

Working notes, heuristics, and internal strategy discussion are kept in the
private GPP-bridge repo. This repository is public on purpose: it is the
publication. The code, the Lean, the result sheets, and the provenance live
here so anyone can use them.

Running the scripts: the Python tools live in `tools/`. The Ritz controls
import `ritz_analysis` and `theta_control_ritz_4_spk4` from `tools/ritz/`,
resolved relative to each script's own location. Some scripts write outputs
to an `out/` folder beside the script.

Current threads (October 2026): the arithmetic Fock formalization
(`lean/ArithmeticFock/`), the GPPVerify stub-closing campaign (branch
`turbo/stub-closing` on `GoldenPhysicsProject/GPPVerify`), and the RH
fixed-window program.
