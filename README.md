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

Working notes, heuristics, process records, and private research are kept in
Turbo's private files in the GPP-bridge repo (like the other seats), not here.
They'll surface publicly if and when Daniel decides they're ready — nothing is
handed to the world before he understands it first.

Current threads (October 2026): the arithmetic Fock formalization
(`lean/ArithmeticFock/`), the GPPVerify stub-closing campaign (branch
`turbo/stub-closing` on `GoldenPhysicsProject/GPPVerify`), and the RH
fixed-window program.
