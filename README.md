# GPPDiscovery4

Turbo's workbench for the Golden Physics Project — the seat Daniel assigned to
Turbo, alongside GPPDiscovery (Claude), GPPDiscovery2 (Codex), GPPDiscovery3 (Grok/Astra).

Follows the Codex/GPPDiscovery2 process conventions so every seat can find and
use everything here. Read `HEURISTICS.md` first, then
`discovery/ACTIVE_FRONTS_MAP_2026-10-12.md`.

## Layout

- `research/` — dated notes by kind: `FORMALIZATION_UPDATE_`, `TURBO_UPDATE_`
  (relay to Daniel), `DRIVE_MINING_DELTA_`, `NOTE_`. Epistemic label up front;
  only independently checkable mathematics is promoted.
- `discovery/` — one directory per active front (`FRONT.md` each) plus the
  integrated `ACTIVE_FRONTS_MAP_<date>.md` carry-forward program.
- `experiments/` — runnable scripts, never pasted outputs.
- `lean/` — kernel-checked Lean formalizations with build provenance.
- `results.jsonl` — append-only machine-readable record of every result.
- `HEURISTICS.md` — how this seat works; process changelog included.

## Notes on notes

Detailed research notes are ALSO deposited in `MuseResearch/` on the private
`GPP-bridge` repo and in Daniel's Drive (`Newer rh physics stuff`), per standing
instruction. This repo is the seat-visible workbench layer.

## Related branches

- `turbo/stub-closing` on `GoldenPhysicsProject/GPPVerify` — the campaign to make
  the project's big claims literally true in Lean and close `open_` stubs with
  real proofs. Workbench-only; nothing merges to `main` without Daniel's review.
