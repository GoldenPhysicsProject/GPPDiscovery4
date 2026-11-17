# ATLAS.md — GPPDiscovery4 (Discovery 4 — Turbo, D4)

*This repo's atlas — one line per holding, with what it connects to. The full inventory lives in GPPVerify's `GppVerify/Atlas.lean` (atlas.json at that repo's root). Keep this file current: when you add something a future seat must find, add its line here.*

Turbo's own discovery repo, the smallest of the five: 56 files, tools-heavy. Turbo's working rigs and staging live here; his audited results bank in GPPVerify from turbo/stub-closing, and his operations state lives in the bridge and ~/workspace/office.

- `tools/` — 44 files: the standing tooling (rigs, helpers, checkers) → the first place to look before writing a new script.
- `experiments/` — 3 files: numeric runs → receipts for Turbo's filings.
- `lean/` — 2 files: Lean staging → migrates to GPPVerify workbench branches, never banked from here.
- `RESULTS.md`, `results.jsonl` — the result ledger → labelled claims.
- `GUIDE.md`, `README.md` — what this repo is and how it is worked.

## Strings out of this repo

- → GPPVerify workbench `turbo/stub-closing`; Turbo also keeps `turbo/*` task branches
- → GPP-bridge: LANE-STATE-turbo.md, the conversation, TEAM-SOPS.md; ~/workspace/office (DECISIONS.md, FOLLOWUPS.md)
- → GPPVerify Atlas.lean: the full estate inventory
