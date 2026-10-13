# tools/ — proven tooling, organized by job

Seeded 2026-10-13 (Turbo seat). Each directory holds scripts that already
produced banked results in the fixed-window RH program, with a short readme
naming the wave/source and the epistemic status of what they check. The
`experiments/` directory remains the home of dated one-off reproductions;
`tools/` is the reusable harness layer.

| Directory | Contents | Status of contents |
|---|---|---|
| `ritz/` | Fixed-window Ritz machinery: prime-side vs zero-side forms, F_θ shifted-product control | NUMERICAL, banked 2026-10-08 |
| `dyadic/` | Dyadic prime-sum harness: exact tiling identities + interval-certified windowed remainder bound | Identities EXACT; bound PROVED (interval-certified, wave 26) |
| `symbolic/` | Symbolic/numeric cross-checks: 2026-09-30 check trio (Codex) + wave-70 sympy antiderivative verification | Computer-algebra PROVED / NUMERICAL as marked per file |

Code and numeric outputs only. No credentials, no unpublished prose —
research notes live in the GPP-bridge (`MuseResearch/`).
