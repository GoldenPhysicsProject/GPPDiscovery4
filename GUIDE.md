# Guide to GPPDiscovery4

GPPDiscovery4 is Turbo's seat workbench for the Golden Physics Project, alongside GPPDiscovery (Claude), GPPDiscovery2 (Codex), and GPPDiscovery3 (Astra). Its standing threads are the arithmetic Fock formalization, the GPPVerify stub-closing campaign, and the RH fixed-window program.

This guide is the map: where things are, how results are recorded, and how
the repository relates to the rest of the Golden Physics Project. It
describes the layout as of 2026-11-09. The workbench moves; the map is
revised when the layout changes.

## Start here

- `README.md`: the repository's own introduction.
- `RESULTS.md`: the human-readable results ledger.
- `results.jsonl`: the raw result rows behind the ledger.

## Layout

- `lean/`: snapshots of kernel-checked Lean formalizations, each with build provenance (what was proved and how it was verified). Live development happens on Turbo's workbench branch of GPPVerify; this repository holds the shareable snapshots.
- `experiments/`: runnable numeric and symbolic experiment scripts.
- `tools/`: the workbench tooling; the Ritz controls live in `tools/ritz/`.

## Conventions

Results carry epistemic labels: THEOREM, PROVED LEMMA, CONDITIONAL,
NUMERICAL, OPEN, KILLED. Numerical work is evidence, not proof; what
hardens is promoted to
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), the
kernel-checked Lean 4 library, through the build seat's audit. Candidate
positivity mechanisms are tested against the controls F_θ = ξ(s−θ)ξ(s+θ)
and the Davenport–Heilbronn function before they are taken seriously.
Working notes, heuristics, and internal strategy live in the private
GPP-bridge repository; this repository is public by design and is the
publication of the workbench's math, code, and results.

## The family

- [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify): the formalization library.
- [GPPDiscovery](https://github.com/GoldenPhysicsProject/GPPDiscovery): Claude's workbench.
- [GPPDiscovery2](https://github.com/GoldenPhysicsProject/GPPDiscovery2): Codex's workbench.
- [GPPDiscovery3](https://github.com/GoldenPhysicsProject/GPPDiscovery3): Grok's (Astra's) workbench.
- [GPPDiscovery4](https://github.com/GoldenPhysicsProject/GPPDiscovery4): Turbo's workbench.
- [GPPDiscovery5](https://github.com/GoldenPhysicsProject/GPPDiscovery5): Grok's lab (guest seats).

License: Apache-2.0. Citation: see `CITATION.cff`. Author: Daniel Toupin,
Golden Physics Project.
