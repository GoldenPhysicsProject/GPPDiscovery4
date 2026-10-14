# Provenance — Yakaboylu reflection-pairing check (2026-10-14)

Turbo, at Daniel's request ("take his ideas and apply what we know"), after the audit of
Yakaboylu arXiv:2408.15135v17 (bridge MuseResearch/turbo-note-2026-10-14-yakaboylu-v17-audit.md).

- `yakaboylu_pair_check.py`: runnable script (~13 s), built on the existing rigs in
  `~/workspace/ym-k-experiments/` (weil_window_scan / ftheta_sector_split conventions,
  `zeta_ordinates_140.npy`). Verified: per-pair identity pair = 2(|m|^2 - |d|^2) to 5.4e-20;
  antisymmetric mass A = alpha0 * theta^2 (alpha0 = 2.5938e-02 carrier, 7.0423e-03 control);
  carrier A/S = 131.6 at theta = 0.25 (99.2% of negative pair mass is mirror-difference mass);
  2(S-A) reproduces the rig zero-side matrix c'Zc to ~2e-11; banked c2 = -0.0520 = -2*alpha0 to 3 digits.
- `yakaboylu_pair_results.md`: full result sheet with setup, tables, and caveats. Label: NUMERICAL.

Communicated to the seats: bridge CONVERSATION.md (commit 98c94a8) and PR #11 comment 6171970517.
