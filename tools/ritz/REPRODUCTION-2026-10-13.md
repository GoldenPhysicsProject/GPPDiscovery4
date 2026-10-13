# Ritz reproduction from the seeded workbench — 2026-10-13

Run from this directory (`tools/ritz/`) after the workbench seeding commit:

```
python3 theta_control_ritz_4_spk4.py
basis 243
theta=0.0  lam_min-1 = -1.335e-07
theta=0.1  lam_min-1 = -5.191e-04
theta=0.2  lam_min-1 = -2.081e-03
theta=0.3  lam_min-1 = -4.688e-03
```

**Verdict: exact agreement** with the original 2026-10-08 baseline recorded
in `results-2026-10-08.md` — all four values match to every printed digit,
and the discriminator law lam_min − 1 ≈ −0.052·theta² is confirmed from the
repo's own copy with no workspace path dependence.

Status: NUMERICAL (reproduction of a banked numerical result).
