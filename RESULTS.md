# Results — GPPDiscovery4

Headline results of the RH fixed-window / sewn-sector programme, at full precision,
each with its status label and provenance. The machine-readable companion is
[`results.jsonl`](results.jsonl) (one JSON object per result; headline values are
strings, to preserve precision). Status labels: **CERTIFIED** = sign certified
against a measured error floor; **NUMERICAL** = high-precision numerical value,
not interval-certified; **THEOREM** = proved (constants at the digits shown);
**EXACT** = exact rational arithmetic, verified.

## 1. Sector certification ladder

Method: the S1 sewn sector form $M_Q$ is assembled from segmented prime sums
(compensated accumulation) plus the archimedean term; the sign of
$\lambda_{\min}(M_Q)$ is certified against the measured Platt zero-data floor
(`err_platt`); the certification ratio is $\mathrm{CR} = |\lambda_{\min}|/\mathrm{err}_{\mathrm{Platt}}$.
Run pairs are bit-identical, and the K ≤ 17 block reproduces the previously banked
sectors with difference exactly 0.0. $\zeta$ is feasible at every $K = 2..18$;
the shifted-product control $F_\theta$, $\theta = 0.25$, is certified infeasible at every $K$.

### ζ

| K | λ_min(M_Q) | λ_min(P⁺, even) | λ_min(P⁻, odd) | Platt floor | CR | status |
|---|-----------|-----------------|----------------|-------------|-----|--------|
| 2 | 3.135007620997027e-05 | 3.203545149599611e-05 | 3.135007620997017e-05 | 1.980557712232173e-14 | 1582891324.8196843 | CERTIFIED |
| 3 | 5.523398602359483e-07 | 1.1209396861496043e-06 | 5.523398602359619e-07 | 2.4520327745994535e-14 | 22525794.3514305 | CERTIFIED |
| 4 | 3.5541713498842336e-07 | 9.043850660408256e-07 | 3.5541713498821965e-07 | 2.851715659766975e-14 | 12463273.951283975 | CERTIFIED |
| 5 | 2.2850005342193665e-07 | 3.586917040792169e-07 | 2.2850005342192458e-07 | 3.251193267235668e-14 | 7028190.4101080755 | CERTIFIED |
| 6 | 2.2153536255041886e-07 | 2.2907971959293336e-07 | 2.2153536255052924e-07 | 3.966059905194464e-14 | 5585779.535509979 | CERTIFIED |
| 7 | 1.8551739132738434e-07 | 1.8551739132738993e-07 | 2.036388345827925e-07 | 6.72471124259206e-14 | 2758741.373940037 | CERTIFIED |
| 8 | 1.5482395763442788e-07 | 1.548239576344435e-07 | 1.5678730181630875e-07 | 1.5863837849674722e-13 | 975955.2455183626 | CERTIFIED |
| 9 | 9.827232568840367e-08 | 1.5061739697215995e-07 | 9.827232568844773e-08 | 2.501388423753918e-13 | 392871.1141187864 | CERTIFIED |
| 10 | 5.2699083975685126e-08 | 1.419255388754848e-07 | 5.269908397563126e-08 | 4.074688184221134e-13 | 129332.8019055731 | CERTIFIED |
| 11 | 2.6631037414073695e-08 | 9.291169033918171e-08 | 2.6631037414087254e-08 | 5.96695558524409e-13 | 44630.86247856545 | CERTIFIED |
| 12 | 2.466516010500742e-08 | 2.8329188921942932e-08 | 2.4665160104963675e-08 | 2.3159171059254033e-12 | 10650.27761222551 | CERTIFIED |
| 13 | 2.0958250498801098e-08 | 2.095825049890295e-08 | 2.4638698499626916e-08 | 1.035775866891893e-11 | 2023.4349118107589 | CERTIFIED |
| 14 | 2.0466410833361367e-08 | 2.0466410833390328e-08 | 2.340356779856925e-08 | 1.3958016602617028e-11 | 1466.2836000298253 | CERTIFIED |
| 15 | 2.005820979948799e-08 | 2.0058209799392096e-08 | 2.0848713925475105e-08 | 4.77640072544222e-11 | 419.9440321798149 | CERTIFIED |
| 16 | 1.9023374767511586e-08 | 1.9023374767472504e-08 | 2.0738793548746713e-08 | 6.670184799400268e-11 | 285.2001157332556 | CERTIFIED |
| 17 | 1.3967321533728195e-08 | 1.8341235636725215e-08 | 1.396732153347196e-08 | 1.3298654749797777e-10 | 105.02807837718011 | CERTIFIED |
| 18 | 1.2566674164415027e-08 | 1.2566674164520848e-08 | 1.368694925297216e-08 | 1.7101111315700061e-10 | 73.48454689536995 | CERTIFIED |

### F_θ control (θ = 0.25)

| K | λ_min(M_Q) | Platt floor | CR | status |
|---|-----------|-------------|-----|--------|
| 2 | -6.943990335619734e-05 | 3.84614113918569e-14 | 1805443452.1063688 | CERTIFIED (infeasible) |
| 3 | -0.0005309193537663755 | 4.707071861914666e-14 | 11279185220.478378 | CERTIFIED (infeasible) |
| 4 | -0.0007732913189579408 | 5.547716054820208e-14 | 13938913082.728094 | CERTIFIED (infeasible) |
| 5 | -0.001723683790313082 | 6.534282488109601e-14 | 26379082836.557194 | CERTIFIED (infeasible) |
| 6 | -0.003809083419344692 | 2.6938139427440924e-13 | 14140113238.349766 | CERTIFIED (infeasible) |
| 7 | -0.005577815614506711 | 1.1407276371161678e-12 | 4889699725.8765335 | CERTIFIED (infeasible) |
| 8 | -0.006101296776480582 | 4.425734323612625e-12 | 1378595353.9796383 | CERTIFIED (infeasible) |
| 9 | -0.01248607749304036 | 9.585617468507295e-12 | 1302584578.8299265 | CERTIFIED (infeasible) |
| 10 | -0.022445545324711474 | 1.5379121123308346e-11 | 1459481666.3933656 | CERTIFIED (infeasible) |
| 11 | -0.02615614170326902 | 7.107186345953288e-11 | 368023862.4693144 | CERTIFIED (infeasible) |
| 12 | -0.040485580664008876 | 2.57481402205912e-10 | 157236912.3251547 | CERTIFIED (infeasible) |
| 13 | -0.06711297748561397 | 2.082724171704174e-09 | 32223651.31082109 | CERTIFIED (infeasible) |
| 14 | -0.08684411270998282 | 2.964359841604424e-09 | 29296076.50567129 | CERTIFIED (infeasible) |
| 15 | -0.09631078035658319 | 1.051325793828184e-08 | 9160888.177763388 | CERTIFIED (infeasible) |
| 16 | -0.20100284027963258 | 3.193081810659347e-08 | 6294948.021958982 | CERTIFIED (infeasible) |
| 17 | -0.28946525122063593 | 4.9839411746439807e-08 | 5807958.823697661 | CERTIFIED (infeasible) |
| 18 | -0.3068065288414963 | 4.3953774800907724e-07 | 698020.8872416577 | CERTIFIED (infeasible) |

## 2. Step scalars of the sector chain

Method: bordering (Haynsworth) step scalars $s_K = c - b^{T} P_K^{-1} b$ for the
chain $P_K \to P_{K+1}$ on the lattice $a_j = 0.35 + 0.7j$; positivity of every
step scalar is equivalent to positive-definiteness of the whole chain. Each value
is cross-checked by four independent evaluations (Schur solve, determinant ratio,
inverse entry, variational) and by finer-grid repricing. Steps 2→6 use the original
builders (bit-exact reproduction gate); later steps use validated extended paths.
ζ's scalar rises to a peak at 5→6, falls to 13→14, and rises again at 14→15.

| step | s_K | status |
|------|-----|--------|
| 2→3 | 1.5057282291990905e-06 | NUMERICAL |
| 3→4 | 2.566885918385976e-06 | NUMERICAL |
| 4→5 | 3.815941469720329e-06 | NUMERICAL |
| 5→6 | 3.926810897066416e-06 | NUMERICAL |
| 6→7 | 2.309087889418814e-06 | NUMERICAL |
| 7→8 | 1.161529303714639e-06 | NUMERICAL |
| 8→9 | 7.144337548918136e-07 | NUMERICAL |
| 9→10 | 5.894734752083714e-07 | NUMERICAL |
| 10→11 | 4.6237030067197184e-07 | NUMERICAL |
| 11→12 | 3.7749869473175e-07 | NUMERICAL |
| 12→13 | 3.748054148115928e-07 | NUMERICAL |
| 13→14 | 2.578222329735285e-07 | NUMERICAL |
| 14→15 | 2.698572135717797e-07 | NUMERICAL |
| 15→16 | 2.0022361232997852e-07 | NUMERICAL |

## 3. Normalised Schur scalars (wave 234)

Method: the sector sewing written as a Schur complement in the $Q$-metric; $s_K$ is
the normalised margin of the sewing at sector size $K$.

| K | s_K | status |
|---|-----|--------|
| 2 | 0.3193423279306039 | NUMERICAL |
| 3 | 0.019060512784618844 | NUMERICAL |
| 4 | 0.016709982562574943 | NUMERICAL |
| 5 | 0.01261425681856232 | NUMERICAL |
| 6 | 0.007117093020035159 | NUMERICAL |
| 7 | 0.005473366483028763 | NUMERICAL |

## 4. The prime-free form Δ: uniform floors

Method: Δ is an explicit sum-of-squares (Toeplitz) integral with closed-form symbol
$\Psi$; its minimum is $\Psi(\pi) = 0.00886271$, giving the uniform floor
$\Delta \succeq 2\min\Psi \cdot I_K$ for every $K$ and both parities (THEOREM; the
symbol infimum is $f(\pi) = 0.0177254124818$). Independently, interval arithmetic
(mpmath, 6000 subintervals) certifies the weaker floor 0.0163 with enclosure bound 0.016371.

| statement | floor | status |
|-----------|-------|--------|
| Δ ⪰ c·I_K, all K, both parities (SOS symbol) | 0.0177254 | THEOREM |
| Δ_K ⪰ c·I_K, all K ≥ 1 (interval-certified) | 0.0163 | CERTIFIED |

## 5. F_θ controls: chain death steps

Method: per-parity Haynsworth step tables for $F_\theta$ ($\theta = 0.25$), rooted at
$K = 1$. The first negative step scalar in each parity is where the control's chain
dies; the even death scalar is identical to the first (ROW+) row-ladder failure,
reproduced by an independent route. (The sector-level infeasibility ladder is §1.)

| parity | death step | s | status |
|--------|-----------|---|--------|
| odd (P⁻) | 1→2 | -7.40047e-05 | NUMERICAL |
| even (P⁺) | 2→3 | -6.01482e-05 | NUMERICAL |

## 6. Exact factorisation at K = 6

Statement (EXACT, rational arithmetic): the $K = 6$ sector form factorises over $\mathbb{Q}$ as

$$\hat P_6 = C^{T} C, \qquad C \in \mathrm{Mat}_{24 \times 6}(\mathbb{Q}),$$

with $C^{T}C$ equal to the target matrix exactly and the LDL reconstruction exact;
each of the 6 pivots is a sum of 4 rational squares. Formalisation-ready form:
`∃ C : Matrix (Fin 24) (Fin 6) ℚ, Phat6 = Cᵀ * C`.
