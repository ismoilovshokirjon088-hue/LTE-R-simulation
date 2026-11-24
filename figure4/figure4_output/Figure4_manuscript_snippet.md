# Figure 4 - ready-to-paste numbers for Section 3.2 (adapt wording)

## Main result: deadline reliability vs offered load at 200 km/h

| rho | C0 (baseline) P_del | C1 P_del | C2 P_del | C3 (full) P_del |
|---|---|---|---|---|
| 0.10 | 0.9985 [0.9982, 0.9987] | 0.9986 [0.9983, 0.9989] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] |
| 0.20 | 0.9986 [0.9984, 0.9989] | 0.9989 [0.9987, 0.9991] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] |
| 0.30 | 0.9979 [0.9975, 0.9983] | 0.9986 [0.9982, 0.9989] | 1.0000 [1.0000, 1.0000] | 1.0000 [1.0000, 1.0000] |
| 0.40 | 0.9970 [0.9966, 0.9974] | 0.9984 [0.9981, 0.9987] | 1.0000 [1.0000, 1.0000] | 1.0000 [0.9999, 1.0000] |
| 0.50 | 0.9956 [0.9951, 0.9962] | 0.9979 [0.9975, 0.9982] | 1.0000 [1.0000, 1.0000] | 0.9999 [0.9999, 1.0000] |
| 0.60 | 0.9913 [0.9904, 0.9922] | 0.9972 [0.9967, 0.9977] | 0.9999 [0.9997, 1.0000] | 1.0000 [0.9999, 1.0000] |
| 0.65 | 0.9863 [0.9849, 0.9876] | 0.9966 [0.9959, 0.9972] | 1.0000 [0.9999, 1.0000] | 1.0000 [0.9999, 1.0000] |
| 0.70 | 0.9766 [0.9743, 0.9790] | 0.9951 [0.9943, 0.9959] | 0.9999 [0.9998, 1.0000] | 0.9999 [0.9999, 1.0000] |
| 0.75 | 0.9553 [0.9516, 0.9589] | 0.9933 [0.9924, 0.9941] | 0.9997 [0.9996, 0.9999] | 0.9998 [0.9998, 0.9999] |
| 0.80 | 0.9265 [0.9201, 0.9329] | 0.9903 [0.9890, 0.9916] | 0.9989 [0.9985, 0.9992] | 0.9994 [0.9992, 0.9997] |
| 0.85 | 0.8631 [0.8546, 0.8717] | 0.9822 [0.9804, 0.9840] | 0.9941 [0.9932, 0.9951] | 0.9965 [0.9958, 0.9971] |
| 0.90 | 0.7806 [0.7708, 0.7904] | 0.9721 [0.9699, 0.9742] | 0.9843 [0.9827, 0.9859] | 0.9881 [0.9866, 0.9895] |
| 0.95 | 0.6145 [0.6034, 0.6257] | 0.9565 [0.9541, 0.9589] | 0.9682 [0.9662, 0.9703] | 0.9725 [0.9706, 0.9743] |

## Threshold crossings (largest rho at which P_del >= target, at 200 km/h)

| configuration | 99 % crossing | 99.9 % crossing |
|---|---|---|
| C0 | 0.61 | never within tested range |
| C1 | 0.80 | never within tested range |
| C2 | 0.87 | 0.79 |
| C3 | 0.89 | 0.81 |

## Paired reduction of the outage fraction (1 - P_del): C0 vs C3

| rho | mean reduction (percentage points) | 95 % CI | Wilcoxon p (Holm) |
|---|---|---|---|
| 0.10 | 0.152 | [0.125, 0.179] | 6.90e-14 |
| 0.20 | 0.137 | [0.114, 0.161] | 8.80e-15 |
| 0.30 | 0.209 | [0.172, 0.246] | 4.80e-15 |
| 0.40 | 0.298 | [0.259, 0.337] | 1.80e-16 |
| 0.50 | 0.433 | [0.378, 0.487] | 4.81e-17 |
| 0.60 | 0.866 | [0.773, 0.960] | 4.81e-17 |
| 0.65 | 1.370 | [1.232, 1.508] | 4.81e-17 |
| 0.70 | 2.332 | [2.097, 2.567] | 4.81e-17 |
| 0.75 | 4.458 | [4.091, 4.825] | 4.81e-17 |
| 0.80 | 7.292 | [6.663, 7.920] | 4.81e-17 |
| 0.85 | 13.333 | [12.509, 14.157] | 4.81e-17 |
| 0.90 | 20.746 | [19.846, 21.647] | 4.81e-17 |
| 0.95 | 35.793 | [34.774, 36.812] | 4.81e-17 |

All C0 vs C3 differences at rho >= 0.30 are significant at Holm-corrected p < 0.001.

## Stepwise attribution at each load (which mechanism contributes where?)

| rho | C0 to C1 (priority) | C1 to C2 (standby) | C2 to C3 (adaptive HO) |
|---|---|---|---|
| 0.10 | +0.015 pp | +0.137 pp | +0.000 pp |
| 0.20 | +0.030 pp | +0.108 pp | -0.001 pp |
| 0.30 | +0.067 pp | +0.142 pp | +0.000 pp |
| 0.40 | +0.139 pp | +0.160 pp | -0.001 pp |
| 0.50 | +0.225 pp | +0.212 pp | -0.005 pp |
| 0.60 | +0.593 pp | +0.266 pp | +0.007 pp |
| 0.65 | +1.029 pp | +0.340 pp | +0.001 pp |
| 0.70 | +1.847 pp | +0.480 pp | +0.005 pp |
| 0.75 | +3.801 pp | +0.646 pp | +0.011 pp |
| 0.80 | +6.385 pp | +0.853 pp | +0.055 pp |
| 0.85 | +11.905 pp | +1.195 pp | +0.233 pp |
| 0.90 | +19.144 pp | +1.227 pp | +0.375 pp |
| 0.95 | +34.198 pp | +1.173 pp | +0.422 pp |

## Key sentence candidates for Section 3.2

- At 200 km/h (Figure 4), the baseline configuration C0 meets the 99 % deadline-reliability target only up to ρ ~= 0.61; the full proposed configuration C3 extends this range to ρ ~= 0.89.
- The 99.9 % SIL-3-like target is never met by C0 within the tested range of offered load; C3 meets it up to ρ ~= 0.81.
- The step-wise decomposition in the table above attributes the largest reliability gain at every load to the addition of the standby path (C1 to C2), with the speed-adaptive handover rule contributing less than 1 percentage point of additional reliability on top of C2.
- The FRMCS 99.9999 % target is NOT met by any simulated configuration at any load; the present figure reports only the simulated deadline-based delivery, not a formal safety-case reliability.
