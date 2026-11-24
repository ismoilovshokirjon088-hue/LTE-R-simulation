# Figure 5 - ready-to-paste numbers for Section 3.3 (adapt wording)

All numbers below are produced by 100 independent replications per (train speed, configuration) at rho = 0.70 and are reported with 95 % CIs over the between-run Student-t distribution. The main figure plots only the baseline C0 and the full proposed C3; the full ablation with the intermediate configurations C1 and C2 is given in Figure5_ablation_outage.csv.

## Main result: outage probability vs train speed (rho = 0.70)

| v (km/h) | C0 baseline P_out | C3 full P_out | C0 to C3 reduction |
|---|---|---|---|
| 160 | 1.357% [1.258, 1.457] | 0.014% [0.003, 0.024] | 99.0% |
| 200 | 1.655% [1.540, 1.771] | 0.025% [0.011, 0.039] | 98.5% |
| 250 | 1.904% [1.777, 2.030] | 0.020% [0.009, 0.031] | 98.9% |
| 300 | 2.340% [2.193, 2.487] | 0.041% [0.020, 0.062] | 98.3% |

All C0 vs C3 differences are significant at Holm-corrected p < 0.001.

## Full ablation (addresses R2 comment 10 ablation precision)

| v (km/h) | C0 baseline | C1 +priority | C2 +standby | C3 full | C3 vs C2 delta |
|---|---|---|---|---|---|
| 160 | 1.357% | 1.357% | 0.014% | 0.014% | +0.000 pp |
| 200 | 1.655% | 1.655% | 0.021% | 0.025% | +0.004 pp |
| 250 | 1.904% | 1.904% | 0.021% | 0.020% | -0.001 pp |
| 300 | 2.340% | 2.340% | 0.022% | 0.041% | +0.019 pp |

Two structural observations should be reported honestly in the manuscript narrative:

  1. **C0 = C1 by construction**. Priority scheduling only re-orders the queue and does not affect the radio-layer availability sampled by Eq. (14), so the simulated P_out for C0 and C1 is numerically identical at every speed. This is a sanity check that the engine distinguishes correctly between scheduling and radio outage.
  2. **C3 vs C2 trade-off at high speed**. The speed-adaptive handover rule triggers more handovers at 300 km/h (approximately double the baseline count), each contributing a nominal 100 ms blackout. The simulated handover-related outage of C3 therefore slightly EXCEEDS that of C2 at 300 km/h (delta = +0.019 pp). Because the standby path still absorbs almost all of this blackout, the net C3 P_out remains more than 50x below C0; but the manuscript narrative should NOT claim that adaptive handover strictly improves outage. The adaptive rule's benefit is in Figure 3 (delay) and Figure 4 (deadline reliability), with a small outage-exposure cost quantified here.

## Per-cause attribution (addresses R2 comment 9)

For C0 and C3 only, the fraction of P_out at each speed attributable to each cause category. Categories partition P_out (each outage window attributed to its first bad-millisecond cause, no double-counting).

| v (km/h) | config | handover-related | radio-link failure | low-SINR (fading) | total P_out |
|---|---|---|---|---|---|
| 160 | C0 | 1.267% | 0.000% | 0.090% | 1.357% |
| 160 | C3 | 0.014% | 0.000% | 0.000% | 0.014% |
| 200 | C0 | 1.474% | 0.000% | 0.182% | 1.655% |
| 200 | C3 | 0.022% | 0.000% | 0.003% | 0.025% |
| 250 | C0 | 1.578% | 0.000% | 0.325% | 1.904% |
| 250 | C3 | 0.019% | 0.000% | 0.001% | 0.020% |
| 300 | C0 | 1.846% | 0.000% | 0.494% | 2.340% |
| 300 | C3 | 0.041% | 0.000% | 0.000% | 0.041% |

## Key sentence candidates for Section 3.3

- Figure 5 reports the simulated communication outage probability at rho = 0.70 for the baseline configuration (C0) and the full proposed configuration (C3) across the four train speeds {160, 200, 250, 300} km/h, with 95 % confidence intervals obtained from 100 independent replications per point.
- At 300 km/h the baseline configuration C0 shows P_out = 2.340%, whereas the full proposed configuration C3 reduces this to 0.041%, a relative reduction of 98.3%.
- The per-cause attribution in Figure 5 (supplementary) partitions P_out into handover-related, radio-link failure and low-SINR (fading / interference) contributions. The handover-related category (nominal 100 ms blackout plus handover command and target failures) dominates at every simulated speed; radio-link failure (T310 timer triggered after 1 s continuous out-of-sync) contributes negligibly, which is consistent with the stable simulated channel used here.
- The speed-adaptive hysteresis and time-to-trigger of configuration C3 preserve the trigger-decision travel distance at the baseline speed of 160 km/h; the per-speed values and derivation are given in Figure5_adaptive_parameter_table.csv (manuscript Section 2.7 revision addressing reviewer R2 comment 8).
- The outage reductions reported in Figure 5 come from a stochastic simulator (not a closed-form formula); the 100 ms hand-over execution time in Table 4 is an assumed procedure duration, not a measured interruption duration (reviewer R3 clarification).
- The configuration C1 (priority scheduling only) is not plotted because its outage is numerically identical to the baseline C0 by construction; the full ablation (including C1 and C2) is reported in Figure5_ablation_outage.csv.
