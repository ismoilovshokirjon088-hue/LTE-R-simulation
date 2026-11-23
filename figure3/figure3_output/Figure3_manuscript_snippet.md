# Figure 3 - ready-to-paste numbers for Section 3.1 (adapt wording)

## Main result (mean end-to-end delay of control messages at rho = 0.70)

| v (km/h) | C0 (ms) | C1 (ms) | C2 (ms) | C3 (ms) | C0->C3 reduction |
|---|---|---|---|---|---|
| 160 | 18.58 [18.20, 18.95] | 14.94 | 13.81 | 13.81 [13.76, 13.85] | 25.7 % [24.4, 27.0] |
| 200 | 20.65 [19.98, 21.32] | 15.82 | 13.94 | 13.81 [13.77, 13.85] | 33.1 % [31.1, 35.2] |
| 250 | 21.07 [20.35, 21.78] | 16.01 | 13.99 | 13.81 [13.77, 13.84] | 34.5 % [32.4, 36.5] |
| 300 | 22.63 [21.79, 23.46] | 16.63 | 14.05 | 13.84 [13.81, 13.88] | 38.8 % [36.7, 40.9] |

All C0 vs C3 differences are significant (Wilcoxon p < 1e-15, Holm-corrected).

## Mechanism attribution (stepwise paired differences, ms)

| v (km/h) | C0->C1 priority | C1->C2 standby | C2->C3 adaptive HO |
|---|---|---|---|
| 160 | 3.64 [3.42, 3.86] | 1.13 [0.98, 1.28] | 0.00 [0.00, 0.00] |
| 200 | 4.82 [4.44, 5.21] | 1.88 [1.61, 2.15] | 0.13 [0.10, 0.16] |
| 250 | 5.05 [4.70, 5.41] | 2.03 [1.63, 2.42] | 0.18 [0.15, 0.21] |
| 300 | 6.00 [5.54, 6.47] | 2.58 [2.23, 2.93] | 0.20 [0.16, 0.24] |

Note: at v = 160 km/h the adaptive rule H(v) = 480/v = 3 dB and TTT(v) = 160 ms reduce to the baseline fixed-handover parameters, so C2 and C3 are identical by construction.

## Deadline reliability (P_del, fraction of control messages delivered within 100 ms)

| v (km/h) | C0 | C1 | C2 | C3 |
|---|---|---|---|---|
| 160 | 0.9821 | 0.9969 | 1.0000 | 1.0000 |
| 200 | 0.9742 | 0.9942 | 1.0000 | 0.9999 |
| 250 | 0.9723 | 0.9942 | 0.9999 | 0.9999 |
| 300 | 0.9652 | 0.9922 | 0.9999 | 1.0000 |

C0 does not reach 99 % deadline reliability above 200 km/h; C3 reaches at least 99.9 % at every speed.

## Handover cost (mean per 120-s run)

| v (km/h) | C0 handovers | C3 handovers | C3 ping-pongs | C3 cmd failures |
|---|---|---|---|---|
| 160 | 9.0 | 9.0 | 1.3 | 0.07 |
| 200 | 10.9 | 13.3 | 3.2 | 0.11 |
| 250 | 11.5 | 19.0 | 7.1 | 0.08 |
| 300 | 12.6 | 23.8 | 10.2 | 0.15 |

The speed-adaptive policy roughly doubles the handover count at 300 km/h and raises the ping-pong count from about 1 to about 12 per 120 s; this signalling cost is the trade-off behind the small additional delay reduction seen only above 200 km/h.
