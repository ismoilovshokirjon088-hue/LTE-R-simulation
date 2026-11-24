# Figure 5 simulation model (adapt for manuscript Section 2.4 / 2.7 / 3.3)

## Scenario
Equally spaced eNodeBs along a 1-D track (ISD = 1700 m),
E-UTRA Band 31 downlink carrier f_c = 465 MHz, channel bandwidth B = 5 MHz.
Figure 5 fixes **rho = 0.70** (same as Figure 3) and sweeps the train speed over
v in [160, 200, 250, 300] km/h. The radio and link layers are identical to Figure 3 and Figure 4;
see those model summaries for the WINNER II D2a path loss, log-normal shadowing, Rician
fading, information-outage decoding threshold and chase-combining HARQ.

## Definition of the outage probability (manuscript Eq. 14)
The outage probability P_out is calculated as the proportion of non-overlapping
observation intervals of duration T_obs = 100 ms, taken from the
post-warm-up portion of each replication, that contain at least one communication
interruption. At every simulated millisecond the engine knows whether the primary and
standby paths are both ineligible; a 100 ms window is classified as an outage window
if any of its milliseconds has no eligible path. Multiple interruption events within
the same observation interval are counted as one outage window (as stated in the
manuscript's Eq. 14 definition).

Across the replication of 120 s after a 10 s warm-up, the
number of observation windows is N_obs = 1099.
P_out is reported as the sample mean of N_out / N_obs over the 100 independent
replications, with 95 % CI from the between-run Student-t distribution.

## Per-cause attribution of the outage intervals (addresses reviewer R2 comment 9)
Each outage window is attributed to its FIRST bad-millisecond cause code, so the
breakdown is partitioned (no double-counting). The five engine cause codes are:

  1 = handover execution blackout  (nominal 100 ms data-plane interruption
                                     following a successful A3 trigger)
  2 = handover command failure      (initial access to the target cell failed; the
                                     scheduler re-establishes after 500 ms)
  3 = handover target failure       (target becomes unusable after hand-over; re-establish)
  4 = radio-link failure            (the serving SINR has stayed below gamma_th for T310 =
                                     1.0 s, triggering a 500 ms
                                     re-establishment)
  5 = low SINR                      (both primary and standby below the decoding threshold
                                     because of deep fading or interference; this is a
                                     RADIO-layer effect only, not an L2 congestion drop
                                     while no handover or RLF is in progress)

In the manuscript narrative the five codes are grouped into three reviewer-named
categories (reviewer R2 comment 9):
  - **Handover-related outage**: codes 1, 2, 3.
  - **Radio-link failure + recovery**: code 4.
  - **Low-SINR (fading / interference)**: code 5.
The per-category contribution to P_out at each speed is reported in
Figure5_cause_breakdown.csv and plotted in Figure5_cause_attribution.png. This
breakdown resolves the earlier reviewer concern that the reported 10 % and 4 %
outage probabilities could not be reconciled with the nominal 100 ms
hand-over execution time alone (manuscript Section 3.3 revision).

## Derivation of the speed-adaptive hysteresis and time-to-trigger (R2 comment 8)
The baseline (C0) uses the 3GPP-default A3 parameters H = 3.0 dB,
TTT = 160 ms. The proposed speed-adaptive rule preserves the
trigger-decision distance across the operating-speed range: with TTT scaled as
TTT(v) = 25600 / v [ms, v in km/h], a train travels the same
7.11 m during TTT at every speed (160 km/h / 3.6 * 0.160 s = 7.11 m, 300 km/h / 3.6 *
0.085 s = 7.11 m). The hysteresis uses the same proportionality,
H(v) = 480 / v [dB], so that H and TTT reduce to the baseline values
at v = 160 km/h and become progressively more aggressive at higher speeds, in line with
the standard mobility-adaptive-handover arguments of the cited references. The values
at the four simulated speeds are listed in Figure5_adaptive_parameter_table.csv.

The timing-feasibility condition T_HO + T_margin <= T_ov (Eq. 17) is satisfied at every
simulated speed: the overlap crossing time T_ov = L_ov / v ranges from 3.60
s at 300 km/h to 6.75 s at 160 km/h, both
substantially above T_HO + T_margin = 150 ms. Eq. 17 is
therefore nominal timing feasibility, and the simulated outage probability is governed
by the stochastic HARQ, SINR and hand-over outcomes of the simulator, not by Eq. 17 alone.

## Common random numbers and run-to-run variation
Within each replication the radio environment, traffic realisation and fading draws are
shared across all four configurations, so the paired differences in Figure5_paired_stats.csv
have very low variance. The reported CIs combine 100 replications per point.

## Why C0 is identical to C1 at the radio layer (sanity check)
Priority scheduling (the only mechanism that differentiates C0 from C1) re-orders the
scheduler queue. It has no effect on which path is eligible at a given millisecond, and
the outage probability in Eq. (14) is defined over the eligibility mask, not the queue
state. The simulated P_out for C0 and C1 is therefore numerically identical at every
speed; Figure5_ablation_outage.csv reports both and includes an explicit identity
check. The main Figure 5 plots only C0 and C3 for clarity.

## Why radio-link failure (RLF, code 4) contributes near-zero outage
The T310 timer only fires after 1.0 s of continuous out-of-sync, which
requires sustained low SINR. The simulated channel (WINNER II D2a path loss with
sigma = 3.73 dB log-normal shadowing and Rician K = 5.1
dB fading) rarely produces sustained deep-fade intervals at the simulated rho, so RLF is
an exceptional event and its contribution to Figure5_cause_breakdown.csv is at most
~0 % at every speed. At larger shadow sigma, higher rho, or more aggressive traffic
bursts, this cause would become visible; Figure 3's sensitivity sweep on shadow sigma
illustrates this.

## Why Figure5_interruption_durations.csv uses fixed engine constants
The interruption durations reported in Figure5_interruption_durations.csv are derived
from the engine's handover and radio-link events using the configured T_HO,
T_reest and T_HO + T_reest values. These are ASSUMED procedure durations (Table 4),
not empirical measurements of each individual blackout interval. The reviewer R3
clarification in Section 3.3 should state this explicitly. A measurement-based
interruption duration study would require per-millisecond traces, which are not
retained by the current script to keep the archived outputs compact.

## C3 vs C2: honest trade-off at high speed
The adaptive handover rule triggers more handovers at higher speeds (at 300 km/h TTT
reduces from 160 ms to 85 ms, so a handover is triggered after a smaller RSRP margin).
The nominal 100 ms blackout therefore occurs more often for C3 than for C2, and the
standby path does not absorb every blackout millisecond (particularly when the standby
path is itself unavailable, which has about a 1 % duty cycle in the ASSUMED backup
availability model). The net effect visible in Figure5_ablation_outage.csv is that C3
has marginally HIGHER P_out than C2 at 300 km/h, by about 0.02 percentage points.
Figure 3 and Figure 4 show that the adaptive rule still improves end-to-end delay and
deadline reliability at 300 km/h; the outage-exposure cost quantified here is a real
engineering trade-off and should be reported in Section 3.3 narrative rather than
obscured.

## Analytic cross-check
Figure5_P_out_analytic.csv compares the simulated handover-related P_out to a nominal
budget estimate (total blackout seconds divided by total measurement seconds). The
simulated value is equal to or larger than the budget because 100 ms observation
windows typically overlap two blackouts when misaligned; the gap quantifies this
window-alignment overhead.

## Engine validation
The engine passed 7/7
self-tests against analytical M/M/1, M/D/1, accumulating-priority (Eq. 15a) and HARQ
truncated-geometric references at rho = 0.70; the tests run on every script invocation
and the script aborts on failure.
