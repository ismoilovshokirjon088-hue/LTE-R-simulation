# Figure 3 simulation model (ready to adapt for manuscript Section 2.4)

## Scenario
Equally spaced eNodeBs along a 1-D track, inter-site distance ISD = 2R - L_ov = 1700 m,
E-UTRA Band 31 downlink carrier f_c = 465 MHz, channel bandwidth B = 5 MHz,
train speeds v in {160, 200, 250, 300} km/h, offered load rho = 0.70.

## Radio layer
Path loss follows the WINNER II D2a open scenario [ref.]:
    PL(d) [dB] = 21.5 log10(d) + 44.2 + 20.0 log10(f_c / 5 GHz)
where d is the slant distance to the eNodeB. Shadow fading is log-normal with sigma = 3.73 dB
and exponential spatial correlation with decorrelation distance 50 m. Small-scale fading
is Rician with K = 5.12 dB (independent per transmission attempt). Receiver noise power follows
kT_0 B + NF with NF = 9.0 dB. Co-channel interference from the six neighbouring cells is
included, each weighted by the offered-load activity factor.

## Link abstraction
The 1 Mbit/s service decodes if and only if the chase-combined post-fading SINR exceeds
    gamma_dec = (2^(R/B) - 1) * 10^(Delta/10),  Delta = 3.0 dB
(the information-outage bound with a 3.0 dB implementation gap). Failed transmissions
re-enter the scheduler after T_attempt = 8 ms and participate in chase combining;
at most 4 transmissions per packet.

## Scheduler
C0 is FIFO. C1-C3 apply accumulating priority: at each scheduling instant the message with the
largest priority-weighted waiting time w_i (t - t_i) is dispatched, with w_ctrl : w_other = 4 : 1.

## Mobility
Fixed-H baseline: hysteresis H = 3.0 dB, time-to-trigger TTT = 160 ms.
Speed-adaptive (H flag): H(v) = 480/v dB, TTT(v) = 25600/v ms.
Handover timing feasibility: T_HO + T_margin <= T_ov (Eq. 17).
Handover execution occupies T_HO = 100 ms with no data-plane service.

## Standby path (B flag)
Second radio link toward the strongest non-serving cell, eligible when it is available
(two-state Markov availability A_b(t), long-run unavailability 0.010, mean down time
0.5 s) and its post-fading SINR >= gamma_th = -8.0 dB. Additional one-way delay
when a message uses the standby path: 10 ms.

## Measurement window and statistics
Simulated horizon = 120 s, warm-up = 10 s, drain = 5 s.
Independent replications per (speed, configuration) = 100, with common random numbers across
configurations of the same replication. All reported differences use paired bootstrap confidence
intervals, Wilcoxon signed-rank tests, and Holm correction for multiple testing.

## Engine validation
The engine passed 7/7
self-tests against the analytical M/M/1, M/D/1, accumulating-priority (Eq. 15a) and HARQ-truncated-geometric
references at rho = 0.70 (deviation < 3 % for every test).

## Parameter-source tagging
T4 = manuscript Table 4; REF = cited reference; ASSUMED = modelling assumption not in Table 4.
The full list with values is in Figure3_manifest.json. The sensitivity of the main conclusion
to the ASSUMED parameters is reported in Figure3_sensitivity.{csv,png}.
