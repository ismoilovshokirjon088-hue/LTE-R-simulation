# Figure 4 simulation model (adapt for manuscript Section 2.4)

## Scenario
Equally spaced eNodeBs along a 1-D track; inter-site distance ISD = 2R - L_ov = 1700 m;
E-UTRA Band 31 downlink carrier f_c = 465 MHz; channel bandwidth B = 5 MHz.
Figure 4 fixes the train speed at **v = 200 km/h** (3GPP Rel-13 LTE-R reference speed for revenue
operation) and sweeps the offered load rho in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95].

## Radio and link layers (identical to Figure 3)
WINNER II D2a path loss PL(d) = 21.5 log10(d) + 44.2 + 20.0 log10(f_c/5 GHz),
log-normal shadowing (sigma = 3.73 dB, exponential spatial correlation with
decorrelation distance 50 m), Rician small-scale fading
(K = 5.12 dB, independent per transmission attempt). Receiver noise is
kT_0 B + NF with NF = 9.0 dB. The 1 Mbit/s service decodes if and only if
the chase-combined post-fading SINR exceeds the information-outage bound of
(2^(R/B) - 1) * 10^(3.0/10). Failed transmissions re-enter the scheduler
after T_attempt = 8 ms and participate in chase combining; at most
4 transmissions per packet.

## Baseline (C0) is a 3GPP-conformant configuration
The baseline uses the 3GPP-default A3 handover parameters: hysteresis **H = 3.0 dB**,
time-to-trigger **TTT = 160 ms**. It is scheduled first-come-first-served and
does not use a standby path. This is the configuration against which the three proposed
mechanisms are compared.

## Deadline-reliability metric (Eq. 13 of the manuscript)
P_del is the sample fraction of control messages whose end-to-end delay
(processing + queueing + transmission + HARQ retransmissions + handover interruption + standby-path
switching) does not exceed T_dl = 100 ms. At each (rho, configuration) the
fraction is computed over 100 independent replications of 120 s of train operation;
the reported 95 % confidence interval is the between-run Student-t interval. Lost messages
(exceeding 4 HARQ attempts) are counted against P_del. The analytic upper bound
(FIFO with M/M/1 waits, from manuscript Eq. 7a) is reported in Figure4_reliability_analytic.csv
and is compared to the simulated C0 at four rho values; the gap quantifies the contribution of
handover interruptions and HARQ retransmissions, neither of which appears in Eq. (7a).

## Heavy-tail mechanism (addresses reviewer R1 comment 2 and R2 comment 3)
The empirical delay CDF (Figure4_delay_CDF.png, raw data in Figure4_cdf_samples.csv) shows that
the delay distribution is NOT exponential. At high load the right tail is dominated by three
distinct mechanisms, all explicit in the simulator and absent from the mean-only Eq. (9):

  (i) **Handover interruption**. During the T_HO = 100 ms data-plane blackout
      that follows every successful A3 trigger, control messages continue to arrive and queue.
      At rho = 0.7 the serving link carries 10 control + 331 background messages per second and
      the queue builds at the full arrival rate during the blackout; the first message to
      transmit after the blackout has waited approximately 100 ms before its own
      service time is added. For C0 (FIFO + single link) a control message that arrived at the
      start of the blackout therefore accumulates roughly 100 ms of waiting time even at
      moderate load; this is the tail of its delay distribution. C2 breaks this chain by
      sending the message over the standby path (switching delay ~10 ms).

 (ii) **HARQ retransmission chain**. A block in deep fading can require up to 4
      transmissions separated by T_attempt = 8 ms, so the retransmission
      contribution alone is bounded by 32 ms before the
      packet is declared lost. For messages with an initially low SINR this is the second main
      source of the right-tail mass.

(iii) **Radio-link failure and re-establishment**. When the serving link stays below the
      out-of-sync threshold for T310 = 1.0 s, an RLF triggers a
      500 ms re-establishment interval; any message queued during the interval
      inherits its full duration as waiting time. This is rare at the simulated rho and shows
      up as the very long end of the C0 tail in Figure4_delay_CDF.png.

The three mechanisms also explain why P_del collapses at high rho while mean delay only grows
to a few tens of ms: a single 100 ms interruption during a handover converts a small fraction of
messages into deadline misses, so small changes in the right tail move P_del much more than they
move the mean.

## Threshold-crossing analysis
For each configuration the script reports the largest rho at which the simulated P_del remains
>= 99 % and >= 99.9 %, using linear interpolation between adjacent rho points. These crossings
are reported in Figure4_threshold_crossings.csv and quoted in the manuscript's Section 3.2.

## Relationship to FRMCS 99.9999 % target (addresses reviewer R3)
The 99.9 % threshold shown on the figure is a generic SIL-3-style benchmark, not the FRMCS
99.9999 % target. The proposed configuration C3 does NOT meet 99.9999 % at any simulated load;
the manuscript's Section 3.6 and Conclusions state this explicitly and identify safety-case
measurement-based reliability growth analysis as future work.

## Engine validation and reproducibility
The engine passed 7/7
self-tests against the analytical M/M/1, M/D/1, accumulating-priority (Eq. 15a) and HARQ
truncated-geometric references at rho = 0.70 (deviation < 3 %). The tests run on every
script invocation and the script aborts on failure. Common random numbers are used across
configurations within a replication, so paired comparisons (Figure4_paired_stats.csv) have
very low variance. The engine SHA-256 and all package versions are recorded in
Figure4_manifest.json.
