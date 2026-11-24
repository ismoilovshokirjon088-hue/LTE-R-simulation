# Figure 4 - internal review-point mapping (NOT a reviewer-response letter)

> This file is a DEV MAPPING: it links every reviewer comment that concerns Figure 4 to
> the exact artifacts and numbers this script produces, so that the actual response letter
> can be drafted later (once Figures 3-5 are finalised and the manuscript body is revised).

This document maps each reviewer comment that concerns Figure 4 to the specific files, numbers and sections of the revised manuscript that address it.

## Reviewer R1, comment 2

> "Reliability results need stronger justification. Average-delay equations cannot reproduce the deadline-based reliability curves in Figure 4. Define the full delay distribution and how the 100 ms deadline probability is calculated."

**Addressed by**:

- Figure4_model_summary.md -> section 'Deadline-reliability metric' (definition of P_del as the sample fraction over 100 independent replications, with the analytic upper bound from Eq. 7a reported in Figure4_reliability_analytic.csv).
- Figure4_model_summary.md -> section 'Heavy-tail mechanism' (three distinct simulator mechanisms responsible for the right tail of the delay distribution, each named and quantified: handover interruption, HARQ retransmission chain, radio-link failure).
- Figure4_delay_CDF.png and Figure4_cdf_samples.csv (empirical delay CDF at rho = 0.70 and rho = 0.90, from which the 100 ms deadline probability can be read directly).

## Reviewer R2, comment 3 (the most critical comment for Figure 4)

> "Figure 4 cannot be reconciled with the rest of the model as presented. ... Please show the delay distribution (a CDF or percentile table at two or three utilization values) and explain which mechanism in the simulator produces the heavy tail."

**Addressed by**:

- The earlier (MATLAB-based) Figure 4 was based on an inverse calibration of the mean retransmission delay through Eq. 9 and Pe (see R2 comment 1). The revised Section 2.4 retires that calibration; the new engine derives Figure 4 from a physically consistent simulator and therefore no longer contradicts the mean-delay figures.
- Figure4_delay_CDF.png shows the empirical delay CDF at rho = 0.70 and rho = 0.90 for every configuration (reviewer-requested 'two or three utilization values').
- Figure4_model_summary.md -> section 'Heavy-tail mechanism' names the three mechanisms that produce the right tail, and explains why P_del can collapse at high load while the mean delay stays at a few tens of milliseconds.
- The present simulated numbers (new model) are: at rho = 0.70, 200 km/h, baseline C0 mean delay = 20.1 ms and P_del = 0.977; at rho = 0.95, proposed C3 mean delay = 34.4 ms and P_del = 0.972. These replace the earlier numbers that could not be reconciled with each other.

## Reviewer R2, comment 5

> "Table 4 states 100 independent runs with 95 percent confidence intervals for the proposed configuration, but no interval appears in any figure or table, and the number of runs for the baseline is not given. Please add error bars to Figures 3 and 5 and confidence bands to Figure 4 for both configurations, and state the random seed handling. If the run-to-run variation is negligible, say so and show the numbers."

**Addressed by**:

- The main Figure 4 (Figure4_reliability_vs_load.png) shows 95 % CI bands as shaded areas around every curve, for all four configurations (C0 through C3). The baseline C0 uses the same 100 independent replications as the proposed configurations.
- Figure4_summary.csv lists the per-(rho, config) mean, 95 % CI, and the between-run standard deviation.
- Figure4_manifest.json records the base random seed, the per-replication seed derivation (numpy SeedSequence.spawn), and the common-random-numbers scheme used across configurations within a replication.

## Reviewer R2, comment 10 (practical recommendations)

> "the recommendation 'rho ~= 0.5-0.6' is taken from Figure 4 at 200 km/h only, so please clarify as 'at 200 km/h (Figure 4)'."

**Addressed by**:

- Figure4_manuscript_snippet.md -> 'Key sentence candidates' -- every sentence that quotes a rho crossing now reads 'At 200 km/h (Figure 4), ...' so the scope limitation is explicit. The exact crossings from the full simulation are: C0 at 0.61; C1 at 0.80; C2 at 0.87; C3 at 0.89.

## Reviewer R3 (FRMCS 99.9999 % target)

> "The reliability results in Figure 4 do not demonstrate compliance with the 99.9999% reliability target considered for FRMCS."

**Addressed by**:

- Figure4_model_summary.md -> section 'Relationship to FRMCS 99.9999 % target'. The 99.9 % threshold shown on the figure is explicitly identified as a generic SIL-3 benchmark, not the FRMCS 99.9999 % target. The proposed configuration C3 does not meet the FRMCS target at any simulated load and the manuscript Section 3.6 and Conclusions state this explicitly. Field measurement-based reliability growth analysis remains future work (manuscript Section 3.7).

## General point on the simulation model (R1 comment 1, R2 comment 1)

> "Provide pseudocode or source code for stochastic channel generation, packet transmission, and random simulation procedures." / "Pe is back-computed from the simulation; please reword the contribution."

**Addressed by**:

- The full Python source code of the simulator is included inside this standalone file (section '===== ENGINE ====='), and the entire run is reproducible from Figure4_manifest.json (seed, parameters with source tags, package versions, SHA-256 of this script).
- The new engine does NOT use the inverse Pe calibration. Block errors are produced by an information-outage link abstraction with chase-combining HARQ, driven by the simulated SINR; the analytical Eq. 7a appears only as a cross-check lower bound for the baseline reliability (Figure4_reliability_analytic.csv), not as a calibration target.
