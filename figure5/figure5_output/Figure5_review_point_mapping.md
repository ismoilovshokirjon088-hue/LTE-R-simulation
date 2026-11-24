# Figure 5 - internal review-point mapping (NOT a reviewer-response letter)

> This file is a DEV MAPPING: it links every reviewer comment that concerns Figure 5 to
> the exact artifacts and numbers this script produces, so that the actual response letter
> can be drafted later (once Figures 3-5 are finalised and the manuscript body is revised).

## Reviewer R2, comment 5 (confidence bands / seed handling)

> "Please add error bars to Figures 3 and 5 and confidence bands to Figure 4 for both configurations, and state the random seed handling."

**Addressed by**:

- Figure5_outage_vs_speed.png shows 95 % CI error bars on every bar for all four
  configurations (C0 through C3). The baseline C0 uses the same 100 independent
  replications as the proposed configurations.
- Figure5_summary.csv lists the per-(speed, config) mean, 95 % CI and between-run
  standard deviation of P_out.
- Figure5_manifest.json records the base random seed, the per-replication seed
  derivation (numpy SeedSequence.spawn) and the common-random-numbers scheme across
  configurations within a replication.

## Reviewer R2, comment 8 (speed-adaptive parameter derivation)

> "The speed-dependent handover parameters H(v) = 480/v dB and tau(v) = 25600/v ms (Table 4) are introduced without justification... a short derivation or a sensitivity sweep over the constants is needed."

**Addressed by**:

- Figure5_model_summary.md -> section 'Derivation of the speed-adaptive hysteresis
  and time-to-trigger'. The derivation preserves the trigger-decision travel distance
  (7.11 m, matching the baseline 160 ms TTT at 160 km/h) and is physically consistent
  with standard mobility-adaptive-handover arguments.
- Figure5_adaptive_parameter_table.csv tabulates H(v) and TTT(v) at each simulated
  speed and the corresponding trigger-decision distance, making the inverse-proportional
  form transparent.
- The sensitivity of the Figure 3 and Figure 4 conclusions to the six ASSUMED engine
  parameters (which include the hand-over timing budget) is reported in Figure 3's
  sensitivity sweep (same engine, same parameters).

## Reviewer R2, comment 9 (per-cause attribution of outage; the most critical Fig-5 comment)

> "The outage probabilities of 10 percent (baseline) and 4 percent (proposed) at 300 km/h are far larger than nominal handover execution can explain... the counts by cause are available and should be reported."

**Addressed by**:

- Figure5_cause_attribution.png (supplementary figure) and
  Figure5_cause_breakdown.csv (raw numbers) partition P_out at every (speed, config)
  into the three reviewer-named categories: handover-related, radio-link failure +
  recovery, low-SINR. The engine assigns each outage window its first bad-millisecond
  cause code, so the breakdown is partitioned without double-counting.
- Figure5_model_summary.md -> section 'Per-cause attribution of the outage intervals'
  names the five engine cause codes and the three reviewer categories, with the
  physical interpretation of each.
- At 300 km/h the attribution is (for the simulated model):
  - C0: handover-related 1.846%, radio-link failure 0.000%, low-SINR (fading) 0.494% -> total 2.340%.
  - C3: handover-related 0.041%, radio-link failure 0.000%, low-SINR (fading) 0.000% -> total 0.041%.
- The manuscript narrative about 'reducing interruptions during cell transitions'
  (Section 3.3) can now be either quantitatively supported by the handover-related
  reduction from C0 to C3, or equivalently reworded to refer to the measured total
  outage reduction; either option is now backed by the simulator.

## Reviewer R3 (handover success rate and interruption duration are distinct KPIs)

> "Handover success probability and communication interruption duration are distinct performance indicators that cannot be inferred directly from the outage probability."

**Addressed by**:

- Figure5_handover_kpis.csv reports the per-(speed, config) total handover count,
  successful handover count, command-failure count, target-failure count, ping-pong
  count and derived handover success rate - all as distinct KPIs independent of P_out.
- Figure5_interruption_durations.csv reports the mean, P95 and maximum interruption
  duration derived from the handover and radio-link events, separately from the
  outage-window count.

## Reviewer R1, comment 5 (simulation transparency)

> "Provide pseudocode or source code for stochastic channel generation, packet transmission and random simulation procedures."

**Addressed by**:

- The full Python source code of the simulator is inlined inside this standalone
  file (section '====== ENGINE ======'), and the entire run is reproducible from
  Figure5_manifest.json (seed, parameters with source tags, package versions, SHA-256
  of this script).
