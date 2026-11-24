# Figure 5 — Communication outage probability vs train speed

Reproduces Figure 5 of the manuscript and its supplementary cause-attribution plot.

## Quick start

```bash
pip install -r ../requirements.txt
python figure5.py
```

Runs in about 1-2 minutes on 2-4 CPU cores. The script:

1. Validates the discrete-event engine against seven analytical queueing references
   (M/M/1, M/D/1, accumulating-priority Eq. 15a, HARQ truncated-geometric). Aborts on
   failure; the check is mandatory and cannot be skipped.
2. Simulates 100 independent replications of 120 s of operation at `rho = 0.70` for
   each train speed in `{160, 200, 250, 300} km/h` and each of the four configurations
   C0, C1, C2, C3. Common random numbers are shared across configurations within a
   replication, so paired comparisons have very low variance.
3. Writes the main figure (outage probability vs train speed for the baseline C0 and
   the full proposed C3 only, with 95 % CI error bars; the full ablation is given in
   `Figure5_ablation_outage.csv`), one supplementary figure (per-cause attribution of
   the 100 ms outage intervals, reviewer R2 comment 9), and the full set of
   reproducibility artifacts.

## CLI options

| Flag | Default | Meaning |
|---|---|---|
| `--outdir` | `figure5_output` | directory for the outputs (created if missing) |
| `--workers` | CPU count | parallel worker processes |
| `--seed` | `20260101` | base random seed (same seed + same package versions = bit-reproducible) |

No `--quick` or `--sensitivity`: the script always runs the full study exactly once.

## Why the main figure plots only C0 and C3

Priority scheduling (the only mechanism that distinguishes C0 from C1) re-orders the
scheduler queue and has no effect on radio-layer availability, so the simulated P_out
of C0 is numerically identical to that of C1 at every speed (a consistency check that
the engine correctly separates the scheduling layer from the radio layer). The full
ablation including C1 and C2 is reported in `Figure5_ablation_outage.csv`; the main
figure plots only the baseline C0 and the full proposed C3 to match the original
manuscript's Figure 5 and to avoid plotting the two outage-indistinguishable
configurations as separate bars.

## Output files (archived in `figure5_output/`)

| File | Content |
|---|---|
| `Figure5_outage_vs_speed.{png,pdf}` | **Main Figure 5 for the manuscript** (C0 vs C3) |
| `Figure5_cause_attribution.{png,pdf}` | Supplementary: per-cause decomposition of P_out; addresses R2 comment 9 |
| `Figure5_PLOT_DATA.csv` | Exact numbers drawn in the main figure |
| `Figure5_summary.csv` | Per-(speed, config) outage/handover/link-failure metrics with 95 % CI |
| `Figure5_run_level.csv` | Every per-run per-config result |
| `Figure5_paired_stats.csv` | Paired bootstrap + Wilcoxon + Holm correction (C0 vs C3 and stepwise at every speed) |
| `Figure5_cause_breakdown.csv` | Per-cause outage interval counts and ratios per (speed, config); addresses R2 comment 9 |
| `Figure5_ablation_outage.csv` | Full 4-configuration ablation with the C0 ≡ C1 identity check and the C3 vs C2 trade-off at 300 km/h |
| `Figure5_P_out_analytic.csv` | Analytic nominal blackout budget vs simulated handover-related P_out (reviewer cross-check) |
| `Figure5_handover_kpis.csv` | Handover KPIs (success rate, failures, ping-pongs) distinct from P_out; addresses R3 |
| `Figure5_interruption_durations.csv` | Interruption-duration summary per (speed, config); derived from engine constants (NOT empirical per-blackout traces) |
| `Figure5_adaptive_parameter_table.csv` | Derivation and values of H(v) = 480/v and TTT(v) = 25600/v; addresses R2 comment 8 |
| `Figure5_model_summary.md` | Full Eq.(14) definition, cause taxonomy, H(v)/TTT(v) derivation, C0=C1 identity, RLF=0 explanation, C3 vs C2 trade-off; ready to adapt for Section 2.4 / 2.7 / 3.3 |
| `Figure5_manuscript_snippet.md` | Ready-to-paste numbers for Section 3.3, including honest trade-off narrative |
| `Figure5_review_point_mapping.md` | INTERNAL DEV MAPPING ONLY: links each reviewer comment that concerns Figure 5 to the exact artifacts here. NOT a reviewer-response letter. |
| `Figure5_manifest.json` | Reproducibility manifest: parameters with source tags (T4 / REF / ASSUMED), base seed, engine self-test, SHA-256 of this script, package versions, wall clock |

## Key result (100 runs per point, rho = 0.70)

Communication outage probability (%, 95 % CI from between-run Student-t):

| v (km/h) | C0 baseline | C3 full proposal | reduction |
|---|---|---|---|
| 160 | 1.357 % | 0.014 % | -99.0 % |
| 200 | 1.655 % | 0.025 % | -98.5 % |
| 250 | 1.904 % | 0.020 % | -98.9 % |
| 300 | 2.340 % | 0.041 % | -98.3 % |

C3 reduces outage probability by 98-99 % across the entire tested speed range, from
1.4-2.3 % to below 0.05 %. Per-cause attribution (Figure5_cause_attribution.png) shows
that the handover-related category (nominal 100 ms blackout plus command/target
failures) dominates the baseline outage at every speed; radio-link failure contributes
negligibly at the simulated rho.

Honest trade-off (manuscript Section 3.3): at 300 km/h the speed-adaptive handover rule
triggers roughly twice as many handovers as the baseline, slightly increasing the
residual outage of C3 (0.041 %) compared to C2 (0.022 %). The adaptive rule's benefit
is in Figure 3 (delay) and Figure 4 (deadline reliability), with a small outage-exposure
cost quantified here.

## Status

Re-implementation prepared in response to the reviewer comments. Not a restored copy of
the original MATLAB code (lost in a computer-virus incident). The engine is byte-for-byte
identical to the ones used for Figures 3 and 4 (same self-test, same cause codes). The
earlier (MATLAB-based) Figure 5 reported 10 % / 4 % at 300 km/h based on a non-physical
calibration; the new physically consistent numbers are 2.34 % / 0.04 % and are
accompanied by the per-cause attribution reviewer R2 explicitly requested.
