# Figure 3 — End-to-end delay of LTE-R control messages vs train speed

Reproduces Figure 3 of the manuscript and its supplementary sensitivity analysis.

## Quick start

```bash
pip install -r ../requirements.txt
python figure3.py                   # full study, ~3 min on 2-4 cores
```

With the sensitivity analysis (used in the reviewer response):

```bash
python figure3.py --sensitivity     # adds ~1 min, produces Figure3_sensitivity.{png,pdf,csv}
```

Options:

| Flag | Meaning |
|---|---|
| `--runs 100` | independent runs per (speed, config); default 100 |
| `--seed 20260101` | base seed; same seed + same package versions = bit-reproducible |
| `--workers 4` | parallel processes; default = CPU count |
| `--outdir figure3_output` | where to write the outputs |
| `--quick` | smoke test (8 runs, 40 s per run; result is tagged "NOT FOR PUBLICATION") |
| `--sensitivity` | also run the 6-parameter sensitivity sweep |
| `--skip-selftest` | skip engine validation against analytic queueing (not recommended) |

## What the script does

1. Validates the discrete-event engine against 7 analytic references
   (M/M/1, M/D/1, accumulating-priority Eq. 15a, HARQ truncated-geometric).
   The script aborts if any test fails.
2. Simulates the 2^3 = 8 on/off combinations of the three proposed mechanisms
   (priority scheduling, standby path, speed-adaptive handover) at
   speeds {160, 200, 250, 300} km/h, offered load rho = 0.70, 100 runs each.
3. Reports mean delay with 95 % CI, paired bootstrap CI for the C0 -> C3
   reduction, Wilcoxon signed-rank p-values with Holm correction, Cohen's dz.
4. Writes Figure 3 (single-panel grouped bar chart) plus every intermediate
   CSV and a reproducibility manifest.

## Output files (archived in `figure3_output/` after a full run)

| File | Content |
|---|---|
| `Figure3_grouped_bar.{png,pdf}` | **This is Figure 3 of the manuscript.** |
| `Figure3_PLOT_DATA.csv` | Exact numbers drawn in the figure |
| `Figure3_summary.csv` | Per-(speed, config) means, CIs, every secondary metric |
| `Figure3_run_level.csv` | Per-run per-config results (100 x 4 x 8 rows) |
| `Figure3_paired_C0_minus_C3.csv` | Paired bootstrap + Wilcoxon + Holm |
| `Figure3_stepwise_paired.csv` | Attribution to each mechanism (C0->C1->C2->C3) |
| `Figure3_ablation_main_effects.csv` | 2^3 factorial main effects |
| `Figure3_delay_decomposition.csv` | Break-down of mean delay (proc / tx / wait / retx / switch) |
| `Figure3_reliability.csv` | Deadline reliability P_del per (speed, config) |
| `Figure3_handover_summary.csv` | Mean handovers, ping-pongs, failures, RLFs |
| `Figure3_handover_events.csv` | Timestamped handover events (C0 and C3) |
| `Figure3_cross_checks.csv` | Analytic M/M/1 lower bound vs simulated C0, independence |
| `Figure3_sensitivity.{csv,png,pdf}` | Sensitivity of C0 and C3 to the six ASSUMED parameters |
| `Figure3_model_summary.md` | Model description (source for manuscript Section 2.4) |
| `Figure3_manuscript_snippet.md` | Ready-to-paste numbers for manuscript Section 3.1 |
| `Figure3_manifest.json` | Full reproducibility manifest: parameters, seeds, hashes, env |

## Key result

With rho = 0.70 and 100 independent replications per point, the full
configuration C3 reduces the mean end-to-end delay of control messages by
25.7 %, 33.1 %, 34.5 % and 38.8 % at 160, 200, 250 and 300 km/h respectively
(all Wilcoxon p < 1e-15, Holm-corrected). The baseline C0 does not reach
99 % deadline reliability above 200 km/h; C3 reaches at least 99.99 % at
every speed. See `Figure3_manuscript_snippet.md` for the full table.

## Status

Re-implementation prepared in response to the reviewer comments.
Not a restored copy of the original MATLAB code (lost in a computer-virus
incident; reviewers were informed by e-mail and granted a 10-day extension).
