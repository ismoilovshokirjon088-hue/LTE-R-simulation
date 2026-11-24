# Figure 4 — Deadline reliability of LTE-R control messages vs offered load

Reproduces Figure 4 of the manuscript and its supplementary delay CDF.

## Quick start

```bash
pip install -r ../requirements.txt
python figure4.py
```

The run takes about 3-4 minutes on 2-4 CPU cores. It:

1. Validates the discrete-event engine against seven analytical queueing references
   (M/M/1, M/D/1, accumulating-priority Eq. 15a, HARQ truncated-geometric). The run
   aborts if any test fails; this check is mandatory and cannot be skipped.
2. Simulates 100 independent replications per `(rho, config)` at `v = 200 km/h` for
   13 offered-load values `rho in {0.10, 0.20, ..., 0.95}` and 4 configurations
   (C0, C1, C2, C3). Common random numbers are used across configurations within a
   replication, so paired comparisons have very low variance.
3. Writes the main figure (deadline reliability vs offered load with 95 % CI bands),
   one supplementary figure (empirical delay CDF at `rho = 0.70` and `rho = 0.90`,
   reviewer-requested to explain the heavy tail), and the full set of reproducibility
   artifacts.

## CLI options

| Flag | Default | Meaning |
|---|---|---|
| `--outdir` | `figure4_output` | directory for the outputs (created if missing) |
| `--workers` | CPU count | parallel worker processes |
| `--seed` | `20260101` | base random seed (same seed + same package versions = bit-reproducible) |

There is no `--quick` or `--sensitivity` flag. The script always runs the full study
exactly once; this is the configuration reviewers are expected to run.

## Output files (archived in `figure4_output/`)

| File | Content |
|---|---|
| `Figure4_reliability_vs_load.{png,pdf}` | **Main Figure 4 for the manuscript** |
| `Figure4_delay_CDF.{png,pdf}` | Supplementary: delay CDF at two load values (addresses reviewer R2 comment 3) |
| `Figure4_PLOT_DATA.csv` | Exact numbers drawn in the main figure |
| `Figure4_summary.csv` | Per-(rho, config) mean, 95 % CI, between-run standard deviation, every secondary metric |
| `Figure4_run_level.csv` | Every per-run per-config result (5200 rows) |
| `Figure4_paired_stats.csv` | Paired bootstrap + Wilcoxon + Holm correction (C0 vs C3, and C0->C1, C1->C2, C2->C3 stepwise at every rho) |
| `Figure4_threshold_crossings.csv` | Largest rho at which each configuration meets the 99 % and 99.9 % reliability targets |
| `Figure4_reliability_analytic.csv` | Simulated C0 reliability vs the M/M/1 upper bound from manuscript Eq. 7a (reviewer cross-check of the engine) |
| `Figure4_independence_check.csv` | Lag-1 inter-run autocorrelation per (rho, config) with 99 % test interval |
| `Figure4_cdf_samples.csv` | Message-level delay samples for the CDF |
| `Figure4_model_summary.md` | Model description including the three heavy-tail mechanisms; ready to adapt for manuscript Section 2.4 |
| `Figure4_manuscript_snippet.md` | Ready-to-paste numbers for Section 3.2 narrative |
| `Figure4_review_point_mapping.md` | **Internal dev mapping only**: links each reviewer comment that concerns Figure 4 to the exact artifacts produced here. NOT a reviewer-response letter. |
| `Figure4_manifest.json` | Full reproducibility manifest: every parameter with source tag (T4 / REF / ASSUMED), base seed, engine self-test report, SHA-256 of this script, package versions, wall clock |

## Key result (100 runs per point, rho in 0.10-0.95, 200 km/h)

Largest offered load at which each configuration meets the reliability target:

| configuration | 99 % crossing | 99.9 % crossing |
|---|---|---|
| C0 baseline | rho = 0.61 | never within tested range |
| C1 +priority | rho = 0.80 | never within tested range |
| C2 +standby | rho = 0.87 | rho = 0.79 |
| C3 full | rho = 0.89 | rho = 0.81 |

All C0 vs C3 differences at `rho >= 0.30` are significant at Holm-corrected `p < 0.001`.
At 200 km/h the proposed full configuration C3 extends the 99 %-reliability-compliant
load range by about 0.28 compared to the baseline. The 99.9 % SIL-3-style target is
never met by C0; the FRMCS 99.9999 % target is not met by any configuration at any
simulated load and this is stated explicitly in the manuscript Sections 3.6 and 6.

## Status

Re-implementation prepared in response to the reviewer comments. Not a restored copy of
the original MATLAB code (lost in a computer-virus incident; reviewers were informed by
e-mail). The engine is identical to the one used for Figure 3.
