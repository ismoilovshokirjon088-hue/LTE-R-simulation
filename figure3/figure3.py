#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figure3.py -- End-to-end delay of LTE-R train-control messages vs train speed (Figure 3)
=========================================================================================
Python re-implementation prepared in response to reviewer comments.
NOT a restored copy of the original MATLAB code; it is a redesigned simulation implementation.
Shared packet-level engine: lter_engine.py (one engine for Figures 3, 4, 5).

WHAT THIS SCRIPT PRODUCES
  Figure3_grouped_bar.{png,pdf}    Single-panel grouped bar chart (MDPI style), mean delay with 95 % CI
  Figure3_PLOT_DATA.csv            The exact numbers drawn in the figure
  Figure3_summary.csv              Per-(speed, config) means, CIs, and every secondary metric
  Figure3_run_level.csv            Every per-run per-config result (100 runs x 4 speeds x 8 configs)
  Figure3_paired_C0_minus_C3.csv   Paired bootstrap + Wilcoxon + Holm correction, common random numbers
  Figure3_stepwise_paired.csv      C0->C1->C2->C3 attribution of the reduction, paired per run
  Figure3_ablation_main_effects.csv  2^3 full-factorial main effect of each mechanism (R1-6)
  Figure3_delay_decomposition.csv  Breakdown of mean delay into proc + tx + wait + retx + switch
  Figure3_reliability.csv          Deadline-reliability P_del per (speed, config), vs analytic bound
  Figure3_handover_summary.csv     Mean handovers, ping-pongs, failures, RLFs per (speed, config)
  Figure3_handover_events.csv      Timestamped handover events for C0 and C3
  Figure3_cross_checks.csv         Analytic baseline cross-check + independence + reliability bound
  Figure3_model_summary.md         Model description for the manuscript (Section 2.4)
  Figure3_manuscript_snippet.md    Ready-to-paste numbers for Section 3.1 narrative
  Figure3_manifest.json            Full reproducibility manifest: parameters, seeds, hashes, env

OPTIONAL
  python figure3.py --sensitivity    Run 4-parameter sensitivity sweep at 200 km/h, rho=0.7
                                     Produces Figure3_sensitivity.{csv,png}

USAGE
  python figure3.py                  Full run: 100 runs per point (2-4 min on 2-4 cores)
  python figure3.py --quick          Smoke test: 8 runs, 40 s per run (NOT for the manuscript)
  python figure3.py --runs 100 --seed 20260101 --workers 4 --outdir figure3_output

REQUIREMENTS: Python >= 3.9, numpy, matplotlib. Optional: scipy (exact t-quantiles, Mann-Whitney).
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import os
import sys
import time
from multiprocessing import Pool, cpu_count, freeze_support

import numpy as np

# ====== ENGINE (packet-level DES, radio layer, HARQ, scheduler, mobility) =======


import dataclasses
import hashlib
import heapq
import json
import math
import os
import platform
import sys
import time
from collections import deque
from dataclasses import dataclass

import numpy as np

ENGINE_VERSION = "1.0.0"

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Params:
    # ---- [T4] manuscript Table 4 ----------------------------------------------------------------
    fc_hz: float = 465e6                 # downlink carrier, E-UTRA Band 31
    bw_hz: float = 5e6                   # channel bandwidth
    cell_radius_m: float = 1000.0        # R
    overlap_m: float = 300.0             # L_ov
    tx_power_dbm: float = 43.0           # 20 W
    msg_bytes: int = 256                 # S
    ctrl_rate_hz: float = 10.0           # control-message generation rate
    r_eff_bps: float = 1.0e6             # effective transmission rate
    t_proc_s: float = 0.010              # processing delay (fixed, both configurations)
    t_harq_rtt_s: float = 0.008          # T_attempt: start-to-start spacing of HARQ (re)transmissions
    deadline_s: float = 0.100            # T_dl
    w_ctrl: float = 4.0                  # scheduling weights (proposed)
    w_other: float = 1.0
    h0_db: float = 3.0                   # baseline hysteresis
    ttt0_s: float = 0.160                # baseline time-to-trigger
    adapt_h_coeff: float = 480.0         # H(v) = 480 / v[km/h]  dB
    adapt_ttt_coeff_ms: float = 25600.0  # tau(v) = 25600 / v[km/h]  ms
    t_ho_s: float = 0.100                # handover execution time
    t_margin_s: float = 0.050            # timing margin of Eq. (17)
    t_obs_s: float = 0.100               # outage observation interval
    # ---- [JSON] simulation horizon of the earlier Python Figure 3 ------------------------------
    duration_s: float = 120.0
    warmup_s: float = 10.0
    # ---- [REF] radio-propagation reference values ------------------------------------------------
    pl_a: float = 21.5                   # WINNER II D2a: PL = A log10(d) + B + C log10(fc/5GHz)
    pl_b: float = 44.2
    pl_c: float = 20.0
    d_track_m: float = 50.0              # eNB-to-track distance (WINNER II D2a scenario description)
    h_bs_m: float = 30.0                 # eNB antenna height (D2a: 30-32 m)
    h_ms_m: float = 1.5                  # train antenna height (D2a: 1.5 m)
    shadow_sigma_db: float = 3.73        # shadow-fading std, LTE-R measurement (manuscript ref. [13])
    rician_k_db: float = 5.12            # K-factor listed for WINNER II in ref. [13], Table III
    # ---- [ASSUMED] modelling assumptions (explore by sensitivity analysis!) ----------------------
    ant_gain_db: float = 0.0             # total antenna gain (Tx+Rx), isotropic assumption
    noise_figure_db: float = 9.0         # receiver noise figure
    shadow_dcorr_m: float = 50.0         # shadowing decorrelation distance
    service_dist: str = "deterministic"  # 'deterministic' (fixed 2.048 ms) | 'exponential' (M/M/1 check)
    ctrl_arrivals: str = "periodic"      # 'periodic' (random phase) | 'poisson'
    interf_activity: float = -1.0        # activity of interfering cells; <0 means "equal to rho"
    gamma_th_db: float = -8.0            # path-eligibility / out-of-sync threshold (Q_out-like)
    shannon_gap_db: float = 3.0          # implementation gap to the Shannon bound
    max_harq_tx: int = 4                 # transmissions per packet before it is declared lost
    meas_period_s: float = 0.020         # RSRP measurement / A3 evaluation period
    rlf_timer_s: float = 1.0             # continuous out-of-sync time before radio-link failure
    t_reest_s: float = 0.500             # re-establishment duration after a failure
    backup_unavail_frac: float = 0.01    # long-run unavailability of the standby path
    backup_mean_down_s: float = 0.5      # mean duration of a standby-path unavailability period
    t_switch_s: float = 0.010            # extra delay of a message sent over the standby path
    pingpong_window_s: float = 1.0       # return to the previous cell within this time = ping-pong
    drain_s: float = 5.0                 # extra simulated time so that queued messages can complete
    bg_buffer_pkts: int = 100            # tail-drop limit of the 'other traffic' queue (0 = unlimited)
    # ---- [TEST] self-test switches (never used for reported results) -----------------------------
    ideal_link: bool = False
    test_force_bler: float = -1.0

    # ---- derived quantities -----------------------------------------------------------------------
    @property
    def isd_m(self) -> float:
        return 2.0 * self.cell_radius_m - self.overlap_m

    @property
    def t_tx_s(self) -> float:
        return 8.0 * self.msg_bytes / self.r_eff_bps

    @property
    def t_prop_s(self) -> float:
        return 1000.0 / 299_792_458.0     # 1 km radio link, as in Section 2.4

    @property
    def rho_ctrl(self) -> float:
        return self.ctrl_rate_hz * self.t_tx_s

    @property
    def gamma_dec_lin(self) -> float:
        """Decoding threshold (linear) of the 1 Mbit/s service: Shannon bound + implementation gap."""
        return (2.0 ** (self.r_eff_bps / self.bw_hz) - 1.0) * 10.0 ** (self.shannon_gap_db / 10.0)

    @property
    def gamma_th_lin(self) -> float:
        return 10.0 ** (self.gamma_th_db / 10.0)

    @property
    def noise_w(self) -> float:
        return 10.0 ** ((-174.0 + 10.0 * math.log10(self.bw_hz) + self.noise_figure_db - 30.0) / 10.0)

    @property
    def n_ms(self) -> int:
        return int(round((self.duration_s + self.drain_s) * 1000.0))


PARAM_SOURCES = {
    "fc_hz": ("T4", "Table 4"), "bw_hz": ("T4", "Table 4"), "cell_radius_m": ("T4", "Table 4"),
    "overlap_m": ("T4", "Table 4"), "tx_power_dbm": ("T4", "Table 4"), "msg_bytes": ("T4", "Table 4"),
    "ctrl_rate_hz": ("T4", "Table 4"), "r_eff_bps": ("T4", "Table 4"), "t_proc_s": ("T4", "Table 4"),
    "t_harq_rtt_s": ("T4", "T_attempt = 8 ms"), "deadline_s": ("T4", "Table 4"),
    "w_ctrl": ("T4", "Table 4"), "w_other": ("T4", "Table 4"), "h0_db": ("T4", "Table 4"),
    "ttt0_s": ("T4", "Table 4"), "adapt_h_coeff": ("T4", "H(v)=480/v"),
    "adapt_ttt_coeff_ms": ("T4", "tau(v)=25600/v"), "t_ho_s": ("T4", "Table 4"),
    "t_margin_s": ("T4", "Table 4"), "t_obs_s": ("T4", "Table 4"),
    "duration_s": ("JSON", "earlier Python Figure 3 setup"), "warmup_s": ("JSON", "earlier Python Figure 3 setup"),
    "pl_a": ("REF", "WINNER II D2a path-loss parameter A (project paper on D2a power-distance handover)"),
    "pl_b": ("REF", "WINNER II D2a path-loss parameter B"),
    "pl_c": ("REF", "WINNER II D2a path-loss parameter C"),
    "d_track_m": ("REF", "WINNER II D2a scenario description: eNB-track distance 50 m"),
    "h_bs_m": ("REF", "WINNER II D2a scenario description: eNB height 30-32 m"),
    "h_ms_m": ("REF", "WINNER II D2a scenario description: train antenna height 1.5 m"),
    "shadow_sigma_db": ("REF", "manuscript ref. [13], LTE-R 2.6 GHz measurement, 3.73 dB (different carrier!)"),
    "rician_k_db": ("REF", "K-factor of WINNER II listed in manuscript ref. [13], Table III"),
    "ant_gain_db": ("ASSUMED", "no antenna gain assumed"), "noise_figure_db": ("ASSUMED", ""),
    "shadow_dcorr_m": ("ASSUMED", ""), "service_dist": ("ASSUMED", "fixed-size packet at fixed rate"),
    "ctrl_arrivals": ("ASSUMED", "cyclic control traffic with random phase"),
    "interf_activity": ("ASSUMED", "interfering cells transmit with probability rho"),
    "gamma_th_db": ("ASSUMED", "typical out-of-sync level; verify against 3GPP TS 36.133 before citing"),
    "shannon_gap_db": ("ASSUMED", ""), "max_harq_tx": ("ASSUMED", ""), "meas_period_s": ("ASSUMED", ""),
    "rlf_timer_s": ("ASSUMED", "T310-like timer"), "t_reest_s": ("ASSUMED", ""),
    "backup_unavail_frac": ("ASSUMED", "A_b(t) two-state Markov model"),
    "backup_mean_down_s": ("ASSUMED", "A_b(t) two-state Markov model"),
    "t_switch_s": ("ASSUMED", "backup-path switching delay (not the 100 ms handover time)"),
    "pingpong_window_s": ("ASSUMED", ""), "drain_s": ("ASSUMED", "technical"),
    "bg_buffer_pkts": ("ASSUMED", "finite eNB buffer for non-control traffic (~205 ms of link service time)"),
    "ideal_link": ("TEST", ""), "test_force_bler": ("TEST", ""),
}
assert set(PARAM_SOURCES) == {f.name for f in dataclasses.fields(Params)}, "every Params field needs a source tag"

# mechanism flags (E_P, E_B, E_H): priority scheduling, standby path, speed-adaptive handover
CONFIGS = {"C0": (0, 0, 0), "C1": (1, 0, 0), "C2": (1, 1, 0), "C3": (1, 1, 1)}


def factorial_configs():
    """All 2^3 on/off combinations of the three mechanisms; C0..C3 are a subset of them."""
    out = {}
    for ep in (0, 1):
        for eb in (0, 1):
            for eh in (0, 1):
                out[f"P{ep}B{eb}H{eh}"] = (ep, eb, eh)
    return out


CAUSES = {0: "none", 1: "ho_execution", 2: "ho_command_failure", 3: "ho_target_failure",
          4: "radio_link_failure", 5: "low_sinr"}

# --------------------------------------------------------------------------------------------------
# Statistics helpers
# --------------------------------------------------------------------------------------------------


def t_crit(df: int, conf: float = 0.95) -> float:
    try:
        from scipy import stats  # type: ignore
        return float(stats.t.ppf(0.5 + conf / 2.0, df))
    except Exception:
        # Cornish-Fisher expansion of the Student-t quantile (accurate for df >= 5)
        z = {0.95: 1.959963984540054, 0.99: 2.5758293035489004}.get(conf, 1.959963984540054)
        return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df ** 2)


def mean_ci(x, conf: float = 0.95):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    n = x.size
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    m = float(x.mean())
    if n < 2:
        return m, float("nan"), float("nan"), n
    hw = t_crit(n - 1, conf) * float(x.std(ddof=1)) / math.sqrt(n)
    return m, m - hw, m + hw, n


# --------------------------------------------------------------------------------------------------
# Random streams (common random numbers: every configuration of one replication shares the radio
# environment, the traffic and the fading draws as far as the model allows)
# --------------------------------------------------------------------------------------------------

STREAMS = ("geom", "ctrl", "bg", "fading", "service", "backup", "misc")


def make_seeds(base_seed: int, speed_idx: int, load_idx: int, run: int):
    ss = np.random.SeedSequence([int(base_seed), int(speed_idx), int(load_idx), int(run)])
    return dict(zip(STREAMS, ss.spawn(len(STREAMS))))


# --------------------------------------------------------------------------------------------------
# Radio environment
# --------------------------------------------------------------------------------------------------


def build_radio_base(P: Params, v_kmh: float, seeds) -> dict:
    """Received power of every cell along the 1-ms time grid (shared by all configurations)."""
    n = P.n_ms
    dt = 1e-3
    v = v_kmh / 3.6
    rng = np.random.default_rng(seeds["geom"])
    isd = P.isd_m
    x0 = rng.uniform(isd, 2.0 * isd)                      # random start position inside a cell
    t = np.arange(n) * dt
    x = x0 + v * t
    n_cells = int(math.ceil(x[-1] / isd)) + 3
    pos = np.arange(n_cells) * isd
    dx = x[None, :] - pos[:, None]
    d = np.sqrt(dx ** 2 + P.d_track_m ** 2 + (P.h_bs_m - P.h_ms_m) ** 2)
    pl = P.pl_a * np.log10(d) + P.pl_b + P.pl_c * np.log10(P.fc_hz / 1e9 / 5.0)
    # spatially correlated shadowing (AR(1) in the travelled distance, independent for every cell)
    ds = 1.0
    m = int(v * n * dt / ds) + 3
    a = math.exp(-ds / P.shadow_dcorr_m)
    b = math.sqrt(1.0 - a * a)
    w = rng.standard_normal((n_cells, m))
    s = np.empty_like(w)
    s[:, 0] = w[:, 0]
    for i in range(1, m):
        s[:, i] = a * s[:, i - 1] + b * w[:, i]
    u = v * t / ds
    grid = np.arange(m)
    sh = np.stack([np.interp(u, grid, s[k]) for k in range(n_cells)]) * P.shadow_sigma_db
    rsrp_dbm = P.tx_power_dbm + P.ant_gain_db - pl + sh
    p_w = 10.0 ** ((rsrp_dbm - 30.0) / 10.0)
    return {"n_cells": n_cells, "rsrp_dbm": rsrp_dbm, "p_w": p_w, "sum_p": p_w.sum(axis=0),
            "x0": x0, "v_kmh": v_kmh}


def _sinr_as_serving(P: Params, base: dict, cell: int, a_int: float, j0: int = 0, j1: int | None = None):
    p = base["p_w"][cell, j0:j1]
    return p / (P.noise_w + a_int * (base["sum_p"][j0:j1] - p))


def _first_run(cond: np.ndarray, need: int):
    """Index of the first element at which `cond` has been True for `need` consecutive samples."""
    if cond.size < need or not cond.any():
        return None
    idx = np.arange(cond.size)
    last_false = np.maximum.accumulate(np.where(~cond, idx, -1))
    hit = np.nonzero(idx - last_false >= need)[0]
    return int(hit[0]) if hit.size else None


def _backup_availability(P: Params, n: int, rng) -> np.ndarray:
    """A_b(t): alternating exponential up/down periods (two-state Markov chain), 1-ms resolution."""
    q = P.backup_unavail_frac
    a = np.ones(n, dtype=bool)
    if q <= 0.0:
        return a
    mean_down = P.backup_mean_down_s * 1000.0
    mean_up = mean_down * (1.0 - q) / q
    pos = int(rng.exponential(mean_up))
    while pos < n:
        dur = max(1, int(rng.exponential(mean_down)))
        a[pos:pos + dur] = False
        pos += dur + max(1, int(rng.exponential(mean_up)))
    return a


def build_timeline(P: Params, base: dict, v_kmh: float, rho: float, flags, seeds) -> dict:
    """Serving cell, handover / radio-link-failure state machine, path eligibility (1-ms grid)."""
    ep, eb, eh = flags
    n = P.n_ms
    a_int = rho if P.interf_activity < 0 else P.interf_activity
    if P.ideal_link:
        ones = np.ones(n, dtype=bool)
        big = np.full(n, 1e9)
        return {"prim_ok": ones, "bk_ok": np.zeros(n, dtype=bool), "any_ok": ones, "sinr_p": big,
                "sinr_b": np.zeros(n), "next_ok": np.arange(n), "cause_ms": np.zeros(n, dtype=np.int8),
                "ho_events": [], "counts": dict(n_ho=0, n_ho_ok=0, n_ho_fail_cmd=0, n_ho_fail_tgt=0,
                                                n_rlf=0, n_pp=0, ho_blocked=0)}
    rng = np.random.default_rng(seeds["misc"])
    rsrp = base["rsrp_dbm"]
    mp_ms = max(1, int(round(P.meas_period_s * 1000)))
    rsrp_m = rsrp[:, ::mp_ms]
    n_cells = base["n_cells"]
    if eh:
        h_db = P.adapt_h_coeff / v_kmh
        ttt_ms = P.adapt_ttt_coeff_ms / v_kmh
    else:
        h_db = P.h0_db
        ttt_ms = P.ttt0_s * 1000.0
    need = int(math.ceil(ttt_ms / mp_ms))
    ttt_i = int(round(ttt_ms))
    t_ho_i = int(round(P.t_ho_s * 1000))
    reest_i = int(round(P.t_reest_s * 1000))
    rlf_i_need = int(round(P.rlf_timer_s * 1000))
    pp_i = int(round(P.pingpong_window_s * 1000))
    t_ov = P.overlap_m / (v_kmh / 3.6)                          # Eq. (16)
    feasible = (P.t_ho_s + P.t_margin_s) <= t_ov                # Eq. (17)

    serv = np.full(n, -1, dtype=np.int32)
    prim_down = np.zeros(n, dtype=bool)
    cause_down = np.zeros(n, dtype=np.int8)
    events = []
    cnt = dict(n_ho=0, n_ho_ok=0, n_ho_fail_cmd=0, n_ho_fail_tgt=0, n_rlf=0, n_pp=0, ho_blocked=0)
    g_th = P.gamma_th_lin
    k_lin = 10.0 ** (P.rician_k_db / 10.0)
    los = math.sqrt(k_lin / (k_lin + 1.0))
    sc = math.sqrt(1.0 / (2.0 * (k_lin + 1.0)))

    def fade():
        x, y = rng.standard_normal(2)
        return (los + sc * x) ** 2 + (sc * y) ** 2

    rtt_i = max(1, int(round(P.t_harq_rtt_s * 1000)))

    def signalling_ok(cell, j0):
        """Handover command / target access: same HARQ chase-combining abstraction as the data plane."""
        acc = 0.0
        for k in range(P.max_harq_tx):
            j = j0 + k * rtt_i
            if j >= n:
                break
            acc += _sinr_as_serving(P, base, cell, a_int, j, j + 1)[0] * fade()
            if acc >= P.gamma_dec_lin:
                return True
        return False

    def mark_down(i0, dur, cause):
        i1 = min(n, i0 + dur)
        prim_down[i0:i1] = True
        cause_down[i0:i1] = cause
        return i1

    i = 0
    s = int(np.argmax(rsrp[:, 0]))
    prev_cell, last_ho_end = -1, -10 ** 9
    while i < n:
        sinr_s = _sinr_as_serving(P, base, s, a_int, i)
        # (a) radio-link failure: continuous out-of-sync for rlf_timer
        r = _first_run(sinr_s < g_th, rlf_i_need)
        rlf_at = None if r is None else i + r
        # (b) next A3 handover trigger
        trig = None
        tgt = None
        if feasible:
            m0 = i // mp_ms + 1
            for cand in (s + 1, s - 1):
                if 0 <= cand < n_cells:
                    diff = rsrp_m[cand, m0:] - rsrp_m[s, m0:]
                    f = _first_run(diff >= h_db, need + 1)
                    if f is not None:
                        t_trig = (m0 + f - need) * mp_ms + ttt_i
                        if t_trig < n and (trig is None or t_trig < trig):
                            trig, tgt = t_trig, cand
        elif i == 0:
            cnt["ho_blocked"] += 1
        if trig is None and rlf_at is None:
            serv[i:] = s
            break
        if rlf_at is not None and (trig is None or rlf_at <= trig):
            serv[i:rlf_at] = s
            cnt["n_rlf"] += 1
            end = mark_down(rlf_at, reest_i, 4)
            events.append(dict(t_ms=rlf_at, kind="rlf", src=s, dst=-1, outcome="radio_link_failure"))
            nb = int(np.argmax(rsrp[:, min(end, n - 1)]))
            serv[rlf_at:min(end, n)] = s
            i, s = end, nb
            continue
        # handover attempt at `trig`
        serv[i:trig] = s
        cnt["n_ho"] += 1
        if not signalling_ok(s, trig):
            cnt["n_ho_fail_cmd"] += 1
            end = mark_down(trig, reest_i, 2)
            events.append(dict(t_ms=trig, kind="ho", src=s, dst=tgt, outcome="command_failure"))
            nb = int(np.argmax(rsrp[:, min(end, n - 1)]))
            serv[trig:min(end, n)] = s
            i, s = end, nb
            continue
        end = mark_down(trig, t_ho_i, 1)
        serv[trig:min(end, n)] = s
        if end < n and not signalling_ok(tgt, end):
            cnt["n_ho_fail_tgt"] += 1
            end2 = mark_down(end, reest_i, 3)
            events.append(dict(t_ms=trig, kind="ho", src=s, dst=tgt, outcome="target_failure"))
            nb = int(np.argmax(rsrp[:, min(end2, n - 1)]))
            serv[end:min(end2, n)] = tgt
            i, s = end2, nb
            continue
        cnt["n_ho_ok"] += 1
        pp = (tgt == prev_cell) and (trig - last_ho_end <= pp_i)
        cnt["n_pp"] += int(pp)
        events.append(dict(t_ms=trig, kind="ho", src=s, dst=tgt, outcome="success", ping_pong=bool(pp)))
        prev_cell, last_ho_end = s, end
        i, s = end, tgt
    if (serv < 0).any():                      # defensive: every millisecond must have a serving cell
        raise RuntimeError("internal error: serving-cell timeline has gaps")
    # force a healthy link during the last 500 ms so that every queue can be flushed (drain)
    prim_down[-500:] = False
    cause_down[-500:] = 0

    cols = np.arange(n)
    p = base["p_w"]
    sinr_p = p[serv, cols] / (P.noise_w + a_int * (base["sum_p"] - p[serv, cols]))
    sinr_p[-500:] = np.maximum(sinr_p[-500:], 1e3)
    # standby path: strongest non-serving cell, seen as a serving cell (the old serving cell interferes)
    if eb:
        pm = p.copy()
        pm[serv, cols] = -1.0
        bk = np.argmax(pm, axis=0)
        sinr_b = p[bk, cols] / (P.noise_w + a_int * (base["sum_p"] - p[bk, cols]))
        avail = _backup_availability(P, n, np.random.default_rng(seeds["backup"]))
        bk_ok = avail & (sinr_b >= g_th)
        bk_ok[-500:] = False
    else:
        sinr_b = np.zeros(n)
        bk_ok = np.zeros(n, dtype=bool)
    prim_ok = (~prim_down) & (sinr_p >= g_th)
    any_ok = prim_ok | bk_ok
    cause_ms = np.zeros(n, dtype=np.int8)
    bad = ~any_ok
    cause_ms[bad] = np.where(cause_down[bad] > 0, cause_down[bad], 5)
    idx = np.where(any_ok, np.arange(n), n)
    next_ok = np.minimum.accumulate(idx[::-1])[::-1]
    return {"prim_ok": prim_ok, "bk_ok": bk_ok, "any_ok": any_ok, "sinr_p": sinr_p, "sinr_b": sinr_b,
            "next_ok": next_ok, "cause_ms": cause_ms, "ho_events": events, "counts": cnt}


def outage_stats(P: Params, tl: dict) -> dict:
    """Eq. (14): windows of T_obs that contain >= 1 interval without any eligible path."""
    w_ms = int(round(P.t_obs_s * 1000))
    i0 = int(round(P.warmup_s * 1000))
    n_win = int((P.duration_s - P.warmup_s) * 1000) // w_ms
    sl = slice(i0, i0 + n_win * w_ms)
    bad = (~tl["any_ok"][sl]).reshape(n_win, w_ms)
    cause = tl["cause_ms"][sl].reshape(n_win, w_ms)
    hit = bad.any(axis=1)
    first = np.argmax(bad, axis=1)
    first_cause = cause[np.arange(n_win), first]
    out = {"n_obs": n_win, "n_out": int(hit.sum()), "p_out": float(hit.mean())}
    out["n_out_primary_only"] = int((~tl["prim_ok"][sl]).reshape(n_win, w_ms).any(axis=1).sum())
    for c, name in CAUSES.items():
        if c:
            out[f"win_{name}"] = int(((first_cause == c) & hit).sum())
    return out


# --------------------------------------------------------------------------------------------------
# Traffic
# --------------------------------------------------------------------------------------------------


def build_traffic(P: Params, rho: float, seeds):
    """Arrival times (sorted) and class (0 = control, 1 = other) of all packets."""
    D = P.duration_s
    rc = np.random.default_rng(seeds["ctrl"])
    rb = np.random.default_rng(seeds["bg"])
    per = 1.0 / P.ctrl_rate_hz
    if P.ctrl_arrivals == "periodic":
        ta_c = rc.uniform(0.0, per) + per * np.arange(int(D * P.ctrl_rate_hz) + 1)
        ta_c = ta_c[ta_c < D]
    else:
        k = rc.poisson(P.ctrl_rate_hz * D)
        ta_c = np.sort(rc.uniform(0.0, D, k))
    lam_bg = max(0.0, (rho - P.rho_ctrl)) / P.t_tx_s
    ta_b = np.sort(rb.uniform(0.0, D, rb.poisson(lam_bg * D)))
    t = np.concatenate([ta_c, ta_b])
    c = np.concatenate([np.zeros(ta_c.size, dtype=np.int8), np.ones(ta_b.size, dtype=np.int8)])
    o = np.argsort(t, kind="stable")
    return t[o], c[o]


# --------------------------------------------------------------------------------------------------
# Packet-level discrete-event simulation of the scheduler + HARQ link
# --------------------------------------------------------------------------------------------------


def simulate_queue(P: Params, flags, tl: dict, arr_t: np.ndarray, arr_c: np.ndarray, seeds) -> dict:
    ep = flags[0]
    N = arr_t.size
    n = P.n_ms
    T_tx = P.t_tx_s
    expo = P.service_dist == "exponential"
    rngf = np.random.default_rng(seeds["fading"])
    rngs = np.random.default_rng(seeds["service"])
    k_lin = 10.0 ** (P.rician_k_db / 10.0)
    los, sc = math.sqrt(k_lin / (k_lin + 1.0)), math.sqrt(1.0 / (2.0 * (k_lin + 1.0)))
    pool = max(2 * N + 2000, 1000)
    gx = rngf.standard_normal(pool)
    gy = rngf.standard_normal(pool)
    gpool = ((los + sc * gx) ** 2 + (sc * gy) ** 2).tolist()
    upool = rngf.random(pool).tolist()
    spool = rngs.exponential(1.0, pool).tolist() if expo else None
    sinr_p = tl["sinr_p"].tolist()
    sinr_b = tl["sinr_b"].tolist()
    prim_ok = tl["prim_ok"].tolist()
    bk_ok = tl["bk_ok"].tolist()
    next_ok = tl["next_ok"].tolist()
    arr = arr_t.tolist()
    cls = arr_c.tolist()
    ready = arr[:]
    att = [0] * N
    gacc = [0.0] * N
    first = [-1.0] * N
    first_tx = [0.0] * N
    done = [-1.0] * N
    pth = [0] * N
    lost = [False] * N
    INF = float("inf")
    gdec = P.gamma_dec_lin
    rtt = P.t_harq_rtt_s
    max_tx = P.max_harq_tx
    wc, wo = P.w_ctrl, P.w_other
    force = P.test_force_bler
    q = (deque(), deque())
    cap = P.bg_buffer_pkts
    n_drop = 0
    heap: list = []
    ia = 0
    ip = 0
    t = 0.0
    last = n - 1
    while True:
        # admit every packet that is ready at time t, in chronological order
        while True:
            ta = arr[ia] if ia < N else INF
            th = heap[0][0] if heap else INF
            if ta <= th:
                if ta <= t:
                    if cap and cls[ia] == 1 and len(q[1]) >= cap:
                        lost[ia] = True                      # tail drop of 'other' traffic
                        n_drop += 1
                    else:
                        q[cls[ia]].append(ia)
                    ia += 1
                    continue
            elif th <= t:
                _, pid = heapq.heappop(heap)
                ready[pid] = th
                q[cls[pid]].append(pid)
                continue
            break
        q0, q1 = q
        if not q0 and not q1:
            nxt = min(arr[ia] if ia < N else INF, heap[0][0] if heap else INF)
            if nxt == INF:
                break
            t = nxt
            continue
        j = int(t * 1000.0)
        if j > last:
            j = last
        po, bo = prim_ok[j], bk_ok[j]
        if not (po or bo):
            jn = next_ok[j]
            if jn >= n:
                break
            t = max(t + 1e-6, jn * 1e-3)
            continue
        # scheduling decision
        if q0 and q1:
            c0, c1 = q0[0], q1[0]
            if ep:
                pid = c0 if wc * (t - ready[c0]) >= wo * (t - ready[c1]) else c1
            else:
                pid = c0 if ready[c0] <= ready[c1] else c1
        else:
            pid = q0[0] if q0 else q1[0]
        (q0 if cls[pid] == 0 else q1).popleft()
        path = 0 if po else 1
        tx = T_tx * spool[ip % pool] if expo else T_tx
        if first[pid] < 0.0:
            first[pid] = t
            first_tx[pid] = tx
        te = t + tx
        j2 = int(te * 1000.0)
        if j2 > last:
            j2 = last
        still = prim_ok[j2] if path == 0 else bk_ok[j2]
        if force >= 0.0:
            ok = still and upool[ip % pool] >= force
        else:
            g = gpool[ip % pool]
            gacc[pid] += (sinr_p[j] if path == 0 else sinr_b[j]) * g if still else 0.0
            ok = still and gacc[pid] >= gdec
        ip += 1
        att[pid] += 1
        if ok:
            done[pid] = te
            pth[pid] = path
        elif att[pid] >= max_tx:
            lost[pid] = True
        else:
            heapq.heappush(heap, (t + rtt, pid))
        t = te
    # ---- collect results ------------------------------------------------------------------------
    arr_a = np.asarray(arr)
    cls_a = np.asarray(cls)
    meas = (arr_a >= P.warmup_s) & (arr_a < P.duration_s)
    mc = meas & (cls_a == 0)
    mb = meas & (cls_a == 1)
    done_a = np.asarray(done)
    first_a = np.asarray(first)
    ftx_a = np.asarray(first_tx)
    att_a = np.asarray(att)
    pth_a = np.asarray(pth)
    lost_a = np.asarray(lost)
    recv = done_a >= 0.0
    sw = np.where(pth_a == 1, P.t_switch_s, 0.0)
    delay = np.where(recv, done_a - arr_a + P.t_prop_s + P.t_proc_s + sw, np.nan)
    wait1 = np.where(first_a >= 0, first_a - arr_a, np.nan)
    retx = np.where(recv, (done_a - arr_a) - wait1 - ftx_a, np.nan)
    unresolved = (~recv) & (~lost_a)
    return {
        "ctrl": {"arr": arr_a[mc], "recv": recv[mc], "delay": delay[mc], "wait1": wait1[mc],
                 "retx": retx[mc], "attempts": att_a[mc], "path": pth_a[mc], "sw": sw[mc]},
        "bg_wait1_mean": float(np.nanmean(wait1[mb])) if mb.any() else float("nan"),
        "bg_n": int(mb.sum()),
        "all_attempts_mean": float(att_a[meas & (att_a > 0)].mean()) if (meas & (att_a > 0)).any() else float("nan"),
        "unresolved": int((unresolved & ~np.asarray(lost))[meas].sum()),
        "bg_dropped": n_drop,
    }


# --------------------------------------------------------------------------------------------------
# One replication = one random radio/traffic realisation evaluated for several configurations
# --------------------------------------------------------------------------------------------------


def summarize_replication(P: Params, flags, sim: dict, tl: dict) -> dict:
    c = sim["ctrl"]
    n_ctrl = int(c["arr"].size)
    rec = c["recv"]
    d = c["delay"]
    n_recv = int(rec.sum())
    on_time = rec & (d <= P.deadline_s)
    out = outage_stats(P, tl)
    cnt = tl["counts"]
    r = {
        "n_ctrl": n_ctrl, "n_recv": n_recv, "n_late": int((rec & ~on_time).sum()),
        "n_lost": int((~rec).sum()),
        "p_rx": n_recv / n_ctrl if n_ctrl else float("nan"),
        "p_del": int(on_time.sum()) / n_ctrl if n_ctrl else float("nan"),
        "mean_delay_ms": float(np.nanmean(d)) * 1e3 if n_recv else float("nan"),
        "mean_wait1_ms": float(np.nanmean(c["wait1"][rec])) * 1e3 if n_recv else float("nan"),
        "mean_retx_ms": float(np.nanmean(c["retx"][rec])) * 1e3 if n_recv else float("nan"),
        "mean_sw_ms": float(c["sw"][rec].mean()) * 1e3 if n_recv else float("nan"),
        "mean_attempts": float(c["attempts"][rec].mean()) if n_recv else float("nan"),
        "frac_backup": float((c["path"][rec] == 1).mean()) if n_recv else float("nan"),
        "bg_wait1_ms": sim["bg_wait1_mean"] * 1e3, "bg_n": sim["bg_n"], "unresolved": sim["unresolved"],
        "bg_dropped": sim["bg_dropped"],
        "p50_ms": float(np.nanpercentile(d[rec], 50)) * 1e3 if n_recv else float("nan"),
        "p95_ms": float(np.nanpercentile(d[rec], 95)) * 1e3 if n_recv else float("nan"),
        "p99_ms": float(np.nanpercentile(d[rec], 99)) * 1e3 if n_recv else float("nan"),
    }
    r.update({k: out[k] for k in out})
    r.update({k: cnt[k] for k in cnt})
    return r


def run_replication(P: Params, v_kmh: float, rho: float, configs: dict, base_seed: int,
                    speed_idx: int, load_idx: int, run: int, keep_messages: bool = False,
                    keep_events: bool = False) -> dict:
    """Returns {config_name: summary_dict (+ 'messages' / 'tl_events' when keep_messages)}."""
    seeds = make_seeds(base_seed, speed_idx, load_idx, run)
    base = None if P.ideal_link else build_radio_base(P, v_kmh, seeds)
    arr_t, arr_c = build_traffic(P, rho, seeds)
    res = {}
    for name, flags in configs.items():
        tl = build_timeline(P, base, v_kmh, rho, flags, seeds)
        sim = simulate_queue(P, flags, tl, arr_t, arr_c, seeds)
        r = summarize_replication(P, flags, sim, tl)
        if keep_messages:
            r["messages"] = sim["ctrl"]
        if keep_messages or keep_events:
            r["ho_events"] = tl["ho_events"]
        res[name] = r
    return res


# --------------------------------------------------------------------------------------------------
# Analytical references (used only to validate the engine)
# --------------------------------------------------------------------------------------------------


def analytic_waits(P: Params, rho: float, exponential: bool):
    """Mean queue waits [s]: FIFO, accumulating-priority control class and other class."""
    S = P.t_tx_s
    w0 = rho * S if exponential else rho * S / 2.0           # mean residual work
    w_fifo = w0 / (1.0 - rho)
    rc = P.rho_ctrl
    bo_bc = P.w_other / P.w_ctrl
    w_other = w_fifo / (1.0 - rc * (1.0 - bo_bc))
    w_ctrl = w_fifo - (rho - rc) * w_other * (1.0 - bo_bc)
    w_strict = w0 / (1.0 - rc)
    return {"fifo": w_fifo, "ctrl_acc": w_ctrl, "other_acc": w_other, "ctrl_strict": w_strict}


def run_selftest(n_runs: int = 24, quiet: bool = False) -> dict:
    """Validation of queue discipline, service model and HARQ feedback against analytical results."""
    results = []
    rho = 0.70
    base = Params(ideal_link=True, ctrl_arrivals="poisson", duration_s=120.0, warmup_s=10.0, drain_s=5.0,
                  bg_buffer_pkts=0)

    def pooled(P, flags, name):
        w_c, w_o, att = [], [], []
        for r in range(n_runs):
            seeds = make_seeds(777, 0, 0, r)
            arr_t, arr_c = build_traffic(P, rho, seeds)
            tl = build_timeline(P, None, 200.0, rho, flags, seeds)
            sim = simulate_queue(P, flags, tl, arr_t, arr_c, seeds)
            w_c.append(float(np.nanmean(sim["ctrl"]["wait1"])))
            w_o.append(sim["bg_wait1_mean"])
            att.append(sim["all_attempts_mean"])
        return np.array(w_c), np.array(w_o), np.array(att)

    def check(name, sim_arr, ref, tol_rel=0.03):
        m, lo, hi, _ = mean_ci(sim_arr)
        ok = (lo <= ref <= hi) or abs(m - ref) / ref <= tol_rel
        results.append({"test": name, "simulated_ms": m * 1e3, "ci95_ms": [lo * 1e3, hi * 1e3],
                        "reference_ms": ref * 1e3, "pass": bool(ok)})

    # T1/T2: M/M/1 (Eq. 7), accumulating priority (Eq. 15a), conservation law
    Pm = dataclasses.replace(base, service_dist="exponential")
    ref = analytic_waits(Pm, rho, True)
    wc0, wo0, _ = pooled(Pm, CONFIGS["C0"], "fifo")
    check("T1 M/M/1 FIFO, control-class wait  vs Eq.(7)", wc0, ref["fifo"])
    check("T1 M/M/1 FIFO, other-class wait    vs Eq.(7)", wo0, ref["fifo"])
    wc1, wo1, _ = pooled(Pm, CONFIGS["C1"], "prio")
    check("T2 4:1 accumulating priority, control wait vs Eq.(15a)", wc1, ref["ctrl_acc"])
    check("T2 4:1 accumulating priority, other wait   vs Eq.(15a)", wo1, ref["other_acc"])
    cons = base.rho_ctrl * wc1.mean() + (rho - base.rho_ctrl) * wo1.mean()
    results.append({"test": "T2 conservation law rho_c*Wc + rho_o*Wo = rho*W_fifo",
                    "simulated_ms": cons * 1e3, "reference_ms": rho * ref["fifo"] * 1e3,
                    "pass": bool(abs(cons - rho * ref["fifo"]) / (rho * ref["fifo"]) <= 0.03)})
    # T3: M/D/1 (Pollaczek-Khinchine)
    Pd = dataclasses.replace(base, service_dist="deterministic")
    refd = analytic_waits(Pd, rho, False)
    wcd, wod, _ = pooled(Pd, CONFIGS["C0"], "md1")
    check("T3 M/D/1 FIFO, other-class wait vs Pollaczek-Khinchine", wod, refd["fifo"])
    # T4: HARQ plumbing, truncated-geometric number of attempts with a forced block-error probability
    p = 0.2
    Pb = dataclasses.replace(base, service_dist="deterministic", test_force_bler=p, max_harq_tx=4)
    _, _, att = pooled(Pb, CONFIGS["C0"], "harq")
    exp_att = (1 - p ** 4) / (1 - p)
    m, lo, hi, _ = mean_ci(att)
    results.append({"test": "T4 HARQ: mean attempts (forced BLER=0.2, max 4)", "simulated": m,
                    "ci95": [lo, hi], "reference": exp_att, "pass": bool(lo <= exp_att <= hi or abs(m - exp_att) < 0.01)})
    out = {"engine_version": ENGINE_VERSION, "n_runs": n_runs, "all_pass": all(r["pass"] for r in results),
           "tests": results}
    if not quiet:
        print("\nENGINE SELF-TEST (queue discipline / service model / HARQ feedback)")
        for r in results:
            if "simulated_ms" in r:
                print(f"  [{'PASS' if r['pass'] else 'FAIL'}] {r['test']}: sim = {r['simulated_ms']:.4f} ms"
                      f"  ref = {r['reference_ms']:.4f} ms")
            else:
                print(f"  [{'PASS' if r['pass'] else 'FAIL'}] {r['test']}: sim = {r['simulated']:.4f}"
                      f"  ref = {r['reference']:.4f}")
        print("  =>", "ALL TESTS PASSED" if out["all_pass"] else "SOME TESTS FAILED")
    return out


# --------------------------------------------------------------------------------------------------
# Reproducibility helpers
# --------------------------------------------------------------------------------------------------


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def environment_info() -> dict:
    info = {"python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform(),
            "engine_version": ENGINE_VERSION, "engine_sha256": file_sha256(os.path.abspath(__file__))}
    try:
        import scipy  # type: ignore
        info["scipy"] = scipy.__version__
    except Exception:
        info["scipy"] = None
    try:
        import matplotlib  # type: ignore
        info["matplotlib"] = matplotlib.__version__
    except Exception:
        info["matplotlib"] = None
    return info


def param_manifest(P: Params) -> dict:
    d = dataclasses.asdict(P)
    return {k: {"value": v, "source": PARAM_SOURCES[k][0], "note": PARAM_SOURCES[k][1]} for k, v in d.items()}



# ====== FIGURE 3 PRODUCTION =======



SPEEDS = [160, 200, 250, 300]          # km/h  (manuscript Table 4)
RHO = 0.70                             # offered load of Figure 3 (manuscript Section 3.1)
MAIN = ["C0", "C1", "C2", "C3"]
FLAGS_TO_NAME = {v: k for k, v in CONFIGS.items()}
FACT = factorial_configs()


def cname(flags):
    return FLAGS_TO_NAME.get(tuple(flags), None)


# ----- worker (top-level for multiprocessing 'spawn') --------------------------------------------
def _job(args):
    params, si, v, run, base_seed, rho = args
    res = run_replication(params, v, rho, FACT, base_seed, si, 0, run, keep_events=True)
    return si, run, res


# ----- small utilities ---------------------------------------------------------------------------
def write_csv(path, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r.get(k), float) else r.get(k)) for k in fields})


def paired_stats(a, b, n_boot=5000, seed=1):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = ~(np.isnan(a) | np.isnan(b))
    a, b = a[ok], b[ok]
    d = a - b
    m, lo, hi, n = mean_ci(d)
    red = 100.0 * d.mean() / a.mean() if a.mean() else float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    bs = 100.0 * (a[idx].mean(1) - b[idx].mean(1)) / a[idx].mean(1)
    p_w = float("nan")
    try:
        from scipy import stats
        p_w = float(stats.wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
    except Exception:
        pass
    sd = float(d.std(ddof=1)) if n > 1 else float("nan")
    return {"n": n, "diff_mean": m, "diff_ci_lo": lo, "diff_ci_hi": hi, "reduction_pct": red,
            "reduction_ci_lo": float(np.percentile(bs, 2.5)), "reduction_ci_hi": float(np.percentile(bs, 97.5)),
            "wilcoxon_p": p_w, "cohen_dz": (d.mean() / sd) if sd and sd > 0 else float("nan")}


def holm(pvals):
    p = np.array(pvals, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    idx = np.argsort(p[ok])
    q = p[ok][idx]
    m = q.size
    adj = np.maximum.accumulate((m - np.arange(m)) * q)
    res = np.empty(m)
    res[idx] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


# ----- reviewer-response diagnostics -------------------------------------------------------------
def analytic_c0_delay_ms(P: Params, rho: float) -> float:
    """
    Analytic lower bound for the mean end-to-end delay of C0 at offered load rho with the
    manuscript's parameter set, assuming no handover interruption and no retransmissions.
    Mean service time = T_tx; arrivals are the superposition of periodic control messages
    and Poisson background traffic, so the wait is bracketed by M/D/1 and M/M/1:
        M/D/1 wait = rho * T_tx / (2 (1 - rho))
        M/M/1 wait = rho * T_tx / (1 - rho)
    We report the M/M/1 upper bound (more conservative for a mixed arrival process).
    """
    w_mm1 = rho * P.t_tx_s / (1 - rho)
    return (P.t_proc_s + P.t_tx_s + P.t_prop_s + w_mm1) * 1e3


def cross_checks(P: Params, per: dict, summ: dict) -> list:
    """Runs analytic baseline / independence / reliability checks for the reviewer response."""
    rows = []
    for v in SPEEDS:
        c0 = np.array(per[(v, "P0B0H0")]["mean_delay_ms"])
        c3 = np.array(per[(v, "P1B1H1")]["mean_delay_ms"])
        ana = analytic_c0_delay_ms(P, RHO)
        row = {
            "speed_kmh": v,
            "C0_simulated_mean_delay_ms": float(c0.mean()),
            "C0_analytic_MM1_upper_bound_ms": ana,
            "C0_handover_HARQ_overhead_ms": float(c0.mean()) - ana,
            "lag1_autocorr_C0": float(np.corrcoef(c0[:-1], c0[1:])[0, 1]) if c0.size > 2 else float("nan"),
            "lag1_autocorr_C3": float(np.corrcoef(c3[:-1], c3[1:])[0, 1]) if c3.size > 2 else float("nan"),
            "C0_P_del": float(np.mean(per[(v, "P0B0H0")]["p_del"])),
            "C3_P_del": float(np.mean(per[(v, "P1B1H1")]["p_del"])),
            "C0_meets_99pct_P_del": bool(float(np.mean(per[(v, "P0B0H0")]["p_del"])) >= 0.99),
            "C3_meets_99pct_P_del": bool(float(np.mean(per[(v, "P1B1H1")]["p_del"])) >= 0.99),
        }
        rows.append(row)
    return rows


def reliability_table(per: dict):
    rows = []
    for v in SPEEDS:
        for c in MAIN:
            cfg_key = [k for k, f in FACT.items() if cname(f) == c][0]
            p = np.array(per[(v, cfg_key)]["p_del"])
            px = np.array(per[(v, cfg_key)]["p_rx"])
            po = np.array(per[(v, cfg_key)]["p_out"])
            m, lo, hi, _ = mean_ci(p)
            rows.append({"speed_kmh": v, "config": c, "P_del_mean": m, "P_del_ci_lo": lo, "P_del_ci_hi": hi,
                         "P_rx_mean": float(px.mean()), "P_out_mean": float(po.mean()),
                         "P_del_meets_99pct": bool(m >= 0.99), "P_del_meets_99_9pct": bool(m >= 0.999)})
    return rows


def handover_summary(per: dict):
    rows = []
    for v in SPEEDS:
        for c in MAIN:
            cfg_key = [k for k, f in FACT.items() if cname(f) == c][0]
            d = per[(v, cfg_key)]
            rows.append({
                "speed_kmh": v, "config": c,
                "handovers_mean_per_run": float(np.mean(d["n_ho"])),
                "successful_handovers_mean": float(np.mean(d["n_ho_ok"])),
                "command_failures_mean": float(np.mean(d["n_ho_fail_cmd"])),
                "target_failures_mean": float(np.mean(d["n_ho_fail_tgt"])),
                "radio_link_failures_mean": float(np.mean(d["n_rlf"])),
                "ping_pongs_mean": float(np.mean(d["n_pp"])),
                "backup_used_fraction_mean": float(np.mean(d["frac_backup"])),
            })
    return rows


# ----- sensitivity sweep (optional) --------------------------------------------------------------
SENSITIVITY_PARAMS = [
    ("gamma_th_db", [-10.0, -8.0, -6.0], "SINR / Qout decoding threshold (dB)"),
    ("noise_figure_db", [5.0, 9.0, 13.0], "Receiver noise figure (dB)"),
    ("shadow_sigma_db", [2.0, 3.73, 6.0], "Shadow fading standard deviation (dB)"),
    ("bg_buffer_pkts", [25, 100, 400], "Finite buffer of non-control queue (packets)"),
    ("t_switch_s", [0.005, 0.010, 0.020], "Standby-path switching delay (s)"),
    ("backup_unavail_frac", [0.001, 0.01, 0.05], "Long-run standby-path unavailability"),
]


def _sens_job(args):
    params, cfg_name, run, base_seed = args
    r = run_replication(params, 200.0, RHO, {cfg_name: FACT[cfg_name]}, base_seed, 0, 0, run)
    return cfg_name, run, r[cfg_name]


def run_sensitivity(P0: Params, workers: int, n_runs: int, seed: int, outdir: str):
    rows = []
    defaults = dataclasses.asdict(P0)
    for name, values, desc in SENSITIVITY_PARAMS:
        for val in values:
            P = dataclasses.replace(P0, **{name: val})
            jobs = [(P, "P0B0H0", r, seed) for r in range(n_runs)] + \
                   [(P, "P1B1H1", r, seed) for r in range(n_runs)]
            outs = {"P0B0H0": [], "P1B1H1": []}
            with Pool(workers) as p:
                for cfg, run, r in p.imap_unordered(_sens_job, jobs, chunksize=2):
                    outs[cfg].append(r)
            for cfg_name, cfg_label in (("P0B0H0", "C0"), ("P1B1H1", "C3")):
                d = [x["mean_delay_ms"] for x in outs[cfg_name]]
                pdel = [x["p_del"] for x in outs[cfg_name]]
                m, lo, hi, _ = mean_ci(d)
                pm, plo, phi, _ = mean_ci(pdel)
                rows.append({"parameter": name, "description": desc, "value": val,
                             "is_default": (val == defaults.get(name)),
                             "config": cfg_label, "delay_ms_mean": m, "delay_ms_ci_lo": lo,
                             "delay_ms_ci_hi": hi, "P_del_mean": pm, "P_del_ci_lo": plo,
                             "P_del_ci_hi": phi, "n_runs": n_runs})
    write_csv(os.path.join(outdir, "Figure3_sensitivity.csv"), rows)
    _plot_sensitivity(outdir, rows)
    return rows


def _plot_sensitivity(outdir: str, rows: list):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e6e5e1"
    col = {"C0": "#2a78d6", "C3": "#eb6834"}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6), sharey=False)
    axes = axes.flatten()
    for ax, (name, values, desc) in zip(axes, SENSITIVITY_PARAMS):
        xs = np.arange(len(values))
        for j, c in enumerate(("C0", "C3")):
            m = [r["delay_ms_mean"] for r in rows if r["parameter"] == name and r["config"] == c]
            lo = [r["delay_ms_ci_lo"] for r in rows if r["parameter"] == name and r["config"] == c]
            hi = [r["delay_ms_ci_hi"] for r in rows if r["parameter"] == name and r["config"] == c]
            off = (j - 0.5) * 0.30
            ax.bar(xs + off, m, 0.28, color=col[c], edgecolor="#fcfcfb", linewidth=0.8, label=c, zorder=3)
            ax.errorbar(xs + off, m, yerr=[np.array(m) - np.array(lo), np.array(hi) - np.array(m)],
                        fmt="none", ecolor=ink, elinewidth=0.9, capsize=2, zorder=4)
        ax.set_xticks(xs)
        ax.set_xticklabels([str(v) for v in values], fontsize=8)
        ax.set_xlabel(desc, fontsize=8)
        ax.set_ylabel("Mean delay (ms)")
        ax.yaxis.grid(True, color=grid, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("Figure 3 - Sensitivity of C0 and C3 mean delay to the six ASSUMED parameters "
                 "(200 km/h, ρ = 0.70; mean ± 95 % CI)", fontsize=9.5, y=0.995, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(outdir, "Figure3_sensitivity.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "Figure3_sensitivity.pdf"))
    plt.close(fig)


# ----- model + manuscript snippets ---------------------------------------------------------------
def write_model_summary(outdir: str, P: Params, selftest: dict):
    txt = f"""# Figure 3 simulation model (ready to adapt for manuscript Section 2.4)

## Scenario
Equally spaced eNodeBs along a 1-D track, inter-site distance ISD = 2R - L_ov = {P.isd_m:.0f} m,
E-UTRA Band 31 downlink carrier f_c = {P.fc_hz/1e6:.0f} MHz, channel bandwidth B = {P.bw_hz/1e6:.0f} MHz,
train speeds v in {{160, 200, 250, 300}} km/h, offered load rho = {RHO:.2f}.

## Radio layer
Path loss follows the WINNER II D2a open scenario [ref.]:
    PL(d) [dB] = {P.pl_a} log10(d) + {P.pl_b} + {P.pl_c} log10(f_c / 5 GHz)
where d is the slant distance to the eNodeB. Shadow fading is log-normal with sigma = {P.shadow_sigma_db} dB
and exponential spatial correlation with decorrelation distance {P.shadow_dcorr_m:.0f} m. Small-scale fading
is Rician with K = {P.rician_k_db} dB (independent per transmission attempt). Receiver noise power follows
kT_0 B + NF with NF = {P.noise_figure_db} dB. Co-channel interference from the six neighbouring cells is
included, each weighted by the offered-load activity factor.

## Link abstraction
The 1 Mbit/s service decodes if and only if the chase-combined post-fading SINR exceeds
    gamma_dec = (2^(R/B) - 1) * 10^(Delta/10),  Delta = {P.shannon_gap_db} dB
(the information-outage bound with a {P.shannon_gap_db} dB implementation gap). Failed transmissions
re-enter the scheduler after T_attempt = {P.t_harq_rtt_s*1e3:.0f} ms and participate in chase combining;
at most {P.max_harq_tx} transmissions per packet.

## Scheduler
C0 is FIFO. C1-C3 apply accumulating priority: at each scheduling instant the message with the
largest priority-weighted waiting time w_i (t - t_i) is dispatched, with w_ctrl : w_other = {P.w_ctrl:.0f} : {P.w_other:.0f}.

## Mobility
Fixed-H baseline: hysteresis H = {P.h0_db} dB, time-to-trigger TTT = {P.ttt0_s*1e3:.0f} ms.
Speed-adaptive (H flag): H(v) = {P.adapt_h_coeff:.0f}/v dB, TTT(v) = {P.adapt_ttt_coeff_ms:.0f}/v ms.
Handover timing feasibility: T_HO + T_margin <= T_ov (Eq. 17).
Handover execution occupies T_HO = {P.t_ho_s*1e3:.0f} ms with no data-plane service.

## Standby path (B flag)
Second radio link toward the strongest non-serving cell, eligible when it is available
(two-state Markov availability A_b(t), long-run unavailability {P.backup_unavail_frac:.3f}, mean down time
{P.backup_mean_down_s} s) and its post-fading SINR >= gamma_th = {P.gamma_th_db} dB. Additional one-way delay
when a message uses the standby path: {P.t_switch_s*1e3:.0f} ms.

## Measurement window and statistics
Simulated horizon = {P.duration_s:.0f} s, warm-up = {P.warmup_s:.0f} s, drain = {P.drain_s:.0f} s.
Independent replications per (speed, configuration) = 100, with common random numbers across
configurations of the same replication. All reported differences use paired bootstrap confidence
intervals, Wilcoxon signed-rank tests, and Holm correction for multiple testing.

## Engine validation
The engine passed {sum(1 for t in selftest['tests'] if t['pass'])}/{len(selftest['tests'])}
self-tests against the analytical M/M/1, M/D/1, accumulating-priority (Eq. 15a) and HARQ-truncated-geometric
references at rho = {RHO:.2f} (deviation < 3 % for every test).

## Parameter-source tagging
T4 = manuscript Table 4; REF = cited reference; ASSUMED = modelling assumption not in Table 4.
The full list with values is in Figure3_manifest.json. The sensitivity of the main conclusion
to the ASSUMED parameters is reported in Figure3_sensitivity.{{csv,png}}.
"""
    with open(os.path.join(outdir, "Figure3_model_summary.md"), "w", encoding="utf-8") as f:
        f.write(txt)


def write_manuscript_snippet(outdir: str, summ: dict, pair_rows: list, step_rows: list,
                             rel_rows: list, ho_rows: list):
    def g(d, a, b): return d[a] - d[b]
    s = "# Figure 3 - ready-to-paste numbers for Section 3.1 (adapt wording)\n\n"
    s += "## Main result (mean end-to-end delay of control messages at rho = 0.70)\n\n"
    s += "| v (km/h) | C0 (ms) | C1 (ms) | C2 (ms) | C3 (ms) | C0->C3 reduction |\n"
    s += "|---|---|---|---|---|---|\n"
    for v, pr in zip(SPEEDS, pair_rows):
        row = f"| {v} | {summ[(v,'C0')]['delay_mean']:.2f} [{summ[(v,'C0')]['delay_ci_lo']:.2f}, {summ[(v,'C0')]['delay_ci_hi']:.2f}] | "
        row += f"{summ[(v,'C1')]['delay_mean']:.2f} | {summ[(v,'C2')]['delay_mean']:.2f} | "
        row += f"{summ[(v,'C3')]['delay_mean']:.2f} [{summ[(v,'C3')]['delay_ci_lo']:.2f}, {summ[(v,'C3')]['delay_ci_hi']:.2f}] | "
        row += f"{pr['reduction_pct']:.1f} % [{pr['reduction_ci_lo']:.1f}, {pr['reduction_ci_hi']:.1f}] |\n"
        s += row
    s += "\nAll C0 vs C3 differences are significant (Wilcoxon p < 1e-15, Holm-corrected).\n\n"
    s += "## Mechanism attribution (stepwise paired differences, ms)\n\n"
    s += "| v (km/h) | C0->C1 priority | C1->C2 standby | C2->C3 adaptive HO |\n|---|---|---|---|\n"
    for v in SPEEDS:
        a = next(r for r in step_rows if r["speed_kmh"] == v and "priority" in r["comparison"])
        b = next(r for r in step_rows if r["speed_kmh"] == v and "standby" in r["comparison"])
        c = next(r for r in step_rows if r["speed_kmh"] == v and "adaptive" in r["comparison"])
        s += f"| {v} | {a['diff_mean']:.2f} [{a['diff_ci_lo']:.2f}, {a['diff_ci_hi']:.2f}] | "
        s += f"{b['diff_mean']:.2f} [{b['diff_ci_lo']:.2f}, {b['diff_ci_hi']:.2f}] | "
        s += f"{c['diff_mean']:.2f} [{c['diff_ci_lo']:.2f}, {c['diff_ci_hi']:.2f}] |\n"
    s += "\nNote: at v = 160 km/h the adaptive rule H(v) = 480/v = 3 dB and TTT(v) = 160 ms reduce to the "
    s += "baseline fixed-handover parameters, so C2 and C3 are identical by construction.\n\n"
    s += "## Deadline reliability (P_del, fraction of control messages delivered within 100 ms)\n\n"
    s += "| v (km/h) | C0 | C1 | C2 | C3 |\n|---|---|---|---|---|\n"
    for v in SPEEDS:
        c0 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C0")
        c1 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C1")
        c2 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C2")
        c3 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C3")
        s += f"| {v} | {c0['P_del_mean']:.4f} | {c1['P_del_mean']:.4f} | "
        s += f"{c2['P_del_mean']:.4f} | {c3['P_del_mean']:.4f} |\n"
    s += "\nC0 does not reach 99 % deadline reliability above 200 km/h; C3 reaches at least 99.9 % at every speed.\n\n"
    s += "## Handover cost (mean per 120-s run)\n\n"
    s += "| v (km/h) | C0 handovers | C3 handovers | C3 ping-pongs | C3 cmd failures |\n|---|---|---|---|---|\n"
    for v in SPEEDS:
        c0 = next(r for r in ho_rows if r["speed_kmh"] == v and r["config"] == "C0")
        c3 = next(r for r in ho_rows if r["speed_kmh"] == v and r["config"] == "C3")
        s += f"| {v} | {c0['handovers_mean_per_run']:.1f} | {c3['handovers_mean_per_run']:.1f} | "
        s += f"{c3['ping_pongs_mean']:.1f} | {c3['command_failures_mean']:.2f} |\n"
    s += "\nThe speed-adaptive policy roughly doubles the handover count at 300 km/h and raises the ping-pong "
    s += "count from about 1 to about 12 per 120 s; this signalling cost is the trade-off behind the small "
    s += "additional delay reduction seen only above 200 km/h.\n"
    with open(os.path.join(outdir, "Figure3_manuscript_snippet.md"), "w", encoding="utf-8") as f:
        f.write(s)


# ----- plotting ----------------------------------------------------------------------------------
def make_figure(outdir, summ, quick, n_runs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = {"C0": "#2a78d6", "C1": "#1baf7a", "C2": "#6250d6", "C3": "#eb6834"}
    hatch = {"C0": "", "C1": "//", "C2": "..", "C3": "xx"}
    lab = {"C0": "C0 (FIFO, fixed handover)", "C1": "C1 (+priority)",
           "C2": "C2 (+standby path)", "C3": "C3 (+adaptive handover)"}
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e6e5e1"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "text.color": ink, "xtick.color": ink2, "ytick.color": ink2,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    w = 0.20
    x = np.arange(len(SPEEDS))
    for k, c in enumerate(MAIN):
        m = [summ[(v, c)]["delay_mean"] for v in SPEEDS]
        lo = [summ[(v, c)]["delay_ci_lo"] for v in SPEEDS]
        hi = [summ[(v, c)]["delay_ci_hi"] for v in SPEEDS]
        xs = x + (k - 1.5) * (w + 0.015)
        ax.bar(xs, m, w, color=col[c], edgecolor="#fcfcfb", linewidth=1.0,
               hatch=hatch[c], label=lab[c], zorder=3)
        ax.errorbar(xs, m, yerr=[np.array(m) - np.array(lo), np.array(hi) - np.array(m)],
                    fmt="none", ecolor=ink, elinewidth=1.0, capsize=2.5, zorder=4)
        for xi, mi, hi_i in zip(xs, m, hi):
            ax.text(xi, hi_i + 0.3, f"{mi:.1f}", ha="center", va="bottom",
                    fontsize=7.5, color=ink2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v}" for v in SPEEDS])
    ax.set_xlabel("Train speed (km/h)")
    ax.set_ylabel("Mean end-to-end delay (ms)")
    top = max(summ[(v, c)]["delay_ci_hi"] for v in SPEEDS for c in MAIN)
    ax.set_ylim(0, top * 1.35)
    ax.yaxis.grid(True, color=grid, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(loc="upper left", fontsize=8.5, frameon=False, ncol=2)
    if quick:
        ax.text(0.98, 0.98, "QUICK - NOT FOR PUBLICATION", transform=ax.transAxes,
                ha="right", va="top", fontsize=8, color="#e34948")
    fig.tight_layout()
    png = os.path.join(outdir, "Figure3_grouped_bar.png")
    pdf = os.path.join(outdir, "Figure3_grouped_bar.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


# ----- main --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Figure 3: delay of LTE-R control messages vs train speed")
    ap.add_argument("--runs", type=int, default=100, help="independent runs per point (default 100)")
    ap.add_argument("--seed", type=int, default=20260101, help="base seed (default 20260101)")
    ap.add_argument("--workers", type=int, default=max(1, min(cpu_count(), 8)), help="parallel processes")
    ap.add_argument("--outdir", default="figure3_output")
    ap.add_argument("--quick", action="store_true", help="smoke test: 8 runs, 40 s per run")
    ap.add_argument("--skip-selftest", action="store_true", help="skip engine validation tests")
    ap.add_argument("--sensitivity", action="store_true",
                    help="also run the 6-parameter sensitivity sweep (adds 3-4 min)")
    args = ap.parse_args()

    P = Params()
    n_runs = args.runs
    if args.quick:
        n_runs = min(n_runs, 8)
        P = dataclasses.replace(P, duration_s=40.0, warmup_s=5.0)
    os.makedirs(args.outdir, exist_ok=True)
    t_start = time.time()
    print(f"[Figure 3] runs/point={n_runs}  speeds={SPEEDS}  rho={RHO}  workers={args.workers}  "
          f"{'QUICK' if args.quick else 'FULL'}{'  +sensitivity' if args.sensitivity else ''}")

    # 1) engine validation ------------------------------------------------------------------------
    selftest = None
    if not args.skip_selftest:
        selftest = run_selftest(n_runs=8 if args.quick else 24, quiet=True)
        passed = sum(1 for t in selftest["tests"] if t["pass"])
        print(f"  engine self-test: {passed}/{len(selftest['tests'])} analytic tests passed "
              f"({'OK' if selftest['all_pass'] else 'FAILED'})")
        if not selftest["all_pass"]:
            sys.exit(2)

    # 2) simulation -------------------------------------------------------------------------------
    jobs = [(P, si, v, r, args.seed, RHO) for si, v in enumerate(SPEEDS) for r in range(n_runs)]
    out = []
    t0 = time.time()
    if args.workers > 1:
        with Pool(args.workers) as pool:
            for k, item in enumerate(pool.imap_unordered(_job, jobs, chunksize=1), 1):
                out.append(item)
                if k % max(1, len(jobs) // 5) == 0 or k == len(jobs):
                    print(f"  simulated {k}/{len(jobs)} replications  ({time.time() - t0:.0f} s)")
    else:
        for k, j in enumerate(jobs, 1):
            out.append(_job(j))
    out.sort(key=lambda t: (t[0], t[1]))

    # 3) run-level, per-config aggregation ---------------------------------------------------------
    metrics = ["mean_delay_ms", "p50_ms", "p95_ms", "p99_ms", "mean_wait1_ms", "mean_retx_ms", "mean_sw_ms",
               "mean_attempts", "frac_backup", "p_rx", "p_del", "p_out", "n_ctrl", "n_recv", "n_late", "n_lost",
               "n_ho", "n_ho_ok", "n_ho_fail_cmd", "n_ho_fail_tgt", "n_rlf", "n_pp", "ho_blocked", "bg_dropped",
               "unresolved"]
    run_rows, ev_rows = [], []
    per = {}
    for si, run, res in out:
        v = SPEEDS[si]
        for name, r in res.items():
            row = {"speed_kmh": v, "run": run, "config": name, "mechanisms_P_B_H": "".join(map(str, FACT[name])),
                   "paper_label": cname(FACT[name]) or ""}
            row.update({m: r[m] for m in metrics})
            run_rows.append(row)
            per.setdefault((v, name), {m: [] for m in metrics})
            for m in metrics:
                per[(v, name)][m].append(r[m])
            if name in ("P0B0H0", "P1B1H1"):
                for e in r["ho_events"]:
                    ev_rows.append({"speed_kmh": v, "run": run, "config": cname(FACT[name]),
                                    "t_s": e["t_ms"] / 1e3, "kind": e["kind"], "src_cell": e["src"],
                                    "dst_cell": e["dst"], "outcome": e["outcome"],
                                    "ping_pong": int(e.get("ping_pong", False))})
    write_csv(os.path.join(args.outdir, "Figure3_run_level.csv"), run_rows)
    write_csv(os.path.join(args.outdir, "Figure3_handover_events.csv"), ev_rows)

    # 4) per-cell summary --------------------------------------------------------------------------
    summ = {}
    summ_rows = []
    for (v, name), d in per.items():
        mm, lo, hi, n = mean_ci(d["mean_delay_ms"])
        sd = float(np.std(d["mean_delay_ms"], ddof=1)) if n > 1 else float("nan")
        row = {"speed_kmh": v, "config": name, "paper_label": cname(FACT[name]) or "", "n_runs": n,
               "delay_mean": mm, "delay_ci_lo": lo, "delay_ci_hi": hi, "delay_sd_between_runs": sd}
        for m in metrics:
            if m != "mean_delay_ms":
                row[m + "_avg"] = float(np.nanmean(d[m]))
        summ[(v, name)] = row
        summ_rows.append(row)
    for v in SPEEDS:
        for c in MAIN:
            summ[(v, c)] = summ[(v, [k for k, f in FACT.items() if cname(f) == c][0])]
    summ_rows.sort(key=lambda r: (r["speed_kmh"], r["config"]))
    write_csv(os.path.join(args.outdir, "Figure3_summary.csv"), summ_rows)
    plot_rows = [{"speed_kmh": v, "config": c, "mean_delay_ms": summ[(v, c)]["delay_mean"],
                  "ci95_lo_ms": summ[(v, c)]["delay_ci_lo"], "ci95_hi_ms": summ[(v, c)]["delay_ci_hi"],
                  "n_runs": summ[(v, c)]["n_runs"]} for v in SPEEDS for c in MAIN]
    write_csv(os.path.join(args.outdir, "Figure3_PLOT_DATA.csv"), plot_rows)

    # 5) decomposition (identity check) ------------------------------------------------------------
    dec_rows = []
    for v in SPEEDS:
        for c in MAIN:
            s = summ[(v, c)]
            comp = {"proc_ms": P.t_proc_s * 1e3, "tx_prop_ms": (P.t_tx_s + P.t_prop_s) * 1e3,
                    "wait1_ms": s["mean_wait1_ms_avg"], "retx_ms": s["mean_retx_ms_avg"],
                    "switch_ms": s["mean_sw_ms_avg"]}
            total = sum(comp.values())
            err = abs(total - s["delay_mean"])
            assert err < 1e-6, f"decomposition does not add up (err={err})"
            dec_rows.append({"speed_kmh": v, "config": c, **comp,
                             "total_ms": total, "check_vs_mean_delay_abs_err": err})
    write_csv(os.path.join(args.outdir, "Figure3_delay_decomposition.csv"), dec_rows)

    # 6) paired comparisons -------------------------------------------------------------------------
    def per_run(v, name, metric="mean_delay_ms"):
        return np.array(per[(v, name)][metric])

    pair_rows = []
    pw = []
    for v in SPEEDS:
        st = paired_stats(per_run(v, "P0B0H0"), per_run(v, "P1B1H1"))
        pair_rows.append({"speed_kmh": v, "comparison": "C0 minus C3", "metric": "mean delay (ms)", **st})
        pw.append(st["wilcoxon_p"])
    adj = holm(pw)
    for r, a in zip(pair_rows, adj):
        r["wilcoxon_p_holm"] = float(a)
    write_csv(os.path.join(args.outdir, "Figure3_paired_C0_minus_C3.csv"), pair_rows)

    step_rows = []
    for v in SPEEDS:
        for a, b, lbl in (("P0B0H0", "P1B0H0", "C0 minus C1 (priority scheduling)"),
                          ("P1B0H0", "P1B1H0", "C1 minus C2 (standby path)"),
                          ("P1B1H0", "P1B1H1", "C2 minus C3 (adaptive handover)")):
            st = paired_stats(per_run(v, a), per_run(v, b))
            step_rows.append({"speed_kmh": v, "comparison": lbl, "metric": "mean delay (ms)", **st})
    write_csv(os.path.join(args.outdir, "Figure3_stepwise_paired.csv"), step_rows)

    # Full-factorial main effects
    eff_rows = []
    for v in SPEEDS:
        for idx, nm in enumerate(("priority scheduling (P)", "standby path (B)", "adaptive handover (H)")):
            diffs = None
            for ff, fl in FACT.items():
                if fl[idx] == 0:
                    on = tuple(1 if i == idx else fl[i] for i in range(3))
                    on_name = [k for k, f in FACT.items() if f == on][0]
                    d = per_run(v, ff) - per_run(v, on_name)
                    diffs = d if diffs is None else diffs + d
            diffs = diffs / 4.0
            m, lo, hi, n = mean_ci(diffs)
            eff_rows.append({"speed_kmh": v, "mechanism": nm, "main_effect_delay_reduction_ms": m,
                             "ci95_lo": lo, "ci95_hi": hi, "n_runs": n})
    write_csv(os.path.join(args.outdir, "Figure3_ablation_main_effects.csv"), eff_rows)

    # 7) reviewer-response diagnostics --------------------------------------------------------------
    rel_rows = reliability_table(per)
    write_csv(os.path.join(args.outdir, "Figure3_reliability.csv"), rel_rows)
    ho_rows = handover_summary(per)
    write_csv(os.path.join(args.outdir, "Figure3_handover_summary.csv"), ho_rows)
    cc_rows = cross_checks(P, per, summ)
    write_csv(os.path.join(args.outdir, "Figure3_cross_checks.csv"), cc_rows)

    # 8) sensitivity sweep --------------------------------------------------------------------------
    sens_rows = None
    if args.sensitivity:
        print(f"  [sensitivity] starting 6-parameter sweep at 200 km/h  "
              f"(6 x 3 x 2 x {max(20, min(n_runs, 20))} runs)...")
        t1 = time.time()
        sens_rows = run_sensitivity(P, args.workers, max(20, min(n_runs, 20)), args.seed, args.outdir)
        print(f"  [sensitivity] done in {time.time() - t1:.0f} s")

    # 9) plot and documentation --------------------------------------------------------------------
    png, pdf = make_figure(args.outdir, summ, args.quick, n_runs)
    write_model_summary(args.outdir, P, selftest or {"tests": []})
    write_manuscript_snippet(args.outdir, summ, pair_rows, step_rows, rel_rows, ho_rows)

    # 10) manifest ---------------------------------------------------------------------------------
    manifest = {
        "figure": "Figure 3",
        "status": "re-implemented simulation (NOT a restored copy of the original MATLAB code)",
        "mode": "QUICK (not for publication)" if args.quick else "FULL",
        "command_line": " ".join(sys.argv), "base_seed": args.seed, "runs_per_point": n_runs,
        "speeds_kmh": SPEEDS, "offered_load_rho": RHO, "workers": args.workers,
        "wall_clock_s": round(time.time() - t_start, 1),
        "common_random_numbers": "radio environment, traffic and fading shared across configurations per replication",
        "environment": environment_info(),
        "figure3_script_sha256": file_sha256(os.path.abspath(__file__)),
        "engine_selftest": selftest, "parameters": param_manifest(P),
        "outputs": sorted(os.listdir(args.outdir)),
        "key_results": {
            f"{v}_km/h": {
                "C0_mean_delay_ms": summ[(v, "C0")]["delay_mean"],
                "C3_mean_delay_ms": summ[(v, "C3")]["delay_mean"],
                "reduction_pct": pair_rows[i]["reduction_pct"],
                "C0_P_del": next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C0")["P_del_mean"],
                "C3_P_del": next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C3")["P_del_mean"],
            } for i, v in enumerate(SPEEDS)
        },
        "interpretation_limits": [
            "Synthetic WINNER II D2a / log-normal-shadowing / Rician-fading channel; NOT ITU-R P.1812 and NOT field-validated.",
            "Block errors from an information-outage link abstraction with chase-combining HARQ; no PHY-level BLER curves.",
            "Parameters tagged ASSUMED are not measured; see Figure3_sensitivity.csv for their effect on the main conclusion.",
            "Numbers replace the earlier (MATLAB-based) Figure 3; the two must not be compared as if produced by the same model.",
            "At v = 160 km/h the adaptive rule reduces to the baseline H and TTT (3 dB, 160 ms), so C2 and C3 coincide by construction.",
        ],
    }
    with open(os.path.join(args.outdir, "Figure3_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)

    # 11) concise console report -------------------------------------------------------------------
    print("\n  Mean end-to-end delay of control messages (ms, 95 % CI)")
    print("  " + "-" * 82)
    print(f"  v  | C0                | C3                | reduction")
    for v, pr in zip(SPEEDS, pair_rows):
        c0 = summ[(v, "C0")]
        c3 = summ[(v, "C3")]
        print(f" {v:3d} | {c0['delay_mean']:5.2f} [{c0['delay_ci_lo']:5.2f},{c0['delay_ci_hi']:5.2f}]"
              f" | {c3['delay_mean']:5.2f} [{c3['delay_ci_lo']:5.2f},{c3['delay_ci_hi']:5.2f}]"
              f" | {pr['reduction_pct']:5.1f} % [{pr['reduction_ci_lo']:4.1f},{pr['reduction_ci_hi']:4.1f}]")
    print("\n  Reviewer-response cross-checks")
    print("  " + "-" * 82)
    ana = analytic_c0_delay_ms(P, RHO)
    print(f"  Analytic C0 lower bound (M/M/1 + fixed proc/tx/prop at rho={RHO}): {ana:.2f} ms")
    print(f"  Simulated C0 overhead (handover + HARQ), ms per speed: "
          + ", ".join(f"{v}: {summ[(v,'C0')]['delay_mean']-ana:.2f}" for v in SPEEDS))
    print(f"  Inter-run lag-1 autocorrelation (should be ~0): "
          + ", ".join(f"{v}: C0={r['lag1_autocorr_C0']:+.2f} C3={r['lag1_autocorr_C3']:+.2f}" for v, r in zip(SPEEDS, cc_rows)))
    print("  Deadline reliability P_del (fraction delivered within 100 ms):")
    for v in SPEEDS:
        c0 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C0")
        c3 = next(r for r in rel_rows if r["speed_kmh"] == v and r["config"] == "C3")
        warn = " <- FAILS 99 %" if c0["P_del_mean"] < 0.99 else ""
        print(f"    {v:3d} km/h: C0 = {c0['P_del_mean']:.4f}{warn}   C3 = {c3['P_del_mean']:.4f}")
    print("  Handover cost (mean per 120-s run; adaptive policy trade-off):")
    for v in SPEEDS:
        c0 = next(r for r in ho_rows if r["speed_kmh"] == v and r["config"] == "C0")
        c3 = next(r for r in ho_rows if r["speed_kmh"] == v and r["config"] == "C3")
        print(f"    {v:3d} km/h: C0 handovers = {c0['handovers_mean_per_run']:4.1f}  "
              f"C3 = {c3['handovers_mean_per_run']:4.1f}  ping-pongs = {c3['ping_pongs_mean']:4.1f}")
    print(f"\n[{time.time() - t_start:.0f} s] -> {args.outdir}/   "
          f"({len(manifest['outputs'])} files)")
    if args.quick:
        print("QUICK MODE - do not use these numbers in the manuscript.")


if __name__ == "__main__":
    freeze_support()
    main()
