# LTE-R Train-Control Communication — Simulation code

Supporting code for the manuscript
**"Latency and Reliability of LTE-R Train Control Communication on High-Speed Railway Corridors"**
(submitted to MDPI *Network*).

## Status

**This code is a re-implementation prepared in response to the reviewer comments.**
It is NOT a restored copy of the original MATLAB simulation (the original records were
lost after a computer-virus incident; the reviewers were informed by e-mail). Numbers
produced here differ from the earlier MATLAB figures and replace them.

Each figure of the manuscript has its own **self-contained Python script** inside its
own folder. No shared module is imported across folders, so every figure can be
reproduced independently by a reviewer.

```
.
├── figure3/
│   ├── figure3.py                 # standalone: engine + Fig 3 production in one file
│   ├── figure3_output/            # outputs of the archived full run (for reference)
│   └── README.md                  # how to reproduce Figure 3
├── figure4/
│   ├── figure4.py                 # standalone: engine + Fig 4 production in one file
│   ├── figure4_output/            # outputs of the archived full run (for reference)
│   └── README.md                  # how to reproduce Figure 4
├── figure5/
│   ├── figure5.py                 # standalone  (to be added)
│   ├── figure5_output/
│   └── README.md
├── requirements.txt               # Python package versions used in the archived runs
├── LICENSE
└── README.md                      # this file
```

## How to reproduce (any figure)

```bash
cd figure3                     # or figure4 / figure5
pip install -r ../requirements.txt
python figure3.py              # full study; no --quick or --sensitivity flags
```

Each figure folder's README describes the exact CLI options and lists every output file.

## What the output guarantees

Each run writes a `manifest.json` containing:
- base random seed and the set of per-replication seeds derived from it,
- every parameter with its source tag (**T4** = manuscript Table 4, **REF** = cited
  reference, **ASSUMED** = modelling assumption not in Table 4),
- the SHA-256 of the script that produced the outputs and of the Python packages used,
- the engine self-test report (M/M/1, M/D/1, accumulating priority Eq. 15a, HARQ
  truncated-geometric; the full set must pass or the script aborts).

With the same base seed and the same package versions the run is bit-reproducible.

## Known limitations (also documented in the manuscript, Section 5)

- The radio layer is a synthetic WINNER II D2a / log-normal-shadowing / Rician-fading
  channel. It is **not** ITU-R P.1812 and the figures are **not** field-validated.
- Block errors come from an information-outage link abstraction with chase-combining
  HARQ, not from PHY-level BLER curves.
- Parameters tagged ASSUMED are not measured values; their effect on the main
  conclusion is reported in the sensitivity sweep of Figure 3.
- The 99.9 % reliability target shown on Figure 4 is a generic SIL-3-style benchmark;
  the FRMCS 99.9999 % target is NOT met by any simulated configuration at any load.

## License

Released under the MIT License. See `LICENSE`.
