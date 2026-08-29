# LunaNMR v1.0

> ⚠️ Experimental and under active development. Expect instability, breaking API changes, and incomplete features. No support or warranty. Use at your own risk.

**Advanced NMR Peak Analysis and Integration Suite**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)](https://github.com/mas-gu/lunaNMR)

PySide6/Qt6 suite for NMR peak detection, Voigt profile fitting, and integration — extended with binding-affinity (Kd / titration) fitting and protein dynamics (T₁/T₂, methyl T₂, spectral density, model-free). Focus on automation: from automatic peak detection and overlap-aware multi-peak fitting through end-to-end batch processing of multi-spectrum datasets.

🌐 **[Landing page →](https://mas-gu.github.io/lunaNMR/)**

---

## Screenshots

<div align="center">
 LunaNMR interface <br>
  <img width="776" height="494" alt="luna_interface" src="https://github.com/user-attachments/assets/21e991e3-e622-4e60-8d32-c3b6ccd6d365" />
</div>

<div align="center">
 Series integration <br>
<img width="405" height="494" alt="serie_integration" src="https://github.com/user-attachments/assets/1b71dcea-e2a6-46d8-ad68-e29f8049ef71" />
</div>

<div align="center">
 Individual-peak fitting vs volume integration <br>
<img width="776" height="494" alt="individualpeak" src="https://github.com/user-attachments/assets/899ef2f6-3f1c-42aa-87c4-3b4b3702de66" />
</div>

<div align="center">
 Multi-peak fitting vs volume integration <br>
<img width="776" height="494" alt="multipeaks" src="https://github.com/user-attachments/assets/163800b4-081c-42e3-988f-ad015242c306" />
</div>

---

## Features

- **Voigt fitting**: Gaussian + Lorentzian lineshapes, ArPLS baseline correction, multi-peak deconvolution with simultaneous fitting of overlapping peaks.
- **2D simultaneous multi-peak fitting (PS2D)**: closely-spaced peaks auto-routed to 2D Levenberg-Marquardt deconvolution.
- **Peak detection**: network-based clustering with prominence analysis, no manual seeding required.
- **Quality assessment**: R², uncertainties (bootstrap / covariance), color-coded markers.
- **Series workflows**: relaxation (T₁, T₂, hetNOE), titration, and multi-spectrum (temperature/pH) processing with automated peak tracking.
- **Kd / titration**: CSP (quadratic 1:1 isotherm) and intensity-decay binding fits.
- **Dynamics (DynamiXs)**: T₁/T₂, methyl T₂ (bi-exponential), T₁ρ, spectral density, model-free.
- **Parallel processing**: two-pass cluster-based, ~2.7× speedup.
- **ML/statistics**: dual-path fit-parameter prediction with statistical fallback (optional CNN classifier).
- **Headless CLI**: `python -m lunaNMR` — every analysis is scriptable with no display. Subcommand reference: [docs/CLI.md](docs/CLI.md); machine contract and gotchas for agents: [docs/CLI_AGENT.md](docs/CLI_AGENT.md).

**File formats**: Bruker TopSpin (`.2ii`, `.2rr`), Varian/Agilent (`.fid`, `.ft`), NMRPipe (`.ft`, `.pipe`), SPARKY (`.ucsf`). Export: CSV, JSON, PNG, PDF, SVG, EPS.

---

## Quick Start

```bash
git clone https://github.com/mas-gu/lunaNMR.git
cd lunaNMR/lunaNMR_v1o0

pip install -r requirements.txt                          # dependencies
python3 lunaNMR/validation/verify_installation.py        # verify
python3 launch_lunaNMR.py                                 # GUI launcher (LunaNMR / DynamiXs selector)
python -m lunaNMR --help                                  # headless CLI
```

Requires **numba** (2D multi-peak fitting) and **PySide6** (GUI). `torch`/`torchvision` optional (CNN classifier).

---

## Package Structure

```
lunaNMR_v1o0/
├── launch_lunaNMR.py          # GUI launcher (tkinter selector → Qt app)
├── lunaNMR/
│   ├── cli.py / __main__.py   # headless CLI (python -m lunaNMR)
│   ├── core/                  # peak picker, Voigt fitter, integrator, PS2D, parallel processor
│   ├── gui/                   # PySide6/Qt6 interface
│   ├── processors/            # single/multi-spectrum workflows
│   ├── batch_processing/      # CLI batch automation
│   ├── ml/                    # ML/statistics parameter prediction
│   ├── integrators/           # specialized integrators
│   ├── utils/                 # config, project bundles, output routing
│   └── validation/            # installation verification
├── modules/dynamiXs_v2o0/     # DynamiXs (dynamics) + Kd / titration
└── docs/                      # CLI.md, ARCHITECTURE.md, ALGORITHMS.md, guide
```

---

## Scientific Background

- **Voigt function**: convolution of Gaussian (instrumental) and Lorentzian (natural linewidth) components; multi-peak deconvolution with constraints.
- **Peak detection**: network-based clustering with prominence analysis.
- **Optimization**: Levenberg-Marquardt with parameter normalization.
- **Uncertainty**: bootstrap and covariance-based error estimation.

---

## Contributing

```bash
python3 -m venv dev_env && source dev_env/bin/activate   # or dev_env\Scripts\activate on Windows
pip install -r requirements.txt
make test        # all three test roots, headless (QT_QPA_PLATFORM=offscreen)
```

Fork → feature branch → add tests → PR. PEP 8, NumPy-style docstrings.

---

## Citation

```bibtex
@software{lunaNMR2025,
  title = {LunaNMR: Advanced NMR Peak Analysis and Integration Suite},
  author = {LunaNMR Contributors},
  year = {2025},
  version = {1.0},
  url = {https://github.com/mas-gu/lunaNMR},
  license = {MIT}
}
```

## License

MIT — see [LICENSE](LICENSE).
</content>
</invoke>
