- 📫 How to reach me : khairy.guillem@gmail.com | guillem.khairy@uni.lu

[VoFBoiling — VoF boiling source model for OpenFOAM]

Volume-of-Fluid (VoF) phase-change model (Lee-type) for OpenFOAM. Other models currently being added.
Adds **mass transfer** in the α-equation and **latent heat** in the energy equation through the `fvModels` framework. Compatible with OpenFOAM **v12–v13 (Foundation)**.

## Features

- Boiling source (Lee, Schrage, etc):
  * bounded α source with interface masking
  * latent heat source $Q=\dot m\,L$ in the energy equation
- Works with modular solvers `compressibleVoF`. incompressibleVoF, incompressibleMultiPhaseVoF and compressibleMultiPhaseVoF in progress
- Run-time selection via `constant/fvModels`
- Example cases (pool boiling, heated wall)

## Requirements

* OpenFOAM Foundation v12 or v13
* C++14 compiler (GCC recommended)
* Linux (tested on Ubuntu/Debian-based)

---

## How to cite

If you use this code in scientific work, please cite:
```
@software{VoFBoiling,
  author  = {Guillem KHAIRY},
  title   = {VoFBoiling: modular boiling fvModel for VoF in OpenFOAM Foundation},
  year    = {2025},
  url     = {https://github.com/<your-user>/VoFBoiling},
  version = {v0.1}
}
```
You can also add a `CITATION.cff` file to the repository so GitHub shows a “Cite this repository” button.

## License

This project is licensed under the **GNU General Public License v3.0** (GPL-3.0).
You are free to use, modify, and distribute it under the terms of the GPL-3.0.

* SPDX: `GPL-3.0-or-later`
* See [LICENSE](./LICENSE) for full text.

> Note: OpenFOAM itself is GPL-3.0; this project remains GPL-compatible.

---

## Acknowledgments

* Built on OpenFOAM Foundation (GPL-3.0).
* Thanks to contributors and the OpenFOAM community.

---

## Contact

* Author: Guillem KHAIRY - guillem.khairy@uni.lu | khairy.guillem@gmail.com
* Repo: [https://github.com/](https://github.com/RadioactiveMoose/RadioactiveMoose.git)

