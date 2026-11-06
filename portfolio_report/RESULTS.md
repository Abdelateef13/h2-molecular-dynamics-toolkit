# H₂ molecular dynamics — campaign overview

This is a comparative *simulation* summary, not an experimental validation of self-diffusion.
Confidence intervals represent **between-replica statistical uncertainty** only.

| H₂ molecules | T (K) | P (MPa) | Replicas | Density (kg/m³) | EOS error | D (10⁻⁶ m²/s) | Warnings |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 256 | 323.15 | 10 | 5 | 7.1087 ± 0.0289 | +0.122% | 1.7417 ± 0.1472 | 0 |
| 500 | 323.15 | 10 | 3 | 7.0885 ± 0.0502 | -0.163% | 1.7386 ± 0.2139 | 0 |
| 864 | 323.15 | 10 | 3 | 7.0955 ± 0.0395 | -0.064% | 1.8264 ± 0.2263 | 0 |

## Reproducibility and interpretation

- Density is from the equilibrated portion of the NPT stage; H₂ self-diffusion is from molecular COM positions during NVT production.
- Diffusion CI does not include force-field bias, finite-size corrections, or fitting-window systematic error.
- Compare sizes only at the same thermodynamic state and protocol; overlapping CIs do not prove size independence.
- Model: two-site LJ site parameters from Adam et al., *International Journal of Hydrogen Energy* (2024), DOI:10.1016/j.ijhydene.2024.06.319. Rigid H–H bond is a separate implementation assumption.
- Source files and exact per-run diagnostics are in `summary.json`, `replicas.csv` and each analysis `REPORT.md`.
