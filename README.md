# Rotation-Curve-Profile-Against-the-Full-SPARC-Catalogue
Data and code for 10.5281/zenodo.18912199
# φ-Core Dark Halo Profile — SPARC Fitting Code and Results

**Associated paper:**
> Berry, J. (2026). *The Golden-Ratio Dark Halo: Testing a Geometry-Motivated
> Rotation Curve Profile Against the Full SPARC Catalogue.*
> MNRAS (submitted). [10.5281/zenodo.18912199]

---

## What this deposit contains

```
phi_core_sparc/
├── README.md                       ← this file
├── code/
│   └── phi_core_fit.py             ← main fitting script (all three profiles)
└── data/
    └── phi_core_chi2_results.csv   ← per-galaxy χ² results, 172 galaxies
```

---

## The phi-core profile

The φ-core profile is:

```
v(r) = v₀ · r / sqrt(φ² · rc² + r²)
```

where **φ = (1+√5)/2 = 1.6180...** is the golden ratio, fixed by the
theoretical framework. It is **not a free parameter**. Each galaxy has
two free parameters: `v0` (asymptotic velocity scale) and `rc` (core radius).

This is compared against:
- **ISO**: pseudo-isothermal sphere (2 free parameters)
- **NFW c=10**: Navarro-Frenk-White with fixed concentration (2 free parameters)
- **NFW c free**: NFW with freely fitted concentration (3 free parameters)

---

## Results file: `phi_core_chi2_results.csv`

| Column | Description |
|--------|-------------|
| `galaxy` | SPARC galaxy name |
| `n_halo_pts` | Number of valid dark halo data points after baryonic subtraction |
| `chi2_phi_core` | Reduced χ² for φ-core profile fit |
| `chi2_ISO` | Reduced χ² for ISO profile fit |
| `chi2_NFW_c10` | Reduced χ² for NFW (c=10 fixed, 2 params) |
| `chi2_NFW_cfree` | Reduced χ² for NFW (c free, 3 params) |
| `hp_ratio` | v(φ·rc)/v₀ from best-fit parameters (half-power consistency check; blank if φ·rc outside radial coverage) |
| `hp_pull` | Per-galaxy pull: (v_obs − v_pred) / σ_v at r = φ·rc |

**Baryonic subtraction:** Y★ = 0.5 M☉/L☉ at 3.6 μm (fixed; stellar population
theory, Schombert et al. 2014). Signed SPARC velocities are preserved.

**N=172 fitted; 3 galaxies excluded** (CamB, NGC4389, UGC02455) for
insufficient dark halo points after baryonic subtraction (minimum 4 required).

---

## Reproducing the results

**Requirements:** Python ≥ 3.8, numpy, scipy, matplotlib

**Input data:** SPARC rotation curve files (`*_rotmod.dat`) from
http://astroweb.cwru.edu/SPARC/

```bash
# Install dependencies
pip install numpy scipy matplotlib

# Download SPARC data (175 galaxies)
# from http://astroweb.cwru.edu/SPARC/ → "Rotation Curves"

# Run primary fit (Y★ = 0.5)
python code/phi_core_fit.py --data_dir /path/to/rotmod/ --output my_results.csv

# Y★ robustness sweep (reproduces Table 3)
python code/phi_core_fit.py --data_dir /path/to/rotmod/ --ystar 0.3 0.4 0.5 0.6 0.7
```

The output CSV should match `data/phi_core_chi2_results.csv` to within
floating-point rounding.

---

## Key results summary (Y★ = 0.5)

| Model | Free params | Median χ² | Mean χ² |
|-------|------------|-----------|---------|
| φ-core | 2 (φ fixed) | 1.221 | 4.725 |
| ISO | 2 | 1.293 | 4.792 |
| NFW (c free) | 3 | 2.439 | 8.580 |

φ-core beats ISO on **108/172 galaxies** (sign test p = 0.0005).
φ-core beats NFW (c free, 3 params) on **131/172 galaxies**.
Y★ robustness: φ-core advantage holds at all tested values Y★ ∈ {0.3, 0.4, 0.5, 0.6, 0.7}.

---

## Citation

If you use this code or data, please cite the associated paper:

```
Berry, J. (2026). The Golden-Ratio Dark Halo: Testing a Geometry-Motivated
Rotation Curve Profile Against the Full SPARC Catalogue.
MNRAS (submitted). doi:10.5281/zenodo.[THIS_RECORD]
```

Please also cite the SPARC catalogue:

```
Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016).
SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry
and Accurate Rotation Curves. AJ, 152, 157.
doi:10.3847/0004-6256/152/6/157
```

---

## Licence

Code: MIT License
Data: CC BY 4.0

© 2026 Jennifer Berry, Fibonacci Research Institute
