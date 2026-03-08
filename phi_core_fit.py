#!/usr/bin/env python3
"""
phi_core_fit.py
===============
Fit three dark halo profiles to SPARC rotation curves with baryonic
subtraction, reproducing the results of:

    Berry, J. (2026). "The Golden-Ratio Dark Halo: Testing a
    Geometry-Motivated Rotation Curve Profile Against the Full
    SPARC Catalogue." MNRAS (submitted).

PROFILES
--------
phi-core :  v(r) = v0 * r / sqrt(phi^2 * rc^2 + r^2)
            phi = (1+sqrt(5))/2 = 1.6180...  FIXED, not fitted
            Free parameters: v0, rc

ISO      :  pseudo-isothermal sphere
            v(r) = v0 * sqrt(1 - (rc/r) * arctan(r/rc))
            Free parameters: v0, rc

NFW      :  Navarro-Frenk-White (two variants)
            (a) c = 10 fixed  (2 free parameters: v200, rs)
            (b) c fitted freely  (3 free parameters: v200, rs, c)

BARYONIC SUBTRACTION
--------------------
Following Li et al. (2020) and standard SPARC practice:

    Vhalo^2 = Vobs^2 - Y* * Vdisk^2 - Y* * Vbul^2 - Vgas^2

where Y* = 0.5 M_sun/L_sun at 3.6 um (fixed; stellar population
theory, Schombert et al. 2014). Signs are preserved for signed
SPARC velocities.

INPUT DATA
----------
SPARC rotmod files (*_rotmod.dat), available at:
    http://astroweb.cwru.edu/SPARC/

Columns: r [kpc], Vobs, errV, Vgas, Vdisk, Vbul [km/s]

USAGE
-----
    python phi_core_fit.py --data_dir /path/to/rotmod_files/
                           --ystar 0.5
                           --output results.csv

    # Y* robustness sweep:
    python phi_core_fit.py --data_dir ./rotmod/ --ystar 0.3 0.4 0.5 0.6 0.7

DEPENDENCIES
------------
    numpy, scipy, matplotlib  (standard scientific Python stack)
    Python >= 3.8
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import wilcoxon, binom_test
import glob, os, csv, argparse, warnings
warnings.filterwarnings('ignore')

# ── Constants ─────────────────────────────────────────────────────────────

PHI = (1.0 + np.sqrt(5.0)) / 2.0   # golden ratio — NEVER FITTED


# ── Profile functions ─────────────────────────────────────────────────────

def v_phi_core(r, v0, rc):
    """
    phi-core dark halo profile.
    phi is FIXED at (1+sqrt(5))/2. It is NOT a free parameter.
    The profile shape in the transition region r ~ phi*rc is a
    prediction of the geometric constraint, not a fit.
    """
    return v0 * r / np.sqrt(PHI**2 * rc**2 + r**2)


def v_iso(r, v0, rc):
    """Pseudo-isothermal sphere."""
    x = np.clip(r / rc, 1e-10, None)
    return v0 * np.sqrt(np.clip(1.0 - (1.0 / x) * np.arctan(x), 0.0, None))


def v_nfw(r, v200, rs, c=10.0):
    """
    NFW profile. c=10 is the fixed-concentration variant (2 free params).
    Pass c as a free parameter for the 3-param variant.
    """
    x    = np.clip(r / rs, 1e-10, None)
    c    = np.clip(c, 1.0, 100.0)
    norm = np.log(1.0 + c) - c / (1.0 + c)
    if norm <= 0:
        return np.zeros_like(r)
    val  = np.clip((np.log(1.0 + x) - x / (1.0 + x)) / (x * norm), 0.0, None)
    return v200 * np.sqrt(val)


def v_nfw_free(r, v200, rs, c):
    """NFW with freely fitted concentration (3 free parameters)."""
    return v_nfw(r, v200, rs, c=c)


# ── Fitting utilities ─────────────────────────────────────────────────────

def reduced_chi2(v_obs, v_err, v_mod, n_par):
    """Reduced chi-squared. Returns nan if dof <= 0."""
    dof = len(v_obs) - n_par
    if dof <= 0:
        return np.nan
    return float(np.sum(((v_obs - v_mod) / v_err) ** 2) / dof)


def fit_profile(func, r, v, e, p0, bounds, n_par):
    """Fit profile, return (params, chi2). Returns (None, nan) on failure."""
    try:
        popt, _ = curve_fit(
            func, r, v, p0=p0, bounds=bounds,
            sigma=e, absolute_sigma=True, maxfev=20000
        )
        chi2 = reduced_chi2(v, e, func(r, *popt), n_par)
        return popt, chi2
    except Exception:
        return None, np.nan


def initial_guess_rc(r, v):
    """Estimate core radius from 70% of peak velocity."""
    vmax = v.max()
    idx  = np.argmin(np.abs(v - 0.7 * vmax))
    return max(r[idx], 0.1)


# ── Baryonic subtraction ──────────────────────────────────────────────────

def compute_halo_velocity(vobs, vgas, vdisk, vbul, ystar):
    """
    Subtract baryonic contribution from observed rotation curve.

    Uses signed velocities as stored in SPARC rotmod files:
        Vhalo^2 = Vobs^2 - Y* * sign(Vdisk)*Vdisk^2
                         - Y* * sign(Vbul)*Vbul^2
                         - sign(Vgas)*Vgas^2

    Returns Vhalo^2 (may contain negative values at poorly-measured points).
    """
    vbar2 = (ystar * np.sign(vdisk) * vdisk**2 +
             ystar * np.sign(vbul)  * vbul**2  +
             np.sign(vgas)          * vgas**2)
    return vobs**2 - vbar2


# ── Per-galaxy fitting ────────────────────────────────────────────────────

def fit_galaxy(fpath, ystar=0.5):
    """
    Load one SPARC rotmod file, subtract baryons, fit all profiles.

    Returns dict with chi2 values and half-power test results,
    or None if insufficient data.
    """
    gname = os.path.basename(fpath).replace('_rotmod.dat', '')

    try:
        d = np.loadtxt(fpath, comments='#')
    except Exception:
        return None

    if d.ndim < 2 or d.shape[1] < 6 or d.shape[0] < 4:
        return None

    r, vobs, errv = d[:, 0], d[:, 1], d[:, 2]
    vgas, vdisk, vbul = d[:, 3], d[:, 4], d[:, 5]

    # Baryonic subtraction
    vhalo2 = compute_halo_velocity(vobs, vgas, vdisk, vbul, ystar)

    # Quality mask: positive halo residual, valid errors, positive radius
    mask = ((errv > 0) & np.isfinite(vobs) & np.isfinite(errv) &
            (r > 0) & (vhalo2 > 0))

    if mask.sum() < 4:
        return None

    r_m = r[mask]
    v_m = np.sqrt(vhalo2[mask])
    e_m = errv[mask]

    vmax = v_m.max()
    rc0  = initial_guess_rc(r_m, v_m)
    bds2 = ([1.0, 0.01], [2000.0, 300.0])

    # Fit phi-core (2 params)
    popt_phi, chi2_phi = fit_profile(
        v_phi_core, r_m, v_m, e_m, [vmax, rc0], bds2, n_par=2
    )

    # Fit ISO (2 params)
    _, chi2_iso = fit_profile(
        v_iso, r_m, v_m, e_m, [vmax, rc0], bds2, n_par=2
    )

    # Fit NFW c=10 fixed (2 params)
    _, chi2_nfw_c10 = fit_profile(
        v_nfw, r_m, v_m, e_m,
        [vmax, rc0 * 3.0], ([1.0, 0.01], [2000.0, 2000.0]), n_par=2
    )

    # Fit NFW c free (3 params)
    _, chi2_nfw_cfree = fit_profile(
        v_nfw_free, r_m, v_m, e_m,
        [vmax, rc0 * 3.0, 10.0],
        ([1.0, 0.01, 1.0], [2000.0, 2000.0, 100.0]), n_par=3
    )

    result = dict(
        galaxy=gname,
        n_halo_pts=int(mask.sum()),
        chi2_phi_core=chi2_phi,
        chi2_ISO=chi2_iso,
        chi2_NFW_c10=chi2_nfw_c10,
        chi2_NFW_cfree=chi2_nfw_cfree,
        hp_ratio=np.nan,
        hp_pull=np.nan,
    )

    # Half-power consistency check
    if popt_phi is not None:
        v0f, rcf = popt_phi
        r_test   = PHI * rcf
        if r_m[0] < r_test < r_m[-1]:
            v_at_test = np.interp(r_test, r_m, v_m)
            e_at_test = np.interp(r_test, r_m, e_m)
            result['hp_ratio'] = float(v_at_test / v0f)
            result['hp_pull']  = float((v_at_test - v0f / np.sqrt(2.0)) / e_at_test)

    return result


# ── Summary statistics ────────────────────────────────────────────────────

def summarise(results):
    """Print summary statistics for a list of per-galaxy result dicts."""
    def arr(key):
        return np.array([r[key] for r in results], dtype=float)

    chi2_phi  = arr('chi2_phi_core')
    chi2_iso  = arr('chi2_ISO')
    chi2_nfwf = arr('chi2_NFW_cfree')
    hp_ratio  = arr('hp_ratio')
    hp_pull   = arr('hp_pull')

    ok_pi = np.isfinite(chi2_phi) & np.isfinite(chi2_iso)
    ok_pn = np.isfinite(chi2_phi) & np.isfinite(chi2_nfwf)
    n_valid = ok_pi.sum()

    n_phi_iso = int((chi2_phi[ok_pi] < chi2_iso[ok_pi]).sum())
    n_phi_nfw = int((chi2_phi[ok_pn] < chi2_nfwf[ok_pn]).sum())

    # Sign test
    from scipy.stats import binom_test as bt
    p_sign = bt(n_phi_iso, n_valid, 0.5, alternative='greater')

    # Wilcoxon signed-rank
    diffs = (chi2_phi - chi2_iso)[ok_pi]
    try:
        _, p_wilcox = wilcoxon(diffs, alternative='less')
    except Exception:
        p_wilcox = np.nan

    hp_ok = np.isfinite(hp_pull)
    n_hp  = hp_ok.sum()
    n_2sigma = int((np.abs(hp_pull[hp_ok]) < 2.0).sum())
    pred  = 1.0 / np.sqrt(2.0)

    print(f"\nResults summary  (N={len(results)} fitted)")
    print(f"{'Model':<16}  {'Mean chi2':>9}  {'Median chi2':>11}")
    print(f"{'phi-core':<16}  {np.nanmean(chi2_phi):9.3f}  {np.nanmedian(chi2_phi):11.3f}")
    print(f"{'ISO':<16}  {np.nanmean(chi2_iso):9.3f}  {np.nanmedian(chi2_iso):11.3f}")
    print(f"{'NFW (c free)':<16}  {np.nanmean(chi2_nfwf):9.3f}  {np.nanmedian(chi2_nfwf):11.3f}")
    print(f"\nphi-core beats ISO:       {n_phi_iso}/{n_valid}")
    print(f"phi-core beats NFW(free): {n_phi_nfw}/{n_valid}")
    print(f"Sign test p-value:        {p_sign:.4f}")
    print(f"Wilcoxon p-value:         {p_wilcox:.4f}")
    print(f"\nHalf-power check (N={n_hp}):")
    print(f"  Predicted:  {pred:.4f}")
    if n_hp > 0:
        pull_mean = (hp_ratio[hp_ok].mean() - pred) / (hp_ratio[hp_ok].std() / np.sqrt(n_hp))
        print(f"  Obs. mean:  {hp_ratio[hp_ok].mean():.4f} +/- {hp_ratio[hp_ok].std():.4f}")
        print(f"  Mean pull:  {pull_mean:+.3f} sigma")
        print(f"  |pull_i|<2: {n_2sigma}/{n_hp}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Fit phi-core, ISO, NFW profiles to SPARC rotation curves.'
    )
    parser.add_argument('--data_dir', default='.',
                        help='Directory containing *_rotmod.dat files')
    parser.add_argument('--ystar', type=float, nargs='+', default=[0.5],
                        help='Stellar M/L value(s) at 3.6 um (default: 0.5)')
    parser.add_argument('--output', default='phi_core_results.csv',
                        help='Output CSV filename')
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir, '*_rotmod.dat')))
    if not files:
        raise FileNotFoundError(f"No *_rotmod.dat files found in {args.data_dir}")

    print(f"Found {len(files)} SPARC rotmod files.")
    print(f"PHI = {PHI:.10f}  (fixed; not fitted)")

    all_results = {}

    for ystar in args.ystar:
        print(f"\n{'='*60}")
        print(f"Y* = {ystar}")
        print(f"{'='*60}")

        results = []
        skipped = 0
        for fpath in files:
            r = fit_galaxy(fpath, ystar=ystar)
            if r is None:
                skipped += 1
            else:
                results.append(r)

        print(f"Fitted: {len(results)}  Skipped: {skipped}")
        summarise(results)
        all_results[ystar] = results

    # Write CSV for the primary Y*=0.5 run
    primary = all_results[0.5] if 0.5 in all_results else all_results[args.ystar[0]]
    fieldnames = ['galaxy', 'n_halo_pts',
                  'chi2_phi_core', 'chi2_ISO', 'chi2_NFW_c10', 'chi2_NFW_cfree',
                  'hp_ratio', 'hp_pull']
    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in primary:
            w.writerow({k: ('' if (isinstance(v, float) and np.isnan(v)) else
                           round(v, 6) if isinstance(v, float) else v)
                       for k, v in row.items()})

    print(f"\nResults written to: {args.output}")


if __name__ == '__main__':
    main()
