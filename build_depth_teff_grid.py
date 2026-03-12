#!/usr/bin/env python3
"""
build_depth_teff_grid.py
========================
Build a grid of intrinsic line depths as a function of effective temperature
using MARCS model atmospheres through iSpec.

Scientific motivation
---------------------
Each spectral line is sensitive to a different depth in the photosphere and
therefore probes a different temperature.  By synthesising a grid of spectra
over a range of Teff (with no rotational or macroturbulent broadening) we
obtain, for each transition, a calibration curve:

    depth(Teff)   —  fit with a 2nd-order polynomial

This calibration is then used downstream to:
  1. Estimate a "local Teff" per line from its observed (broadening-corrected) depth.
  2. Look for systematic trends of vsini / vmacro vs that local Teff as a
     proxy for photospheric formation depth.

Algorithm
---------
For each Teff in the grid:
  1. Interpolate MARCS atmosphere at (Teff, logg, [M/H], alpha).
  2. Synthesise the spectrum with:
       vsini = 0,  vmacro = 0  → pure intrinsic depths
       vmic  = empirical estimate from Teff / logg / [M/H]
       R     = 300 000  (resolve all thermal line widths)
  3. For every line in the target list, measure
       depth = 1 - min(flux in ±window nm around line centre)

Outputs
-------
  <output_stem>_raw.csv          one row per line, columns: wave_nm, element,
                                  loggf, EP, depth_T3500, depth_T3750, ...
  <output_stem>_poly.csv         2nd-order polynomial coefficients per line
  <output_stem>_raw.pkl          pandas DataFrame (binary, for fast loading)

Usage
-----
  python build_depth_teff_grid.py \\
        --ispec-dir /path/to/iSpec \\
        --teff-min 3500 --teff-max 5500 --teff-step 250 \\
        --logg 4.5 --mh 0.0 \\
        --wave-min 500 --wave-max 700 \\
        --linelist VALD.300_1100nm \\
        --depth-threshold 0.02 \\
        --output depth_grid

The script can also be imported and the main functions called directly.
"""

import argparse
import logging
import sys
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# iSpec bootstrap
# ---------------------------------------------------------------------------

def _bootstrap_ispec(ispec_dir: Path) -> None:
    """Add iSpec to sys.path if it is not already importable."""
    ispec_dir = str(ispec_dir)
    if ispec_dir not in sys.path:
        sys.path.insert(0, ispec_dir)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_atmosphere_grid(ispec_dir: Path, atmosphere_name: str = "MARCS.GES"):
    """
    Load a packed model-atmosphere grid from iSpec's input directory.

    Parameters
    ----------
    ispec_dir : Path
    atmosphere_name : str
        Sub-directory name under input/atmospheres/.
        Typical choices: "MARCS.GES", "MARCS", "MARCS.APOGEE",
        "ATLAS9.Castelli", "ATLAS9.Kurucz".

    Returns
    -------
    modeled_layers_pack : iSpec atmosphere object
    is_atlas : bool   True when the grid is an ATLAS9 model.
    """
    import ispec  # noqa: F401  (only importable after bootstrap)

    model_dir = ispec_dir / "input" / "atmospheres" / atmosphere_name
    if not model_dir.exists():
        # Try plain MARCS as fallback
        fallback = ispec_dir / "input" / "atmospheres" / "MARCS"
        if fallback.exists():
            log.warning(
                "%s not found; falling back to %s", atmosphere_name, fallback
            )
            model_dir = fallback
        else:
            raise FileNotFoundError(
                f"Atmosphere grid not found.\n"
                f"Checked: {model_dir}\n"
                f"Install the iSpec input data from https://www.blancocuaresma.com/s/iSpec"
            )

    log.info("Loading atmosphere grid from %s …", model_dir)
    modeled_layers_pack = ispec.load_modeled_layers_pack(str(model_dir))
    is_atlas = "ATLAS" in atmosphere_name.upper()
    return modeled_layers_pack, is_atlas


def load_linelist(ispec_dir: Path, linelist_name: str,
                  wave_min: float, wave_max: float,
                  depth_threshold: float = 0.01):
    """
    Load an atomic line list from iSpec's input directory.

    Parameters
    ----------
    ispec_dir : Path
    linelist_name : str
        Sub-directory name under input/linelists/transitions/.
        E.g. "VALD.300_1100nm", "GESv6_atom_hfs_iso.420_920nm".
    wave_min, wave_max : float
        Wavelength range in nm.
    depth_threshold : float
        Minimum theoretical depth (solar) to include a line.

    Returns
    -------
    numpy recarray  (iSpec atomic linelist format)
    """
    import ispec

    candidates = [
        ispec_dir / "input" / "linelists" / "transitions" / linelist_name / "atomic_lines.tsv",
        ispec_dir / "input" / "linelists" / "transitions" / linelist_name / "linelist.tsv",
    ]
    linelist_file = None
    for c in candidates:
        if c.exists():
            linelist_file = c
            break
    if linelist_file is None:
        raise FileNotFoundError(
            f"Linelist '{linelist_name}' not found under "
            f"{ispec_dir / 'input' / 'linelists' / 'transitions'}"
        )

    log.info("Loading linelist from %s …", linelist_file.name)
    linelist = ispec.read_atomic_linelist(
        str(linelist_file), wave_base=wave_min, wave_top=wave_max
    )
    n_before = len(linelist)
    linelist = linelist[linelist["theoretical_depth"] >= depth_threshold]
    log.info(
        "  %d lines in range, %d kept after depth threshold %.3f",
        n_before, len(linelist), depth_threshold,
    )
    return linelist


def load_support_data(ispec_dir: Path, is_atlas: bool):
    """
    Load solar abundances and isotope data.

    Returns
    -------
    solar_abundances, isotopes  (iSpec format)
    """
    import ispec

    if is_atlas:
        abund_file = ispec_dir / "input" / "abundances" / "Grevesse.1998" / "stdatom.dat"
    else:
        abund_file = ispec_dir / "input" / "abundances" / "Grevesse.2007" / "stdatom.dat"

    if not abund_file.exists():
        # try any abundance file we can find
        for subdir in (ispec_dir / "input" / "abundances").iterdir():
            candidate = subdir / "stdatom.dat"
            if candidate.exists():
                abund_file = candidate
                log.warning("Preferred abundances not found; using %s", abund_file)
                break
        else:
            raise FileNotFoundError(
                f"Solar abundances not found under {ispec_dir / 'input' / 'abundances'}"
            )

    isotope_file = ispec_dir / "input" / "isotopes" / "SPECTRUM.lst"
    if not isotope_file.exists():
        raise FileNotFoundError(f"Isotope file not found: {isotope_file}")

    solar_abundances = ispec.read_solar_abundances(str(abund_file))
    isotopes = ispec.read_isotope_data(str(isotope_file))
    log.info("Loaded abundances from %s", abund_file.parent.name)
    return solar_abundances, isotopes


# ---------------------------------------------------------------------------
# Synthesis helpers
# ---------------------------------------------------------------------------

def _make_segments_for_linelist(linelist, window_nm: float = 0.4):
    """
    Build a minimal set of synthesis segments covering every line in the
    linelist, merging overlapping windows.

    Parameters
    ----------
    linelist   : iSpec atomic linelist recarray
    window_nm  : half-width of the window around each line centre (nm)

    Returns
    -------
    numpy recarray with dtype [('wave_base', float), ('wave_top', float)]
    """
    if len(linelist) == 0:
        return np.recarray(0, dtype=[("wave_base", float), ("wave_top", float)])

    waves = np.sort(linelist["wave_nm"])
    bases = waves - window_nm
    tops  = waves + window_nm

    # Merge overlapping intervals
    merged = []
    b, t = bases[0], tops[0]
    for bi, ti in zip(bases[1:], tops[1:]):
        if bi <= t:
            t = max(t, ti)
        else:
            merged.append((b, t))
            b, t = bi, ti
    merged.append((b, t))

    segs = np.recarray(len(merged), dtype=[("wave_base", float), ("wave_top", float)])
    for i, (b, t) in enumerate(merged):
        segs["wave_base"][i] = b
        segs["wave_top"][i]  = t
    return segs


def synthesise_at_teff(
    teff: float,
    logg: float,
    MH: float,
    modeled_layers_pack,
    linelist,
    solar_abundances,
    isotopes,
    segments,
    *,
    code: str = "spectrum",
    resolution: int = 300_000,
    wave_step_nm: float = 0.001,
    verbose: int = 0,
):
    """
    Synthesise a spectrum at the given Teff with zero rotational and
    macroturbulent broadening to obtain intrinsic line depths.

    Parameters
    ----------
    teff, logg, MH : float
        Stellar parameters.  alpha is estimated from MH.
    modeled_layers_pack : iSpec atmosphere pack
    linelist, solar_abundances, isotopes : iSpec data structures
    segments : recarray with wave_base / wave_top columns
        Wavelength windows to synthesise (use _make_segments_for_linelist).
    code : str
        Radiative transfer code: 'spectrum', 'turbospectrum', 'moog'.
    resolution : int
        Instrumental resolution R = λ/δλ.  Use ≥ 300 000 to avoid smearing.
    wave_step_nm : float
        Wavelength sampling step (nm).
    verbose : int
        iSpec verbosity level (0 = silent).

    Returns
    -------
    waveobs : ndarray  (nm)
    flux    : ndarray  (normalised, 0–1)
    """
    import ispec

    alpha = ispec.determine_abundance_enchancements(MH)
    vmic  = ispec.estimate_vmic(teff, logg, MH)

    # Validate that the parameters fall within the grid
    if not ispec.valid_atmosphere_target(
        modeled_layers_pack, {"teff": teff, "logg": logg, "MH": MH, "alpha": alpha}
    ):
        raise ValueError(
            f"Teff={teff}, logg={logg}, [M/H]={MH} is outside the atmosphere grid."
        )

    atmosphere_layers = ispec.interpolate_atmosphere_layers(
        modeled_layers_pack,
        {"teff": teff, "logg": logg, "MH": MH, "alpha": alpha},
        code=code,
    )

    # Build the waveobs grid covering all synthesis segments
    wave_ranges = []
    for seg in segments:
        w = np.arange(seg["wave_base"], seg["wave_top"] + wave_step_nm, wave_step_nm)
        wave_ranges.append(w)
    waveobs = np.concatenate(wave_ranges)
    waveobs = np.unique(waveobs)  # deduplicate where segments touch

    flux = ispec.generate_spectrum(
        waveobs,
        atmosphere_layers,
        teff, logg, MH, alpha,
        linelist, isotopes, solar_abundances,
        fixed_abundances=None,
        microturbulence_vel=vmic,
        macroturbulence=0.0,   # ← no macroturbulence: intrinsic depth
        vsini=0.0,             # ← no rotation: intrinsic depth
        limb_darkening_coeff=0.6,
        R=resolution,
        regions=None,          # regions arg filters the linelist, segments already do that
        verbose=verbose,
        code=code,
    )

    return waveobs, flux


# ---------------------------------------------------------------------------
# Line-depth measurement
# ---------------------------------------------------------------------------

def measure_line_depths(waveobs, flux, line_waves, window_nm: float = 0.15):
    """
    Measure depth = 1 - min(flux) for each line in *line_waves*.

    A line is considered "covered" by the synthesis if at least one waveobs
    point falls within ±window_nm of its centre.  Lines not covered receive
    NaN.

    Parameters
    ----------
    waveobs    : 1-D array, nm
    flux       : 1-D array, same length as waveobs
    line_waves : 1-D array, line centre wavelengths in nm
    window_nm  : half-width of the measurement window (nm)

    Returns
    -------
    depths : 1-D array, same length as line_waves
    """
    depths = np.full(len(line_waves), np.nan)
    for i, lw in enumerate(line_waves):
        mask = (waveobs >= lw - window_nm) & (waveobs <= lw + window_nm)
        if mask.sum() < 3:
            continue
        depths[i] = 1.0 - np.min(flux[mask])
    return depths


# ---------------------------------------------------------------------------
# Polynomial fitting
# ---------------------------------------------------------------------------

def fit_depth_vs_teff(teff_grid, depth_array, degree: int = 2):
    """
    Fit a polynomial to depth(Teff) for each line.

    Parameters
    ----------
    teff_grid   : 1-D array of Teff values (K), shape (N_T,)
    depth_array : 2-D array, shape (N_lines, N_T)
    degree      : polynomial degree (default 2)

    Returns
    -------
    coeffs : 2-D array, shape (N_lines, degree+1)
             Each row is [a0, a1, ..., a_degree] in np.polyfit order
             (highest power first).  NaN rows where < 3 valid points.
    """
    n_lines = depth_array.shape[0]
    coeffs  = np.full((n_lines, degree + 1), np.nan)
    for i in range(n_lines):
        d = depth_array[i]
        valid = np.isfinite(d)
        if valid.sum() < degree + 1:
            continue
        try:
            coeffs[i] = np.polyfit(teff_grid[valid], d[valid], degree)
        except np.linalg.LinAlgError:
            pass
    return coeffs


# ---------------------------------------------------------------------------
# Main grid-building function
# ---------------------------------------------------------------------------

def build_depth_grid(
    ispec_dir: Path,
    teff_grid: np.ndarray,
    logg: float,
    MH: float,
    wave_min: float,
    wave_max: float,
    *,
    linelist_name: str = "VALD.300_1100nm",
    atmosphere_name: str = "MARCS.GES",
    code: str = "spectrum",
    depth_threshold: float = 0.02,
    synthesis_window_nm: float = 0.4,
    measure_window_nm: float = 0.15,
    resolution: int = 300_000,
    wave_step_nm: float = 0.001,
    poly_degree: int = 2,
    verbose: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the depth-vs-Teff grid for all lines in the specified range.

    Returns
    -------
    df_raw : DataFrame
        One row per line.  Columns: wave_nm, species, loggf, EP_eV,
        plus depth_T<teff> for each grid temperature.
    df_poly : DataFrame
        One row per line.  Same identifier columns plus polynomial
        coefficients poly_c0 … poly_c<degree> (highest power first)
        and poly_valid (bool: True if enough valid depths for a fit).
    """
    _bootstrap_ispec(ispec_dir)
    import ispec  # noqa: E402

    # --- load data ----------------------------------------------------------
    modeled_layers_pack, is_atlas = load_atmosphere_grid(ispec_dir, atmosphere_name)
    linelist = load_linelist(ispec_dir, linelist_name, wave_min, wave_max, depth_threshold)
    solar_abundances, isotopes = load_support_data(ispec_dir, is_atlas)

    if len(linelist) == 0:
        raise ValueError("No lines found in the specified wavelength range with the given depth threshold.")

    # --- build synthesis segments -------------------------------------------
    segments = _make_segments_for_linelist(linelist, window_nm=synthesis_window_nm)
    log.info(
        "%d lines → %d synthesis segments covering %.1f–%.1f nm",
        len(linelist), len(segments), wave_min, wave_max,
    )

    line_waves = linelist["wave_nm"]

    # --- identifier columns -------------------------------------------------
    # iSpec recarray field names vary slightly by linelist; handle gracefully
    def _get_field(arr, *names, default=np.nan):
        for n in names:
            if n in arr.dtype.names:
                return arr[n]
        return np.full(len(arr), default)

    species_col = _get_field(linelist, "element", "species", "name")
    loggf_col   = _get_field(linelist, "loggf", "log_gf")
    ep_col      = _get_field(linelist, "lower_state_eV", "excitation_potential",
                              "EP", "lower_energy_eV")

    # --- iterate over Teff grid ---------------------------------------------
    depth_matrix = np.full((len(linelist), len(teff_grid)), np.nan)

    for j, teff in enumerate(teff_grid):
        log.info("Synthesising at Teff = %d K  (logg=%.2f, [M/H]=%.2f) …",
                 int(teff), logg, MH)
        try:
            waveobs, flux = synthesise_at_teff(
                teff, logg, MH,
                modeled_layers_pack, linelist,
                solar_abundances, isotopes, segments,
                code=code, resolution=resolution,
                wave_step_nm=wave_step_nm, verbose=verbose,
            )
        except Exception as exc:
            log.error("  Synthesis failed at Teff=%d: %s", int(teff), exc)
            continue

        depths = measure_line_depths(waveobs, flux, line_waves,
                                     window_nm=measure_window_nm)
        depth_matrix[:, j] = depths

        n_measured = np.sum(np.isfinite(depths))
        log.info("  Measured %d / %d line depths.", n_measured, len(linelist))

    # --- assemble raw DataFrame ---------------------------------------------
    teff_cols = {f"depth_T{int(t)}": depth_matrix[:, j]
                 for j, t in enumerate(teff_grid)}

    df_raw = pd.DataFrame({
        "wave_nm":  line_waves,
        "species":  species_col,
        "loggf":    loggf_col,
        "EP_eV":    ep_col,
        **teff_cols,
    })

    # --- polynomial fits ----------------------------------------------------
    coeffs = fit_depth_vs_teff(teff_grid, depth_matrix, degree=poly_degree)

    coeff_cols = {f"poly_c{k}": coeffs[:, k] for k in range(poly_degree + 1)}
    n_valid = np.sum(np.isfinite(depth_matrix), axis=1)
    poly_valid = n_valid >= (poly_degree + 1)

    df_poly = pd.DataFrame({
        "wave_nm":    line_waves,
        "species":    species_col,
        "loggf":      loggf_col,
        "EP_eV":      ep_col,
        **coeff_cols,
        "n_teff_pts": n_valid,
        "poly_valid": poly_valid,
    })

    log.info(
        "Grid complete.  %d / %d lines have valid polynomial fits.",
        poly_valid.sum(), len(linelist),
    )
    return df_raw, df_poly


# ---------------------------------------------------------------------------
# Convenience: evaluate polynomial at any Teff
# ---------------------------------------------------------------------------

def evaluate_depth_at_teff(df_poly: pd.DataFrame, teff: float,
                            poly_degree: int = 2) -> pd.Series:
    """
    Evaluate the fitted polynomial for every line at a given Teff.

    Parameters
    ----------
    df_poly     : DataFrame returned by build_depth_grid
    teff        : effective temperature (K)
    poly_degree : degree used when building the grid

    Returns
    -------
    pd.Series of predicted depths, indexed by wave_nm
    """
    coeff_cols = [f"poly_c{k}" for k in range(poly_degree + 1)]
    coeffs = df_poly[coeff_cols].values      # shape (N_lines, degree+1)

    # np.poly1d / Horner evaluation for all lines at once
    depths = np.zeros(len(df_poly))
    for k, col in enumerate(coeff_cols):
        power = poly_degree - k
        depths += df_poly[col].values * teff ** power

    result = pd.Series(depths, index=df_poly["wave_nm"].values, name=f"depth_T{int(teff)}")
    result[~df_poly["poly_valid"].values] = np.nan
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--ispec-dir", required=True, type=Path,
                   help="Root directory of your iSpec installation.")
    p.add_argument("--teff-min", type=float, default=3500,
                   help="Minimum Teff in grid (K). Default: 3500")
    p.add_argument("--teff-max", type=float, default=5500,
                   help="Maximum Teff in grid (K). Default: 5500")
    p.add_argument("--teff-step", type=float, default=250,
                   help="Step size between Teff grid points (K). Default: 250")
    p.add_argument("--logg", type=float, default=4.5,
                   help="Surface gravity log g (dex). Default: 4.5")
    p.add_argument("--mh", type=float, default=0.0,
                   help="Metallicity [M/H] (dex). Default: 0.0")
    p.add_argument("--wave-min", type=float, required=True,
                   help="Minimum wavelength (nm).")
    p.add_argument("--wave-max", type=float, required=True,
                   help="Maximum wavelength (nm).")
    p.add_argument("--linelist", default="VALD.300_1100nm",
                   help="Linelist sub-directory under input/linelists/transitions/. "
                        "Default: VALD.300_1100nm")
    p.add_argument("--atmosphere", default="MARCS.GES",
                   help="Atmosphere grid sub-directory under input/atmospheres/. "
                        "Default: MARCS.GES")
    p.add_argument("--code", default="spectrum",
                   choices=["spectrum", "turbospectrum", "moog"],
                   help="Radiative transfer code. Default: spectrum")
    p.add_argument("--depth-threshold", type=float, default=0.02,
                   help="Minimum theoretical depth (solar) to keep a line. "
                        "Default: 0.02")
    p.add_argument("--resolution", type=int, default=300_000,
                   help="Synthesis resolution R = lambda/delta_lambda. Default: 300000")
    p.add_argument("--wave-step", type=float, default=0.001,
                   help="Synthesis wavelength step (nm). Default: 0.001")
    p.add_argument("--poly-degree", type=int, default=2,
                   help="Degree of polynomial fit to depth(Teff). Default: 2")
    p.add_argument("--output", default="depth_grid",
                   help="Output file stem (no extension). Default: depth_grid")
    p.add_argument("--verbose", action="store_true",
                   help="Pass verbose=1 to iSpec synthesis (very chatty).")
    return p.parse_args()


def main():
    args = _parse_args()

    teff_grid = np.arange(args.teff_min, args.teff_max + 0.1, args.teff_step)
    log.info("Teff grid: %s", teff_grid.tolist())

    df_raw, df_poly = build_depth_grid(
        ispec_dir       = args.ispec_dir,
        teff_grid       = teff_grid,
        logg            = args.logg,
        MH              = args.mh,
        wave_min        = args.wave_min,
        wave_max        = args.wave_max,
        linelist_name   = args.linelist,
        atmosphere_name = args.atmosphere,
        code            = args.code,
        depth_threshold = args.depth_threshold,
        resolution      = args.resolution,
        wave_step_nm    = args.wave_step,
        poly_degree     = args.poly_degree,
        verbose         = int(args.verbose),
    )

    stem = args.output
    raw_csv  = Path(stem + "_raw.csv")
    poly_csv = Path(stem + "_poly.csv")
    raw_pkl  = Path(stem + "_raw.pkl")

    df_raw.to_csv(raw_csv,   index=False)
    df_poly.to_csv(poly_csv, index=False)
    df_raw.to_pickle(raw_pkl)

    log.info("Wrote %s", raw_csv)
    log.info("Wrote %s", poly_csv)
    log.info("Wrote %s  (binary, fastest to reload)", raw_pkl)

    # Quick sanity print
    print("\n--- Depth grid summary ---")
    print(f"  Lines in grid:          {len(df_raw)}")
    print(f"  Teff range:             {teff_grid[0]:.0f} – {teff_grid[-1]:.0f} K "
          f"({len(teff_grid)} points)")
    print(f"  Polynomial degree:      {args.poly_degree}")
    print(f"  Lines with valid fit:   {df_poly['poly_valid'].sum()}")
    print(f"\nOutputs saved with stem '{stem}'")


if __name__ == "__main__":
    main()
