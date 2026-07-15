#!/usr/bin/env python3
"""
fit_harps_spectrum.py
=====================
Fit stellar parameters from a HARPS spectrum using iSpec.

Fitted parameters
-----------------
  Teff          effective temperature
  log g         surface gravity
  [M/H]         global metallicity
  vmic          microturbulence velocity
  vsini         projected rotational velocity
  vmac          macroturbulence velocity (optional; default from empirical GES relation)

  Individual element abundances (second-pass, optional):
    Fe, Ca, Si, Mg, Ti, Na, Al, Cr, Ni, ...

Prior / initial guess (from stepar-syn + isochrone fitting for ~10 Myr, ~K-type, ~1 M_sun)
-------------------------------------------------------------------------------------------
  Teff  ~ 4400 K
  log g ~ 4.3
  [M/H] ~ 0.0   (near-solar assumed for a young disk star)

HARPS input formats
-------------------
  --format s1d   : merged 1-D FITS       (*_s1d_A.fits)
  --format e2ds  : 2-D per-order FITS    (*_e2ds_A.fits)
  --format auto  : try to detect from filename / FITS structure (default)

Instrument notes
----------------
  HARPS (La Silla 3.6 m):  R ≈ 115 000,  wavelength range 378–691 nm
  HARPS-N (TNG 3.58 m):    R ≈ 115 000,  wavelength range 383–693 nm
  Both produce identical DRS output formats.

  The e2ds format stores 72 echelle orders as a 2-D FITS image [72 × 4096].
  Wavelengths are derived from per-order polynomial coefficients stored in
  the FITS header under the keywords
      HIERARCH ESO DRS CAL TH COEFF LLD00 … LLD287
  (72 orders × 4 polynomial terms; wavelength in Ångström).
  If those keywords are absent the script falls back to a companion
  *_wave_A.fits file if present alongside the e2ds file.

Usage examples
--------------
  python fit_harps_spectrum.py --input spectrum_s1d_A.fits
  python fit_harps_spectrum.py --input spectrum_e2ds_A.fits --format e2ds
  python fit_harps_spectrum.py --input spectrum.fits \\
      --teff 4400 --logg 4.3 --mh 0.0 --vsini 15 \\
      --free-vmac --abundances Fe Ca Si Mg Ti \\
      --output-dir results/ --code spectrum
"""

import os
import sys
import argparse
import logging
import json
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# iSpec setup – script must live inside (or alongside) the iSpec directory
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ISPEC_DIR = SCRIPT_DIR          # script is placed at the iSpec root

sys.path.insert(0, str(ISPEC_DIR))

try:
    import ispec
except ImportError as e:
    sys.exit(
        f"ERROR: Cannot import ispec from {ISPEC_DIR}.\n"
        "Make sure this script is in the iSpec root directory and "
        "that iSpec dependencies (numpy, scipy, astropy, …) are installed.\n"
        f"Details: {e}"
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ===========================================================================
# Constants
# ===========================================================================

HARPS_R = 115_000             # resolving power (HARPS / HARPS-N standard mode)

# Default initial stellar parameters (from stepar-syn + isochrone fitting
# for a ~10 Myr, ~K-type, ~1 M_sun star)
INITIAL_TEFF  = 4400.0        # K
INITIAL_LOGG  = 4.3           # dex (cgs)
INITIAL_MH    = 0.00          # [M/H]; near-solar assumed for a young disk star
INITIAL_VSINI = 10.0          # km/s; young stars often rotate faster

# HARPS wavelength coverage (DRS default, nm)
HARPS_WMIN = 378.0
HARPS_WMAX = 691.0

ISPEC_INPUT = ISPEC_DIR / "input"


def _data_path(*parts: str) -> Path:
    """Return an absolute path inside iSpec's input/ directory."""
    return ISPEC_INPUT.joinpath(*parts)


def _check_data():
    """Verify that the iSpec input data directory is present."""
    if not ISPEC_INPUT.is_dir():
        sys.exit(
            "ERROR: iSpec input data directory not found at:\n"
            f"  {ISPEC_INPUT}\n\n"
            "Please download the iSpec data files and place them in the\n"
            "iSpec root directory.  The data are available from:\n"
            "  https://www.blancocuaresma.com/s/iSpec/manual\n"
            "(look for the 'Download input data' section)"
        )


# ===========================================================================
# HARPS FITS readers
# ===========================================================================

def _read_harps_s1d(fits_path: str) -> np.recarray:
    """
    Read a HARPS S1D merged 1-D spectrum (*_s1d_A.fits).

    The HARPS DRS writes a single 1-D flux array into the PRIMARY HDU with
    linear WCS keywords:
        CRVAL1   – wavelength of pixel 1 (Ångström)
        CDELT1   – dispersion (Å per pixel)
        CRPIX1   – reference pixel (usually 1)

    There is no separate error extension in standard HARPS S1D products;
    the error array is left as zeros and is not used by the fitter when
    all values are zero (continuum normalisation handles that case).
    """
    from astropy.io import fits

    log.info("Reading HARPS S1D file: %s", fits_path)

    with fits.open(fits_path) as hdul:
        hdul.info()
        primary = hdul[0]
        hdr  = primary.header
        flux = primary.data.astype(float).flatten()

        npix  = len(flux)
        crval = float(hdr.get("CRVAL1", 0.0))    # Å
        cdelt = float(hdr.get("CDELT1", hdr.get("CD1_1", 1.0)))  # Å/pix
        crpix = float(hdr.get("CRPIX1", 1.0))

        wave = crval + (np.arange(npix) - (crpix - 1)) * cdelt   # Å

        # Error: look for a STAT or ERR extension; otherwise zeros
        err = np.zeros_like(flux)
        for hdu in hdul[1:]:
            from astropy.io.fits import ImageHDU
            if isinstance(hdu, ImageHDU) and hdu.data is not None:
                name = hdu.name.upper()
                if any(k in name for k in ("ERR", "SIGMA", "NOISE", "STAT")):
                    err = hdu.data.astype(float).flatten()
                    break

    return _build_spectrum(wave, flux, err, angstrom=True)


def _harps_wave_from_header(hdr, n_orders: int, npix: int) -> np.ndarray:
    """
    Reconstruct the HARPS per-order wavelength solution from the DRS header
    polynomial coefficients.

    The HARPS DRS stores the wavelength solution as a flat array of polynomial
    coefficients under keys:
        HIERARCH ESO DRS CAL TH COEFF LLD00, LLD01, …, LLD{N-1}

    Typically 4 coefficients per order (cubic polynomial), so N = n_orders × 4.
    For order i and pixel p (0-indexed):
        λ(i, p) = c[i*deg + 0] + c[i*deg + 1]*p + c[i*deg + 2]*p² + …   (Å)

    Returns
    -------
    wave : ndarray, shape (n_orders, npix), in Ångström
    """
    # Collect all LLDxx coefficients in order
    coeffs = []
    idx = 0
    while True:
        key = f"HIERARCH ESO DRS CAL TH COEFF LLD{idx:02d}"
        if key not in hdr:
            # Try without leading zero
            key = f"HIERARCH ESO DRS CAL TH COEFF LLD{idx}"
            if key not in hdr:
                break
        coeffs.append(float(hdr[key]))
        idx += 1

    if not coeffs:
        return None

    n_coeff = len(coeffs)
    deg = n_coeff // n_orders          # polynomial terms per order (usually 4)
    pixels = np.arange(npix, dtype=float)

    wave = np.zeros((n_orders, npix), dtype=float)
    for i in range(n_orders):
        c = coeffs[i * deg : (i + 1) * deg]
        # Evaluate as a polynomial (lowest power first)
        w = np.zeros(npix, dtype=float)
        for k, ck in enumerate(c):
            w += ck * pixels ** k
        wave[i] = w

    return wave   # Å


def _read_harps_e2ds(fits_path: str,
                     wmin_nm: float = 378.0,
                     wmax_nm: float = 691.0) -> np.recarray:
    """
    Read and order-merge a HARPS E2DS per-order spectrum (*_e2ds_A.fits).

    HARPS E2DS files store 72 echelle orders as a [72 × 4096] 2-D image
    in the PRIMARY HDU.  The wavelength solution is recovered from polynomial
    coefficients embedded in the header (``HIERARCH ESO DRS CAL TH COEFF LLDxx``).

    Fallback chain
    --------------
    1. Header polynomial coefficients (preferred; always present in DRS output).
    2. Companion ``*_wave_A.fits`` file in the same directory.
    3. Linear WCS from CRVAL1 / CDELT1 per order (approximate; last resort).

    Parameters
    ----------
    fits_path : str
        Path to the e2ds FITS file.
    wmin_nm, wmax_nm : float
        Wavelength range to keep after merging orders [nm].

    Returns
    -------
    spectrum : np.recarray with fields (waveobs [nm], flux, err)
    """
    from astropy.io import fits as pyfits

    log.info("Reading HARPS E2DS file: %s", fits_path)

    with pyfits.open(fits_path) as hdul:
        hdul.info()
        primary = hdul[0]
        hdr   = primary.header
        f_arr = primary.data.astype(float)   # [n_orders, npix]

        if f_arr.ndim != 2:
            raise ValueError(
                f"Expected a 2-D [orders × pixels] PRIMARY array in {fits_path}; "
                f"got shape {f_arr.shape}."
            )

        n_orders, npix = f_arr.shape
        log.info("HARPS E2DS: %d orders × %d pixels", n_orders, npix)

        # ---- Wavelength: strategy 1 – header polynomial -------------------
        w_arr = _harps_wave_from_header(hdr, n_orders, npix)

        # ---- Wavelength: strategy 2 – companion _wave_A.fits --------------
        if w_arr is None:
            wave_path = fits_path.replace("_e2ds_A.fits", "_wave_A.fits")
            if not Path(wave_path).exists():
                wave_path = fits_path.replace("e2ds", "wave")
            if Path(wave_path).exists():
                log.info("Using companion wavelength file: %s", wave_path)
                with pyfits.open(wave_path) as whdul:
                    w_arr = whdul[0].data.astype(float)
            else:
                log.warning(
                    "No header polynomial coefficients and no companion wave file found. "
                    "Falling back to linear WCS per order (approximate)."
                )

        # ---- Wavelength: strategy 3 – linear WCS per order ----------------
        if w_arr is None:
            crval = float(hdr.get("CRVAL1", 0.0))
            cdelt = float(hdr.get("CDELT1", hdr.get("CD1_1", 1.0)))
            crpix = float(hdr.get("CRPIX1", 1.0))
            pix   = np.arange(npix, dtype=float)
            w_row = crval + (pix - (crpix - 1)) * cdelt
            w_arr = np.tile(w_row, (n_orders, 1))

        # Convert Å → nm if clearly in Ångström (median > 1000)
        if np.nanmedian(w_arr) > 1000:
            w_arr = w_arr / 10.0

        # ---- Error array ---------------------------------------------------
        # Standard HARPS E2DS has no dedicated error extension.
        # Look for a companion *_e2ds_A_blaze.fits or *_e2ds_A_ccf*.fits;
        # absent that, leave errors as zero (continuum normalisation will cope).
        e_arr = np.zeros_like(f_arr)
        for hdu in hdul[1:]:
            from astropy.io.fits import ImageHDU
            if isinstance(hdu, ImageHDU) and hdu.data is not None:
                name = hdu.name.upper()
                if any(k in name for k in ("ERR", "SIGMA", "NOISE")):
                    e_arr = hdu.data.astype(float)
                    break

    # ---- Merge orders -------------------------------------------------------
    all_wave, all_flux, all_err = [], [], []
    for i in range(n_orders):
        _append_order(all_wave, all_flux, all_err,
                      w_arr[i], f_arr[i], e_arr[i])

    if not all_wave:
        raise ValueError(
            f"Could not read any spectral orders from {fits_path}.\n"
            "Check the FITS structure with astropy.io.fits.open(file).info()."
        )

    return _merge_orders(all_wave, all_flux, all_err, wmin_nm, wmax_nm)


def _first_match(names, candidates):
    """Return the first name from *candidates* found in *names* (case-insensitive)."""
    names_upper = [n.upper() for n in names]
    for c in candidates:
        if c.upper() in names_upper:
            return names[names_upper.index(c.upper())]
    return None


def _append_order(all_wave, all_flux, all_err, w, f, e,
                  trim_fraction: float = 0.15):
    """
    Trim the blaze-edge pixels (trim_fraction on each side) and append
    a single order's arrays to the running lists.
    """
    # Remove NaN / non-positive wavelengths
    good = np.isfinite(w) & np.isfinite(f) & (w > 0)
    w, f, e = w[good], f[good], e[good]
    if len(w) == 0:
        return

    # Trim blaze edges
    n = len(w)
    lo = int(n * trim_fraction)
    hi = n - lo
    if lo >= hi:
        return

    all_wave.append(w[lo:hi])
    all_flux.append(f[lo:hi])
    all_err.append(e[lo:hi])


def _merge_orders(all_wave, all_flux, all_err,
                  wmin_nm: float, wmax_nm: float) -> np.recarray:
    """
    Concatenate all orders, keep only wavelengths in [wmin_nm, wmax_nm],
    sort by wavelength, and in overlap regions keep the pixel with higher SNR.
    """
    wave = np.concatenate(all_wave)
    flux = np.concatenate(all_flux)
    err  = np.concatenate(all_err)

    # Wavelength range filter
    mask = (wave >= wmin_nm) & (wave <= wmax_nm)
    wave, flux, err = wave[mask], flux[mask], err[mask]

    # Sort by wavelength
    idx = np.argsort(wave)
    wave, flux, err = wave[idx], flux[idx], err[idx]

    # In overlap regions keep higher-SNR pixel.
    # We use a simple approach: round to nearest 0.001 nm grid and keep max-SNR.
    snr = np.where(err > 0, flux / err, flux / (np.nanmedian(np.abs(flux)) + 1e-30))
    _, keep = np.unique(np.round(wave, 3), return_index=True)

    # For each unique wavelength bin keep the pixel with the highest SNR
    wave_r = np.round(wave, 3)
    unique_w = np.unique(wave_r)
    out_w, out_f, out_e = [], [], []
    for uw in unique_w:
        sel = wave_r == uw
        if sel.sum() == 1:
            out_w.append(wave[sel][0])
            out_f.append(flux[sel][0])
            out_e.append(err[sel][0])
        else:
            best = np.argmax(snr[sel])
            out_w.append(wave[sel][best])
            out_f.append(flux[sel][best])
            out_e.append(err[sel][best])

    return _build_spectrum(np.array(out_w), np.array(out_f), np.array(out_e),
                           angstrom=False)


def _build_spectrum(wave, flux, err, angstrom: bool = True) -> np.recarray:
    """Build an iSpec spectrum recarray from raw arrays."""
    if angstrom:
        wave = wave / 10.0        # Å → nm

    # Remove invalid pixels
    good = np.isfinite(wave) & np.isfinite(flux) & (wave > 0) & (flux > 0)
    wave, flux, err = wave[good], flux[good], err[good]

    if len(wave) == 0:
        raise ValueError("Spectrum is empty after filtering invalid pixels.")

    spectrum = ispec.create_spectrum_structure(wave, flux, err)
    spectrum.sort(order="waveobs")
    return spectrum


def detect_harps_format(fits_path: str) -> str:
    """
    Attempt to auto-detect whether a FITS file is a HARPS S1D or E2DS product.

    Heuristics (in order):
      1. Filename contains 's1d'  → 's1d'
      2. Filename contains 'e2ds' → 'e2ds'
      3. PRIMARY NAXIS == 2 and NAXIS2 in [50..100] (order count) → 'e2ds'
      4. PRIMARY NAXIS == 1 → 's1d'
      5. Default → 's1d'
    """
    from astropy.io import fits

    name = Path(fits_path).name.lower()
    if "s1d" in name:
        return "s1d"
    if "e2ds" in name:
        return "e2ds"

    with fits.open(fits_path) as hdul:
        hdr   = hdul[0].header
        naxis  = hdr.get("NAXIS", 0)
        naxis2 = hdr.get("NAXIS2", 0)
        naxis1 = hdr.get("NAXIS1", 0)

        # E2DS: 2-D array where NAXIS2 = number of echelle orders (50–100)
        if naxis == 2 and 50 <= naxis2 <= 100 and naxis1 > 100:
            return "e2ds"
        if naxis == 1 and naxis1 > 100:
            return "s1d"

    log.warning("Could not auto-detect HARPS format – assuming S1D.")
    return "s1d"


# ===========================================================================
# Radial velocity measurement
# ===========================================================================

def measure_and_correct_rv(spectrum, ispec_dir: Path,
                            initial_teff: float = 4400.0,
                            output_dir: Path = None) -> tuple:
    """
    Measure the radial velocity via CCF with a HARPS/SOPHIE binary mask and
    return the RV-corrected spectrum plus the measured RV [km/s].
    """
    if initial_teff < 3800:
        mask_name = "HARPS_SOPHIE.M5.400_687nm"
    elif initial_teff < 4500:
        mask_name = "HARPS_SOPHIE.K5.378_680nm"
    elif initial_teff < 5200:
        mask_name = "HARPS_SOPHIE.K0.378_679nm"
    elif initial_teff < 6000:
        mask_name = "HARPS_SOPHIE.G2.375_679nm"
    else:
        mask_name = "HARPS_SOPHIE.F0.360_698nm"

    mask_file = ispec_dir / f"input/linelists/CCF/{mask_name}/mask.lst"
    if not mask_file.exists():
        mask_file = ispec_dir / "input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst"
    if not mask_file.exists():
        log.warning("No CCF mask available. Skipping RV measurement; assuming RV = 0 km/s.")
        return spectrum, 0.0, 0.0

    ccf_mask = ispec.read_cross_correlation_mask(str(mask_file))
    mask_label = mask_file.name
    log.info("Measuring radial velocity with mask: %s", mask_label)
    models, ccf = ispec.cross_correlate_with_mask(
        spectrum, ccf_mask,
        lower_velocity_limit=-200,
        upper_velocity_limit=200,
        velocity_step=0.5,
        mask_depth=0.01,
        fourier=False,
        only_one_peak=True,
    )

    if not models:
        log.warning("RV CCF fit failed; assuming RV = 0 km/s.")
        return spectrum, 0.0, 0.0

    g = models[0]
    rv      = float(g.mu())
    rv_err  = float(g.emu())
    fwhm    = float(g.sig()) * 2.3548  # σ → FWHM in km/s
    depth   = float(g.A())
    baseline = float(g.baseline())

    log.info("CCF fit results:")
    log.info("  RV       = %+.4f ± %.4f km/s", rv, rv_err)
    log.info("  FWHM     = %.4f km/s  (σ = %.4f km/s)", fwhm, float(g.sig()))
    log.info("  depth    = %.4f  (relative to baseline %.4f)", depth, baseline)
    log.info("  mask     = %s", mask_label)

    # Save the CCF profile so it can be compared against known RVs externally.
    # The CCF recarray returned by iSpec uses fields 'x' (velocity) and 'y' (CCF).
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        ccf_file = output_dir / "ccf_profile.txt"
        with open(ccf_file, "w") as f:
            f.write("# CCF profile written by fit_harps_spectrum.py\n")
            f.write(f"# mask: {mask_label}\n")
            f.write(f"# fitted RV = {rv:+.6f} km/s  (err = {rv_err:.6f} km/s)\n")
            f.write(f"# FWHM = {fwhm:.6f} km/s\n")
            f.write("# velocity_km_s\tccf\n")
            for vel, val in zip(ccf['x'], ccf['y']):
                f.write(f"{vel:.6f}\t{val:.8f}\n")
        log.info("CCF profile saved to %s", ccf_file)

    corrected = ispec.correct_velocity(spectrum, rv)
    return corrected, round(rv, 4), round(rv_err, 4)


# ===========================================================================
# Continuum normalisation
# ===========================================================================

def normalize_spectrum(spectrum, resolution: int):
    """
    Fit and normalise the continuum using iSpec's spline fitting.

    For a cool K-type star with a dense absorption spectrum, we use a
    conservative spline with wide median/max windows to avoid over-fitting.
    """
    # Diagnose the spectrum before continuum fitting so problems are visible.
    n_pix = len(spectrum)
    wmin  = float(spectrum["waveobs"].min())
    wmax  = float(spectrum["waveobs"].max())
    fmed  = float(np.nanmedian(spectrum["flux"]))
    emed  = float(np.nanmedian(spectrum["err"]))
    log.info(
        "Spectrum going into continuum fit: %d pixels, %.3f–%.3f nm, "
        "median flux=%.4g, median err=%.4g",
        n_pix, wmin, wmax, fmed, emed,
    )

    if n_pix == 0:
        raise ValueError(
            "Spectrum is empty before continuum fitting. "
            "Check that --wmin/--wmax match your data's wavelength range "
            "and that the file is in the expected format (waveobs in nm, "
            "tab-separated, header: waveobs\\tflux\\terr)."
        )

    # If all errors are zero (common for reduced spectra without error arrays),
    # error-weighted continuum fitting hits 1/0 and produces NaN in the
    # continuum model, which then makes every normalised flux value NaN.
    use_errors = bool(np.any(spectrum["err"] > 0))
    if not use_errors:
        log.warning(
            "All flux errors are zero — fitting continuum without error "
            "weighting (use_errors_for_fitting=False)."
        )

    log.info("Fitting continuum...")
    continuum_model = ispec.fit_continuum(
        spectrum,
        from_resolution=resolution,
        nknots=None,          # automatic: 1 knot per 5 nm
        degree=2,
        median_wave_range=0.10,   # wider for cool star blanketing
        max_wave_range=1.0,
        model="Splines",
        order="median+max",
        automatic_strong_line_detection=True,
        strong_line_probability=0.5,
        use_errors_for_fitting=use_errors,
    )
    normalized = ispec.normalize_spectrum(
        spectrum, continuum_model, consider_continuum_errors=False
    )
    # Replace the continuum model with a fixed '1.0' for the fitter
    flat_continuum = ispec.fit_continuum(
        normalized, fixed_value=1.0, model="Fixed value"
    )
    return normalized, flat_continuum


# ===========================================================================
# Stellar-parameter fit (synthetic-spectrum method)
# ===========================================================================

def load_ispec_data(ispec_dir: Path, wmin_nm: float, wmax_nm: float):
    """
    Load model atmospheres, linelist, solar abundances, and isotopes.

    Uses MARCS.GES atmospheres, which are well-suited for cool K-type stars
    (Teff 2500–8000 K, log g 0–5, [M/H] –5 to +1).
    """
    # Model atmospheres – MARCS.GES covers cool stars down to Teff=2500 K
    model_dir = ispec_dir / "input/atmospheres/MARCS.GES/"
    if not model_dir.exists():
        # fall back to plain MARCS
        model_dir = ispec_dir / "input/atmospheres/MARCS/"
    if not model_dir.exists():
        raise FileNotFoundError(
            f"Model atmosphere directory not found.\n"
            f"Checked: {ispec_dir / 'input/atmospheres/MARCS.GES/'}"
        )
    log.info("Loading model atmospheres from %s", model_dir)
    modeled_layers_pack = ispec.load_modeled_layers_pack(str(model_dir))

    # Atomic linelist – GESv6 with hyperfine splitting and isotopes
    linelist_file = ispec_dir / "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    if not linelist_file.exists():
        linelist_file = ispec_dir / "input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    if not linelist_file.exists():
        raise FileNotFoundError(f"Atomic linelist not found near {linelist_file}")
    log.info("Loading linelist from %s", linelist_file.name)
    atomic_linelist = ispec.read_atomic_linelist(
        str(linelist_file),
        wave_base=wmin_nm,
        wave_top=wmax_nm,
    )
    # Keep only lines with meaningful depth in the Sun (proxy for relevance)
    atomic_linelist = atomic_linelist[atomic_linelist["theoretical_depth"] >= 0.01]

    # Solar abundances (Grevesse 2007 for MARCS; 1998 for ATLAS9)
    solar_file = ispec_dir / "input/abundances/Grevesse.2007/stdatom.dat"
    if not solar_file.exists():
        solar_file = ispec_dir / "input/abundances/Grevesse.1998/stdatom.dat"
    if not solar_file.exists():
        raise FileNotFoundError(f"Solar abundances file not found near {solar_file}")
    solar_abundances = ispec.read_solar_abundances(str(solar_file))

    # Isotope data for the SPECTRUM code
    isotope_file = ispec_dir / "input/isotopes/SPECTRUM.lst"
    if not isotope_file.exists():
        raise FileNotFoundError(f"Isotope file not found: {isotope_file}")
    isotopes = ispec.read_isotope_data(str(isotope_file))

    return modeled_layers_pack, atomic_linelist, solar_abundances, isotopes


def load_line_regions(ispec_dir: Path, code: str,
                      spectrum: np.recarray) -> np.recarray:
    """
    Load the spectral line regions (wavelength windows) used during fitting.

    iSpec ships pre-defined regions optimised for R ~ 47,000 spectra observed
    with GES.  These are wavelength windows, not resolution-dependent, so they
    are valid for HARPS at R ~ 115,000 as well.

    For a cool K-type star we use all available regions (including TiI, TiII,
    CaI, FeI lines that are strong in K-type spectra) and trim to the observed
    wavelength range.
    """
    code_name = "moog" if code == "moog-scat" else code

    region_candidates = [
        ispec_dir / f"input/regions/47000_GES/{code_name}_synth_good_for_params_all.txt",
        ispec_dir / f"input/regions/47000_VALD/{code_name}_synth_good_for_params_all.txt",
        ispec_dir / "input/regions/47000_GES/spectrum_synth_good_for_params_all.txt",
        ispec_dir / "input/regions/47000_VALD/spectrum_synth_good_for_params_all.txt",
    ]

    region_file = None
    for candidate in region_candidates:
        if candidate.exists():
            region_file = candidate
            break

    if region_file is None:
        raise FileNotFoundError(
            "Could not find line-region files for parameter fitting.\n"
            f"Searched in {ispec_dir / 'input/regions/'}\n"
            "Make sure the iSpec input data are installed."
        )

    log.info("Loading line regions from %s", region_file.name)
    line_regions = ispec.read_line_regions(str(region_file))

    # Restrict to the observed wavelength range (±0.5 nm margin)
    wmin = spectrum["waveobs"].min() + 0.5
    wmax = spectrum["waveobs"].max() - 0.5
    mask = (line_regions["wave_peak"] >= wmin) & (line_regions["wave_peak"] <= wmax)
    line_regions = line_regions[mask]

    if len(line_regions) == 0:
        raise ValueError(
            "No line regions fall within the spectrum's wavelength range "
            f"({wmin:.1f}–{wmax:.1f} nm).  "
            "Check that the linelist and spectrum wavelength ranges overlap."
        )

    log.info("Using %d spectral line regions for fitting.", len(line_regions))
    return line_regions


def setup_free_abundances(elements: list, solar_abundances, atomic_linelist):
    """
    Build the free-abundance structure for individual-element fitting.

    Parameters
    ----------
    elements : list of str
        Element symbols, e.g. ["Fe", "Ca", "Si"].
    solar_abundances : recarray
        iSpec solar abundances table.
    atomic_linelist : recarray
        Atomic linelist (used to get chemical element data).

    Returns
    -------
    free_abundances : recarray or None
    """
    if not elements:
        return None

    chemical_elements_file = ISPEC_INPUT / "abundances/chemical_elements_symbols.dat"
    if not chemical_elements_file.exists():
        log.warning(
            "Chemical elements file not found (%s). "
            "Skipping individual abundance fitting.",
            chemical_elements_file,
        )
        return None

    chemical_elements = ispec.read_chemical_elements(str(chemical_elements_file))
    free_abundances = ispec.create_free_abundances_structure(
        elements, chemical_elements, solar_abundances
    )
    return free_abundances


def fit_stellar_parameters(
    spectrum,
    flat_continuum,
    modeled_layers_pack,
    atomic_linelist,
    solar_abundances,
    isotopes,
    line_regions,
    initial_teff: float,
    initial_logg: float,
    initial_mh: float,
    initial_vsini: float,
    free_vmac: bool,
    free_abundances,
    code: str,
    resolution: int,
    max_iterations: int,
    limb_darkening: float = 0.6,
    fit_limb_darkening: bool = False,
    vsini_range: tuple = (0.0, 300.0),
):
    """
    Run iSpec's synthetic-spectrum fitting and return the full result tuple.
    """
    initial_alpha = ispec.determine_abundance_enchancements(initial_mh)

    # Microturbulence: empirical GES relation, then allowed to vary freely.
    # For cool K dwarfs (Teff ~ 4400K, log g ~ 4.3) this gives ~ 0.9–1.1 km/s.
    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_mh)
    log.info("Initial vmic (empirical) = %.2f km/s", initial_vmic)

    # Macroturbulence: GES Bergemann relation for cool main-sequence stars.
    # For Teff < 5000 K the relation gives ~ 1–2 km/s.
    initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_mh,
                                       relation="GES")
    log.info("Initial vmac (empirical) = %.2f km/s", initial_vmac)

    initial_limb_darkening = limb_darkening
    initial_vrad = 0.0              # spectrum already RV-corrected; vrad fixed

    free_params = ["teff", "logg", "MH", "vmic", "vsini"]
    if free_vmac:
        free_params.append("vmac")
        log.info("vmac will be fitted freely.")
    else:
        log.info("vmac fixed by GES empirical relation.")
    if fit_limb_darkening:
        free_params.append("limb_darkening_coeff")
        log.info("limb_darkening_coeff will be fitted freely (start=%.3f, bounds=[0,1]).",
                 initial_limb_darkening)
    else:
        log.info("limb_darkening_coeff fixed to %.3f.", initial_limb_darkening)
    log.info("vrad fixed at 0.0 km/s (RV already corrected by CCF step).")

    # Segments: fit windows around each line (+/- 0.25 nm)
    # min_margin preserves the user's window width even in crowded spectral regions
    # where the nearest continuum point would otherwise collapse the window to ~1 pixel.
    line_regions = ispec.adjust_linemasks(spectrum, line_regions, max_margin=0.5, min_margin=0.04)
    segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

    log.info("Starting synthetic-spectrum fitting with %d lines...", len(line_regions))
    log.info("Free parameters: %s", free_params)
    log.info("Synthesis code  : %s", code)
    log.info("Resolution      : R = %d", resolution)

    result = ispec.model_spectrum(
        spectrum,
        flat_continuum,
        modeled_layers_pack,
        atomic_linelist,
        isotopes,
        solar_abundances,
        free_abundances,          # None → scaled-solar; list → individual elements
        None,                     # linelist_free_loggf (not fitting log gf)
        initial_teff,
        initial_logg,
        initial_mh,
        initial_alpha,
        initial_vmic,
        initial_vmac,
        initial_vsini,
        initial_limb_darkening,
        resolution,               # instrumental resolution R
        initial_vrad,
        free_params,
        segments=segments,
        linemasks=line_regions,
        enhance_abundances=True,  # apply alpha-enhancement consistent with [M/H]
        use_errors=True,          # weight fit by spectrum errors (recommended)
        vmic_from_empirical_relation=False,   # vmic is a free parameter
        vmac_from_empirical_relation=(not free_vmac),
        max_iterations=max_iterations,
        code=code,
        use_molecules=False,      # molecular lines not in GESv6 atom list
        tmp_dir=None,
        timeout=3600,
        vsini_range=vsini_range,
    )
    return result


# ===========================================================================
# Abundance fitting (second pass)
# ===========================================================================

def fit_individual_abundances(
    spectrum,
    flat_continuum,
    params: dict,
    modeled_layers_pack,
    atomic_linelist,
    solar_abundances,
    isotopes,
    abundance_elements: list,
    code: str,
    ispec_dir: Path,
):
    """
    Determine individual element abundances from the best-fit stellar parameters.

    Uses the equivalent-width / synthetic method in iSpec for a set of clean
    diagnostic lines per element.
    """
    if not abundance_elements:
        return None

    log.info(
        "Second-pass: fitting individual abundances for: %s",
        ", ".join(abundance_elements),
    )

    chemical_elements_file = ispec_dir / "input/abundances/chemical_elements_symbols.dat"
    if not chemical_elements_file.exists():
        log.warning("chemical_elements_symbols.dat not found; skipping abundance pass.")
        return None

    chemical_elements = ispec.read_chemical_elements(str(chemical_elements_file))

    # Build free-abundance structure
    free_abundances = ispec.create_free_abundances_structure(
        abundance_elements, chemical_elements, solar_abundances
    )

    teff  = params["teff"]
    logg  = params["logg"]
    mh    = params["MH"]
    alpha = params["alpha"]
    vmic  = params["vmic"]

    # Atmosphere layers at best-fit parameters
    if not ispec.valid_atmosphere_target(
        modeled_layers_pack, {"teff": teff, "logg": logg, "MH": mh, "alpha": alpha}
    ):
        log.warning("Best-fit parameters are outside atmosphere grid; skipping abundance pass.")
        return None

    atmosphere_layers = ispec.interpolate_atmosphere_layers(
        modeled_layers_pack,
        {"teff": teff, "logg": logg, "MH": mh, "alpha": alpha},
        code=code,
    )

    # Load per-element line regions
    abund_region_file = ispec_dir / "input/regions/47000_GES/spectrum_synth_good_for_params_all.txt"
    if not abund_region_file.exists():
        log.warning("Abundance line region file not found; skipping abundance pass.")
        return None

    abund_line_regions = ispec.read_line_regions(str(abund_region_file))

    # Filter to lines of the requested elements and within the spectrum range
    wmin = spectrum["waveobs"].min() + 0.5
    wmax = spectrum["waveobs"].max() - 0.5

    # 'note' field in line_regions contains the element label (e.g. 'Fe 1')
    element_filter = np.zeros(len(abund_line_regions), dtype=bool)
    for elem in abundance_elements:
        element_filter |= np.array([
            n.startswith(elem) for n in abund_line_regions["note"]
        ])
    wave_filter = (
        (abund_line_regions["wave_peak"] >= wmin) &
        (abund_line_regions["wave_peak"] <= wmax)
    )
    abund_line_regions = abund_line_regions[element_filter & wave_filter]

    if len(abund_line_regions) == 0:
        log.warning(
            "No abundance-sensitive lines found for elements %s; skipping.",
            abundance_elements,
        )
        return None

    abund_line_regions = ispec.adjust_linemasks(
        spectrum, abund_line_regions, max_margin=0.5
    )

    log.info(
        "Fitting abundances using %d lines for elements: %s",
        len(abund_line_regions),
        ", ".join(abundance_elements),
    )

    found_abundances = ispec.determine_abundances(
        atmosphere_layers,
        teff, logg, mh, alpha,
        abund_line_regions,
        solar_abundances,
        microturbulence_vel=vmic,
        ignore=None,
        verbose=0,
        code=code,
    )

    return found_abundances


# ===========================================================================
# Result reporting and saving
# ===========================================================================

def print_results(params: dict, errors: dict, abundances=None, rv: float = 0.0):
    """Pretty-print the fitted stellar parameters."""
    separator = "─" * 60

    print()
    print(separator)
    print("  iSpec Stellar Parameter Fit (HARPS) — Results")
    print(separator)
    print(f"  Radial velocity (pre-corrected):  {rv:+.2f} km/s")
    print()
    print(f"  Teff   = {params['teff']:7.1f} ± {errors.get('teff', 0):5.1f}  K")
    print(f"  log g  = {params['logg']:7.3f} ± {errors.get('logg', 0):5.3f}  dex (cgs)")
    print(f"  [M/H]  = {params['MH']:+7.3f} ± {errors.get('MH', 0):5.3f}  dex")
    print(f"  vmic   = {params['vmic']:7.2f} ± {errors.get('vmic', 0):5.2f}  km/s")
    print(f"  vmac   = {params['vmac']:7.2f} ± {errors.get('vmac', 0):5.2f}  km/s")
    print(f"  vsini  = {params['vsini']:7.2f} ± {errors.get('vsini', 0):5.2f}  km/s")
    print(f"  alpha  = {params['alpha']:+7.3f}                  dex")
    print()

    if abundances is not None and len(abundances) > 0:
        print("  Individual abundances [X/H]:")
        for row in abundances:
            elem = row["element"] if hasattr(row, "__getitem__") else "?"
            abund = row["Abund"] if hasattr(row, "__getitem__") else 0.0
            print(f"    {elem:<8s}  {abund:+.3f} dex")
        print()

    print(separator)
    print()


def save_results(output_dir: Path, params: dict, errors: dict,
                 modeled_spectrum, obs_spectrum,
                 abundances=None, loggf_found=None,
                 status=None, stats_linemasks=None,
                 rv: float = 0.0, rv_err: float = 0.0):
    """
    Save fitting results to files in *output_dir*:

      stellar_params.json     — best-fit parameters + uncertainties
      modeled_spectrum.fits   — synthetic spectrum at best-fit params
      observed_spectrum.fits  — RV-corrected, normalised observed spectrum
      results.dump            — iSpec binary dump (can be restored later)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # JSON summary
    summary = {
        "rv_kms":    rv,
        "rv_err_kms": rv_err,
        "params":    {k: float(v) for k, v in params.items()},
        "errors":    {k: float(v) for k, v in errors.items()},
    }
    if abundances is not None:
        summary["abundances"] = [
            {"element": str(r["element"]), "Abund": float(r["Abund"])}
            for r in abundances
        ]

    json_path = output_dir / "stellar_params.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    log.info("Results saved to %s", json_path)

    # Synthetic spectrum
    synth_path = output_dir / "modeled_spectrum.fits"
    ispec.write_spectrum(modeled_spectrum, str(synth_path))
    log.info("Synthetic spectrum saved to %s", synth_path)

    # Observed (normalised, RV-corrected) spectrum
    obs_path = output_dir / "observed_spectrum_normalized.fits"
    ispec.write_spectrum(obs_spectrum, str(obs_path))
    log.info("Observed normalised spectrum saved to %s", obs_path)

    # Binary dump (full results for later inspection/restoration)
    dump_path = output_dir / "results.dump"
    ispec.save_results(
        str(dump_path),
        (params, errors, abundances, loggf_found, status, stats_linemasks),
    )
    log.info("Binary dump saved to %s", dump_path)


# ===========================================================================
# CLI argument parsing
# ===========================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- Input/output ---------------------------------------------------
    p.add_argument(
        "--input", "-i", required=True,
        help="HARPS FITS spectrum (*_s1d_A.fits or *_e2ds_A.fits)"
    )
    p.add_argument(
        "--format", "-f",
        choices=["auto", "s1d", "e2ds"],
        default="auto",
        help="HARPS file format: s1d (merged 1-D) or e2ds (per-order 2-D). "
             "Default: auto-detect from filename."
    )
    p.add_argument(
        "--output-dir", "-o", default="ispec_results",
        help="Directory for output files (default: ispec_results/)"
    )

    # ---- Initial stellar parameters (override defaults) -----------------
    p.add_argument("--teff",  type=float, default=INITIAL_TEFF,
                   help=f"Initial Teff [K] (default: {INITIAL_TEFF})")
    p.add_argument("--logg",  type=float, default=INITIAL_LOGG,
                   help=f"Initial log g (default: {INITIAL_LOGG})")
    p.add_argument("--mh",    type=float, default=INITIAL_MH,
                   help=f"Initial [M/H] (default: {INITIAL_MH})")
    p.add_argument("--vsini", type=float, default=INITIAL_VSINI,
                   help=f"Initial vsini [km/s] (default: {INITIAL_VSINI})")

    # ---- Fitting options ------------------------------------------------
    p.add_argument(
        "--vsini-min", type=float, default=0.0,
        metavar="KM_S",
        help="Lower bound on vsini during fitting [km/s] (default: 0)"
    )
    p.add_argument(
        "--vsini-max", type=float, default=300.0,
        metavar="KM_S",
        help="Upper bound on vsini during fitting [km/s] (default: 300)"
    )
    p.add_argument(
        "--free-vmac", action="store_true", default=False,
        help=(
            "Fit vmac as a free parameter instead of using the "
            "GES empirical relation (recommended for young/active stars)"
        )
    )
    p.add_argument(
        "--limb-darkening", type=float, default=0.6,
        metavar="U",
        help=(
            "Linear limb-darkening coefficient u (default: 0.6). "
            "Used as the fixed value when --fit-limb-darkening is not given, "
            "or as the starting point when it is."
        )
    )
    p.add_argument(
        "--fit-limb-darkening", action="store_true", default=False,
        help=(
            "Fit the limb-darkening coefficient as a free parameter "
            "within [0, 1]. Without this flag the coefficient is fixed "
            "to --limb-darkening (equivalent to a delta-function prior). "
            "iSpec uses a least-squares optimiser, not a sampler, so "
            "only box constraints are available; for a tight prior set "
            "--limb-darkening to your best estimate and keep it fixed."
        )
    )
    p.add_argument(
        "--abundances", nargs="+", default=[],
        metavar="ELEMENT",
        help=(
            "Elements for individual abundance fitting (second pass), "
            "e.g.  --abundances Fe Ca Si Mg Ti Na"
        )
    )
    p.add_argument(
        "--code",
        choices=["spectrum", "turbospectrum", "moog", "moog-scat", "synthe"],
        default="spectrum",
        help=(
            "Spectral synthesis code (default: spectrum). "
            "turbospectrum/moog require the respective binaries compiled."
        )
    )
    p.add_argument(
        "--resolution", type=int, default=HARPS_R,
        help=f"Instrumental resolving power R (default: {HARPS_R} for HARPS/HARPS-N)"
    )
    p.add_argument(
        "--max-iterations", type=int, default=20,
        help="Maximum fitting iterations (default: 20)"
    )
    p.add_argument(
        "--no-rv-correction", action="store_true", default=False,
        help="Skip automatic radial velocity correction"
    )
    p.add_argument(
        "--wmin", type=float, default=420.0,
        help="Minimum wavelength to use in fit [nm] (default: 420; "
             "HARPS coverage starts at 378 nm but 420 nm avoids the blue-end blaze drop)"
    )
    p.add_argument(
        "--wmax", type=float, default=691.0,
        help="Maximum wavelength to use in fit [nm] (default: 691, HARPS red limit)"
    )

    return p.parse_args()


# ===========================================================================
# Main
# ===========================================================================

def main():
    args = parse_args()

    # ---- Sanity checks --------------------------------------------------
    _check_data()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    # Check that the requested synthesis code is available
    code = args.code
    code_checks = {
        "turbospectrum": ispec.is_turbospectrum_support_enabled,
        "moog":          ispec.is_moog_support_enabled,
        "moog-scat":     ispec.is_moog_scat_support_enabled,
        "synthe":        ispec.is_synthe_support_enabled,
        "spectrum":      ispec.is_spectrum_support_enabled,
    }
    if code in code_checks and not code_checks[code]():
        log.warning(
            "Synthesis code '%s' is not compiled/available. "
            "Falling back to 'spectrum'.",
            code,
        )
        code = "spectrum"
        if not ispec.is_spectrum_support_enabled():
            sys.exit(
                "ERROR: The 'spectrum' synthesis code is also unavailable.\n"
                "Please compile iSpec's C extension first:\n"
                "  cd <iSpec_dir> && python setup.py build_ext --inplace"
            )

    output_dir = Path(args.output_dir)
    log.info("Output directory: %s", output_dir.resolve())

    # ---- 1. Read spectrum -----------------------------------------------
    input_path = args.input.lower()
    is_fits = any(input_path.endswith(ext)
                  for ext in (".fits", ".fit", ".fits.gz", ".fit.gz"))

    if not is_fits:
        # Plain text file (tab-separated waveobs/flux/err, or NARVAL/ESPaDOnS
        # format) — iSpec's built-in reader handles all text variants.
        log.info("Reading text spectrum: %s", args.input)
        spectrum = ispec.read_spectrum(args.input)
    else:
        fmt = args.format
        if fmt == "auto":
            fmt = detect_harps_format(args.input)
        log.info("HARPS FITS format: %s", fmt.upper())

        if fmt == "e2ds":
            spectrum = _read_harps_e2ds(
                args.input, wmin_nm=HARPS_WMIN, wmax_nm=HARPS_WMAX
            )
        else:
            spectrum = _read_harps_s1d(args.input)

    log.info(
        "Spectrum loaded: %d pixels, %.2f – %.2f nm, median flux = %.2f",
        len(spectrum),
        spectrum["waveobs"].min(),
        spectrum["waveobs"].max(),
        np.nanmedian(spectrum["flux"]),
    )

    # Strip NaN/inf pixels – these silently corrupt the Jacobian in the
    # Levenberg-Marquardt fitter, causing parameters to become NaN.
    finite_mask = np.isfinite(spectrum["flux"]) & np.isfinite(spectrum["waveobs"])
    n_bad = int((~finite_mask).sum())
    if n_bad:
        log.warning("Removing %d NaN/Inf pixels from spectrum before fitting.", n_bad)
    spectrum = spectrum[finite_mask]

    # Trim to the requested wavelength range (GES linelist: 420–920 nm)
    wmin, wmax = args.wmin, args.wmax
    mask = (spectrum["waveobs"] >= wmin) & (spectrum["waveobs"] <= wmax)
    spectrum = spectrum[mask]
    log.info(
        "After wavelength trimming (%g – %g nm): %d pixels",
        wmin, wmax, len(spectrum),
    )

    # ---- 2. Radial velocity correction ----------------------------------
    rv, rv_err = 0.0, 0.0
    if not args.no_rv_correction:
        spectrum, rv, rv_err = measure_and_correct_rv(
            spectrum, ISPEC_DIR, initial_teff=args.teff,
            output_dir=output_dir,
        )

    # ---- 3. Continuum normalisation -------------------------------------
    spectrum, flat_continuum = normalize_spectrum(spectrum, args.resolution)

    # ---- 4. Load iSpec model data ---------------------------------------
    log.info("Loading model atmospheres and linelists...")
    modeled_layers_pack, atomic_linelist, solar_abundances, isotopes = \
        load_ispec_data(ISPEC_DIR, wmin, wmax)

    # ---- 5. Load line regions -------------------------------------------
    # line_regions = load_line_regions(ISPEC_DIR, code, spectrum)

    line_regions = ispec.read_line_regions("/Users/apx061/Desktop/CRIRES/mylinefits_identifications_contrast007to03.txt")
    # Trim to observed wavelength range
    wmin = spectrum["waveobs"].min() + 0.5
    wmax = spectrum["waveobs"].max() - 0.5
    mask = (line_regions["wave_peak"] >= wmin) & (line_regions["wave_peak"] <= wmax)
    line_regions = line_regions[mask]

    if len(line_regions) == 0:
        sys.exit(
            f"ERROR: No line regions fall within the spectrum wavelength range "
            f"({wmin:.2f}–{wmax:.2f} nm).\n"
            "Check that your line mask file wavelengths are in nm and overlap "
            "with the observed spectrum."
        )
    log.info("Using %d line regions for fitting.", len(line_regions))

    # ---- 6. (Optional) free individual abundances -----------------------
    free_abundances = setup_free_abundances(
        args.abundances, solar_abundances, atomic_linelist
    ) if args.abundances else None

    # ---- 7. Fit stellar parameters --------------------------------------
    (obs_spec, modeled_synth_spectrum,
     params, errors,
     abundances_found, loggf_found,
     status, stats_linemasks) = fit_stellar_parameters(
        spectrum=spectrum,
        flat_continuum=flat_continuum,
        modeled_layers_pack=modeled_layers_pack,
        atomic_linelist=atomic_linelist,
        solar_abundances=solar_abundances,
        isotopes=isotopes,
        line_regions=line_regions,
        initial_teff=args.teff,
        initial_logg=args.logg,
        initial_mh=args.mh,
        initial_vsini=args.vsini,
        free_vmac=args.free_vmac,
        free_abundances=free_abundances,
        code=code,
        resolution=args.resolution,
        max_iterations=args.max_iterations,
        limb_darkening=args.limb_darkening,
        fit_limb_darkening=args.fit_limb_darkening,
        vsini_range=(args.vsini_min, args.vsini_max),
    )

    # ---- 8. (Optional) second-pass individual abundances ----------------
    individual_abundances = None
    if args.abundances and free_abundances is None:
        # Fit abundances in a dedicated pass using best-fit stellar parameters
        individual_abundances = fit_individual_abundances(
            spectrum=obs_spec,
            flat_continuum=flat_continuum,
            params=params,
            modeled_layers_pack=modeled_layers_pack,
            atomic_linelist=atomic_linelist,
            solar_abundances=solar_abundances,
            isotopes=isotopes,
            abundance_elements=args.abundances,
            code=code,
            ispec_dir=ISPEC_DIR,
        )
    elif abundances_found is not None and len(abundances_found) > 0:
        individual_abundances = abundances_found

    # ---- 9. Report results ----------------------------------------------
    print_results(params, errors, individual_abundances, rv)

    # ---- 10. Save results -----------------------------------------------
    save_results(
        output_dir=output_dir,
        params=params,
        errors=errors,
        modeled_spectrum=modeled_synth_spectrum,
        obs_spectrum=obs_spec,
        abundances=individual_abundances,
        loggf_found=loggf_found,
        status=status,
        stats_linemasks=stats_linemasks,
        rv=rv,
        rv_err=rv_err,
    )

    log.info("Done.")


if __name__ == "__main__":
    main()
