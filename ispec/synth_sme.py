#!/usr/bin/env python
import os
import sys
import numpy as np
import ctypes
import sys
from common import is_sme_support_enabled
from spectrum import create_spectrum_structure, convolve_spectrum, resample_spectrum

class IDL_STRING(ctypes.Structure):
    _fields_ = [
                ('slen', ctypes.c_int),
                ('stype', ctypes.c_short), # To be ignore
                ('s', ctypes.c_char_p),
            ]

def _sme_rtint(mu,inten):
    #;+
    #;NAME:
    #;    RTINT
    #;
    #;PURPOSE:
    #;    Produces a flux profile by integrating intensity profiles (sampled
    #;      at various mu angles) over the visible stellar surface.
    #;
    #;CALLING SEQUENCE:
    #;    flux = RTINT(mu, inten, deltav, vsini, vrt)
    #;
    #;INPUTS:
    #;    MU: (vector(nmu)) cosine of the angle between the outward normal and
    #;      the line of sight for each intensity spectrum in INTEN.
    #;    INTEN: (array(npts,nmu)) intensity spectra at specified values of MU.
    #;    DELTAV: (scalar) velocity spacing between adjacent spectrum points
    #;      in INTEN (same units as VSINI and VRT).
    #;    VSINI (scalar) maximum radial velocity, due to solid-body rotation.
    #;    VRT (scalar) radial-tangential macroturbulence parameter, i.e.
    #;         sqrt(2) times the standard deviation of a Gaussian distribution
    #;         of turbulent velocities. The same distribution function describes
    #;      the radial motions of one component and the tangential motions of
    #;      a second component. Each component covers half the stellar surface.
    #;         See _The Observation and Analysis of Stellar Photospheres_, Gray.
    #;
    #;INPUT KEYWORDS:
    #;    OSAMP: (scalar) internal oversampling factor for convolutions. By
    #;      default convolutions are done using the input points (OSAMP=1),
    #;      but when OSAMP is set to higher integer values, the input spectra
    #;      are first oversampled by cubic spline interpolation.
    #;
    #;OUTPUTS:
    #;    Function Value: (vector(npts)) Disk integrated flux profile.
    #;
    #;RESTRICTIONS:
    #;    Intensity profiles are weighted by the fraction of the projected
    #;      stellar surface they represent, apportioning the area between
    #;      adjacent MU points equally. Additional weights (such as those
    #;      used in a Gauss-Legendre quadrature) can not meaningfully be
    #;      used in this scheme.  About twice as many points are required
    #;      with this scheme to achieve the precision of Gauss-Legendre
    #;      quadrature.
    #;    DELTAV, VSINI, and VRT must all be in the same units (e.g. km/s).
    #;    If specified, OSAMP should be a positive integer.
    #;
    #;AUTHOR'S REQUEST:
    #;    If you use this algorithm in work that you publish, please cite
    #;      Valenti & Anderson 1996, PASP, currently in preparation.
    #;
    #;MODIFICATION HISTORY:
    #;       Feb-88  GM    Created ANA version.
    #;    13-Oct-92 JAV    Adapted from G. Marcy's ANA routine of the same name.
    #;    03-Nov-93 JAV    Switched to annular convolution technique.
    #;    12-Nov-93 JAV    Fixed bug. Intensity components not added when vsini=0.
    #;    14-Jun-94 JAV    Reformatted for "public" release. Heavily commented.
    #;            Pass deltav instead of 2.998d5/deltav. Added osamp
    #;              keyword. Added rebinning logic at end of routine.
    #;            Changed default osamp from 3 to 1.
    #;    20-Feb-95 JAV    Added mu as an argument to handle arbitrary mu sampling
    #;              and remove ambiguity in intensity profile ordering.
    #;            Interpret VTURB as sqrt(2)*sigma instead of just sigma.
    #;            Replaced call_external with call to spl_{init|interp}.
    #;    03-Apr-95 JAV    Multiply flux by !pi to give observed flux.
    #;    24-Oct-95 JAV    Force "nmk" padding to be at least 3 pixels.
    #;    18-Dec-95 JAV    Renamed from dskint() to rtint(). No longer make local
    #;              copy of intensities. Use radial-tangential instead
    #;              of isotropic Gaussian macroturbulence.
    #;    26-Jan-99 JAV    For NMU=1 and VSINI=0, assume resolved solar surface;
    #;              apply R-T macro, but supress vsini broadening.
    #;    01-Apr-99 GMH    Use annuli weights, rather than assuming equal area.
    #;       07-Mar-12 JAV   Force vsini and vmac to be scalars.
    #;-


    #;Convert input MU to projected radii, R, of annuli for a star of unit radius
    #;  (which is just sine, rather than cosine, of the angle between the outward
    #;  normal and the line of sight).
    rmu = np.sqrt(1.0 - np.power(mu, 2))            # use simple trig identity

    ##;Sort the projected radii and corresponding intensity spectra into ascending
    ##;  order (i.e. from disk center to the limb), which is equivalent to sorting
    ##;  MU in descending order.
    isort = np.argsort(rmu)
    rmu = rmu[isort]                #sorted indicies
    nmu = len(mu)                #number of radii

    ###Calculate projected radii for boundaries of disk integration annuli.  The n+1
    ###  boundaries are selected such that r(i+1) exactly bisects the area between
    ###  rmu(i) and rmu(i+1). The innermost boundary, r(0) is set to 0 (disk center)
    ###  and the outermost boundary, r(nmu) is set to 1 (limb).
    r = np.sqrt(0.5 * ( np.power(rmu[0:-1], 2)  + np.power(rmu[1:], 2))) #area midpoints between rmu
    r = np.hstack((0, r, 1))            #bookend with center and limb

    ###Calculate integration weights for each disk integration annulus.  The weight
    ###  is just given by the relative area of each annulus, normalized such that
    ###  the sum of all weights is unity.  Weights for limb darkening are included
    ###  explicitly in the intensity profiles, so they aren't needed here.
    wt = np.power(r[1:], 2) - np.power(r[0:-1], 2)        #weights = relative areas

    ###Loop through annuli, constructing and convolving with rotation kernels.
    flux = np.zeros(inten.shape[0])

    for imu in xrange(nmu):      #loop thru integration annuli
        ###Add contribution from current annulus to the running total.
        flux = flux + wt[imu] * inten[:, isort[imu]]

    flux = np.pi * flux
    return flux


def _sme_librrayversion(sme):
    sme.SMELibraryVersion.restype = ctypes.POINTER(ctypes.c_char_p)
    # External call:
    ptr_err_string = sme.SMELibraryVersion()
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg

def _sme_inputlinelist(sme, atomic_linelist):
    nlines = atomic_linelist.shape[0]
    atomic = np.zeros((8, nlines))
    # atomic number, ionization state, wavelength (in A), excitation energy of lower level (in eV), log(gf), radiative, Stark, and van der Waals damping parameters
    atomic[0] = np.round(map(float,  atomic_linelist['spectrum_moog_species']))
    atomic[1] = atomic_linelist['ion']
    atomic[2] = atomic_linelist['wave_A']
    atomic[3] = atomic_linelist['lower_state_eV']
    atomic[4] = atomic_linelist['loggf']
    atomic[5] = atomic_linelist['rad']
    atomic[6] = atomic_linelist['stark']
    atomic[7] = atomic_linelist['waals']

    elements = []
    for element_s in atomic_linelist['element']:
        # Ensure extra space ("V  1", but "Ni 1")
        element_v = element_s.split()
        if len(element_v) == 2:
            element_s = "%-2s %s" % (element_v[0], element_v[1])
        # Hack: SME requires all molecules to finish with the number 1 or they are not recognized
        #       The VALD linelist have MgH and TiO without the number, so we add it here
        if "SiH" in element_s:
            element_s = "SiH 1   "
        if "MgH" in element_s:
            element_s = "MgH 1   "
        if "TiO" in element_s:
            element_s = "TiO 1   "
        element_s = "%-8s" % (element_s)
        element = IDL_STRING(len(element_s), 0, ctypes.cast(ctypes.create_string_buffer(element_s + " "), ctypes.c_char_p))
        elements.append(element)

    # Transform to ctypes
    Nlines = ctypes.c_int(nlines)
    Species = (IDL_STRING * nlines)(*elements)

    ptr_Nlines = ctypes.pointer(Nlines)
    ptr_Species = ctypes.pointer(Species)
    ptr_Atomic = atomic.ctypes.data_as(ctypes.POINTER((ctypes.c_double * nlines) * 8)) # np.zeros((8, nlines)))


    # Arguments:
    InputLineList_arguments = (ctypes.c_void_p * 3)(ctypes.cast(ptr_Nlines, ctypes.c_void_p), ctypes.cast(ptr_Species, ctypes.c_void_p), ctypes.cast(ptr_Atomic, ctypes.c_void_p))

    # External call:
    sme.InputLineList.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.InputLineList(ctypes.c_int(3), ctypes.byref(InputLineList_arguments))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg


def _sme_inputmodel(sme, teff, logg, MH, microturbulence_vel, atmosphere_layers, solar_abundances, spherical_model):
    OpFlag = [1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0] # Default flags

    nlayers = atmosphere_layers.shape[0]
    #depth = atmosphere_layers[:,7] # lgTau5000
    #depth = 10**atmosphere_layers[:,7] # Tau5000
    #MoTypeT_s = 'TAU'
    depth = atmosphere_layers[:,0] # RHOX
    if spherical_model:
        MoTypeT_s = 'SPH'
    else:
        MoTypeT_s = 'RHOX'
    T = atmosphere_layers[:,1]
    xNe = atmosphere_layers[:,3]

    k = 1.380658e-23      # boltzmanns constant J/K
    Pgas = atmosphere_layers[:,2] # dyn cm^-3
    Pgas = Pgas * 1e-7 # J cm^-3
    xNa = Pgas / (k * T) # cm^-3

    amass= [ 1.008,  4.003,  6.941,  9.012, 10.811, 12.011, 14.007, 15.999,
            18.998, 20.179, 22.990, 24.305, 26.982, 28.086, 30.974, 32.060,
            35.453, 39.948, 39.102, 40.080, 44.956, 47.900, 50.941, 51.996,
            54.938, 55.847, 58.933, 58.710, 63.546, 65.370, 69.720, 72.590,
            74.922, 78.960, 79.904, 83.800, 85.468, 87.620, 88.906, 91.220,
            92.906, 95.940, 98.906,101.070,102.905,106.400,107.868,112.400,
           114.820,118.690,121.750,127.600,126.905,131.300,132.905,137.340,
           138.906,140.120,140.908,144.240,146.000,150.400,151.960,157.250,
           158.925,162.500,164.930,167.260,168.934,170.040,174.970,178.490,
           180.948,183.850,186.200,190.200,192.200,195.090,196.967,200.590,
           204.370,207.190,208.981,210.000,210.000,222.000,223.000,226.025,
           227.000,232.038,230.040,238.029,237.048,242.000,242.000,245.000,
           248.000,252.000,253.000]

    # Compute mass density
    abund = solar_abundances['Abund'][:92]
    abund = np.hstack((abund, [-20]*7)) # Complete until atomic number 99 with unkown abundances represented by -20 (which will be very small since it is 10^-20)
    abund = np.power(10, abund)
    z = np.arange(len(abund)) > 1 # heavier than helium
    abund[z] *= 10**MH # scale
    weight_mole = np.sum(abund*amass) * 1.660e-24
    Rho = xNa * weight_mole + xNe*9.108e-28 # Mass density
    Vturb = np.ones(nlayers) * microturbulence_vel

    # Transform to ctypes
    nDep = ctypes.c_int(nlayers)
    #nDep = ctypes.c_long(nlayers)
    Teff = ctypes.c_double(teff)
    Grav = ctypes.c_double(logg)
    WLstd = ctypes.c_double(5000.) # Reference wavelength
    MoTypeT = IDL_STRING(len(MoTypeT_s), 0, ctypes.cast(ctypes.create_string_buffer(MoTypeT_s), ctypes.c_char_p))
    c_OpFlag = (ctypes.c_int * len(OpFlag))(*OpFlag)

    ptr_nDep = ctypes.pointer(nDep)
    ptr_Teff = ctypes.pointer(Teff)
    ptr_Grav = ctypes.pointer(Grav)
    ptr_WLstd = ctypes.pointer(WLstd)
    ptr_MoTypeT = ctypes.pointer(MoTypeT)
    ptr_OpFlag = ctypes.pointer(c_OpFlag)
    ptr_depth = depth.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(depth)))
    ptr_T = T.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(T)))
    ptr_xNe = xNe.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(xNe)))
    ptr_xNa = xNa.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(xNa)))
    ptr_Rho = Rho.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(Rho)))
    ptr_Vturb = Vturb.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(Vturb))) # microturbulence (in cm/s) at each depth

    if spherical_model:
        radius = atmosphere_layers[0][-1]
        radius = ctypes.c_double(radius)
        #radius = ctypes.c_long(int(radius))
        height = -1*atmosphere_layers[:,8]
        ptr_radius = ctypes.pointer(radius)
        ptr_height = height.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(height)))

        # Arguments:
        InputModel_arguments = (ctypes.c_void_p * 14)(ctypes.cast(ptr_nDep, ctypes.c_void_p), ctypes.cast(ptr_Teff, ctypes.c_void_p), ctypes.cast(ptr_Grav, ctypes.c_void_p), ctypes.cast(ptr_WLstd, ctypes.c_void_p), ctypes.cast(ptr_MoTypeT, ctypes.c_void_p), ctypes.cast(ptr_radius, ctypes.c_void_p), ctypes.cast(ptr_OpFlag, ctypes.c_void_p), ctypes.cast(ptr_depth, ctypes.c_void_p), ctypes.cast(ptr_T, ctypes.c_void_p), ctypes.cast(ptr_xNe, ctypes.c_void_p), ctypes.cast(ptr_xNa, ctypes.c_void_p), ctypes.cast(ptr_Rho, ctypes.c_void_p), ctypes.cast(ptr_Vturb, ctypes.c_void_p), ctypes.cast(ptr_height, ctypes.c_void_p))

        # External call:
        sme.InputModel.restype = ctypes.POINTER(ctypes.c_char_p)
        ptr_err_string = sme.InputModel(ctypes.c_int(14), ctypes.byref(InputModel_arguments))
    else:
        # Arguments:
        InputModel_arguments = (ctypes.c_void_p * 12)(ctypes.cast(ptr_nDep, ctypes.c_void_p), ctypes.cast(ptr_Teff, ctypes.c_void_p), ctypes.cast(ptr_Grav, ctypes.c_void_p), ctypes.cast(ptr_WLstd, ctypes.c_void_p), ctypes.cast(ptr_MoTypeT, ctypes.c_void_p), ctypes.cast(ptr_OpFlag, ctypes.c_void_p), ctypes.cast(ptr_depth, ctypes.c_void_p), ctypes.cast(ptr_T, ctypes.c_void_p), ctypes.cast(ptr_xNe, ctypes.c_void_p), ctypes.cast(ptr_xNa, ctypes.c_void_p), ctypes.cast(ptr_Rho, ctypes.c_void_p), ctypes.cast(ptr_Vturb, ctypes.c_void_p))

        # External call:
        sme.InputModel.restype = ctypes.POINTER(ctypes.c_char_p)
        ptr_err_string = sme.InputModel(ctypes.c_int(12), ctypes.byref(InputModel_arguments))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg

def _sme_inputabund(sme, solar_abundances, MH):
    abund = solar_abundances['Abund'][:92]
    abund = np.hstack((abund, [-20]*7)) # Complete until atomic number 99 with unkown abundances represented by -20 (which will be very small since it is 10^-20)
    abund = np.power(10, abund)
    z = np.arange(len(abund)) > 1 # heavier than helium
    abund[z] *= 10**MH # scale
    Abund = abund.copy()
    Abund[1:] = np.log10(Abund[1:])

    ptr_Abund = Abund.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(Abund)))

    # Arguments:
    ptr_InputAbund_arguments = ctypes.pointer(ptr_Abund)

    # External call:
    sme.InputAbund.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.InputAbund(ctypes.c_int(1), ctypes.byref(ptr_Abund))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg

def _sme_ionization(sme):
    #;Calculate ionization balance for current atmosphere and abundances.
    #;Ionization state is stored in the external library.
    #;Set adopt_eos bit mask to 7 = 1 + 2 + 4 to:
    #;
    #;  (1) adopt particle number densities from EOS,
    #;  (2) adopt electron number densities from EOS,
    #;  (4) and adopt gas densities (g/cm**3) from EOS,
    #;
    #;instead of using values from model atmosphere. Different abundance patterns
    #;in the model atmosphere (usually scaled solar) and SME (may be non-solar)
    #;can affect line shape, e.g. shape of hydrogen lines.
    adopt_eos = 1 + 2 + 4
    ptr_adopt_eos = ctypes.pointer(ctypes.c_int(adopt_eos))

    # Arguments:
    ptr_Ionization_arguments = ctypes.pointer(ptr_adopt_eos)

    # External call:
    sme.Ionization.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.Ionization(ctypes.c_int(1), ctypes.byref(ptr_adopt_eos))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg


def _sme_setvwscale(sme):
    #;Set collisional broadening damping parameter (gamma_6) enhancement factor.
    gam6 = 1.0
    ptr_gam6 = ctypes.pointer(ctypes.c_double(gam6))

    # Arguments:
    SetVWscale_arguments = (ctypes.c_void_p * 1)(ctypes.cast(ptr_gam6, ctypes.c_void_p))

    # External call:
    sme.SetVWscale.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.SetVWscale(ctypes.c_int(1), ctypes.byref(SetVWscale_arguments))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg


def _sme_inputwaverange(sme, wave_base, wave_top):
    Wfirst = ctypes.c_double(wave_base) # A
    Wlast = ctypes.c_double(wave_top) # A
    ptr_Wfirst = ctypes.pointer(Wfirst)
    ptr_Wlast = ctypes.pointer(Wlast)

    # Arguments:
    InputWaveRange_arguments = (ctypes.c_void_p * 2)(ctypes.cast(ptr_Wfirst, ctypes.c_void_p), ctypes.cast(ptr_Wlast, ctypes.c_void_p))

    # External call:
    sme.InputWaveRange.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.InputWaveRange(ctypes.c_int(2), ctypes.byref(InputWaveRange_arguments))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg

def _sme_opacity(sme):
    #;Call external module to compute continuous opacities for current segment.
    #;The external module has a GetOpacity entry point for reading opacities.

    # External call:
    sme.Opacity.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.Opacity(ctypes.c_int(0), None)
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)
    return msg


def _sme_transf(sme, sme_dir, nwmax, keep_lineop=False):
    #;Set number of mu values where specific intensity will be computed.
    nmu = 7 # Number of "equal-area" mu angles at which to calculate specific intensity.
    # The equal-area midpoints of each equal-area annulus for which specific intensities
    # were calculated:
    mu = np.array(np.sqrt(0.5 * (2. * np.arange(nmu)+1) / nmu)[::-1])
    #mu = np.sqrt(0.5 * (2. * np.arange(nmu)+1) / nmu)[::-1]

    #nwmax = 8000 # Maximum number of wavelenghts/fluxes
    #accrt = 0.0001 # Minimum accuracy for sme.sint at wavelength grid points in sme.wint
    #accwi = 0.002 # Minimum accuracy for linear spectrum interpolation vs. wavelength.
    accrt = 0.001 # Minimum accuracy for sme.sint at wavelength grid points in sme.wint
    accwi = 0.003 # Minimum accuracy for linear spectrum interpolation vs. wavelength.

    #;Set flag that controls calculation of line center opacity for every line
    #;in the line list. Calculate line center opacities for the first segment.
    #;No need to recalculate opacities for later segments because one call
    #;handles all lines for all segments.
    #;  KEEP_LINEOP=0, need to calculate line center opacities (first segment)
    #;  KEEP_LINEOP=1, use existing line center opacities (subsequent segments)
    if keep_lineop:
        keep_lineop = 1
    else:
        keep_lineop = 0

    #;Flag that controls continuum intensities returned by external module. Older
    #;versions of the external module did not have the LONG_CONTINUUM flag and
    #;returned continuum intensities only at the first and last wavelength.
    #;
    #;  LONG_CONTINUUM=0, continuum intensities at first and last wavelength
    #;  LONG_CONTINUUM=1, continuum intensities at every wavelength [default]
    long_continuum = 1

    #;Call external module to calculate intensities for the current segment.
    #;
    #;If wavelengths are being reused (CALLRT=0), then NW and WINT_SEG are inputs
    #;to the external module. NW is the number of valid wavelengths in WINT_SEG.
    #;If new wavelengths are being calculated, then the input value of NW is 0L
    #;and the output value is the number of valid wavelength returned in WINT_SEG.
    wint_seg = np.zeros(nwmax) # wavelength range (internal to the library)
    sint_seg = np.zeros((nwmax, nmu)) # line+continuum intensities
    cint_seg = np.zeros((nwmax, nmu)) # all continuum intensities
    cintr_seg = np.zeros(nmu) # red continuum intensity

    nw = 0 # flag value -> new wavelengths
    path = IDL_STRING(len(sme_dir), 0, ctypes.cast(ctypes.create_string_buffer(sme_dir), ctypes.c_char_p))

    nmu = ctypes.c_int(nmu)
    nwmax = ctypes.c_long(nwmax)
    nw = ctypes.c_long(nw)
    accrt = ctypes.c_double(accrt)
    accwi = ctypes.c_double(accwi)
    keep_lineop = ctypes.c_int(keep_lineop)
    long_continuum = ctypes.c_int(long_continuum)

    ptr_nmu = ctypes.pointer(nmu)
    ptr_mu = mu.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(mu)))
    ptr_cintr_seg = cintr_seg.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(cintr_seg)))
    ptr_nwmax = ctypes.pointer(nwmax)
    ptr_nw = ctypes.pointer(nw)
    ptr_wint_seg = wint_seg.ctypes.data_as(ctypes.POINTER(ctypes.c_double * len(wint_seg)))
    ptr_accrt = ctypes.pointer(accrt)
    ptr_accwi = ctypes.pointer(accwi)
    ptr_keep_lineop = ctypes.pointer(keep_lineop)
    ptr_long_continuum = ctypes.pointer(long_continuum)
    ptr_path = ctypes.pointer(path)
    ptr_sint_seg = sint_seg.ctypes.data_as(ctypes.POINTER((ctypes.c_double * nmu.value) * nwmax.value)) # np.zeros((nwmax, nmu)))
    ptr_cint_seg = cint_seg.ctypes.data_as(ctypes.POINTER((ctypes.c_double * nmu.value) * nwmax.value)) # np.zeros((nwmax, nmu)))

    # Arguments:
    Transf_arguments = (ctypes.c_void_p * 13)(
        ctypes.cast(ptr_nmu, ctypes.c_void_p),
        ctypes.cast(ptr_mu, ctypes.c_void_p),
        ctypes.cast(ptr_cint_seg, ctypes.c_void_p),
        ctypes.cast(ptr_cintr_seg, ctypes.c_void_p),
        ctypes.cast(ptr_nwmax, ctypes.c_void_p),
        ctypes.cast(ptr_nw, ctypes.c_void_p),
        ctypes.cast(ptr_wint_seg, ctypes.c_void_p),
        ctypes.cast(ptr_sint_seg, ctypes.c_void_p),
        ctypes.cast(ptr_accrt, ctypes.c_void_p),
        ctypes.cast(ptr_accwi, ctypes.c_void_p),
        ctypes.cast(ptr_keep_lineop, ctypes.c_void_p),
        ctypes.cast(ptr_long_continuum, ctypes.c_void_p),
        ctypes.cast(ptr_path, ctypes.c_void_p))

    # External call:
    sme.Transf.restype = ctypes.POINTER(ctypes.c_char_p)
    ptr_err_string = sme.Transf(ctypes.c_int(13), ctypes.byref(Transf_arguments))
    msg = repr(ctypes.cast(ptr_err_string, ctypes.c_char_p).value)

    nw = nw.value

    # Transpose and trim the intensity arrays. Only points 0 to NW-1 are valid.
    # External module returns arrays that are (nmu,nw) while SME expects (nw,nmu).
    wint_seg = wint_seg[0:nw-1]
    sint_seg = sint_seg[:nw-1,]
    cint_seg = cint_seg[:nw-1,]

    #Calculate the continuum flux spectrum from the continuum intensities
    #on the adaptive wavelength grid used by the external module.
    cflx_seg = _sme_rtint(mu, cint_seg)
    flx_seg = _sme_rtint(mu, sint_seg)

    return wint_seg, flx_seg/cflx_seg # wavelengths in Armstrongs




