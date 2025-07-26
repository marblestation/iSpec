#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from scipy import spatial
from scipy.signal import fftconvolve
import logging

from ispec.spectrum import create_spectrum_structure, resample_spectrum, convolve_spectrum, resample_spectrum, correct_velocity
from ispec.lines import _sampling_uniform_in_velocity

def apply_post_fundamental_effects(waveobs, fluxes, segments, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.60, R=500000, vrad=(0,), verbose=0):
    """
    Apply macroturbulence, rotation (visini), limb darkening coeff and resolution to already generated fundamental synthetic spectrum.
    """
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10

    spectrum = create_spectrum_structure(waveobs, fluxes)
    spectrum.sort(order=['waveobs'])

    if (macroturbulence is not None and macroturbulence > 0) or (vsini is not None and vsini > 0):
        # Build spectrum with sampling uniform in velocity (required by vmac and vsini broadening):
        wave_base = spectrum['waveobs'][0]
        wave_top = spectrum['waveobs'][-1]
        velocity_step = __determine_velocity_step(spectrum)
        waveobs_uniform_in_velocity = _sampling_uniform_in_velocity(wave_base, wave_top, velocity_step)
        fluxes_uniform_in_velocity = np.interp(waveobs_uniform_in_velocity, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0)

        # Apply broadening
        #fluxes_uniform_in_velocity = __vsini_broadening_limbdarkening2(waveobs_uniform_in_velocity, fluxes_uniform_in_velocity, velocity_step, vsini, limb_darkening_coeff)
        fluxes_uniform_in_velocity = __vsini_broadening_limbdarkening(fluxes_uniform_in_velocity, velocity_step, vsini, limb_darkening_coeff)
        fluxes_uniform_in_velocity = __vmac_broadening(fluxes_uniform_in_velocity, velocity_step, macroturbulence)

        # Resample to origin wavelength grid
        fluxes = np.interp(spectrum['waveobs'], waveobs_uniform_in_velocity, fluxes_uniform_in_velocity, left=0.0, right=0.0)
        spectrum['flux'] = fluxes

    if R is not None and R > 0:
        # Convolve (here it is not needed to be with a sampling uniform in velocity, the function is capable of dealing with that)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']

    # Make sure original zeros are set to 1.0 and not modified by the previous broadening operations
    fluxes[zeros] = 1.0e-10

    if type(vrad) not in (tuple, list, np.ndarray):
        raise Exception("Velocity should be an array")

    if np.any(np.asarray(vrad) != 0):
        if len(vrad) != len(segments):
            raise Exception("Velocity should be an array with as many numbers as segments when segments are provided")
        modified = waveobs < 0 # All to false
        for velocity, segment in zip(vrad, segments):
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
            modified = np.logical_or(modified, wfilter)
            spectrum = create_spectrum_structure(waveobs[wfilter], fluxes[wfilter])
            # Opposite velocity to simplify input and output interpretation for spectroscopic binaries; for simple individual line shifts, the sense is not that important
            spectrum = correct_velocity(spectrum, -velocity)
            spectrum = resample_spectrum(spectrum, waveobs[wfilter], method="linear", zero_edges=True)
            fluxes[wfilter] = spectrum['flux']
        fluxes[~modified] = 1.
    return fluxes



def _filter_linelist(linelist, segments):
    # Provide some margin or near-by deep lines might be omitted
    margin = 2. # 2 nm
    # Build wavelength points from regions
    lfilter = None
    for region in segments:
        wave_base = region['wave_base'] - margin
        wave_top = region['wave_top'] + margin

        if lfilter is None:
            lfilter = np.logical_and(linelist['wave_A'] >= wave_base*10., linelist['wave_A'] <= wave_top*10.)
        else:
            lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave_A'] >= wave_base*10., linelist['wave_A'] <= wave_top*10.))

    if lfilter is not None:
        return linelist[lfilter]
    else:
        return linelist


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions for vsini and vmac broadening
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def __lsf_rotate(deltav,vsini,epsilon=0.6):
    # Based on lsf_rotate.pro:
    #  http://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro
    #
    # Adapted from rotin3.f in the SYNSPEC software of Hubeny & Lanz
    # http://nova.astro.umd.edu/index.html    Also see Eq. 17.12 in
    # "The Observation and Analysis of Stellar Photospheres" by D. Gray (1992)
    e1 = 2.0*(1.0 - epsilon)
    e2 = np.pi*epsilon/2.0
    e3 = np.pi*(1.0 - epsilon/3.0)

    npts = np.ceil(2*vsini/deltav)
    if npts % 2 == 0:
        npts += 1
    nwid = np.floor(npts/2)
    x = np.arange(npts) - nwid
    x = x*deltav/vsini
    x1 = np.abs(1.0 - x**2)

    velgrid = x*vsini
    return velgrid, (e1*np.sqrt(x1) + e2*x1)/e3


def __vmac_broadening(flux, velocity_step, vmac):
    """
        velocity_step: fluxes should correspond to a spectrum homogeneously sampled in velocity space
                    with a fixed velocity step [km/s]
        vmac   : macroturbulence velocity [km/s]

        Based on SME's rtint
        It uses radial-tangential instead of isotropic Gaussian macroturbulence.
    """
    if vmac is not None and vmac > 0:
        # mu represent angles that divide the star into equal area annuli,
        # ordered from disk center (mu=1) to the limb (mu=0).
        # But since we don't have intensity profiles at various viewing (mu) angles
        # at this point, we just take a middle point:
        m = 0.5
        # Calc projected simga for radial and tangential velocity distributions.
        sigma = vmac/np.sqrt(2.0) / velocity_step
        sigr = sigma * m
        sigt = sigma * np.sqrt(1.0 - m**2.)
        # Figure out how many points to use in macroturbulence kernel
        nmk = max(min(round(sigma*10), (len(flux)-3)/2), 3)
        # Construct radial macroturbulence kernel w/ sigma of mu*vmac/sqrt(2)
        if sigr > 0:
            xarg = (np.arange(2*nmk+1)-nmk) / sigr   # exponential arg
            #mrkern = np.exp(max((-0.5*(xarg**2)),-20.0))
            mrkern = np.exp(-0.5*(xarg**2))
            mrkern = mrkern/mrkern.sum()
        else:
            mrkern = np.zeros(2*nmk+1)
            mrkern[nmk] = 1.0    #delta function

        # Construct tangential kernel w/ sigma of sqrt(1-mu**2)*vmac/sqrt(2.)
        if sigt > 0:
            xarg = (np.arange(2*nmk+1)-nmk) /sigt
            mtkern = np.exp(-0.5*(xarg**2))
            mtkern = mtkern/mtkern.sum()
        else:
            mtkern = np.zeros(2*nmk+1)
            mtkern[nmk] = 1.0

        ## Sum the radial and tangential components, weighted by surface area
        area_r = 0.5
        area_t = 0.5
        mkern = area_r*mrkern + area_t*mtkern

        # Convolve the flux with the kernel
        flux_conv = 1 - fftconvolve(1-flux, mkern, mode='same') # Fastest
        #import scipy
        #flux_conv = scipy.convolve(flux, mkern, mode='same') # Equivalent but slower

        return flux_conv
    else:
        return flux

def __vsini_broadening_limbdarkening2(waveobs_uniform_in_velocity, flux, velocity_step, vsini, epsilon):
    """
        waveobs_uniform_in_velocity: wavelength
        velocity_step: fluxes should correspond to a spectrum homogeneously sampled in velocity space
                    with a fixed velocity step [km/s]
        vsini   : rotation velocity [km/s]
        epsilon : numeric scalar giving the limb-darkening coefficient,
               default = 0.6 which is typical for  photospheric lines.

        NOTE: The correction is negligible and the new implementation issues a warning
              due to a NaN in the edge of the kernel, that later on is corrected. Thus,
              by default we use the original lsf_rotate function translated to python.

        Bug in lsf_rotate Corrected by Dr. J.I. Bailey (private communication):
            "The crux of the issue was that for small kernels (e.g. on  a small
            multiple of the resolution element) the fractional effects of
            quantization were not taken into account. I used Grey's equation as
            in the Astro rotbrod function but integrated it and made it work
            properly all the way down to 2 pixels/samples (though at that level
             you really don't have much traction on what the vsini is!)

            I initially discovered this via noticing discontinuities in
            chi-square space as I fit for vsini."

        Based on lsf_rotate.pro:
        http://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro

        Adapted from rotin3.f in the SYNSPEC software of Hubeny & Lanz
        http://nova.astro.umd.edu/index.html    Also see Eq. 17.12 in
        "The Observation and Analysis of Stellar Photospheres" by D. Gray (1992)
    """
    if vsini is not None and vsini > 0:
        if epsilon is None:
            epsilon = 0.

        dl = waveobs_uniform_in_velocity[1]-waveobs_uniform_in_velocity[0]
        l0 = (waveobs_uniform_in_velocity[1]+waveobs_uniform_in_velocity[0])*0.5

        dlL = l0*(vsini/2.99792458e5)

        # Nondimensional grid spacing
        dx = dl/dlL

        # Make sure vsini isn't too small to do anything ~.2 for my data
        if dx/2. >= 1.:
            return flux

        # Go out to the the grid point in which dl/dlL=1 falls
        n = np.ceil((2. - dx)/2./dx)*2. + 1.

        # The wavelength grid
        k = np.abs(np.arange(n)- np.floor(n/2))

        kernel_x = (np.arange(n)-np.floor(n/2))*dx

        # Useful constants
        dx2 = dx**2.
        c1 =2.*(1. -epsilon)/np.pi/dlL/(1. - epsilon/3.)
        c2 = 0.5*epsilon/dlL/(1. - epsilon/3.)

        # Compute bulk of kernel
        kernel_y = c2 - c2*dx2/12. - c2*dx2*k**2. + \
                c1/8. * (     np.sqrt(4. - dx2*(1. - 2.*k)**2.) - \
                         2.*k*np.sqrt(4. - dx2*(1. - 2.*k)**2.) + \
                              np.sqrt(4. - dx2*(1. + 2.*k)**2.) + \
                         2.*k*np.sqrt(4. - dx2*(1. + 2.*k)**2.) - \
                         4.*np.arcsin(dx*(k-0.5))/dx + 4.*np.arcsin(dx*(k+0.5))/dx)

        ## Central point
        kernel_y[np.floor(n/2.)] = c2 - (c2*dx2)/12. + \
                       c1*np.sqrt(4. - dx2)/4. + c1*np.arcsin(dx/2.)/dx

        # Edge points
        kernel_y[0] = 1./24./dx*(3.*c1*dx*np.sqrt(4. - dx2*(1. -2.*k[0])**2.)*(1. -2.*k[0]) + \
                          c2*(2. + dx - 2.*dx*k[0])**2.*(4. + dx*(2.*k[0]-1.)) + \
                          12.*c1*np.arccos(dx*(k[0]-.5)))
        kernel_y[0] *= (1. - (k[0]-0.5)*dx)/dx  # Edge point flux compensation
        kernel_y[n-1] = kernel_y[0] # Mirror last point last

        # Integrals done as the average, compensate
        kernel_y *= dx

        # Normalize
        kernel_y /= kernel_y.sum()

        #-- convolve the flux with the kernel
        flux_conv = 1 - fftconvolve(1-flux, kernel_y, mode='same') # Fastest
        #import scipy
        #flux_conv = 1 - scipy.convolve(1-flux, kernel_y, mode='same') # Equivalent but slower
        return flux_conv
    else:
        return flux

def __vsini_broadening_limbdarkening(flux, velocity_step, vsini, epsilon):
    """
        velocity_step: fluxes should correspond to a spectrum homogeneously sampled in velocity space
                    with a fixed velocity step [km/s]
        vsini   : rotation velocity [km/s]
        epsilon : numeric scalar giving the limb-darkening coefficient,
               default = 0.6 which is typical for  photospheric lines.

        Based on lsf_rotate.pro:
        http://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro

        Adapted from rotin3.f in the SYNSPEC software of Hubeny & Lanz
        http://nova.astro.umd.edu/index.html    Also see Eq. 17.12 in
        "The Observation and Analysis of Stellar Photospheres" by D. Gray (1992)
    """
    if vsini is not None and vsini > 0:
        if epsilon is None:
            epsilon = 0.
        kernel_x, kernel_y = __lsf_rotate(velocity_step, vsini, epsilon=epsilon)
        kernel_y /= kernel_y.sum()

        #-- convolve the flux with the kernel
        flux_conv = 1 - fftconvolve(1-flux, kernel_y, mode='same') # Fastest
        #import scipy
        #flux_conv = 1 - scipy.convolve(1-flux, kernel_y, mode='same') # Equivalent but slower
        return flux_conv
    else:
        return flux

def __determine_velocity_step(spectrum):
    # Determine step size for a new model wavelength scale, which must be uniform
    # in velocity to facilitate convolution with broadening kernels. The uniform
    # step size is the largest of:
    wave_base = spectrum['waveobs'][0]
    wave_top = spectrum['waveobs'][-1]
    wmid = (wave_top + wave_base) / 2. # midpoint
    wspan = wave_top - wave_base # width
    # Light speed in vacuum
    #c = 299792458.0 # m/s
    c = 299792.4580 # km/s


    # [1] smallest wavelength step considering the wavelength sampling
    wave_diff = spectrum['waveobs'][1:] - spectrum['waveobs'][:-1]
    min_wave_step = np.min(wave_diff)
    min_wave_step_index = np.argmin(wave_diff)
    vstep1 = min_wave_step / (spectrum['waveobs'][min_wave_step_index] * c)

    # [2] 10% the mean dispersion
    vstep2 = 0.1 * wspan / len(spectrum) / (wmid * c)

    # [3] 0.05 km/s, which is 1% the width of solar line profiles
    vstep3 = 0.05e0

    # Select the largest between 1, 2 and 3:
    velocity_step = np.max((vstep1, vstep2, vstep3))

    return velocity_step
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
