"""
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
    
    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
from common import *
from interpolate import *
from lines import *
from pymodelfit import GaussianModel


## Determines the radial velocity profile using the lines regions:
## - For each point of each line region, get the spectra flux
## - The result of adding all these spectra fluxes will be a noisy base and a deep gaussian
## Return the radial velocity x coordenates and the normalized fluxes (relative intensities)
def build_radial_velocity_profile(spectra, continuum, lines, rv_lower_limit = -100, rv_upper_limit = 100, rv_step=0.5, frame=None):
    rv = 0 # km/s
    
    rv_lower_limit = int(rv_lower_limit)
    rv_upper_limit = int(rv_upper_limit)
    if rv_lower_limit >= rv_upper_limit:
        raise Exception("Upper radial velocity limit should be greater than lower limit")
    
    if (np.abs(rv_lower_limit) + np.abs(rv_upper_limit)) <= 4*rv_step:
        raise Exception("Radial velocity step too small for the established limits")
    
    min_wave = np.min(spectra['waveobs'])
    max_wave = np.max(spectra['waveobs'])
    
    # Only use the lines that are on the spectra wavelength
    wfilter = (lines['wave_peak'] >= min_wave) & (lines['wave_peak'] <= max_wave)
    lines = lines[wfilter]
    
    if frame != None:
        frame.update_progress(0)
        total_lines = len(lines)
    
    # Speed of light in m/s
    c = 299792458.0
    
    # Build a window with the maximum size as to store the biggest region
    total_window = (np.abs(rv_lower_limit)+np.abs(rv_upper_limit)) / rv_step
    total_window = int(np.ceil(total_window))
    fluxes = np.zeros(total_window)
    
    lower_delta_lambda = ((rv_lower_limit*1000) / c)
    upper_delta_lambda = ((rv_upper_limit*1000) / c)
    line_num = 0
    # For each point of each line region, get the flux of the spectra and sum it all
    for line in lines:
        i = 0
        sindex = 0
        wave_peak = line['wave_peak']
        # Common sized window in km/s for all lines (but different in wavelength)
        window_base = wave_peak + (wave_peak * lower_delta_lambda) # nm
        window_top = wave_peak + (wave_peak * upper_delta_lambda)   # nm
        increment = (window_top - window_base) / total_window
        
        wavelength = window_base
        from_index = 0 # Optimization: discard regions already processed
        while wavelength <= window_top and i < total_window:
            # Outside spectra
            if wavelength > max_wave or wavelength < min_wave:
                fluxes[i] += 0 # Outside spectra
            else:
                flux, sindex = get_flux(spectra[from_index:], wavelength)
                fluxes[i] += flux - continuum(wavelength) # Normalize continuum
            #print wavelength
            if sindex > 4:
                from_index = sindex - 4
            wavelength += increment
            i += 1
        if frame != None:
            current_work_progress = (line_num*1.0 / total_lines) * 100
            frame.update_progress(current_work_progress)
        line_num += 1
    
    #fluxes = fluxes/line_num # Average
    fluxes = fluxes/np.abs(np.sum(fluxes)) # Normalize (relative intensities)
    xcoord = np.linspace( rv_lower_limit, rv_upper_limit, total_window ) # z numbers from x to y
    
    return xcoord, fluxes, len(lines)


## Fits a Gaussian to the radial velocity profile
## Return the fitted model
def model_radial_velocity_profile(xcoord, fluxes, renormalize=False):
    if renormalize:
        # Renormalize continuum
        #renorm_model = UniformKnotSplineModel(nknots=1)
        #renorm_model.fitData(xcoord, fluxes)
        #fluxes = fluxes - renorm_model(xcoord)
        #cont = 2*np.median(fluxes) - 1*np.mean(fluxes)
        cont = np.median(fluxes)
        fluxes = fluxes - cont
    
    # The added flux should have a big gaussian
    model = GaussianModel()
    
    # First estimate for the gaussian mean: the place with the min flux
    model.mu = xcoord[np.argmin(fluxes)]
    model.sig = 20
    model.A = -0.70
    model.fitData(xcoord, fluxes, fixedpars=[])
    
    return model, fluxes # return fluxes just in case they have been renormalized

## Constructs the radial velocity profile, fits a Gaussian and returns the mean (km/s)
def calculate_radial_velocity(spectra, continuum, lines, rv_limit = 200, rv_step=0.5, renormalize=False, frame=None):
    xcoord, fluxes = build_radial_velocity_profile(spectra, continuum, lines, rv_limit=rv_limit, rv_step=rv_step, frame=frame)
    model = model_radial_velocity_profile(xcoord, fluxes, renormalize=renormalize)

    rv = np.round(model.mu, 2)  # km/s
    ##rms = np.mean(model.residuals()) + np.std(model.residuals())
    ## Residual peak should be at least below 1 to be confident about the fit
    #residual_peak = model(model.mu) - np.min(fluxes)
    return rv

# Radial vel in km/s
def correct_radial_velocity(spectra, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    # Correct wavelength scale for radial velocity
    spectra['waveobs'] = spectra['waveobs'] / ((radial_vel / c) + 1)
    return spectra

def correct_radial_velocity_regions(regions, radial_vel, with_peak=False):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    regions['wave_base'] = regions['wave_base'] / ((radial_vel / c) + 1)
    regions['wave_top'] = regions['wave_top'] / ((radial_vel / c) + 1)
    if with_peak:
        regions['wave_peak'] = regions['wave_peak'] / ((radial_vel / c) + 1)
    return regions

# Shift regions considering the known radial velocity of the star
def correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000
    # Oposite sense because we want to correct regions, not the spectra
    radial_vel = radial_vel * -1

    continuum_regions = read_continuum_regions(continuum_regions_filename)
    continuum_regions['wave_base'] = continuum_regions['wave_base'] / ((radial_vel / c) + 1)
    continuum_regions['wave_top'] = continuum_regions['wave_top'] / ((radial_vel / c) + 1)
    write_continuum_regions(continuum_regions, continuum_regions_filename_out)
     

    # Lines
    line_regions = read_line_regions(line_regions_filename)
    line_regions['wave_base'] = line_regions['wave_base'] / ((radial_vel / c) + 1)
    line_regions['wave_top'] = line_regions['wave_top'] / ((radial_vel / c) + 1)
    line_regions['wave_peak'] = line_regions['wave_peak'] / ((radial_vel / c) + 1)
    write_line_regions(line_regions, line_regions_filename_out)

    segment_regions = read_segment_regions(segment_regions_filename)
    segment_regions['wave_base'] = segment_regions['wave_base'] / ((radial_vel / c) + 1)
    segment_regions['wave_top'] = segment_regions['wave_top'] / ((radial_vel / c) + 1)
    write_segment_regions(segment_regions, segment_regions_filename_out)


# Return the precession matrix needed to go from EQUINOX1 (i.e. 1950.0) to EQUINOX2 (i.e. 1975.0).
# Source: http://code.google.com/p/astrolibpy/source/browse/trunk/astrolib/
def premat(equinox1, equinox2, fk4=False):
    """
     INPUTS:
             EQUINOX1 - Original equinox of coordinates, numeric scalar.
             EQUINOX2 - Equinox of precessed coordinates.
     OUTPUT:
            matrix - double precision 3 x 3 precession matrix, used to precess
                        equatorial rectangular coordinates
     OPTIONAL INPUT KEYWORDS:
             fk4    - If this keyword is set, the FK4 (B1950.0) system precession
                        angles are used to compute the precession matrix.    The
                        default is to use FK5 (J2000.0) precession angles
     EXAMPLES:
             Return the precession matrix from 1950.0 to 1975.0 in the FK4 system
                matrix = premat( 1950.0, 1975.0, fk4=True)
     PROCEDURE:
             FK4 constants from "Computational Spherical Astronomy" by Taff (1983),
             p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
             Supplement 1992, page 104 Table 3.211.1.
    """

    deg_to_rad = pi / 180.0e0
    sec_to_rad = deg_to_rad / 3600.e0
    
    t = 0.001e0 * (equinox2 - equinox1)
    
    if not fk4:    
        st = 0.001e0 * (equinox1 - 2000.e0)
        #  Compute 3 rotation angles
        a = sec_to_rad * t * (23062.181e0 + st * (139.656e0 + 0.0139e0 * st) + t * (30.188e0 - 0.344e0 * st + 17.998e0 * t))
        b = sec_to_rad * t * t * (79.280e0 + 0.410e0 * st + 0.205e0 * t) + a
        c = sec_to_rad * t * (20043.109e0 - st * (85.33e0 + 0.217e0 * st) + t * (-42.665e0 - 0.217e0 * st - 41.833e0 * t))
    else:    
        st = 0.001e0 * (equinox1 - 1900.e0)
        #  Compute 3 rotation angles
        a = sec_to_rad * t * (23042.53e0 + st * (139.75e0 + 0.06e0 * st) + t * (30.23e0 - 0.27e0 * st + 18.0e0 * t))
        b = sec_to_rad * t * t * (79.27e0 + 0.66e0 * st + 0.32e0 * t) + a
        c = sec_to_rad * t * (20046.85e0 - st * (85.33e0 + 0.37e0 * st) + t * (-42.67e0 - 0.37e0 * st - 41.8e0 * t))
    
    sina = np.sin(a)
    sinb = np.sin(b)
    sinc = np.sin(c)
    cosa = np.cos(a)
    cosb = np.cos(b)
    cosc = np.cos(c)
    
    r = np.zeros((3, 3))
    r[0,:] = np.array([cosa * cosb * cosc - sina * sinb, sina * cosb + cosa * sinb * cosc, cosa * sinc])
    r[1,:] = np.array([-cosa * sinb - sina * cosb * cosc, cosa * cosb - sina * sinb * cosc, -sina * sinc])
    r[2,:] = np.array([-cosb * sinc, -sinb * sinc, cosc])
    
    return r

# Calculates heliocentric and barycentric velocity components of Earth.
# Source: http://code.google.com/p/astrolibpy/source/browse/trunk/astrolib/
def baryvel(datetime, deq=0):
    """
     EXPLANATION:
             BARYVEL takes into account the Earth-Moon motion, and is useful for
             radial velocity work to an accuracy of  ~1 m/s.
     CALLING SEQUENCE:
             dvel_hel, dvel_bary = baryvel(dje, deq)
     INPUTS:
             datetime - (year, month, day, [hour, minute, second])
             DEQ - (scalar) epoch of mean equinox of dvelh and dvelb. If deq=0
                        then deq is assumed to be equal to dje.
     OUTPUTS:
             DVELH: (vector(3)) heliocentric velocity component. in km/s
             DVELB: (vector(3)) barycentric velocity component. in km/s
    
             The 3-vectors DVELH and DVELB are given in a right-handed coordinate
             system with the +X axis toward the Vernal Equinox, and +Z axis
             toward the celestial pole.
     NOTES:
             Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
             Stumpf claimed an accuracy of 42 cm/s for the velocity.     A
             comparison with the JPL FORTRAN planetary ephemeris program PLEPH
             found agreement to within about 65 cm/s between 1986 and 1994
     EXAMPLE:
             Compute the radial velocity of the Earth toward Altair 
                (19 50 46.999 +08 52 05.96) on 15-Feb-2012
                using the original Stumpf algorithm
            
heliocentric_vel, barycentric_vel = baryvel((2012, 2, 15, 00, 00, 00))

ra = (19.0 + 50.0/60 + 46.999/(60*60)) # hours
ra = ra * 360/24 # degrees
ra = ra * ((2*np.pi) / 360) # radians
dec = (8.0 + 52.0/60 + 5.96/(60*60)) # degrees
dec = dec * ((2*np.pi) / 360) # radians

# Project velocity toward star
vb = barycentric_vel
v = vb[0]*np.cos(dec)*np.cos(ra) + vb[1]*np.cos(dec)*np.sin(ra) + vb[2]*np.sin(dec)

measured_rv = -45.0 # Radial velocity measured from a spectrum 
corrected_rv = measured_rv + v
    """
    dje = obstools.calendar_to_jd(datetime) # Julian ephemeris date.
    
    #Define constants
    dc2pi = 2 * np.pi
    cc2pi = 2 * np.pi
    dc1 = 1.0e0
    dcto = 2415020.0e0
    dcjul = 36525.0e0                            #days in Julian year
    dcbes = 0.313e0
    dctrop = 365.24219572e0                    #days in tropical year (...572 insig)
    dc1900 = 1900.0e0
    au = 1.4959787e8
    
    #Constants dcfel(i,k) of fast changing elements.
    dcfel = np.array([1.7400353e00, 6.2833195099091e02, 5.2796e-6, 6.2565836e00, 6.2830194572674e02, -2.6180e-6, 4.7199666e00, 8.3997091449254e03, -1.9780e-5, 1.9636505e-1, 8.4334662911720e03, -5.6044e-5, 4.1547339e00, 5.2993466764997e01, 5.8845e-6, 4.6524223e00, 2.1354275911213e01, 5.6797e-6, 4.2620486e00, 7.5025342197656e00, 5.5317e-6, 1.4740694e00, 3.8377331909193e00, 5.6093e-6])
    dcfel = np.reshape(dcfel, (8, 3))
    
    #constants dceps and ccsel(i,k) of slowly changing elements.
    dceps = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
    ccsel = np.array([1.675104e-2, -4.179579e-5, -1.260516e-7, 2.220221e-1, 2.809917e-2, 1.852532e-5, 1.589963e00, 3.418075e-2, 1.430200e-5, 2.994089e00, 2.590824e-2, 4.155840e-6, 8.155457e-1, 2.486352e-2, 6.836840e-6, 1.735614e00, 1.763719e-2, 6.370440e-6, 1.968564e00, 1.524020e-2, -2.517152e-6, 1.282417e00, 8.703393e-3, 2.289292e-5, 2.280820e00, 1.918010e-2, 4.484520e-6, 4.833473e-2, 1.641773e-4, -4.654200e-7, 5.589232e-2, -3.455092e-4, -7.388560e-7, 4.634443e-2, -2.658234e-5, 7.757000e-8, 8.997041e-3, 6.329728e-6, -1.939256e-9, 2.284178e-2, -9.941590e-5, 6.787400e-8, 4.350267e-2, -6.839749e-5, -2.714956e-7, 1.348204e-2, 1.091504e-5, 6.903760e-7, 3.106570e-2, -1.665665e-4, -1.590188e-7])
    ccsel = np.reshape(ccsel, (17, 3))
    
    #Constants of the arguments of the short-period perturbations.
    dcargs = np.array([5.0974222e0, -7.8604195454652e2, 3.9584962e0, -5.7533848094674e2, 1.6338070e0, -1.1506769618935e3, 2.5487111e0, -3.9302097727326e2, 4.9255514e0, -5.8849265665348e2, 1.3363463e0, -5.5076098609303e2, 1.6072053e0, -5.2237501616674e2, 1.3629480e0, -1.1790629318198e3, 5.5657014e0, -1.0977134971135e3, 5.0708205e0, -1.5774000881978e2, 3.9318944e0, 5.2963464780000e1, 4.8989497e0, 3.9809289073258e1, 1.3097446e0, 7.7540959633708e1, 3.5147141e0, 7.9618578146517e1, 3.5413158e0, -5.4868336758022e2])
    dcargs = np.reshape(dcargs, (15, 2))
    
    #Amplitudes ccamps(n,k) of the short-period perturbations.
    ccamps = np.array([-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5, -2.490817e-7, -3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5, -1.823138e-7, 6.593466e-7, 1.322572e-5, 9.258695e-6, -4.674248e-7, -3.646275e-7, 1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7, 9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7, 7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7, -2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6, -1.655307e-7, -3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6, -3.736225e-7, 3.442177e-7, 2.671323e-6, 1.832858e-6, -2.394688e-7, -3.478444e-7, 8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8, -1.488378e-6, -1.251789e-5, 5.226868e-7, -2.049301e-7, 0.e0, -8.043059e-6, -2.991300e-6, 1.473654e-7, -3.154542e-7, 0.e0, 3.699128e-6, -3.316126e-6, 2.901257e-7, 3.407826e-7, 0.e0, 2.550120e-6, -1.241123e-6, 9.901116e-8, 2.210482e-7, 0.e0, -6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.e0])
    ccamps = np.reshape(ccamps, (15, 5))
    
    #Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
    ccsec3 = -7.757020e-8
    ccsec = np.array([1.289600e-6, 5.550147e-1, 2.076942e00, 3.102810e-5, 4.035027e00, 3.525565e-1, 9.124190e-6, 9.990265e-1, 2.622706e00, 9.793240e-7, 5.508259e00, 1.559103e01])
    ccsec = np.reshape(ccsec, (4, 3))
    
    #Sidereal rates.
    dcsld = 1.990987e-7                         #sidereal rate in longitude
    ccsgd = 1.990969e-7                         #sidereal rate in mean anomaly
    
    #Constants used in the calculation of the lunar contribution.
    cckm = 3.122140e-5
    ccmld = 2.661699e-6
    ccfdi = 2.399485e-7
    
    #Constants dcargm(i,k) of the arguments of the perturbations of the motion
    # of the moon.
    dcargm = np.array([5.1679830e0, 8.3286911095275e3, 5.4913150e0, -7.2140632838100e3, 5.9598530e0, 1.5542754389685e4])
    dcargm = np.reshape(dcargm, (3, 2))
    
    #Amplitudes ccampm(n,k) of the perturbations of the moon.
    ccampm = np.array([1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7, -2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8, 1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8])
    ccampm = np.reshape(ccampm, (3, 4))
    
    #ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
    ccpamv = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
    dc1mme = 0.99999696e0
    
    #Time arguments.
    dt = (dje - dcto) / dcjul
    tvec = np.array([1e0, dt, dt * dt])
    
    #Values of all elements for the instant(aneous?) dje.
    temp = (np.transpose(np.dot(np.transpose(tvec), np.transpose(dcfel)))) % dc2pi
    dml = temp[0]
    forbel = temp[1:8]
    g = forbel[0]                                 #old fortran equivalence
    
    deps = (tvec * dceps).sum() % dc2pi
    sorbel = (np.transpose(np.dot(np.transpose(tvec), np.transpose(ccsel)))) % dc2pi
    e = sorbel[0]                                 #old fortran equivalence
    
    #Secular perturbations in longitude.
    dummy = np.cos(2.0)
    sn = np.sin((np.transpose(np.dot(np.transpose(tvec[0:2]), np.transpose(ccsec[:,1:3])))) % cc2pi)
    
    #Periodic perturbations of the emb (earth-moon barycenter).
    pertl = (ccsec[:,0] * sn).sum() + dt * ccsec3 * sn[2]
    pertld = 0.0
    pertr = 0.0
    pertrd = 0.0
    for k in range(0, 15):
        a = (dcargs[k,0] + dt * dcargs[k,1]) % dc2pi
        cosa = np.cos(a)
        sina = np.sin(a)
        pertl = pertl + ccamps[k,0] * cosa + ccamps[k,1] * sina
        pertr = pertr + ccamps[k,2] * cosa + ccamps[k,3] * sina
        if k < 11:    
            pertld = pertld + (ccamps[k,1] * cosa - ccamps[k,0] * sina) * ccamps[k,4]
            pertrd = pertrd + (ccamps[k,3] * cosa - ccamps[k,2] * sina) * ccamps[k,4]
    
    #Elliptic part of the motion of the emb.
    phi = (e * e / 4e0) * (((8e0 / e) - e) * np.sin(g) + 5 * np.sin(2 * g) + (13 / 3e0) * e * np.sin(3 * g))
    f = g + phi
    sinf = np.sin(f)
    cosf = np.cos(f)
    dpsi = (dc1 - e * e) / (dc1 + e * cosf)
    phid = 2 * e * ccsgd * ((1 + 1.5 * e * e) * cosf + e * (1.25 - 0.5 * sinf * sinf))
    psid = ccsgd * e * sinf / np.sqrt(dc1 - e * e)
    
    #Perturbed heliocentric motion of the emb.
    d1pdro = dc1 + pertr
    drd = d1pdro * (psid + dpsi * pertrd)
    drld = d1pdro * dpsi * (dcsld + phid + pertld)
    dtl = (dml + phi + pertl) % dc2pi
    dsinls = np.sin(dtl)
    dcosls = np.cos(dtl)
    dxhd = drd * dcosls - drld * dsinls
    dyhd = drd * dsinls + drld * dcosls
    
    #Influence of eccentricity, evection and variation on the geocentric
    # motion of the moon.
    pertl = 0.0
    pertld = 0.0
    pertp = 0.0
    pertpd = 0.0
    for k in range(0, 3):
        a = (dcargm[k,0] + dt * dcargm[k,1]) % dc2pi
        sina = np.sin(a)
        cosa = np.cos(a)
        pertl = pertl + ccampm[k,0] * sina
        pertld = pertld + ccampm[k,1] * cosa
        pertp = pertp + ccampm[k,2] * cosa
        pertpd = pertpd - ccampm[k,3] * sina
    
    #Heliocentric motion of the earth.
    tl = forbel[1] + pertl
    sinlm = np.sin(tl)
    coslm = np.cos(tl)
    sigma = cckm / (1.0 + pertp)
    a = sigma * (ccmld + pertld)
    b = sigma * pertpd
    dxhd = dxhd + a * sinlm + b * coslm
    dyhd = dyhd - a * coslm + b * sinlm
    dzhd = -sigma * ccfdi * np.cos(forbel[2])
    
    #Barycentric motion of the earth.
    dxbd = dxhd * dc1mme
    dybd = dyhd * dc1mme
    dzbd = dzhd * dc1mme
    for k in range(0, 4):
        plon = forbel[k + 3]
        pomg = sorbel[k + 1]
        pecc = sorbel[k + 9]
        tl = (plon + 2.0 * pecc * np.sin(plon - pomg)) % cc2pi
        dxbd = dxbd + ccpamv[k] * (np.sin(tl) + pecc * np.sin(pomg))
        dybd = dybd - ccpamv[k] * (np.cos(tl) + pecc * np.cos(pomg))
        dzbd = dzbd - ccpamv[k] * sorbel[k + 13] * np.cos(plon - sorbel[k + 5])
        
    
    #Transition to mean equator of date.
    dcosep = np.cos(deps)
    dsinep = np.sin(deps)
    dyahd = dcosep * dyhd - dsinep * dzhd
    dzahd = dsinep * dyhd + dcosep * dzhd
    dyabd = dcosep * dybd - dsinep * dzbd
    dzabd = dsinep * dybd + dcosep * dzbd
    
    #Epoch of mean equinox (deq) of zero implies that we should use
    # Julian ephemeris date (dje) as epoch of mean equinox.
    if deq == 0:    
        dvelh = au * (np.array([dxhd, dyahd, dzahd]))
        dvelb = au * (np.array([dxbd, dyabd, dzabd]))
        return (dvelh,dvelb)
    
    #General precession from epoch dje to deq.
    deqdat = (dje - dcto - dcbes) / dctrop + dc1900
    prema = premat(deqdat, deq, fk4=True)
    
    dvelh = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxhd, dyahd, dzahd])))))
    dvelb = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxbd, dyabd, dzabd])))))
    
    return (dvelh, dvelb)



#~ if __name__ == '__main__':
    #~ continuum_regions_filename = "input/LUMBA/UVES_MRD_sun_cmask.txt"
    #~ line_regions_filename = "input/LUMBA/UVES_MRD_sun_Fe-linelist.txt"
    #~ segment_regions_filename = "input/LUMBA/UVES_MRD_sun_segments.txt"
#~ 
    #~ radial_vel = -97.2 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+cas+a
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 14.03 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+leo
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 0 # The original spectra for Arcturus has been corrected for radial velocity
    #~ #radial_vel = -5.19 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=arcturus
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)











