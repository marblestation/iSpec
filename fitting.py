import numpy as np
from plotting import *
from continuum import *
from common import *
from radial_velocity import *
from interpolate import *
from pymodelfit import GaussianModel
from pymodelfit import QuadraticModel
from pymodelfit import UniformKnotSplineModel
from scipy.ndimage import gaussian_filter1d

# Reduces outliers considering the median value and 3 sigma (stdev), iterating several times 
# From pymodelfit/utils/alg.py
def sigma_clip(data,sig=3,iters=1,varfunc=np.var,maout=False):
    """
    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.
    
    :param data: input data (will be flattened to 1D)
    :type data: array-like
    :param sig: 
        The number of standard deviations to use as the clipping limit, or 
        the square root of the variance limit.
    :type sig: scalar
    :param iters: 
        The number of iterations to perform clipping for, or None to clip until
        convergence is achieved
    :param varfunc: 
        The method to compute the variance about the. Should take a 1D array as
        input and output a scalar. This will be compared to the square of the
        data as if it is the variance.
    :type varfunc: a callable
    :param maout: If True, return a masked array (see return value for details).
    :type maout: bool
    
    :returns: 
        A :class:`numpy.ma.Maskedarray` with the rejected points masked, if
        `maout` is True. If maout is False, a tuple (filtereddata,mask) is
        returned where the mask is False for rejected points (and matches the
        shape of the input).
    
    """
    data = np.array(data,copy=False)
    oldshape = data.shape
    data = data.ravel()
    
    mask = np.ones(data.size,bool)
    if iters is None:
        lastrej = sum(mask)+1
        while(sum(mask)!=lastrej):
            lastrej = sum(mask)
            do = data - np.median(data[mask])
            mask = do*do <= varfunc(data[mask])*sig**2
    else:
        for i in range(iters):
            do = data - np.median(data[mask])
            mask = do*do <= varfunc(data[mask])*sig**2
        
    if maout:
        return np.ma.MaskedArray(data,~mask,copy='maout'=='copy')
    else:
        return data[mask],mask.reshape(oldshape)


# Remove outliers and fits a uniform knot spline model with 100 knots
# Returns the model
def fit_continuum(spectra):
    # Find outliers considering 3 sigma around the median
    flux_selected, filter_outliers = sigma_clip(spectra['flux'], sig=3, iters=10)

    # Find continuum using a spline model with uniform knots and excluding outliers for the fit
    continuum_model = UniformKnotSplineModel(nknots=100)
    
    # Each point has different weigth depending on its error
    weights = 1/spectra[filter_outliers]['err']

    # Fit
    continuum_model.fitData(spectra[filter_outliers]['waveobs'], spectra[filter_outliers]['flux'], weights=weights)
    
    return continuum_model

# For a given point near the peak of a line, selects the window needed
# for fitting a gaussian
def find_line_window(spectra, loc, smooth=2, window=1):
    # Construct a uniform spaced spectra by interpolating
    wave_filter = (spectra['waveobs'] >= loc - window) & (spectra['waveobs'] <= loc + window)
    spectra_window = spectra[wave_filter]
    min = np.min(spectra_window['waveobs'])
    max = np.max(spectra_window['waveobs'])
    step = (max - min) / len(spectra_window['waveobs'])
    xaxis = np.arange(min, max, step)
    resample_spectra_window = resample_spectra(spectra_window, xaxis)
    
    # Index of the center of the window
    xwi = len(resample_spectra_window['waveobs']) / 2

    # Smooth the spectra to avoid noise problem and improve segment acotation
    resample_spectra_window['flux'] = gaussian_filter1d(resample_spectra_window['flux'], smooth)
    
    # Find inflection points
    # - dysg contains changes of sign (from 1 to -1 or viceversa) when the numbers
    #   change their decreasing or increasing tendency
    dysg = np.sign(np.convolve(resample_spectra_window['flux'], [0,-1,1], mode='same'))
    # - trans contains a 0 where there is no change in tendency and 
    #   2 or -2 if it is the last point of that tendency (inflection point)
    trans = np.convolve(dysg, [-1,1,0], mode='same')

    # Select only the inflection points and find their index position
    transi = np.where(trans!=0)[0]
    # Substract the position of the window's center to the position of the inflection points
    # and extract the sorted index position
    transorti = np.argsort(np.abs(transi-xwi))
    # The center of the window is near transorti[0]
    sgn = trans[transorti[0]]
    # Find the final segment to be used by the model
    # - Position lower/upper than the center of the window & different sign
    loweri = transi[transorti][(transi[transorti] < xwi) & (trans[transorti]!=sgn)][0]
    upperi = transi[transorti][(transi[transorti] > xwi) & (trans[transorti]!=sgn)][0]
    spectra_smooth_slice = resample_spectra_window[loweri:upperi]
    
    wave_filter = (spectra_window['waveobs'] >= resample_spectra_window['waveobs'][loweri]) & (spectra_window['waveobs'] <= resample_spectra_window['waveobs'][upperi])
    spectra_slice = spectra_window[wave_filter]

    return spectra_smooth_slice, spectra_slice

# Fits a gaussian at a given wavelength location
# The continuum can be determined by multiplying per 2 the median and 
# substracting the mean (the best aproximation tested)
# or passing a previously calculated model
def fit_line(spectra_slice, loc, continuum_model=21):
    model = GaussianModel()
    #~ cont = np.zeros(len(spectra_slice['waveobs']))
    #~ cont = continuum_model(spectra_slice['waveobs'])
    if continuum_model == 21:
        cont = 2*np.median(spectra_slice['flux']) - 1*np.mean(spectra_slice['flux'])
    else:
        cont = continuum_model(spectra_slice['waveobs'])
    conterr = 0

    model.mu = loc
    model.sig = 0.02
    model.A = -0.025

    #~ model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=1/spectra_slice['err'], fixedpars=['mu','sig','A'])
    model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, fixedpars=[])

    lineflux = model.integrate(spectra_slice['waveobs'][0], spectra_slice['waveobs'][len(spectra_slice) - 1])
    #~ linefluxerr = 0 # TODO: Figure out

    # Equivalent width        
    mcont, mconterr = np.mean(cont), np.mean(conterr) #continuum should be array by here
    ew = lineflux / mcont 
    #~ ewerr = linefluxerr/mcont - lineflux*mconterr*mcont**-2
    
    return model, cont, lineflux, ew

# Calculate flux for a wavelength grid using a given model
def get_spectra_from_model(model, spectra_wave_grid):
    total_points = len(spectra_wave_grid)
    spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    spectra['waveobs'] = spectra_wave_grid
    spectra['flux'] = model(spectra['waveobs'])
    spectra['err'] = np.zeros(total_points)
    return spectra


#############################################################################################
if __name__ == '__main__':
    #~ spectra = read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s")
    spectra = read_spectra("input/L082N03_spec_norm/05oct08/sp2_Normal/004_vesta_001.s")

    continuum_model = fit_continuum(spectra)
    spectra_continuum = get_spectra_from_model(continuum_model, spectra['waveobs'])

    #~ plot_spectra([spectra, spectra_continuum])
    #~ print np.std(continuum_model.residuals()), continuum_model.pardict
    #~ continuum_model.plot() # Requires --pylabs

    smooth = 2
    window = 1 # nm
    #~ loc = 486.11
    #~ loc = 485.984
    #~ loc = 473.151
    #~ loc = 474.162
    loc = 486.14

    spectra_smooth_slice, spectra_slice = find_line_window(spectra, loc, smooth, window)
    line_model, continuum_value, lineflux, ew = fit_line(spectra_smooth_slice, loc, continuum_model=21)

    spectra_line = get_spectra_from_model(line_model, spectra_smooth_slice['waveobs'])
    # Add the continuum base because the line_model has substracted it
    spectra_line['flux'] += continuum_value

    #~ plot_spectra([spectra_slice, spectra_line])
    #~ print np.std(line_model.residuals()), line_model.pardict
    #~ line_model.plot() # Requires --pylabs
