#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
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
from scipy.interpolate import UnivariateSpline
import numpy.lib.recfunctions as rfn # Extra functions
import numpy as np
import calendar
import re
import shutil
import gzip
import tempfile
import os, errno
#import ipdb
import random
import sys
import log
import logging
import cPickle as pickle
import gzip

def is_spectrum_support_enabled():
    try:
        import synthesizer as __synthesizer_ignore__
        return True
    except:
        return False

def is_turbospectrum_support_enabled():
    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    turbospectrum_dir = ispec_dir + "/synthesizer/turbospectrum/"
    turbospectrum_data = turbospectrum_dir + "/DATA/"
    turbospectrum_bsyn_lu = turbospectrum_dir + "bin/bsyn_lu"
    turbospectrum_eqwidt_lu = turbospectrum_dir + "bin/eqwidt_lu"
    turbospectrum_babsma_lu = turbospectrum_dir + "bin/babsma_lu"

    if not os.path.exists(turbospectrum_eqwidt_lu) or \
            not os.path.exists(turbospectrum_bsyn_lu) or \
            not os.path.exists(turbospectrum_babsma_lu) or \
            not os.path.exists(turbospectrum_data):
        return False
    else:
        return True

def is_moog_support_enabled():
    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    moog_dir = ispec_dir + "/synthesizer/moog/"
    moog_executable = moog_dir + "MOOGSILENT"

    if not os.path.exists(moog_executable) or \
            not os.path.exists(moog_dir):
        return False
    else:
        return True

def is_width_support_enabled():
    from sys import platform as _platform
    if "linux" not in _platform:
        return False
    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    atmos_dir = ispec_dir + "/synthesizer/atmos/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        width_executable = atmos_dir + "bin.amd64/width9.exe"
    else:
        width_executable = atmos_dir + "bin.ia32/width9.exe"

    if not os.path.exists(width_executable) or \
            not os.path.exists(atmos_dir):
        return False
    else:
        return True

def is_ares_support_enabled():
    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    ares_dir = ispec_dir + "/synthesizer/ARES/"
    ares_executable = ares_dir + "bin/ARES"

    if not os.path.exists(ares_executable):
        return False
    else:
        return True

def is_synthe_support_enabled():
    from sys import platform as _platform
    if "linux" not in _platform:
        return False
    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    atmos_dir = ispec_dir + "/synthesizer/atmos/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        xnfpelsyn_executable = atmos_dir + "bin.amd64/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.amd64/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.amd64/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.amd64/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.amd64/synthe.exe"
        spectrv_executable = atmos_dir + "bin.amd64/spectrv.exe"
        rotate_executable = atmos_dir + "bin.amd64/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.amd64/syntoascanga.exe"
    else:
        xnfpelsyn_executable = atmos_dir + "bin.ia32/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.ia32/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.ia32/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.ia32/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.ia32/synthe.exe"
        spectrv_executable = atmos_dir + "bin.ia32/spectrv.exe"
        rotate_executable = atmos_dir + "bin.ia32/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.ia32/syntoascanga.exe"

    if not os.path.exists(atmos_dir) \
            or not os.path.exists(xnfpelsyn_executable) \
            or not os.path.exists(synbeg_executable) \
            or not os.path.exists(rgfallinesnew_executable) \
            or not os.path.exists(rmolescasc_executable) \
            or not os.path.exists(synthe_executable) \
            or not os.path.exists(spectrv_executable) \
            or not os.path.exists(rotate_executable) \
            or not os.path.exists(syntoascanga_executable):
        return False
    else:
        return True

def save_results(dump_filename, data):
    pickle.dump(data, gzip.open(dump_filename, "wb", compresslevel=3), protocol=2)

def restore_results(dump_filename):
    return pickle.load(gzip.open(dump_filename, "rb"))


def report_progress(current_work_progress, last_reported_progress):
    """

    :returns:

        True every 10% of progress.

    """
    return (int(current_work_progress) % 10 == 0 and current_work_progress - last_reported_progress > 10) or last_reported_progress < 0 or current_work_progress == 100


def mkdir_p(path):
    """
    Creates a directory. Same behaviour as 'mkdir -p'.
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

def find_duplicates(a, key):
    """
    Find duplicates in a column of a recarray. This is a simplified version of:
    ::

        import numpy.lib.recfunctions as rfn
        rfn.find_duplicates(...)
    """
    a = np.asanyarray(a).ravel()
    # Get the sorting data (by selecting the corresponding field)
    base = a[key]
    # Get the sorting indices and the sorted data
    sortidx = base.argsort()
    sorteddata = base[sortidx]
    # Compare the sorting data
    flag = (sorteddata[:-1] == sorteddata[1:])
    flag = np.concatenate(([False], flag))
    # We need to take the point on the left as well (else we're missing it)
    flag[:-1] = flag[:-1] + flag[1:]
    duplicates = a[sortidx][flag]
    duplicates_index = sortidx[flag]
    return (duplicates, duplicates_index)


def interquartile_range_filtering(data, k=1.5):
    """
    Interquartile range (IQR) is used to find outliers in data. By default, outliers
    are observations that fall below Quartile1 - k*(IQR) or above Quartile3 + k*(IQR).

    * k = 1.5 represents +/-2.698 * sigma (or standard dev) of a gaussian\
    distribution, which includes the 99.3% of the data.


    """
    # First and third quartile (the second is the median)
    q1 = np.percentile(data, 25) # 25% of the data (left to right)
    q3 = np.percentile(data, 75) # 75%
    # Interquartile range
    iqr = q3 - q1
    sfilter = np.logical_and(data > q1 - k * iqr, data < q3 + k * iqr)
    return data[sfilter], sfilter


def sigma_clipping(data, sig=3, meanfunc=np.mean):
    """
    Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value and 3 sigma (3*stdev),
    iterating until convergence.
    """
    last_total = len(data)

    # First iteration
    stdev = np.std(data)
    diff = data - meanfunc(data)
    sfilter = np.abs(diff) < sig*stdev
    current_total = len(data[sfilter])
    # Continue iterating until convergence (no more points are removed)
    while last_total > current_total:
        #print current_total, stdev
        last_total = current_total

        stdev = np.std(data[sfilter])
        diff = data - meanfunc(data[sfilter])
        sfilter = np.abs(diff) < sig*stdev

        current_total = len(data[sfilter])

    return data[sfilter], sfilter


def find_max_win(x, span=3):
    """
    For an array of values, find local maximum values considering a window
    of "span" elements.
    """
    ret = []
    n = len(x)
    dist = (span + 1) / 2;
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-dist+1, 0])
        l_max = i-1
        r_min = i+1
        r_max = np.min([i+dist-1, n-1])
        is_max = 1;
        # left side
        j = l_min
        while j <= l_max:
            if (x[j] > x[i]):
                is_max = 0;
                break
            j += 1

        # right side
        if (is_max == 1):
            j = r_min
            while j <= r_max:
                if (x[j] > x[i]):
                    is_max = 0;
                    break
                j += 1
        if (is_max == 1):
            ret.append(i)
    return np.asarray(ret)

def find_min_win(x, span=3):
    """
    For an array of values, find local minimum values considering a window
    of "span" elements.
    """
    ret = []
    n = len(x)
    dist = (span + 1) / 2;
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-dist+1, 0])
        l_max = i-1
        r_min = i+1
        r_max = np.min([i+dist-1, n-1])
        is_min = 1;
        # left side
        j = l_min
        while j <= l_max:
            if (x[j] < x[i]):
                is_min = 0;
                break
            j += 1

        # right side
        if (is_min == 1):
            j = r_min
            while j <= r_max:
                if (x[j] < x[i]):
                    is_min = 0;
                    break
                j += 1
        if (is_min == 1):
            ret.append(i)
    return np.asarray(ret)


try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs':[np.get_include()]})
    from common_c import find_local_max_values
    from common_c import find_local_min_values
except:
    print "*********************************************************************"
    print "Not optimized version loaded!"
    print "*********************************************************************"

    def find_local_max_values(x):
        """
        For an array of values, find the position of local maximum values considering only
        the next and previous elements, except they have the same value.
        In that case, the next/previous different value is checked. Therefore,
        ::

            find_local_max([1,2,3,3,2,1,4,3])

        would return:
        ::

            [2, 3, 6]
        """
        ret = []
        n = len(x)
        m = 0;
        for i in np.arange(n):
            l_min = np.max([i-1, 0])
            #l_max = i-1
            #r_min = i+1
            #r_max = np.min([i+1, n-1])
            r_min = np.min([i+1, n-1])
            is_max = True

            # left side
            j = l_min
            # If value is equal, search for the last different value
            while j >= 0 and x[j] == x[i]:
                j -= 1

            if (j < 0 or x[j] > x[i]) and i > 0:
                is_max = False

            # right side
            if is_max:
                j = r_min
                # If value is equal, search for the next different value
                while j < n and x[j] == x[i]:
                    j += 1
                if (j >= n or x[j] > x[i]) and i < n-1:
                    is_max = False

            if is_max:
                ret.append(i)
        return np.asarray(ret)

    def find_local_min_values(x):
        """
        For an array of values, find the position of local maximum values considering only
        the next and previous elements, except they have the same value.
        In that case, the next/previous different value is checked. Therefore,
        ::

            find_local_max([10,9,3,3,9,10,4,30])

        would return:
        ::

            [2, 3, 6]
        """
        ret = []
        n = len(x)
        m = 0;
        for i in np.arange(n):
            l_min = np.max([i-1, 0])
            #l_max = i-1
            #r_min = i+1
            #r_max = np.min([i+1, n-1])
            r_min = np.min([i+1, n-1])
            is_min = True
            # left side
            j = l_min
            # If value is equal, search for the last different value
            while j >= 0 and x[j] == x[i]:
                j -= 1

            if j < 0 or x[j] < x[i]:
                is_min = False

            # right side
            if is_min:
                j = r_min
                # If value is equal, search for the next different value
                while j < n and x[j] == x[i]:
                    j += 1

                if j >= n or x[j] < x[i]:
                    is_min = False

            if is_min:
                ret.append(i)
        return np.asarray(ret)

############## [start] Barycentric vel
def __precession_matrix(equinox1, equinox2, fk4=False):
    """
    Return the precession matrix needed to go from EQUINOX1 (i.e. 1950.0) to EQUINOX2 (i.e. 1975.0).
    The code has been copied from: `astrolib <http://code.google.com/p/astrolibpy/source/browse/trunk/astrolib/>`_
    """

    deg_to_rad = np.pi / 180.0e0
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

def __baryvel(datetime, deq=0):
    """
    Calculates heliocentric and barycentric velocity components of Earth.
    The code has been copied from: `astrolib <http://code.google.com/p/astrolibpy/source/browse/astrolib/baryvel.py>`_
    """
    #dje = astropysics.obstools.calendar_to_jd(datetime) # Julian ephemeris date.
    dje = calendar_to_jd(datetime) # Julian ephemeris date

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
    prema = __precession_matrix(deqdat, deq, fk4=True)

    dvelh = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxhd, dyahd, dzahd])))))
    dvelb = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxbd, dyabd, dzabd])))))
    return (dvelh, dvelb)


def calculate_barycentric_velocity_correction(datetime, coordinates, deq=0):
    """
    Calculates barycentric velocity correction for a given star.
    The code is based on: `astrolib <http://code.google.com/p/astrolibpy/source/browse/astrolib/baryvel.py>`_
    """
    dvelh, dvelb = __baryvel(datetime, deq=2000) # J2000.0

    # Calculate velocity toward a star in a given position
    ra_hours, ra_minutes, ra_seconds, dec_degrees,  dec_minutes, dec_seconds = coordinates
    ra = (ra_hours + ra_minutes/60 + ra_seconds/(60*60)) # hours
    ra = ra * 360/24 # degrees
    ra = ra * ((2*np.pi) / 360) # radians
    dec = (dec_degrees + dec_minutes/60 + dec_seconds/(60*60)) # degrees
    dec = dec * ((2*np.pi) / 360) # radians

    # Project velocity toward star
    barycentric_vel = dvelb[0]*np.cos(dec)*np.cos(ra) + dvelb[1]*np.cos(dec)*np.sin(ra) + dvelb[2]*np.sin(dec) # km/s
    barycentric_vel = np.round(barycentric_vel, 2) # km/s

    return barycentric_vel




############## [end] Barycentric vel

################################################################################
#### [start] Copied from astropysics.obsutils (if not, pyinstaller fails)
# http://packages.python.org/Astropysics/
# https://github.com/eteq/astropysics/blob/master/astropysics/obstools.py

def jd_to_calendar(jd,rounding=1000000,output='datetime',gregorian=None,mjd=False):
    """
    Converts a julian date to a calendar date and time.
    This piece of code has been copied from astropysics.obsutils:

    * `Astropysics <http://packages.python.org/Astropysics/>`_
    * `obstools module <https://github.com/eteq/astropysics/blob/master/astropysics/obstools.py>`_
    """
    import datetime
    from dateutil import tz

    if jd is None:
        jd = calendar_to_jd(datetime.datetime.now(tz.tzlocal()))

    jd = np.array(jd,copy=True,dtype=float)
    scalar = jd.shape == ()
    jd = jd.ravel()

    if mjd:
        jd += mjdoffset

    if rounding > 1000000:
        raise ValueError('rounding cannot exceed a second')
    elif rounding <= 0:
        jd += .5
    else:
        rounding = int(rounding)
        roundingfrac = rounding/86400000000
        jd += .5 + roundingfrac

    z = np.floor(jd).astype(int)
    dec = jd - z #fractional piece

    #fix slight floating-point errors if they hapepn TOOD:check
    dgtr1 = dec>=1.0
    dec[dgtr1] -= 1.0
    z[dgtr1] += 1


    if gregorian is None:
        gregorian = 2299161

    if gregorian is True:
        alpha = ((z-1867216.25)/36524.25).astype(int)
        z += 1 + alpha - alpha//4
    elif gregorian is False:
        pass
    else:
        gmask = z >= gregorian
        alpha = ((z[gmask]-1867216.25)/36524.25).astype(int)
        z[gmask] += 1 + alpha - alpha//4

    b = z + 1524
    c = ((b-122.1)/365.25).astype(int)
    d = (365.25*c).astype(int)
    e = ((b-d)/30.6001).astype(int)

    day = b - d - (30.6001*e).astype(int)

    mmask = e<14
    month = e
    month[mmask] -= 1
    month[~mmask] -= 13
    year = c
    year[month>2] -= 4716
    year[month<=2] -= 4715

    if output == 'fracarray':
        dec = dec-roundingfrac
        dec[dec<0]=0
        return np.array((year,month,day+dec)).T

    if rounding == 1000000:
        secdec = dec*86400
        sec = secdec.astype(int)
        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr
        #sec[sec==secdec] -= 1
        msec = None
    else:
        msec = (dec*86400000000.).astype('int64')
        if rounding > 0:
            div = (msec//1000000)*1000000
            toround = (msec - div)<(2*rounding)
            msec[toround] = div + rounding
            msec  -= rounding

        sec = msec//1000000
        msec -= 1000000*sec

        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr

    if output == 'datetime':
        tzi = tz.tzutc()
        if msec is None:
            ts = (year,month,day,hr%24,min%60,sec%60)
        else:
            ts = (year,month,day,hr%24,min%60,sec%60,msec%1000000)
        res = [datetime.datetime(*t,**dict(tzinfo=tzi)) for t in zip(*ts)]
    elif output == 'array':
        msec = np.zeros_like(sec) if msec is None else msec
        res = np.array([year,month,day,hr%24,min%60,sec%60,msec]).T
    else:
        raise ValueError('invlid output form '+str(output))
    if scalar:
        return res[0]
    else:
        return res




def calendar_to_jd(caltime,tz=None,gregorian=True,mjd=False):
    """
    Convert a calendar date and time to julian date.
    This piece of code has been copied from astropysics.obsutils:

    * `Astropysics <http://packages.python.org/Astropysics/>`_
    * `obstools module <https://github.com/eteq/astropysics/blob/master/astropysics/obstools.py>`_
    """

    #Adapted from xidl  jdcnv.pro
    from datetime import datetime,date,tzinfo

    if caltime is None:
        from dateutil.tz import tzlocal
        datetimes = [datetime.now(tzlocal())]
        scalarout = True
    elif isinstance(caltime,datetime) or isinstance(caltime,date):
        datetimes = [caltime]
        scalarout = True
    elif all([isinstance(ct,datetime) or isinstance(ct,date) for ct in caltime]):
        datetimes = caltime
        scalarout = False
    else:
        datetimes = None
        caltime = list(caltime)
        if not (3 <= len(caltime) < 8):
            raise ValueError('caltime input sequence is invalid size')
        while len(caltime) < 7:
            if len(caltime) == 3:
                #make hours 12
                caltime.append(12*np.ones_like(caltime[-1]))
            else:
                caltime.append(np.zeros_like(caltime[-1]))
        yr,month,day,hr,min,sec,msec = caltime
        scalarout = all([np.shape(v) is tuple() for v in caltime])

    #if input objects are datetime objects, generate arrays
    if datetimes is not None:
        yr,month,day,hr,min,sec,msec = [],[],[],[],[],[],[]
        for dt in datetimes:
            if not hasattr(dt,'hour'):
                dt = datetime(dt.year,dt.month,dt.day,12)

            if tz is None:
                off = dt.utcoffset()
                if off is not None:
                    dt = dt - off

            yr.append(dt.year)
            month.append(dt.month)
            day.append(dt.day)
            hr.append(dt.hour)
            min.append(dt.minute)
            sec.append(dt.second)
            msec.append(dt.microsecond)



    yr = np.array(yr,dtype='int64',copy=False).ravel()
    month = np.array(month,dtype='int64',copy=False).ravel()
    day = np.array(day,dtype='int64',copy=False).ravel()
    hr = np.array(hr,dtype=float,copy=False).ravel()
    min = np.array(min,dtype=float,copy=False).ravel()
    sec = np.array(sec,dtype=float,copy=False).ravel()
    msec = np.array(msec,dtype=float,copy=False).ravel()

    #do tz conversion if tz is provided
    if isinstance(tz,basestring) or isinstance(tz,tzinfo):
        if isinstance(tz,basestring):
            from dateutil import tz
            tzi = tz.gettz(tz)
        else:
            tzi = tz

        utcoffset = []
        for t in zip(yr,month,day,hr,min,sec,msec):
            #microsecond from float component of seconds

            dt = datetime(*[int(ti) for ti in t],**dict(tzinfo=tzi))
            utcdt = dt.utcoffset()
            if utcdt is None:
                utcoffset.append(0)
            else:
                utcoffset.append(utcdt.days*24 + (utcdt.seconds + utcdt.microseconds*1e-6)/3600)
    else:
        utcoffset = tz

#    ly = ((month-14)/12).astype(int) #In leap years, -1 for Jan, Feb, else 0
#    jdn = day - 32075l + 1461l*(yr+4800l+ly)//4

#    jdn += 367l*(month - 2-ly*12)//12 - 3*((yr+4900l+ly)//100)//4

#    res = jdn + (hr/24.0) + min/1440.0 + sec/86400.0 - 0.5

    #this algorithm from meeus 2ed
    m3 = month < 3
    yr[m3] -= 1
    month[m3] += 12

    cen = yr//100

    if gregorian is None:
        gregorian = (1582,10,4)
    if gregorian is True:
        gregoffset = 2 - cen + cen//4
    elif gregorian is False:
        gregoffset = 0
    else:
        gregoffset = 2 - cen + cen//4
        gmask = (yr>gregorian[0])&(month>gregorian[1])&(day>gregorian[2])
        gregoffset[~gmask] = 0


    jdn = (365.25*(yr+4716)).astype(int) + \
          (30.6001*(month + 1)).astype(int) + \
               day + gregoffset - 1524.5
    res = jdn + hr/24.0 + min/1440.0 + sec/86400.0

    if mjd:
        res -= mjdoffset

    if np.any(utcoffset):
        res -= np.array(utcoffset)/24.0

    if scalarout:
        return res[0]
    else:
        return res
#### [end]   Astropysics

def estimate_vmic(teff, logg, feh):
    """
    Estimate Microturbulence velocity (Vmic) by using an empirical relation
    considering the effective temperature, surface gravity and metallicity.

    The relation was constructed based on the UVES Gaia ESO Survey iDR1 data,
    results for the benchmark stars (Jofre et al. 2013),
    and globular cluster data from external literature sources.

    Source: http://great.ast.cam.ac.uk/GESwiki/GesWg/GesWg11/Microturbulence
    """
    t0 = 5500
    g0 = 4.0

    if logg >= 3.5:
        if teff >= 5000:
            # main sequence and subgiants (RGB)
            vmic = 1.05 + 2.51e-4*(teff-t0) + 1.5e-7*(teff-t0)**2 - 0.14*(logg-g0) - 0.05e-1*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
        else:
            # main sequence
            vmic = 1.05 + 2.51e-4*(5000-t0) + 1.5e-7*(5000-t0)**2 - 0.14*(logg-g0) - 0.05e-1*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
    else:
        # giants (RGB/AGB)
        vmic = 1.25 + 4.01e-4*(teff-t0) + 3.1e-7*(teff-t0)**2 - 0.14*(logg-g0) - 0.05*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
    vmic = float("%.2f" % vmic)
    return vmic



def estimate_vmac(teff, logg, feh):
    """
    Estimate Microturbulence velocity (Vmic) by using an empirical relation
    considering the effective temperature, surface gravity and metallicity.

    The relation was constructed by Maria Bergemann for the Gaia ESO Survey.
    """
    t0 = 5500
    g0 = 4.0

    if logg >= 3.5:
        if teff >= 5000:
            # main sequence and subgiants (RGB)
            vmac = 3*(1.15 + 7e-4*(teff-t0) + 1.2e-6*(teff-t0)**2 - 0.13*(logg-g0) + 0.13*(logg-g0)**2 - 0.37*feh - 0.07*feh**2)
        else:
            # main sequence
            vmac = 3*(1.15 + 2e-4*(teff-t0) + 3.95e-7*(teff-t0)**2 - 0.13*(logg-g0) + 0.13*(logg-g0)**2)
    else:
        # giants (RGB/AGB)
        vmac = 3*(1.15 + 2.2e-5*(teff-t0) - 0.5e-7*(teff-t0)**2 - 0.1*(logg-g0) + 0.04*(logg-g0)**2 - 0.37*feh - 0.07*feh**2)

    vmac = float("%.2f" % vmac)
    return vmac


