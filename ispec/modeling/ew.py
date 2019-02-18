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
import sys
import time
from datetime import datetime, timedelta
import numpy as np
import logging

from mpfitmodels import MPFitModel
from ispec.abundances import determine_abundances
from ispec.abundances import determine_abundance_enchancements
from ispec.atmospheres import interpolate_atmosphere_layers, model_atmosphere_is_closest_copy
from common import Constants


class EquivalentWidthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, abundances, teff=5000, logg=3.0, MH=0.0, alpha=0.0, vmic=2.0, adjust_model_metalicity=False, enhance_abundances=True, scale=None):
        self.elements = {}
        #self.elements["1"] = "H"
        #self.elements["2"] = "He"
        self.elements["3"] = "Li"
        self.elements["4"] = "Be"
        self.elements["5"] = "B"
        self.elements["6"] = "C"
        self.elements["7"] = "N"
        self.elements["8"] = "O"
        self.elements["9"] = "F"
        self.elements["10"] = "Ne"
        self.elements["11"] = "Na"
        self.elements["12"] = "Mg"
        self.elements["13"] = "Al"
        self.elements["14"] = "Si"
        self.elements["15"] = "P"
        self.elements["16"] = "S"
        self.elements["17"] = "Cl"
        self.elements["18"] = "Ar"
        self.elements["19"] = "K"
        self.elements["20"] = "Ca"
        self.elements["21"] = "Sc"
        self.elements["22"] = "Ti"
        self.elements["23"] = "V"
        self.elements["24"] = "Cr"
        self.elements["25"] = "Mn"
        self.elements["26"] = "Fe"
        self.elements["27"] = "Co"
        self.elements["28"] = "Ni"
        self.elements["29"] = "Cu"
        self.elements["30"] = "Zn"
        self.elements["31"] = "Ga"
        self.elements["32"] = "Ge"
        self.elements["33"] = "As"
        self.elements["34"] = "Se"
        self.elements["35"] = "Br"
        self.elements["36"] = "Kr"
        self.elements["37"] = "Rb"
        self.elements["38"] = "Sr"
        self.elements["39"] = "Y"
        self.elements["40"] = "Zr"
        self.elements["41"] = "Nb"
        self.elements["42"] = "Mo"
        self.elements["43"] = "Tc"
        self.elements["44"] = "Ru"
        self.elements["45"] = "Rh"
        self.elements["46"] = "Pd"
        self.elements["47"] = "Ag"
        self.elements["48"] = "Cd"
        self.elements["49"] = "In"
        self.elements["50"] = "Sn"
        self.elements["51"] = "Sb"
        self.elements["52"] = "Te"
        self.elements["53"] = "I"
        self.elements["54"] = "Xe"
        self.elements["55"] = "Cs"
        self.elements["56"] = "Ba"
        self.elements["57"] = "La"
        self.elements["58"] = "Ce"
        self.elements["59"] = "Pr"
        self.elements["60"] = "Nd"
        self.elements["61"] = "Pm"
        self.elements["62"] = "Sm"
        self.elements["63"] = "Eu"
        self.elements["64"] = "Gd"
        self.elements["65"] = "Tb"
        self.elements["66"] = "Dy"
        self.elements["67"] = "Ho"
        self.elements["68"] = "Er"
        self.elements["69"] = "Tm"
        self.elements["70"] = "Yb"
        self.elements["71"] = "Lu"
        self.elements["72"] = "Hf"
        self.elements["73"] = "Ta"
        self.elements["74"] = "W"
        self.elements["75"] = "Re"
        self.elements["76"] = "Os"
        self.elements["77"] = "Ir"
        self.elements["78"] = "Pt"
        self.elements["79"] = "Au"
        self.elements["80"] = "Hg"
        self.elements["81"] = "Tl"
        self.elements["82"] = "Pb"
        self.elements["83"] = "Bi"
        self.elements["84"] = "Po"
        self.elements["85"] = "At"
        self.elements["86"] = "Rn"
        self.elements["87"] = "Fr"
        self.elements["88"] = "Ra"
        self.elements["89"] = "Ac"
        self.elements["90"] = "Th"
        self.elements["91"] = "Pa"
        self.elements["92"] = "U"
        self.elements["101"] = "Md"
        self.elements["106"] = "Sg"
        self.elements["107"] = "Bh"
        self.elements["108"] = "Hs"
        self.elements["112"] = "Cn"
        self.elements["113"] = "Uut"
        self.elements["114"] = "Uuq"

        self.modeled_layers_pack = modeled_layers_pack
        self.abundances = abundances
        self.enhance_abundances = enhance_abundances
        self.scale = scale
        self.adjust_model_metalicity = adjust_model_metalicity
        self.lines_for_teff = None
        self.lines_for_vmic = None
        #
        self.calculation_time = 0
        self.cache = {}
        self.m1 = None
        self.c1 = None
        self.m2 = None
        self.c2 = None
        self.fe1 = None
        self.fe2 = None
        self.fe1_std = None
        self.fe2_std = None
        self.fe1_filter = None
        self.fe2_filter = None
        p = [teff, logg, vmic]
        self._MH = MH
        self._eMH = 0.0
        self._alpha = alpha
        self._ealpha = 0.0

        ranges = modeled_layers_pack[7]
        MH_range = ranges['MH']
        self.min_MH = np.min(MH_range)
        self.max_MH = np.max(MH_range)
        #
        super(EquivalentWidthModel, self).__init__(p)



    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        key = "%.0f %.2f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.alpha(), self.vmic())
        if self.cache.has_key(key):
            hit_cache = True
            if not self.quiet:
                print "Cache:", key
            self.last_final_values = self.cache[key]
            spec_abund, absolute_abund, x_over_h, x_over_fe = self.cache[key]
        else:
            hit_cache = False
            if not self.quiet:
                print "Generating:", key
            # Optimization to avoid too small changes in parameters or repetition
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, {'teff':self.teff(), 'logg':self.logg(), 'MH':self.MH(), 'alpha':self.alpha()}, code=self.code)
            if self.fe1_filter is None or self.fe2_filter is None:
                ignore = np.ones(len(self.linemasks)) # Do not ignore any line since it's the first execution and it has not been done any selection
            else:
                ignore = np.zeros(len(self.linemasks))
                ignore[np.where(np.logical_or(self.fe1_filter, self.fe2_filter))[0]] = 1.0 # Do not ignore selected fe1/2 lines

            spec_abund, absolute_abund, x_over_h, x_over_fe = determine_abundances(atmosphere_layers, \
                    self.teff(), self.logg(), self.MH(), self.alpha(), self.linemasks, self.abundances, microturbulence_vel = self.vmic(), \
                    ignore=ignore, verbose=0, code=self.code, tmp_dir=self.tmp_dir)


            if 'EW_absolute_abund_median' in self.linemasks.dtype.names:
                # Instead of the literature solar abundance, use the solar abundance determined by iSpec (differencial analysis)
                differential_x_over_h = absolute_abund - self.linemasks['EW_absolute_abund_median']
                differential_feh = np.median(differential_x_over_h[self.linemasks['element'] == "Fe 1"])
                differential_x_over_fe = differential_x_over_h - differential_feh
                x_over_h = differential_x_over_h
                x_over_fe = differential_x_over_fe

            self.cache[key] = (spec_abund, absolute_abund, x_over_h, x_over_fe)

        # First iteration
        if self.fe1_filter is None or self.fe2_filter is None:
            #self.outliers_detection = None # Don't identify and filter outliers
            self.select_good_lines(x_over_h)

        values_to_evaluate = []
        fitted_lines_params = []
        selected_x_over_h = []

        import statsmodels.api as sm
        ### Temperature
        ## y = mx + c
        #x = self.linemasks['lower_state_eV'][self.fe1_filter]
        #y = x_over_h[self.fe1_filter]
        x = self.linemasks['lower_state_eV'][np.logical_or(self.fe1_filter, self.fe2_filter)]
        y = x_over_h[np.logical_or(self.fe1_filter, self.fe2_filter)]
        unknown = np.isnan(y)
        x = x[~unknown]
        y = y[~unknown]
        if len(x) < 2:
            raise Exception("Not enough abundances were calculated")
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m1 = linear_model.params[0]
        self.c1 = linear_model.params[1]
        self.fe1 = np.nanmedian(x_over_h[self.fe1_filter])
        self.fe1_std = np.nanstd(x_over_h[self.fe1_filter])
        ##self.fe1 = np.median(linear_model.fittedvalues)
        ##self.fe1_std = np.std(linear_model.fittedvalues)
        #print "Fe 1", np.median(linear_model.fittedvalues), np.nanmedian(x_over_h[self.fe1_filter])
        #print "    ", np.std(linear_model.fittedvalues), np.nanstd(x_over_h[self.fe1_filter])
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.scatter(x, y)
        #plt.plot(x, self.m1*x + self.c1)
        #plt.xlabel("Lower state (eV)")
        #plt.ylabel("[Fe/H]")
        #plt.grid()
        #plt.show()

        ### Vmic
        ## y = mx + c
        #x = self.linemasks['ewr'][self.fe1_filter]
        #y = x_over_h[self.fe1_filter]
        x = self.linemasks['ewr'][np.logical_or(self.fe1_filter, self.fe2_filter)]
        y = x_over_h[np.logical_or(self.fe1_filter, self.fe2_filter)]
        unknown = np.isnan(y)
        x = x[~unknown]
        y = y[~unknown]
        if len(x) < 2:
            raise Exception("Not enough abundances were calculated")
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m2 = linear_model.params[0]
        self.c2 = linear_model.params[1]
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.scatter(x, y)
        #plt.plot(x, self.m2*x + self.c2)
        #plt.xlabel("Reduced EW")
        #plt.ylabel("[Fe/H]")
        #plt.grid()
        #plt.show()

        ### Fe2
        ## y = mx + c
        x = self.linemasks['ewr'][self.fe2_filter]
        if len(x) > 1:
            y = x_over_h[self.fe2_filter]
            unknown = np.isnan(y)
            x = x[~unknown]
            y = y[~unknown]
            if len(x) < 2:
                raise Exception("Not enough abundances were calculated")
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
            self.fe2 = np.nanmedian(x_over_h[self.fe2_filter])
            self.fe2_std = np.nanstd(x_over_h[self.fe2_filter])
            ##self.fe2 = np.median(linear_model.fittedvalues)
            ##self.fe2_std = np.std(linear_model.fittedvalues)
            #print "Fe 2", np.median(linear_model.fittedvalues), np.nanmedian(x_over_h[self.fe2_filter])
            #print "    ", np.std(linear_model.fittedvalues), np.nanstd(x_over_h[self.fe2_filter])
            #import matplotlib.pyplot as plt
            #plt.scatter(x, y)
            #plt.plot(x, m2*x + c2)
            #plt.show()
        else:
            self.fe2 = np.nanmedian(x_over_h[self.fe2_filter])
            self.fe2_std = np.nanstd(x_over_h[self.fe2_filter])

        ## Gravity
        abundance_diff = self.fe1 - self.fe2
        abundance_diff2 = self.MH() - self.fe1

        # Rounded to 3 and 2 decimals (using string convertion works better than np.round)
        values_to_evaluate.append(float("%.6f" % self.m1))
        values_to_evaluate.append(float("%.6f" % self.m2))
        values_to_evaluate.append(float("%.6f" % abundance_diff))
        if self.adjust_model_metalicity:
            values_to_evaluate.append(float("%.2f" % abundance_diff2))
        residuals = np.asarray(values_to_evaluate) - self.y

        if not hit_cache:
            print " # Element:                   Fe 1 / Fe 2\n",
            print "   Teff/Vmic slopes:            %.6f %.6f" % (self.m1, self.m2)
            print "   Abundances diff:             %.6f" % abundance_diff
            print "   Abundances diff with model:  %.6f" % abundance_diff2
            print "   Abundances stdev:            %.6f %.6f" % (np.std(x_over_h[self.fe1_filter]), np.std(x_over_h[self.fe2_filter]))
            print "   Abundances median:           %.6f %.6f" % (np.median(self.fe1), np.median(self.fe2))
            print " - Chisq:                       %.10g" % np.sum((self.weights*residuals)**2)

        fitted_lines_params.append(self.m1)
        fitted_lines_params.append(self.c1)
        fitted_lines_params.append(self.m2)
        fitted_lines_params.append(self.c2)
        selected_x_over_h.append(self.fe1_filter.copy())
        selected_x_over_h.append(self.fe2_filter.copy())
        values_to_evaluate = np.asarray(values_to_evaluate)
        if model_atmosphere_is_closest_copy(self.modeled_layers_pack, {'teff':self.teff(), 'logg':self.logg(), 'MH':self.MH(), 'alpha':self.alpha(), 'vmic': self.vmic()}):
            # Penalize these cases
            values_to_evaluate *= 100.
        self.last_final_values = (values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params)

        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        return values_to_evaluate.copy()

    # Default procedure to be called every iteration.  It simply prints
    # the parameter values.
    import scipy
    blas_enorm32, = scipy.linalg.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float32))
    blas_enorm64, = scipy.linalg.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float64))
    def defiter(self, fcn, x, iter, fnorm=None, functkw=None,
                       quiet=0, iterstop=None, parinfo=None,
                       format=None, pformat='%.10g', dof=1):

        if quiet:
            return
        if fnorm is None:
            [status, fvec] = fcn(x, fjac=None, **functkw)
            # If the returned fvec has more than four bits I assume that we have
            # double precision
            # It is important that the machar is determined by the precision of
            # the returned value, not by the precision of the input array
            if np.array([fvec]).dtype.itemsize>4:
                self.blas_enorm = mpfit.blas_enorm64
            else:
                self.blas_enorm = mpfit.blas_enorm32
            fnorm = self.enorm(fvec)**2

        # Determine which parameters to print
        nprint = len(x)
        print "*Iter ", ('%6i' % iter),"   CHI-SQUARE = ",('%.10g' % fnorm)," DOF = ", ('%i' % dof)
        for i in range(nprint):
            if (parinfo is not None) and (parinfo[i].has_key('parname')):
                p = '   ' + parinfo[i]['parname'] + ' = '
            else:
                p = '   P' + str(i) + ' = '
            if (parinfo is not None) and (parinfo[i].has_key('mpprint')):
                iprint = parinfo[i]['mpprint']
            else:
                iprint = 1
            if iprint:
                print p + (pformat % x[i]) + '  '

        ##### Metallicity
        #self._MH = self.fe1
        #self._MH = np.min((self._MH, self.max_MH))
        #self._MH = np.max((self._MH, self.min_MH))
        #self._eMH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])


        ##### Lines
        #values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        #self.select_good_lines(x_over_h) # Modifies self.lines_for_teff and self.lines_for_vmic

        return 0


    def select_good_lines(self, x_over_h, strict=True):
        """
            Modifies self.fe1_filter and self.fe2_filter
        """
        # Out of range
        unknown = np.isnan(x_over_h)
        bad = np.logical_or(x_over_h > 1.0, x_over_h < -5)
        bad = np.logical_or(bad, unknown)
        #### Line selection
        fe1_filter = self.linemasks['element'] == "Fe 1"
        fe2_filter = self.linemasks['element'] == "Fe 2"

        if self.outliers_detection not in ['robust', 'sigma_clipping']:
            strict = False
        else:
            scrict = True

        if strict and len(np.where(~bad)[0]) > 1:
            # Outliers
            import statsmodels.api as sm
            # Do not use NaN values but allow the use of out of range values since
            # they could come back to normal values later on
            x = self.linemasks['lower_state_eV'][~unknown]
            y = x_over_h[~unknown]
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            if self.outliers_detection == "robust":
                # RLM (Robust least squares)
                # Huber's T norm with the (default) median absolute deviation scaling
                # - http://en.wikipedia.org/wiki/Huber_loss_function
                # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
                huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
                linear_model = huber_t.fit()
                reject_filter1 = linear_model.weights < self.outliers_weight_limit
                #reject_filter1 = np.logical_or(reject_filter1, bad) # Done later
            elif self.outliers_detection == "sigma_clipping":
                linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
                m1 = linear_model.params[0]
                c1 = linear_model.params[1]
                corrected_y = y - (m1*x + c1)
                sigma = np.std(corrected_y)
                reject_filter1 = np.logical_or(corrected_y > + self.sigma_level*sigma, corrected_y < -self.sigma_level*sigma)
            else:
                logging.warn("Unknown outlier detection technique: %s" % (self.outliers_detection))
                reject_filter1 = np.asarray([False]*len(x))
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['lower_state_eV'], x_over_h)
            #plt.scatter(self.linemasks['lower_state_eV'][reject_filter1], x_over_h[reject_filter1], color="red")
            #m1 = linear_model.params[0]
            #c1 = linear_model.params[1]
            #plt.plot(x, m1*x + c1, color="green")
            #plt.xlabel("Lower excitation energies (eV)")
            #plt.ylabel("[Fe 1/H] (dex)")
            #plt.grid()
            #plt.savefig("excitation_equilibrium.png")
            #plt.savefig("excitation_equilibrium.eps")
            #plt.savefig("excitation_equilibrium.pdf")
            #plt.show()

            # Outliers
            import statsmodels.api as sm
            # Do not use NaN values but allow the use of out of range values since
            # they could come back to normal values later on
            x = self.linemasks['ewr'][~unknown]
            y = x_over_h[~unknown]
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            if self.outliers_detection == "robust":
                # RLM (Robust least squares)
                # Huber's T norm with the (default) median absolute deviation scaling
                # - http://en.wikipedia.org/wiki/Huber_loss_function
                # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
                huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
                linear_model = huber_t.fit()
                reject_filter2 = linear_model.weights < self.outliers_weight_limit
                #reject_filter2 = np.logical_or(reject_filter2, bad) # Done later
            elif self.outliers_detection == "sigma_clipping":
                linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
                m2 = linear_model.params[0]
                c2 = linear_model.params[1]
                corrected_y = y - (m1*x + c1)
                sigma = np.std(corrected_y)
                reject_filter2 = np.logical_or(corrected_y > self.sigma_level*sigma, corrected_y < -self.sigma_level*sigma)
            else:
                logging.warn("Unknown outlier detection technique: %s" % (self.outliers_detection))
                reject_filter2 = np.asarray([False]*len(x))
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['ewr'], x_over_h)
            #plt.scatter(self.linemasks['ewr'][reject_filter2], x_over_h[reject_filter2], color="red")
            #plt.show()

            reject_filter_tmp = np.logical_or(reject_filter1, reject_filter2)
            known_idx = np.where(~unknown)[0]

            # unknown abundances where excluded, recover them and keep the array size
            # coherent
            reject_filter = unknown.copy()
            reject_filter[known_idx] = reject_filter_tmp

            # Discard bad lines and outliers
            clean_fe1_filter = np.logical_and(~reject_filter, fe1_filter)
            # Ensure that there are at least some lines
            if len(np.where(clean_fe1_filter)[0]) <= 1:
                # Discard only bad lines
                clean_fe1_filter = np.logical_and(~bad, fe1_filter)
                if len(np.where(clean_fe1_filter)[0]) <= 1:
                    clean_fe1_filter = fe1_filter
            if len(np.where(clean_fe1_filter)[0]) <= 1:
                raise Exception("Not enought lines for Fe 1 (%i lines)" % len(np.where(clean_fe1_filter)[0]))
            else:
                self.fe1_filter = clean_fe1_filter

            # Discard bad lines and outliers
            clean_fe2_filter = np.logical_and(~reject_filter, fe2_filter)
            # Ensure that there are at least some lines
            if len(np.where(clean_fe2_filter)[0]) <= 1:
                # Discard only bad lines
                clean_fe2_filter = np.logical_and(~bad, fe2_filter)
                if len(np.where(clean_fe2_filter)[0]) <= 1:
                    clean_fe2_filter = fe2_filter

            ## Discard ONLY bad lines for Fe 2
            #clean_fe2_filter = np.logical_and(~bad, fe2_filter)
            ## Ensure that there are at least some lines
            #if len(np.where(clean_fe2_filter)[0]) >= 1:
                #clean_fe2_filter = fe2_filter

            if len(np.where(clean_fe2_filter)[0]) <= 1:
                raise Exception("Not enought lines for Fe 1 (%i lines)" % len(np.where(clean_fe2_filter)[0]))
            else:
                self.fe2_filter = clean_fe2_filter
        else:
            ##### ACCEPT all
            if len(np.where(~bad & fe1_filter)[0]) > 0:
                self.fe1_filter = np.logical_and(fe1_filter, np.logical_not(bad))
            else:
                self.fe1_filter = fe1_filter
            if len(np.where(~bad & fe2_filter)[0]) > 0:
                self.fe2_filter = np.logical_and(fe2_filter, np.logical_not(bad))
            else:
                self.fe2_filter = fe2_filter
        print " > Selected Fe 1 lines for teff:", len(np.where(self.fe1_filter)[0]), "of", len(np.where(fe1_filter)[0])
        print " > Selected Fe 2 lines for vmic:", len(np.where(self.fe2_filter)[0]), "of", len(np.where(fe2_filter)[0])



    def fitData(self, linemasks, outliers_detection='robust', sigma_level=3, outliers_weight_limit=0.90, parinfo=None, max_iterations=20, quiet=True, code="spectrum", tmp_dir=None):
        base = 5
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        code = code.lower()
        if code not in ['spectrum', 'turbospectrum', 'moog', 'width']:
            raise Exception("Unknown radiative transfer code: %s" % (code))


        self.code = code
        self.tmp_dir = tmp_dir
        self.outliers_detection = outliers_detection
        self.sigma_level = sigma_level
        self.outliers_weight_limit = outliers_weight_limit


        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.linemasks = linemasks
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 0.0   # Not active: Residuals are limited between -1.0 and 1.0 (np.tanh(residuals/1.0))
        #chisq_limit = 4.0e-4 # 0.0004 = np.sum(np.asarray([0.01, 0.01, 0.01, 0.01])**2))
        if code == "moog":
            chisq_limit = None
        else:
            chisq_limit = 3 # = np.sum((np.asarray([0.01, 0.01, 0.01])*100)**2) # weight 100)

        _t0 = default_timer()

        #index = np.asarray([0, 1, 2])
        #index = np.arange(len(linemasks))
        if self.adjust_model_metalicity:
            index = np.arange(4) # 4 values: zero slopes and zero difference between element1 and element2, difference with model
        else:
            index = np.arange(3) # 3 values: zero slopes and zero difference between element1 and element2
        target_values = np.zeros(len(index))
        weights = np.ones(len(index)) * 100
        #weights = np.asarray([3,1,2])
        #weights = np.asarray([100,100,1])
        #weights = np.asarray([1000,1,100])
        #weights = np.asarray([100,1,1000])
        #weights = np.asarray([1,1000,1])
        super(EquivalentWidthModel, self).fitData(index, target_values, weights=weights, parinfo=parinfo, chisq_limit=chisq_limit, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet, iterfunct=self.defiter)

        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        residuals = values_to_evaluate - target_values
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
        # Unweighted
        self.chisq = np.sum((residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof
        # Weighted
        self.wchisq = np.sum((weights * residuals)**2)
        self.reduced_wchisq = self.wchisq / self.m.dof

        #self.cache = {}

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self): return self._parinfo[0]['value']
    def logg(self): return self._parinfo[1]['value']
    def vmic(self): return self._parinfo[2]['value']

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def evmic(self): return self.m.perror[2]

    def MH(self): return self._parinfo[3]['value']
    def eMH(self): return self.m.perror[3]

    #def MH(self): return self._MH
    #def eMH(self): return self._eMH

    def alpha(self):
        if self.enhance_abundances:
            alpha_enhancement = determine_abundance_enchancements(self.MH(), scale=self.scale)
        else:
            alpha_enhancement = self._parinfo[4]['value']
        return alpha_enhancement
    def ealpha(self): return self.m.perror[4]

    def print_solution(self):
        # Calculate MH
        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        #MH = np.median(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        #MH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        MH = self.fe1
        eMH = self.fe1_std

        header = "%8s\t%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","alpha","vmic")
        #solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), self.MH(), self.alpha(), self.vmic())
        #errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), self.eMH(), self.ealpha(), self.evmic())
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), MH, self.alpha(), self.vmic())
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), eMH, self.ealpha(), self.evmic())

        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","wchisq","rwchisq","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f" % (self.m.dof, self.m.niter, self.m.nfev, self.wchisq, self.reduced_wchisq, self.chisq, self.reduced_chisq, self.rms)
        if model_atmosphere_is_closest_copy(self.modeled_layers_pack, {'teff':self.teff(), 'logg':self.logg(), 'MH':MH, 'alpha':self.alpha(), 'vmic': self.vmic()}):
            print ""
            print "WARNING: Model atmosphere used for the final solution was not interpolated, it is a copy of the closest model."
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status


def model_spectrum_from_ew(linemasks, modeled_layers_pack, abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, free_params=["teff", "logg", "vmic"], adjust_model_metalicity=False, enhance_abundances=True, scale=None, max_iterations=20, outliers_detection='robust', sigma_level=3, outliers_weight_limit=0.90, code="spectrum", tmp_dir=None):
    """
    - outlier_detection:
        - 'robust': Fit a robust least square linear model, outliers_weight_limit will be use as a threshold. If it is set to zero, no outliers are filtered.
        - 'sigma_clipping': Fit a tradition least square linear model and filter X times the standard deviation (sigma_level)
    - If enhance_abundances is True, alpha elements and CNO abundances will be scaled
      depending on the metallicity.
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'width']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    ranges = modeled_layers_pack[7]
    teff_range = ranges['teff']
    logg_range = ranges['logg']
    MH_range = ranges['MH']
    alpha_range = ranges.get('alpha', (-1.5, 1.5)) # Default if 'alpha' is not a free parameter for atmosphere interpolation

    # Do not allow users to set free MH in free_params to avoid confusions
    # because metallicity is always free in this method, what we make by including MH in free_params
    # is turning on the adjustment in the metallicity models
    if "MH" in free_params or "mh" in free_params:
        raise Exception("Metallicity cannot be a free parameter!")

    if "alpha" in free_params:
        raise Exception("Alpha enhancement cannot be a free parameter!")

    if adjust_model_metalicity:
        free_params = free_params[:] # copy to avoid modifying user's free_params list
        free_params.append("MH")

    parinfo = __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, teff_range, logg_range, MH_range, alpha_range, free_params, adjust_model_metalicity=adjust_model_metalicity)


    EW_model = EquivalentWidthModel(modeled_layers_pack, abundances, MH=initial_MH, alpha=initial_alpha, adjust_model_metalicity=adjust_model_metalicity, \
                                        enhance_abundances=enhance_abundances, scale=scale)

    lfilter = linemasks['element'] == "Fe 1"
    lfilter = np.logical_or(lfilter, linemasks['element'] == "Fe 2")
    linemasks = linemasks[lfilter]
    EW_model.fitData(linemasks, parinfo=parinfo, max_iterations=max_iterations, quiet=False, outliers_detection=outliers_detection, sigma_level=sigma_level, outliers_weight_limit=outliers_weight_limit, code=code, tmp_dir=tmp_dir)
    print "\n"
    EW_model.print_solution()

    # Calculate MH
    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values
    MH = EW_model.fe1
    eMH = EW_model.fe1_std

    used_linemasks = EW_model.linemasks[np.logical_or(EW_model.fe1_filter, EW_model.fe2_filter)]

    # Collect information to be returned
    params = {}
    params['teff'] = EW_model.teff()
    params['logg'] = EW_model.logg()
    params['MH'] = MH
    params['alpha'] = EW_model.alpha()
    params['vmic'] = EW_model.vmic()

    errors = {}
    errors['teff'] = EW_model.eteff()
    errors['logg'] = EW_model.elogg()
    errors['MH'] = eMH
    errors['alpha'] = EW_model.ealpha()
    errors['vmic'] = EW_model.evmic()

    status = {}
    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values
    # Save parameters (only for Fe, if there are more elements they will not be saved)
    status['slope_excitation_potential'] = values_to_evaluate[0]
    status['slope_ewr'] = values_to_evaluate[1]
    status['abundance_diff'] = values_to_evaluate[2]
    status['fe1_lines'] = len(np.where(selected_x_over_h[0])[0])
    status['fe2_lines'] = len(np.where(selected_x_over_h[1])[0])
    status['model_MH'] = EW_model.MH()

    status['days'] = EW_model.calculation_time.day-1
    status['hours'] = EW_model.calculation_time.hour
    status['minutes'] = EW_model.calculation_time.minute
    status['seconds'] = EW_model.calculation_time.second
    status['dof'] = EW_model.m.dof
    status['error'] = EW_model.m.errmsg
    status['rms'] = EW_model.rms
    status['chisq'] = EW_model.chisq
    status['rchisq'] = EW_model.reduced_chisq
    status['niter'] = EW_model.m.niter
    status['nsynthesis'] = EW_model.m.nfev
    status['status'] = EW_model.m.status

    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks


def __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, teff_range, logg_range, MH_range, alpha_range, free_params, adjust_model_metalicity=False):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 5
    free_params = [param.lower() for param in free_params]
    #parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base)]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(base)]
    ##
    min_teff = np.min(teff_range)
    max_teff = np.max(teff_range)
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = Constants.EW_STEP_TEFF # For auto-derivatives
    #parinfo[0]['mpside'] = 2
    #parinfo[0]['mpmaxstep'] = parinfo[0]['step'] * 1.5
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [min_teff, max_teff]
    if parinfo[0]['value'] > parinfo[0]['limits'][1] or parinfo[0]['value'] < parinfo[0]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[0]['parname'], parinfo[0]['value'], parinfo[0]['limits'][0], parinfo[0]['limits'][1]))
    ##
    min_logg = np.min(logg_range)
    max_logg = np.max(logg_range)
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = Constants.EW_STEP_LOGG # For auto-derivatives
    #parinfo[1]['mpside'] = 2
    #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    #parinfo[1]['mpmaxstep'] = parinfo[1]['step'] * 1.5
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [min_logg, max_logg]
    if parinfo[1]['value'] > parinfo[1]['limits'][1] or parinfo[1]['value'] < parinfo[1]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[1]['parname'], parinfo[1]['value'], parinfo[1]['limits'][0], parinfo[1]['limits'][1]))
    #
    parinfo[2]['parname'] = "Vmic"
    parinfo[2]['value'] = initial_vmic
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = Constants.EW_STEP_VMIC # For auto-derivatives
    #parinfo[2]['mpside'] = 2
    #parinfo[2]['mpmaxstep'] = parinfo[2]['step'] * 2.0
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., 50.0]
    if parinfo[2]['value'] > parinfo[2]['limits'][1] or parinfo[2]['value'] < parinfo[2]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[2]['parname'], parinfo[2]['value'], parinfo[2]['limits'][0], parinfo[2]['limits'][1]))
    ##
    min_MH = np.min(MH_range)
    max_MH = np.max(MH_range)
    #
    parinfo[3]['parname'] = "MH"
    parinfo[3]['value'] = initial_MH
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = Constants.EW_STEP_MH # For auto-derivatives
    #parinfo[3]['mpside'] = 2
    #if not parinfo[3]['fixed']:
        #parinfo[3]['mpmaxstep'] = parinfo[3]['step'] * 1.5
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [min_MH, max_MH]
    if parinfo[3]['value'] > parinfo[3]['limits'][1] or parinfo[3]['value'] < parinfo[3]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[3]['parname'], parinfo[3]['value'], parinfo[3]['limits'][0], parinfo[3]['limits'][1]))
    #
    min_alpha = np.min(alpha_range)
    max_alpha = np.max(alpha_range)
    #
    parinfo[4]['parname'] = "alpha"
    parinfo[4]['value'] = initial_alpha
    parinfo[4]['fixed'] = not parinfo[4]['parname'].lower() in free_params
    parinfo[4]['step'] = Constants.EW_STEP_ALPHA # For auto-derivatives
    #parinfo[4]['mpside'] = 2
    #if not parinfo[4]['fixed']:
        #parinfo[4]['mpmaxstep'] = parinfo[4]['step'] * 1.5
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [min_alpha, max_alpha]
    if parinfo[4]['value'] > parinfo[4]['limits'][1] or parinfo[4]['value'] < parinfo[4]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[4]['parname'], parinfo[4]['value'], parinfo[4]['limits'][0], parinfo[4]['limits'][1]))

    return parinfo

