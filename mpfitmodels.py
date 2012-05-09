#! /usr/bin/env python
import mpfit
import numpy as np
import matplotlib.pyplot as plt

class MPFitModel(object):
    def __init__(self, p):
        # Parinfo: structure where the parameters properties are stored
        self._parinfo = []
        for i in np.arange(len(p)):
            self._parinfo.append({'value':p[i], 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]})
        # Data for the fitting
        self.x = None
        self.y = None
        self.weights = None
    
    def __call__(self, x):
        return self._model_function(x)
    
    def _model_function(self, x, p=None):
        pass
    
    def _model_evaluation_function(self, p, fjac=None):
        # Function that return the weighted deviates
       
        # Parameter values are passed in "p"
        # If fjac==None then partial derivatives should not be
        # computed.  It will always be None if MPFIT is called with default
        # flag.
        model = self._model_function(self.x, p)
        # Non-negative status value means MPFIT should continue, negative means
        # stop the calculation.
        status = 0
        if self.weights != None:
            return([status, (self.y - model)*self.weights])
        else:
            return([status, (self.y - model)])
            
    def fitData(self, x, y, weights=None, parinfo=None):
        self.x = x
        self.y = y
        self.weights = weights
        
        # Parameters' constraints
        if parinfo != None:
            self._parinfo = parinfo
        
        m = mpfit.mpfit(self._model_evaluation_function, parinfo=self._parinfo, quiet=True)
        
        if (m.status <= 0): 
           raise Exception(m.errmsg)
        else:
            for i in np.arange(len(m.params)):
                self._parinfo[i]['value'] = m.params[i]
            # Num iterations: m.niter
            # Uncertainties: m.perror
    
    def residuals(self):
        model = self._model_function(self.x)
        return((self.y - model))
    
    def integrate(self):
        raise NotImplementedError()
    

class GaussianModel(MPFitModel):
    # WARNING: Dot not modify attributes A, sig or mu directly from outside the class!
    def __init__(self, A=-0.025, sig=0.25, mu=0):
        p = [A, sig, mu]
        super(GaussianModel, self).__init__(p)
        # Update public attributes:
        self.A = p[0]
        self.sig = p[1]
        self.mu = p[2]
    
    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p != None:
            # Update public attributes:
            self.A = p[0]
            self.sig = p[1]
            self.mu = p[2]
            # Update internal structure for fitting:
            self._parinfo[0]['value'] = self.A
            self._parinfo[1]['value'] = self.sig
            self._parinfo[2]['value'] = self.mu
        return ((self.A*1.)/np.sqrt(2*np.pi*self.sig**2))*np.exp(-(x-self.mu)**2/(2*self.sig**2))
    
    def fitData(self, x, y, weights=None, parinfo=None):
        super(GaussianModel, self).fitData(x, y, weights, parinfo)
        # Update public attributes with the results:
        self.A = self._parinfo[0]['value']
        self.sig = self._parinfo[1]['value']
        self.mu = self._parinfo[2]['value']
    
    def _make_gauss(self):
        k = self.A / (self.sig * np.sqrt(2*np.pi))
        s = -1.0 / (2 * self.sig * self.sig)
        def f(x):
            return k * np.exp(s * (x - self.mu)*(x - self.mu))
        return f
        
    def integrate(self):
        if self.x == None:
            return 0
        else:
            from scipy.integrate import quad
            # Include 99.97% of the gaussian area
            from_x = self.mu - 3*self.sig
            to_x = self.mu + 3*self.sig
            integral, estimated_error = quad(self._make_gauss(), from_x, to_x)
            return integral

        


class VoigtModel(MPFitModel):
    # WARNING: Dot not modify attributes A, sig, mu or gamma directly from outside the class!
    def __init__(self, A=-0.025, sig=0.25, mu=0, gamma=0.025):
        p = [A, sig, mu, gamma]
        super(VoigtModel, self).__init__(p)
        # Update public attributes:
        self.A = p[0]
        self.sig = p[1]
        self.mu = p[2]
        self.gamma = p[3]
    
    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p != None:
            # Update public attributes:
            self.A = p[0]
            self.sig = p[1]
            self.mu = p[2]
            self.gamma = p[3]
            # Update internal structure for fitting:
            self._parinfo[0]['value'] = self.A
            self._parinfo[1]['value'] = self.sig
            self._parinfo[2]['value'] = self.mu
            self._parinfo[3]['value'] = self.gamma
        if self.sig == 0:
            # Equivalent to a Lorentzian model
            voigt_result = self.A*self.gamma/np.pi/(x*x - 2*x*self.mu+self.mu*self.mu+self.gamma*self.gamma)
        else:
            # Voigt model (Gaussian and Lorentzian)
            from scipy.special import wofz
            w = wofz(((x - self.mu) + 1j*self.gamma)* 2**-0.5/self.sig)
            voigt_result = self.A * w.real*(2*np.pi)**-0.5/self.sig
        return voigt_result
    
    def fitData(self, x, y, weights=None, parinfo=None):
        super(VoigtModel, self).fitData(x, y, weights, parinfo)
        # Update public attributes with the results:
        self.A = self._parinfo[0]['value']
        self.sig = self._parinfo[1]['value']
        self.mu = self._parinfo[2]['value']
        self.gamma = self._parinfo[3]['value']
    
    def _make_voigt(self):
        if self.sig == 0:
            # Equivalent to a Lorentzian model
            k = self.A*self.gamma/np.pi
            s = self.mu*self.mu+self.gamma*self.gamma
            def f(x):
                return k/(x*x - 2*x*self.mu+s)
        else:
            # Voigt model (Gaussian and Lorentzian)
            k = self.A * (2*np.pi)**-0.5/self.sig
            s = 2**-0.5/self.sig
            def f(x):
                from scipy.special import wofz
                return k * wofz(((x - self.mu) + 1j*self.gamma)*s).real
        
        return f
        
    def integrate(self):
        if self.x == None:
            return 0
        else:
            from scipy.integrate import quad
            # Include 99.97% of the gaussian area
            from_x = self.mu - 3*self.sig
            to_x = self.mu + 3*self.sig
            integral, estimated_error = quad(self._make_voigt(), from_x, to_x)
            return integral
    

if __name__ == '__main__':
    ### Full range
    x = np.arange(-10.,10., 1)
    ### Reduced range
    ##x = np.arange(-3.,3., 0.5)

    x_fine = np.arange(-10.,10., 20./1000)

    ############### GAUSSIAN
    ## Generate model data for a Gaussian with param mu and sigma and add noise
    A=-1
    sig=2
    mu=0.5
    gaussian_model = GaussianModel(A, sig, mu)
    y_true = gaussian_model(x)

    #preal = [-1, -2, .5]
    #y_true = gaussian_model(x, preal)

    mu, sigma = 0, 0.7
    y      = y_true + 0.01 * np.random.normal(mu, sigma, len(x) )
    ##err    = 1.0 + 0.01 * np.random.normal(mu, sigma, len(x) )
    err = np.ones(len(x))

    gaussian_model.fitData(x, y)

    print "Fitted pars: "
    print "\tA:\t", gaussian_model.A
    print "\tsig:\t", gaussian_model.sig
    print "\tmu:\t", gaussian_model.mu
    print "RMS: ", np.sqrt(np.sum(np.power(gaussian_model.residuals(), 2)) / len(gaussian_model.residuals()))

    plt.clf()
    plt.plot(x,y,'r', label="Noisy data")
    plt.plot( x_fine, gaussian_model(x_fine), label="Fit" )
    plt.plot( x, y_true, 'g', label="True data" )
    plt.xlabel( "X" )
    plt.ylabel( "Measurement data" )
    plt.title( "Least-squares fit to noisy data using MPFIT" )
    plt.legend()
    plt.show()

    ############### VOIGT
    ## Generate model data for a Gaussian with param mu and sigma and add noise
    A=-1
    sig=2
    mu=0.5
    gamma=0.5
    voigt_model = VoigtModel(A, sig, mu, gamma)
    y_true = voigt_model(x)

    #preal = [-1, -2, .5]
    #y_true = gaussian_model(x, preal)

    mu, sigma = 0, 0.7
    y      = y_true + 0.01 * np.random.normal(mu, sigma, len(x) )
    ##err    = 1.0 + 0.01 * np.random.normal(mu, sigma, len(x) )
    err = np.ones(len(x))

    voigt_model.fitData(x, y)

    print "Fitted pars: "
    print "\tA:\t", voigt_model.A
    print "\tsig:\t", voigt_model.sig
    print "\tmu:\t", voigt_model.mu
    print "\tgamma:\t", voigt_model.gamma
    print "RMS: ", np.sqrt(np.sum(np.power(voigt_model.residuals(), 2)) / len(voigt_model.residuals()))

    plt.clf()
    plt.plot(x,y,'r', label="Noisy data")
    plt.plot( x_fine, voigt_model(x_fine), label="Fit" )
    plt.plot( x, y_true, 'g', label="True data" )
    plt.xlabel( "X" )
    plt.ylabel( "Measurement data" )
    plt.title( "Least-squares fit to noisy data using MPFIT" )
    plt.legend()
    plt.show()


