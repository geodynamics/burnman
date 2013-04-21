# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import seismic
import warnings

def voigt_reuss_hill(phase_volume,X):
    
    """ Compute Voigt-Reuss-Hill average.
    # This is based on Matas 2007, Appendix C
    Calculate the voigt, reuss and hill averages of the phase assemblage for 
    a certain property X.
    
    Input: phase_volume: array of n volumes, X: array of n properties
    Returns: mixture of property X
    
    Source: Matas 2007, Appendix D """

    X_vrh = (voigt_average(phase_volume,X)+reuss_average(phase_volume,X))/2.
    return X_vrh

def voigt_average(phase_volume,X):
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    X_voigt = sum( V_i[i]/V_tot * X[i] for i in it)
    return X_voigt

def reuss_average(phase_volume,X):
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    if (min(X)<=0.0):
        X_reuss = 0.0
        warnings.warn("Oops, called reuss_average with Xi<0!")
    else:
        X_reuss = 1./sum(  V_i[i]/V_tot* 1./X[i] for i in it)
    return X_reuss

