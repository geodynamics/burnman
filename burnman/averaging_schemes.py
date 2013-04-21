# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings

class averaging_scheme:
    """  base class for different averaging schemes """
    def average(self, fractions, volumes, X):
        """
        Average the quantities X of the phases with given fractions 
        (adding up to 1) and volumes. Each parameter is an array
        of the same size.
        """
        raise NotImplementedError("")



def voigt_average_function(phase_volume,X):
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    X_voigt = sum( V_i[i]/V_tot * X[i] for i in it)
    return X_voigt

def reuss_average_function(phase_volume,X):
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    if (min(X)<=0.0):
        X_reuss = 0.0
        warnings.warn("Oops, called reuss_average with Xi<0!")
    else:
        X_reuss = 1./sum(  V_i[i]/V_tot* 1./X[i] for i in it)
    return X_reuss

def voigt_reuss_hill_function(phase_volume,X):
    X_vrh = (voigt_average_function(phase_volume,X) + reuss_average_function(phase_volume,X))/2.
    return X_vrh


class voigt_reuss_hill(averaging_scheme):
    """ Compute Voigt-Reuss-Hill average.
    # This is based on Matas 2007, Appendix C
    Calculate the voigt, reuss and hill averages of the phase assemblage for 
    a certain property X.
    
    Input: phase_volume: array of n volumes, X: array of n properties
    Returns: mixture of property X
    
    Source: Matas 2007, Appendix D """
    
    def average(self, fractions, volumes, X):
        return voigt_reuss_hill_function(volumes, X)

class linear(averaging_scheme):
    """ simple linear averaging """
    def average(self, fractions, volumes, X):
        return sum( [ fractions[i]*X[i] for i in range(len(fractions)) ] )

class voigt(averaging_scheme):
    """ voigt averaging """
    def average(self, fractions, volumes, X):
        return voigt_average_function(volumes,X)

class reuss(averaging_scheme):
    """ reuss averaging """
    def average(self, fractions, volumes, X):
        return reuss_average_function(volumes,X)
    


