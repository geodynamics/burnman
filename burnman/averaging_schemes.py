# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings

class averaging_scheme:
    """ 
    Base class for different averaging schemes
    New averagimg schemes should define the functions
    average_bulk_moduli and average_shear_moduli, as
    specified here.
    """
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the bulk moduli of an assemblage, given
        a list of volume fractions, bulk moduli, and shear moduli
        Returns: a single modulus
        """
        raise NotImplementedError("")
    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        """
        Average the shear moduli of an assemblage, given
        a list of volume fractions, bulk moduli, and shear moduli
        Returns: a single modulus
        """
        raise NotImplementedError("")
    def average_density(self, volumes, densities):
        """
        Average the densities of the rock, given a list of volume
        fractions and densitites.  This is the only
        one that is implemented in the base class, as it should
        not be controvsersial... :)
        Returns: a single density
        """
        total_mass = np.sum(np.array(densities)*np.array(volumes))
        total_vol = np.sum(np.array(volumes)) #should sum to one
        assert(total_vol < 1.0001 and total_vol > 0.9999)
        density = total_mass/total_vol
        return density
         
        

class voigt_reuss_hill(averaging_scheme):
    """ Compute Voigt-Reuss-Hill average.
    # This is based on Matas 2007, Appendix C
    Calculate the voigt, reuss and hill averages of the phase assemblage for 
    a certain property X.
    
    Input: phase_volume: array of n volumes, X: array of n properties
    Returns: mixture of property X
    
    Source: Matas 2007, Appendix D """
    
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        return voigt_reuss_hill_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        return voigt_reuss_hill_function(volumes, shear_moduli)


class voigt(averaging_scheme):
    """ Compute Voigt (iso-strain) bound. """
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        return voigt_average_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        return voigt_average_function(volumes, shear_moduli)


class reuss(averaging_scheme):
    """ Compute Reuss (iso-stress) bound."""
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        return reuss_average_function(volumes, bulk_moduli)

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        return reuss_average_function(volumes, shear_moduli)
    

class hashin_shtrikman_upper(averaging_scheme):
    """
    Lower of the two Hashin-Shtrikman bounds.  
    Implements Formulas from Watt et al (1976)
    """
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):

      K_n = max(bulk_moduli)
      G_n = max(shear_moduli)

      vol_frac = volumes/sum(volumes)
 
      alpha_n = -3. / (3.*K_n+4.*G_n)
      A_n = 0
      for i in range(len(vol_frac)):
          if  bulk_moduli[i] != K_n:
              A_n = A_n + vol_frac[i]/(1./(bulk_moduli[i] - K_n) - alpha_n)

      K_upper = K_n + A_n/(1. + alpha_n*A_n)
      return K_upper

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):

      K_n = max(bulk_moduli)
      G_n = max(shear_moduli)

      vol_frac = volumes/sum(volumes)
 
      beta_n = -3. * (K_n + 2.*G_n)  / (5.*G_n * (3.*K_n+4.*G_n))
      B_n = 0
      for i in range(len(vol_frac)):
          if  shear_moduli[i] != G_n:
              B_n = B_n + vol_frac[i]/(1./(2.*(shear_moduli[i] - G_n)) - beta_n)

      G_upper = G_n + (0.5)*B_n/(1. + beta_n*B_n)
      return G_upper

class hashin_shtrikman_lower(averaging_scheme):
    """
    Lower of the two Hashin-Shtrikman bounds.  
    Implements Formulas from Watt et al (1976)
    """
    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):

      K_1 = min(bulk_moduli)
      G_1 = min(shear_moduli)

      vol_frac = volumes/sum(volumes)
 
      alpha_1 = -3. / (3.*K_1+4.*G_1)
      A_1 = 0
      for i in range(len(vol_frac)):
          if  bulk_moduli[i] != K_1:
              A_1 = A_1 + vol_frac[i]/(1./(bulk_moduli[i] - K_1) - alpha_1)

      K_lower = K_1 + A_1/(1. + alpha_1*A_1)
      return K_lower

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):

      K_1 = min(bulk_moduli)
      G_1 = min(shear_moduli)

      vol_frac = volumes/sum(volumes)
 
      beta_1 = -3. * (K_1 + 2.*G_1)  / (5.*G_1 * (3.*K_1+4.*G_1))
      B_1 = 0
      for i in range(len(vol_frac)):
          if  shear_moduli[i] != G_1:
              B_1 = B_1 + vol_frac[i]/(1./(2.*(shear_moduli[i] - G_1)) - beta_1)

      G_lower = G_1 + (0.5)*B_1/(1. + beta_1*B_1)
      return G_lower

class hashin_shtrikman_average(averaging_scheme):
    """
    Arithmetic mean of the upper and lower
    Hashin-Shtrikman bounds
    """
    def __init__(self):
        self.upper = hashin_shtrikman_upper()
        self.lower = hashin_shtrikman_lower()

    def average_bulk_moduli(self, volumes, bulk_moduli, shear_moduli):
        return (  self.upper.average_bulk_moduli(volumes, bulk_moduli, shear_moduli)\
                + self.lower.average_bulk_moduli(volumes, bulk_moduli, shear_moduli))/2.0

    def average_shear_moduli(self, volumes, bulk_moduli, shear_moduli):
        return (   self.upper.average_shear_moduli(volumes, bulk_moduli, shear_moduli)\
                 + self.lower.average_shear_moduli(volumes, bulk_moduli, shear_moduli))/2.0
        
      

def voigt_average_function(phase_volume,X):
    """
    Do Voigt (iso-strain) average.  Rather like
    resistors in series.  Called by voigt and
    voigt_reuss_hill classes, takes a list of
    volumes and moduli, returns a modulus.
    """
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    X_voigt = sum( V_i[i]/V_tot * X[i] for i in it)
    return X_voigt

def reuss_average_function(phase_volume,X):
    """
    Do Reuss (iso-stress) average.  Rather like
    resistors in parallel.  Called by reuss and
    voigt_reuss_hill classes, takes a list of
    volumes and moduli, returns a modulus.
    """
    it = range(len(phase_volume))
    V_i = phase_volume
    V_tot = sum(V_i)
    if (min(X)<=0.0):
        X_reuss = 0.0
        warnings.warn("Oops, called reuss_average with Xi<=0!")
    else:
        X_reuss = 1./sum(  V_i[i]/V_tot* 1./X[i] for i in it)
    return X_reuss

def voigt_reuss_hill_function(phase_volume,X):
    """
    Do Voigt-Reuss-Hill average (arithmetic mean
    of Voigt and Reuss bounds).  Called by
    voigt_reuss_hill class, takes a list of
    volumes and moduli, returns a modulus.
    """
    X_vrh = (voigt_average_function(phase_volume,X) + reuss_average_function(phase_volume,X))/2.
    return X_vrh

