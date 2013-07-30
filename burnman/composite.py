# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
from minerals import *
import warnings
from collections import namedtuple
import averaging_schemes as avg


phase = namedtuple('phase', ['mineral', 'fraction'])


class composite:
    """
    Base class for a composite material.  The constructor takes a tuple of tuples, 
    where the inner tuple is a mineral/molar-fraction pair.  This can then be passed
    to averaging schemes, the adiabatic geotherm function, or anything else that
    expects a composite material
    """
    def __init__(self, phase_tuples):
        total = 0
        self.phases = [] # make the rock composition an immutable tuple
        for ph in phase_tuples:
            total += ph[1]
        if total != 1.0:
            warnings.warn('Warning: list of molar fractions does not add up to one. Normalizing')
        for ph in phase_tuples:
            self.phases.append( phase(ph[0], ph[1]/total) )

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for ph in self.phases:
            ph.mineral.set_method(method) 

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        """ Update the material to the given pressure [Pa] and temperature [K].
        """
        self.pressure = pressure
        self.temperature = temperature
        for ph in self.phases:
            ph.mineral.set_state(pressure, temperature) 

    def density(self):
        """ Compute the density of the composite based on the molar volumes and masses """
        densities = np.array([ph.mineral.density() for ph in self.phases])
        volumes = np.array([ph.mineral.molar_volume() for ph in self.phases])
        return np.sum(densities*volumes)/np.sum(volumes)
        
                



if __name__ == "__main__":
    pyrolite = composite( [ (SLB2005.mg_fe_perovskite(0.2), 0.8), (SLB2005.ferropericlase(0.4), 0.8) ] )
    pyrolite.set_method('slb3')
    pyrolite.set_state(40.e9, 2000)
    print pyrolite.phases[0].mineral.density(), pyrolite.phases[1].mineral.density(), pyrolite.density()

