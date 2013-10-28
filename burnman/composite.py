# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings
from collections import namedtuple

class abstract_material:
    def set_method(self, method):
        raise NotImplementedError("need to implement this in derived class!")

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        self.pressure = pressure
        self.temperature = temperature

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        raise NotImplementedError("need to implement this in derived class!")
        return ()

    def density(self):
        raise NotImplementedError("need to implement this in derived class!")
        return inf            
        

class composite_base(abstract_material):
    """
    base class for writing your own composites that need to dynamically pick
    the fractions and or minerals. The only function that needs to be implemented
    is unroll, which returns (fractions,minerals), both arrays. unroll may depend
    on temperature and pressure, and the fractions are molar fractions, rather 
    than volume. 
    """
    def set_state(self, pressure, temperature):
        abstract_material.set_state(self, pressure, temperature)




phase = namedtuple('phase', ['mineral', 'fraction'])

    
def check_pairs(fractions, minerals):
        if len(fractions)<1:
            raise Exception('ERROR: we need at least one mineral')

        if len(fractions) != len(minerals):
            raise Exception('ERROR: different array lengths');

        total = sum(fractions)
        if abs(total-1.0)>1e-10:
            raise Exception('ERROR: list of molar fractions does not add up to one')
        for p in minerals:
            if not isinstance(p,abstract_material):
                print type(p)
                raise Exception('ERROR: object is not of type abstract_material')


# static composite of minerals/composites
class composite(composite_base):
    """
    Base class for a composite material.  The constructor takes a tuple of tuples, 
    where the inner tuple is a mineral/molar-fraction pair.  This can then be passed
    to averaging schemes, the adiabatic geotherm function, or anything else that
    expects a composite material
    """
    def __init__(self, phase_tuples):
        total = 0
        self.staticphases = [] # make the rock composition an immutable tuple
        for ph in phase_tuples:
            total += ph[1]
        if total != 1.0:
            warnings.warn('Warning: list of molar fractions does not add up to one. Normalizing')
        for ph in phase_tuples:
            self.staticphases.append( phase(ph[0], ph[1]/total) )

    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for ph in self.staticphases:
            ph.mineral.set_method(method) 

    def unroll(self):
        fractions = []
        minerals = []

        for p in self.staticphases:
            p_fr,p_min = p.mineral.unroll()
            check_pairs(p_fr, p_min)
            fractions.extend([i*p.fraction for i in p_fr])
            minerals.extend(p_min)
        return (fractions, minerals)

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    def debug_print(self):
        (fr,mins) = self.unroll()
        for (fr,mi) in zip(fr,mins):
            print "%g of phase %s"%(fr,mi.to_string())

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        self.pressure = pressure
        self.temperature = temperature
        for ph in self.staticphases:
            ph.mineral.set_state(pressure, temperature) 

    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        densities = np.array([ph.mineral.density() for ph in self.staticphases])
        volumes = np.array([ph.mineral.molar_volume()*ph.fraction for ph in self.staticphases])
        return np.sum(densities*volumes)/np.sum(volumes)

        
                



if __name__ == "__main__":
    import minerals

    pyrolite = composite( [ (minerals.SLB_2005.mg_fe_perovskite(0.2), 0.8), (minerals.SLB_2005.ferropericlase(0.4), 0.8) ] )
    pyrolite.set_method('slb3')
    pyrolite.set_state(40.e9, 2000)
    print pyrolite.staticphases[0].mineral.density(), pyrolite.staticphases[1].mineral.density(), pyrolite.density()

