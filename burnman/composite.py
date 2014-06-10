# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings
from collections import namedtuple

from burnman.material import Material
from burnman.mineral import Mineral

def check_pairs(fractions, minerals):
        if len(fractions)<1:
            raise Exception('ERROR: we need at least one mineral')

        if len(fractions) != len(minerals):
            raise Exception('ERROR: different array lengths')

        total = sum(fractions)
        if abs(total-1.0)>1e-10:
            raise Exception('ERROR: list of molar fractions does not add up to one')
        for p in minerals:
            if not isinstance(p,Mineral):
                raise Exception('ERROR: object of type ''%s'' is not of type material' % (type(p)))


# static composite of minerals/composites
class Composite(Material):
    """
    Base class for a static composite material with fixed molar fractions. The
    elements can be minerals or materials, meaning composite can be nested
    arbitrarily.

    This class is available as ``burnman.Composite``.
    """
    def __init__(self, fractions, phases=None):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        Parameters
        ----------
        fractions: list of floats
            molar fraction for each phase.
        phases: list of :class:`burnman.Material`
            list of phases.
        """

        if phases is None:
            # compatibility hack:
            tmp = fractions
            fractions = [pt[1] for pt in tmp]
            phases = [pt[0] for pt in tmp]

        assert(len(phases)==len(fractions))
        assert(len(phases)>0)
        for f in fractions:
            assert(f>=0)

        self.children = zip(fractions, phases)

        total = sum(fractions)

        if (total-1.0)>1e-8:
            warnings.warn('Warning: list of molar fractions does not add up to one. Normalizing.')
            for (fraction,phase) in self.children:
                fraction = fraction / total


    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for (fraction, phase) in self.children:
            phase.set_method(method)

    def unroll(self):
        fractions = []
        minerals = []

        for (fraction, phase) in self.children:
            p_fr,p_min = phase.unroll()
            check_pairs(p_fr, p_min)
            fractions.extend([i*fraction for i in p_fr])
            minerals.extend(p_min)
        return (fractions, minerals)

    def to_string(self):
        """
        return the name of the composite
        """
        return "'" + self.__class__.__name__ + "'"

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].
        """
        self.pressure = pressure
        self.temperature = temperature
        for (fraction, phase) in self.children:
            phase.set_state(pressure, temperature)

    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        densities = np.array([ph.density() for (_,ph) in self.children])
        volumes = np.array([ph.molar_volume()*fraction for (fraction, ph) in self.children])
        return np.sum(densities*volumes)/np.sum(volumes)

