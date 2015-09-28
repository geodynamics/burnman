# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


import numpy as np
import warnings

from burnman import Material
from burnman import Mineral

def check_pairs(phases, fractions):
        if len(fractions) < 1:
            raise Exception('ERROR: we need at least one phase')

        if len(phases) != len(fractions):
            raise Exception('ERROR: different array lengths for phases and fractions')

        total = sum(fractions)
        if abs(total-1.0)>1e-10:
            raise Exception('ERROR: list of molar fractions does not add up to one')
        for p in phases:
            if not isinstance(p, Mineral):
                raise Exception('ERROR: object of type ''%s'' is not of type Mineral' % (type(p)))


# static composite of minerals/composites
class Composite(Material):
    """
    Base class for a static composite material with fixed molar fractions. The
    elements can be minerals or materials, meaning composite can be nested
    arbitrarily.

    This class is available as ``burnman.Composite``.
    """
    def __init__(self, phases, fractions=None):
        """
        Create a composite using a list of phases and their fractions (adding to 1.0).

        Parameters
        ----------
        phases: list of :class:`burnman.Material`
            list of phases.
        fractions: list of floats
            molar fraction for each phase.
        """

        
        assert(len(phases)>0)
        
        if fractions is not None:
            assert(len(phases)==len(fractions))
            for f in fractions:
                assert (f >= -1e-12)
                fractions = [max(0.0, fraction) for fraction in fractions]  # turn -1e-12 into 0.0
                try:
                    total = sum(fractions)
                except TypeError:
                    raise Exception("Since v0.8, burnman.Composite takes an array of Materials, then an array of fractions")
                if abs(total - 1.0) > 1e-12:
                    warnings.warn("Warning: list of molar fractions does not add up to one but %g. Normalizing." % total)
                    fractions = [fr / total for fr in fractions]

        self.children = zip(phases, fractions)

    def debug_print(self, indent=""):
        print "%sComposite:" % indent
        indent += "  "
        for (phase, fraction) in self.children:
            print "%s%g of" % (indent, fraction)
            phase.debug_print(indent + "  ")


    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the composite
        """
        for (phase, fraction) in self.children:
            phase.set_method(method)

    def unroll(self):
        phases = []
        fractions = []
        for (phase, fraction) in self.children:
            p_mineral, p_fraction = phase.unroll()
            check_pairs(p_mineral, p_fraction)
            fractions.extend([i*fraction for i in p_fraction])
            phases.extend(p_mineral)
        return (phases, fractions)

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
        for (phase, fraction) in self.children:
            phase.set_state(pressure, temperature)

    def density(self):
        """
        Compute the density of the composite based on the molar volumes and masses
        """
        densities = np.array([phase.density() for (phase,_) in self.children])
        volumes = np.array([phase.molar_volume()*fraction for (phase, fraction) in self.children])
        return np.sum(densities*volumes)/np.sum(volumes)

