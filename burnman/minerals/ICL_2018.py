# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Irving et al. 2018
^^^^^^^^^^^^^

Seismically derived EoS parameters for the outer core using normal mode center frequencies.
These parameters are derived parameterizing the outer core velocity and density with an equation-of-state, and inverting for K, Kprime, and molar_volume (molar_mass is fixed to 0.05, as sensitivy of the normal modes is only to molar density). The equations-of-state are assumed to capture the isentropic behaviour across the outer core. Potentially, these parameters can be used to represent cores of other planets.
See:
Irving, Cottaar, Lekic (2018) Seismically determined elastic parameters for the outer core, Science Advances

"""
from __future__ import absolute_import

from ..classes import mineral_helpers as helpers
from ..classes.mineral import Mineral


class EPOC_vinet (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'vinet',
            'V_0': 8.18e-6,
            'K_0': 67.5e9,
            'Kprime_0': 6.12,
            'molar_mass': .05,}
        Mineral.__init__(self)

class EPOC_bm (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'bm3',
            'V_0': 7.63e-6,
            'K_0': 120e9,
            'Kprime_0': 4.6,
            'molar_mass': .05,}
        Mineral.__init__(self)
