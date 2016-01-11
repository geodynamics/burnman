# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


"""
This module provides several helper minerals/materials.

"""
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import warnings

from .material import Material
from .mineral import Mineral
from .composite import Composite


class HelperSpinTransition(Composite):
    """
    Helper class that makes a mineral that switches between two materials
    (for low and high spin) based on some transition pressure [Pa]
    """

    def __init__(self, transition_pressure, ls_mat, hs_mat):
        """
        Takes a transition pressure, and two minerals.  Use the
        thermoelastic parameters for ls_mat below the transition
        pressure, and the thermoelastic parameters for hs_mat
        above the transition pressure
        """
        Material.__init__(self)
        self.transition_pressure = transition_pressure
        self.ls_mat = ls_mat
        self.hs_mat = hs_mat
        Composite.__init__(self, [ls_mat, hs_mat])

    def debug_print(self, indent=""):
        print("%sHelperSpinTransition:" % indent)
        self.ls_mat.debug_print(indent+"  ")
        self.hs_mat.debug_print(indent+"  ")

    def set_state(self, pressure, temperature):
        if (pressure >= self.transition_pressure):
            Composite.set_fractions(self, [1.0, 0.0])
        else:
            Composite.set_fractions(self, [0.0, 1.0])

        Composite.set_state(self, pressure, temperature)
