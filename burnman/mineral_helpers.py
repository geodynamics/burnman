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



class HelperSolidSolution(Mineral):
    """
    This material is deprecated!

    Class for coming up with a new mineral based based on a solid
    solution between two or more end member minerals.  It is not
    completely clear how to do this, or how valid this approximation
    is, but here we just do a weighted arithmetic average of the
    thermoelastic properties of the end members according to their molar fractions
    """
    def __init__(self, endmembers, molar_fractions):
        """
        Takes a list of end member minerals, and a matching list of
        molar fractions of those minerals for mixing them.  Simply
        comes up with a new mineral by doing a weighted arithmetic
        average of the end member minerals
        """
        self.endmembers = endmembers
        self.molar_fractions = molar_fractions
        assert(len(endmembers) == len(molar_fractions))
        assert(sum(molar_fractions) > 0.9999)
        assert(sum(molar_fractions) < 1.0001)

        self.method = endmembers[0].method

        #does not make sense to do a solid solution with different number of
        #atoms per formula unit or different equations of state, at least not simply...
        for m in endmembers:
            m.set_method(self.method)
            if('n' in endmembers[0].params):
                assert(m.params['n'] == endmembers[0].params['n'])
        
        self.params = {}


    def debug_print(self, indent=""):
        print("%sHelperSolidSolution(%s):" % (indent, self.to_string()))
        indent += "  "
        for (fraction, mat) in zip(self.molar_fractions, self.endmembers):
            print("%s%g of" % (indent, fraction))
            mat.debug_print(indent + "  ")

    def set_method(self, method):
        for mat in self.endmembers:
            mat.set_method(method)
        self.method = self.endmembers[0].method

    def set_state(self, pressure, temperature):
        for mat in self.endmembers:
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.endmembers))
        self.params = {}
        for prop in self.endmembers[0].params:
           try:
                self.params[prop] = sum([ self.endmembers[i].params[prop] * self.molar_fractions[i] for i in itrange ])
           except TypeError:
                #if there is a type error, it is probably a string. Just go with the value of the first endmembers.
                self.params[prop] = self.endmembers[0].params[prop]

        Mineral.set_state(self,pressure,temperature)

class HelperSpinTransition(Mineral):
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
        self.transition_pressure = transition_pressure
        self.ls_mat = ls_mat
        self.hs_mat = hs_mat
        self.active_mat = None
        Mineral.__init__(self)

    def debug_print(self, indent=""):
        print("%sHelperSpinTransition:" % indent)
        self.ls_mat.debug_print(indent+"  ")
        self.hs_mat.debug_print(indent+"  ")

    def set_method(self, method):
        self.ls_mat.set_method(method)
        self.hs_mat.set_method(method)

    def set_state(self, pressure, temperature):
        if (pressure >= self.transition_pressure):
            self.active_mat = self.ls_mat
        else:
            self.active_mat = self.hs_mat
        Material.set_state(self, pressure, temperature)
        self.active_mat.set_state(pressure, temperature)

    def evaluate(self,vars_list,pressures, temperatures):
        output = np.empty((len(vars_list),len(pressures)))
        for i in range(len(pressures)):
            if (pressures[i] >= self.transition_pressure):
                self.active_mat = self.ls_mat
            else:
                self.active_mat = self.hs_mat
            self.active_mat.set_state(pressures[i], temperatures[i])
            for j in range(len(vars_list)):
                output[j,i]=getattr(self.active_mat,vars_list[j])
        return output

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        return ([self.active_mat],[1.0])







class HelperFeDependent(Mineral):

    """
    This material is deprecated!

    Helper to implement a rock that does iron exchange (two minerals with
    changing iron content based on P,T).

    Classes deriving from this helper need to implement
    create_inner_material() and provide a function in __init__ that computes
    the iron exchange.
    """

    def __init__(self, iron_number_with_pt, idx):
        self.iron_number_with_pt = iron_number_with_pt
        self.which_index = idx  # take input 0 or 1 from iron_number_with_pt()
        self.method = None

    def debug_print(self, indent=""):
        print("%sHelperFeDependent:" % indent)
        print("%s  %s" % (indent, self.to_string()))

    def create_inner_material(self, iron_number):
        raise NotImplementedError("need to implement create_inner_material() in derived class!")
        return None

    def set_method(self, method):
        self.method = method

    def iron_number(self):
        return self.iron_number_with_pt(self.pressure, self.temperature)[self.which_index]

    def set_state(self, pressure, temperature):
        Material.set_state(self, pressure, temperature)
        self.endmembers = self.create_inner_material(self.iron_number())
        if self.method:
            self.endmembers.set_method(self.method)
        self.endmembers.set_state(pressure, temperature)

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        return ([self.endmembers],[1.0])




