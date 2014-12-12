# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This module provides several helper minerals/materials.

"""


import numpy as np
import warnings

from burnman import Material, Mineral

#import burnman.equation_of_state as eos
#import burnman.birch_murnaghan as bm
#import burnman.slb as slb
#import burnman.mie_grueneisen_debye as mgd


class HelperSolidSolution(Mineral):
    """
    Class for coming up with a new mineral based based on a solid
    solution between two or more end member minerals.  It is not
    completely clear how to do this, or how valid this approximation
    is, but here we just do a weighted arithmetic average of the
    thermoelastic properties of the end members according to their molar fractions
    """
    def __init__(self, base_materials, molar_fraction):
        """
        Takes a list of end member minerals, and a matching list of
        molar fractions of those minerals for mixing them.  Simply
        comes up with a new mineral by doing a weighted arithmetic
        average of the end member minerals
        """
        self.base_materials = base_materials
        self.molar_fraction = molar_fraction
        assert(len(base_materials) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        self.method = base_materials[0].method

        #does not make sense to do a solid solution with different number of
        #atoms per formula unit or different equations of state, at least not simply...
        for m in base_materials:
            m.set_method(self.method)
            if(base_materials[0].params.has_key('n')):
                assert(m.params['n'] == base_materials[0].params['n'])
        
        self.params = {}

    def set_method(self, equation_of_state):
        for m in self.base_materials:
            m.set_method(equation_of_state)

    def debug_print(self, indent=""):
        print "%sHelperSolidSolution(%s):" % (indent, self.to_string())
        indent += "  "
        for (fraction, mat) in zip(self.molar_fraction, self.base_materials):
            print "%s%g of" % (indent, fraction)
            mat.debug_print(indent + "  ")

    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_materials))
        self.params = {}
        for prop in self.base_materials[0].params:
           try:
                self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
                #if there is a type error, it is probably a string. Just go with the value of the first base_material.
                self.params[prop] = self.base_materials[0].params[prop]

        Mineral.set_state(self,pressure,temperature)

class HelperSpinTransition(Material):
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
        self.method=ls_mat.method

    def debug_print(self, indent=""):
        print "%sHelperSpinTransition:" % indent
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

        self.active_mat.set_state(pressure, temperature)

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        return ([1.0],[self.active_mat])

    def density(self):
        return self.active_mat.density()






class HelperFeDependent(Material):

    """
    Helper to implement a rock that does iron exchange (two minerals with
    changing iron content based on P,T).

    Classes deriving from this helper need to implement
    create_inner_material() and provide a function in __init__ that computes
    the iron exchange.
    """

    def __init__(self, iron_number_with_pt, idx):
        self.iron_number_with_pt = iron_number_with_pt
        self.which_index = idx  # take input 0 or 1 from iron_number_with_pt()

    def debug_print(self, indent=""):
        print "%sHelperFeDependent:" % indent
        print "%s  %s" % (indent, self.to_string())

    def create_inner_material(self, iron_number):
        raise NotImplementedError("need to implement create_inner_material() in derived class!")
        return None

    def set_method(self, method):
        self.method = method

    def iron_number(self):
        return self.iron_number_with_pt(self.pressure, self.temperature)[self.which_index]

    def set_state(self, pressure, temperature):
        Material.set_state(self, pressure, temperature)
        self.base_material = self.create_inner_material(self.iron_number())
        self.base_material.set_method(self.method)
        self.base_material.set_state(pressure, temperature)

    def unroll(self):
        """ return (fractions, minerals) where both are arrays. May depend on current state """
        return ([1.0],[self.base_material])

    def density(self):
        return self.base_material.density()


