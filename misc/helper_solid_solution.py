from __future__ import absolute_import
import numpy as np

import burnman
from burnman.minerals import *


class HelperSolidSolution(burnman.Mineral):

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
        burnman.Mineral.__init__(self)
        self.endmembers = endmembers
        self.molar_fractions = molar_fractions
        assert(len(endmembers) == len(molar_fractions))
        assert(sum(molar_fractions) > 0.9999)
        assert(sum(molar_fractions) < 1.0001)

        self.method = endmembers[0].method

        # does not make sense to do a solid solution with different number of
        # atoms per formula unit or different equations of state, at least not
        # simply...
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
                self.params[prop] = sum(
                    [self.endmembers[i].params[prop] * self.molar_fractions[i] for i in itrange])
            except TypeError:
                # if there is a type error, it is probably a string. Just go
                # with the value of the first endmembers.
                self.params[prop] = self.endmembers[0].params[prop]

        burnman.Mineral.set_state(self, pressure, temperature)


class SLB_2005_ferropericlase(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [SLB_2005.periclase(), SLB_2005.wuestite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)


class SLB_2005_mg_fe_perovskite(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [SLB_2005.mg_perovskite(), SLB_2005.fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)


class Murakami_2013_ferropericlase(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [Murakami_2013.periclase(), Murakami_2013.wuestite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)


class Murakami_2013_mg_fe_perovskite(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [
            Murakami_2013.mg_perovskite(), Murakami_2013.fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)

# This is ferropericlase with the deprecated solid solution setup,
# although it is still used in some of the /misc/paper* scripts


class other_ferropericlase(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [other.periclase(), other.wuestite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)

# this is mg_fe_perovskite with the depricated solid solution setup.
# Better not use...


class other_mg_fe_perovskite(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [other.mg_perovskite(), other.fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)

# similar to ferropericlase, using the old solid solution setup. These
# values are based on Zhang, Stixrude and Brodholt 2013


class other_ZSB_2013_mg_fe_perovskite(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [
            other.ZSB_2013_mg_perovskite(), other.ZSB_2013_fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)


class SLB_2011_ZSB_2013_ferropericlase(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [
            SLB_2011_ZSB_2013.periclase(), SLB_2011_ZSB_2013.wuestite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)


class SLB_2011_ZSB_2013_mg_fe_perovskite(HelperSolidSolution):

    def __init__(self, fe_num):
        endmembers = [
            SLB_2011_ZSB_2013.mg_perovskite(), SLB_2011_ZSB_2013.fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num]
            # keep the 0.0 +, otherwise it is an array sometimes
        HelperSolidSolution.__init__(self, endmembers, molar_fractions)
