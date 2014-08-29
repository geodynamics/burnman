# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import warnings

import numpy as np

kd = lambda x,y : 1 if x==y else 0

class SolidSolution(Mineral):
    """
    This is the base class for all solid solutions. 
    States of the solid solution can only be queried after setting 
    the pressure and temperature using set_state(). 
    The method for computing properties of
    the solution is set using set_method(), which should be done
    once after creating the material.

    This class is available as ``burnman.SolidSolution``.

    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T 
    and P derivatives in J/K/mol and m^3/(mol molecule).
    """

    # init sets up matrices to speed up calculations when P, T, X is defined.
    def __init__(self, base_material, interaction_parameter):
        self.base_material = base_material
        self.molar_fraction = molar_fraction
        assert(len(base_material) == len(molar_fraction))
        assert(len(base_material) == len(van_laar_parameter))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        # make sure the number of atoms in each endmember is the same
        # (when we start adding minerals with vacancies
        # this will need to be changed).
        for m in base_material:
            if(base_material[0].params.has_key('n')):
                assert(m.params['n'] == base_material[0].params['n'])

    # Set the state of all the endmembers
    def set_state(self, pressure, temperature):
        for mat in self.base_materials:
            mat.method = self.method
            mat.set_state(pressure, temperature)

        itrange = range(0, len(self.base_material))

        self.params = {}

        # NEEDS CHANGING FROM HERE!!!
        # Description of RTlngamma terms given by the asymmetric formalism
        # is included in Holland and Powell (2003)

        # the volume is calculated in the same way as Gex (p.493),
        # replacing Wij with Wvij
        Vex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*self.interaction_parameter[i][2] / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        V= sum([ self.base_materials[i].molar_volume() * self.molar_fraction[i] for i in itrange ]) + Vex

        # the volume is calculated in the same way as Gex (p.493),
        # replacing Wij with Wvij
        # NEED TO CORRECT INDICES FOR INT PARAM
        # W[i][j] == W[(itrange-1)*i+(j-1)]
        Vex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*self.interaction_parameter[(itrange-1)*i+(j-1)][2] / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        V= sum([ self.base_materials[i].molar_volume() * self.molar_fraction[i] for i in itrange ]) + Vex


        # Same for the Gibbs free energy
        # NEED TO CORRECT INDICES FOR INT PARAM
        Gex=-2.*sum([ sum([ self.molar_fraction[i]*self.van_laar_parameter[i]*self.molar_fraction[j]*self.van_laar_parameter[j]*(pressure*self.interaction_parameter[(itrange-1)*i+(j-1)][0] + temperature*self.interaction_parameter[(itrange-1)*i+(j-1)][1] + pressure*self.interaction_parameter[(itrange-1)*i+(j-1)][2]) / (self.van_laar_parameter[i] + self.van_laar_parameter[j]) for j in range(i, len(self.base_material)) ]) for i in range(0, len(self.base_material)-1) ])/sum([ self.van_laar_parameter[l]*self.molar_fraction[l] for l in itrange ])
        G= sum([ self.base_materials[i].molar_gibbs() * self.molar_fraction[i] for i in itrange ]) + Gex

        

        

        for prop in self.base_materials[0].params:
           try:
               self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
               #if there is a type error, it is probably a string.  Just go with the value of the first base_material.
               self.params[prop] = self.base_materials[0].params[prop]
        Mineral.set_state(self, pressure, temperature)

