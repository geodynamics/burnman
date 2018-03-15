# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
JH_2015 
Solid solutions from Jennings and Holland, 2015 and references therein
(10.1093/petrology/egv020).
The values in this document are all in S.I. units,
unlike those in the original tc file.
"""

from ..mineral import Mineral
from ..solidsolution import SolidSolution
from ..combinedmineral import CombinedMineral
from ..solutionmodel import *
from ..processchemistry import dictionarize_formula, formula_mass
from .HP_2011_ds62 import *

"""
SOLID SOLUTIONS

The parameters in Jennings and Holland (2015) are given in the following units:
[kJ/mol], [kJ/K/mol], [kJ/kbar/mol]

N.B. The excess entropy terms in these solution models have the opposite sign 
to the thermal parameters in Jennings and Holland, 2015. 
This is consistent with its treatment as an excess entropy term (W=W_H-T*W_S+P*W_V), 
rather than a thermal correction to the interaction parameter   (W=W_0+T*W_T+P*W_P).
"""
class ferropericlase(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'ferropericlase (FM)'
        self.endmembers = [[per(), '[Mg]O'],
                           [fper(), '[Fe]O']]
        self.solution_type = 'symmetric'
        self.energy_interaction = [[18.e3]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)
        
class plagioclase(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'plagioclase (NCAS)'
        self.endmembers = [[an(), '[Ca][Al]2Si2O8'],
                           [abh(), '[Na][Al1/2Si1/2]2Si2O8']] # Al-avoidance model
        self.solution_type = 'asymmetric'
        self.alphas = [0.39, 1.]
        self.energy_interaction = [[22.4e3]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)

class clinopyroxene(SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = 'clinopyroxene (NCFMASCrO)'
        self.endmembers = [[di(),   '[Mg][Ca][Si]2O6'],
                           [cfs(),  '[Fe][Fe][Si]2O6'],
                           [cats(), '[Al][Ca][Si1/2Al1/2]2O6'],
                           [crdi(), '[Cr][Ca][Si1/2Al1/2]2O6'],
                           [cess(), '[Fef][Ca][Si1/2Al1/2]2O6'],
                           [jd(),   '[Al][Na][Si]2O6'],
                           [cen(),  '[Mg][Mg][Si]2O6'],
                           [cfm(),  '[Mg][Fe][Si]2O6']] # note cfm ordered endmember 
        self.solution_type = 'asymmetric'
        self.alphas = [1.2, 1.0, 1.9, 1.9, 1.9, 1.2, 1.0, 1.0]
        self.energy_interaction = [[20.e3, 12.3e3, 8.e3, 8.e3, 26.e3, 29.8e3, 18.e3],
                                   [25.e3, 34.e3, 34.e3, 24.e3, 7.e3, 4.e3],
                                   [2.e3, 2.e3, 6.e3, 45.7e3, 27.e3],
                                   [2.e3, 3.e3, 48.e3, 36.e3],
                                   [3.e3, 58.e3, 36.e3],
                                   [40.e3, 40.e3],
                                   [4.e3]]
        self.volume_interaction = [[0., -0.1e-5, 0., 0., 0., -0.03e-5, 0.],
                                   [-0.1e-5, 0., 0., 0., 0., 0.],
                                   [0., 0., 0., -0.29e-5, -0.1e-5],
                                   [0., 0., 0., 0.],
                                   [0., 0., 0.],
                                   [0., 0.],
                                   [0.]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)

class cfs(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'clinoferrosilite',
                                 mineral_list = [fs()],
                                 molar_amounts = [1.],
                                 free_energy_adjustment=[3.8e3, 3., 0.03e-5])  
class crdi(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'chromium diopside',
                                 mineral_list = [cats(), kos(), jd()],
                                 molar_amounts = [1., 1., -1.],
                                 free_energy_adjustment=[-3.e3, 0., 0.])  
class cess(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'ferric diopside',
                                 mineral_list = [cats(), acm(), jd()],
                                 molar_amounts = [1., 1., -1.],
                                 free_energy_adjustment=[-6.e3, 0., 0.])  
class cen(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'clinoenstatite',
                                 mineral_list = [en()],
                                 molar_amounts = [1.],
                                 free_energy_adjustment=[3.5e3, 2., 0.048e-5])

class cfm(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'ordered clinoferroenstatite',
                                 mineral_list = [en(), fs()],
                                 molar_amounts = [0.5, 0.5],
                                 free_energy_adjustment=[-3.e3, 0., 0.])

class olivine(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'olivine (FMS)'
        self.endmembers = [[fo(), '[Mg]2SiO4'],
                           [fa(), '[Fe]2SiO4']]
        self.solution_type = 'symmetric'
        self.energy_interaction = [[9.e3]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)
      

class spinel(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'disordered spinel (CFMASO)'
        self.endmembers = [[sp(),   '[Al2/3Mg1/3]3O4'],
                           [herc(), '[Al2/3Fe1/3]3O4'],
                           [mt(),   '[Fef2/3Fe1/3]3O4'],
                           [picr(), '[Cr2/3Mg1/3]3O4']]
        self.solution_type = 'symmetric'
        self.energy_interaction = [[4.e3, 56.e3, 39.e3],
                                   [32.e3, 27.e3],
                                   [36.e3]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)

        
class garnet(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'garnet (CFMASCrO, low pressure)'
        self.endmembers = [[py(),   '[Mg]3[Al]2Si3O12'],
                           [alm(),  '[Fe]3[Al]2Si3O12'],
                           [gr(),   '[Ca]3[Al]2Si3O12'],
                           [andr(), '[Ca]3[Fef]2Si3O12'],
                           [knor(), '[Mg]3[Cr]2Si3O12']]
        self.solution_type = 'symmetric'
        self.energy_interaction = [[4.e3, 35.e3, 91.e3, 2.e3],
                                   [4.e3, 60.e3, 6.e3],
                                   [2.e3, 47.e3],
                                   [101.e3]]
        self.entropy_interaction = [[0., 0., -1.7, 0.],
                                    [0., -1.7, 0.],
                                    [0., 33.8],
                                    [32.1]] # note huge entropy additions! (and sign change from a + bT + cP format)
        self.volume_interaction = [[0.1e-5, 0.1e-5, 0.032e-5, 0.],
                                   [0.1e-5, 0.032e-5, 0.01e-5],
                                   [0., 0.221e-5],
                                   [0.153e-5]]
        SolidSolution.__init__(self, molar_fractions=molar_fractions)

class orthopyroxene(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'orthopyroxene (CFMASCrO)'
        self.endmembers = [[en(),   '[Mg][Mg][Si]0.5Si1.5O6'],
                           [fs(),   '[Fe][Fe][Si]0.5Si1.5O6'],
                           [fm(),   '[Fe][Mg][Si]0.5Si1.5O6'], 
                           [odi(),  '[Mg][Ca][Si]0.5Si1.5O6'],
                           [mgts(), '[Al][Mg][Si1/2Al1/2]0.5Si1.5O6'],
                           [cren(), '[Cr][Mg][Si1/2Al1/2]0.5Si1.5O6'],
                           [mess(), '[Fef][Mg][Si1/2Al1/2]0.5Si1.5O6']] # fm ordered phase, fake T-site multiplicity
        self.solution_type = 'asymmetric'
        self.alphas = [1., 1., 1., 1.2, 1., 1., 1.]
        self.energy_interaction = [[5.2e3, 4.e3, 32.2e3, 13.e3, 8.e3, 8.e3],
                                   [4.e3, 24.e3, 7.e3, 10.e3, 10.e3],
                                   [18.e3, 2.e3, 12.e3, 12.e3],
                                   [75.4e3, 30.e3, 30.e3],
                                   [2.e3, 2.e3],
                                   [2.e3]]
        self.volume_interaction = [[0., 0., 0.12e-5, -0.15e-5, 0., 0.],
                                   [0., 0., -0.15e-5, 0., 0.],
                                   [0., -0.15e-5, 0., 0.],
                                   [-0.94e-5, 0., 0.],
                                   [0., 0.],
                                   [0.]]
        
        SolidSolution.__init__(self, molar_fractions=molar_fractions)
   


class fm(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'ordered ferroenstatite',
                                 mineral_list = [en(), fs()],
                                 molar_amounts = [0.5, 0.5],
                                 free_energy_adjustment=[-6.e3, 0., 0.])
class odi(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'orthodiopside',
                                 mineral_list = [di()],
                                 molar_amounts = [1.],
                                 free_energy_adjustment=[-0.1e3, -0.211, 0.005e-5]) # note sign of *entropy* change.
class cren(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'chromium enstatite',
                                 mineral_list = [mgts(), kos(), jd()],
                                 molar_amounts = [1., 1., -1.],
                                 free_energy_adjustment=[3.e3, 0., 0.])

class mess(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(self,
                                 name = 'ferrienstatite',
                                 mineral_list = [mgts(), acm(), jd()],
                                 molar_amounts = [1., 1., -1.],
                                 free_energy_adjustment=[-15.e3, 0., 0.15e-5])

        
