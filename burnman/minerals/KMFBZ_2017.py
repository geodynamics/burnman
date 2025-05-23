# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
from Kurnosov et al.  Nature 2017
Kurnosov, A., et al. "Evidence for a Fe3+-rich pyrolitic lower mantle from
(Al, Fe)-bearing bridgmanite elasticity data."
Nature 543.7646 (2017): 543-546. doi:10.1038/nature21390
"""


from ..classes.mineral import Mineral
from ..classes.solution import Solution
from ..classes.solutionmodel import IdealSolution, SymmetricRegularSolution
from ..utils.chemistry import dictionarize_formula, formula_mass

"""
SOLID SOLUTIONS
"""


class bridgmanite(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "bridgmanite/perovskite"
        self.solution_model = IdealSolution(
            endmembers=[
                [mg_si_perovskite(), "[Mg][Si]O3"],
                [fe_si_perovskite(), "[Fe][Si]O3"],
                [al_al_perovskite(), "[Al][Al]O3"],
                [fe_al_perovskite(), "[Fe][Al]O3"],
            ]
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class ferropericlase(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "magnesiowustite/ferropericlase"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[[periclase(), "[Mg]O"], [wuestite(), "[Fe]O"]],
            energy_interaction=[[13.0e3]],
        )

        Solution.__init__(self, molar_fractions=molar_fractions)


"""
ENDMEMBERS
"""


class mg_si_perovskite(Mineral):
    def __init__(self):
        formula = "MgSiO3"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "MgSiO3 perovskite",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -1368000.0,
            "V_0": 2.4445e-05,
            "K_0": 2.51e11,
            "Kprime_0": 4.1,
            "Debye_0": 905.0,
            "grueneisen_0": 1.57,
            "q_0": 1.1,
            "G_0": 1.73e11,
            "Gprime_0": 1.7,
            "eta_s_0": 2.3,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class fe_si_perovskite(Mineral):
    def __init__(self):
        formula = "FeSiO3"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "FeSiO3 perovskite",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -1043000.0,
            "V_0": 2.534e-05,
            "K_0": 2.72e11,
            "Kprime_0": 4.1,
            "Debye_0": 871.0,
            "grueneisen_0": 1.57,
            "q_0": 1.1,
            "G_0": 1.33 + 11,
            "Gprime_0": 1.4,
            "eta_s_0": 2.3,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class fe_al_perovskite(Mineral):
    def __init__(self):
        formula = "FeAlO3"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "FeAlO3 perovskite",
            "formula": formula,
            "equation_of_state": "slb3",
            "V_0": 2.69e-05,
            "K_0": 2.20e11,
            "Kprime_0": 1.3,
            "Debye_0": 886.0,
            "grueneisen_0": 1.57,
            "q_0": 1.1,
            "G_0": 96.0e09,
            "Gprime_0": 3.4,
            "eta_s_0": 2.5,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class al_al_perovskite(Mineral):
    def __init__(self):
        formula = "Al2O3"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "Al2O3 perovskite",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -1533878.0,
            "V_0": 2.494e-05,
            "K_0": 2.58e11,
            "Kprime_0": 4.1,
            "Debye_0": 886.0,
            "grueneisen_0": 1.57,
            "q_0": 1.1,
            "G_0": 1.71e11,
            "Gprime_0": 1.5,
            "eta_s_0": 2.5,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class periclase(Mineral):
    def __init__(self):
        formula = "MgO"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "Periclase",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -569000.0,
            "V_0": 1.124e-05,
            "K_0": 1.602e11,
            "Kprime_0": 3.99,
            "Debye_0": 767.0,
            "grueneisen_0": 1.36,
            "q_0": 1.7,
            "G_0": 1.31e11,
            "Gprime_0": 2.1,
            "eta_s_0": 2.8,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class wuestite(Mineral):
    def __init__(self):
        formula = "FeO"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "Wuestite",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -242000.0,
            "V_0": 1.226e-05,
            "K_0": 1.49e11,
            "Kprime_0": 3.6,
            "Debye_0": 454.0,
            "grueneisen_0": 1.53,
            "q_0": 1.7,
            "G_0": 60.0e9,
            "Gprime_0": 1.8,
            "eta_s_0": 0.6,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class ca_perovskite(Mineral):
    def __init__(self):
        formula = "CaSiO3"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "CaSiO3 perovskite",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -1463358.0,
            "V_0": 2.754e-05,
            "K_0": 2.36e11,
            "Kprime_0": 3.9,
            "Debye_0": 802.0,
            "grueneisen_0": 1.89,
            "q_0": 0.9,
            "G_0": 1.57e11,
            "Gprime_0": 2.2,
            "eta_s_0": 1.3,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


class hcp_iron(Mineral):
    def __init__(self):
        formula = "Fe"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "Hexagonal close packed iron",
            "formula": formula,
            "equation_of_state": "slb3",
            "F_0": -1463358.0,
            "V_0": 0.673e-05,
            "K_0": 1.64e11,
            "Kprime_0": 4.0,
            "Debye_0": 422.0,
            "grueneisen_0": 1.71,
            "q_0": 1.4,
            "G_0": 0.815e11,
            "Gprime_0": 1.9,
            "eta_s_0": 7,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)


"""
Mineral aliases
"""
perovskite = bridgmanite
