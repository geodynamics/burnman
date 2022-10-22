# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
JH_2015
^^^^^^^

Solid solutions from Jennings and Holland, 2015 and references therein
(10.1093/petrology/egv020).
The values in this document are all in S.I. units,
unlike those in the original tc file.
"""
import inspect
import numpy as np
from ..classes.solution import Solution
from ..classes.solutionmodel import SymmetricRegularSolution, AsymmetricRegularSolution
from ..classes.combinedmineral import CombinedMineral
from copy import copy

"""
ENDMEMBERS

Direct import from HP_2011_ds62
"""

from . import HP_2011_ds62

# The next few lines import the classes from HP_2011_ds62,
# so that they can be read as attributes of JH_2015.
# This avoids copying the covariance matrix into the JH_2015,
# as would happen with import HP_2011_ds62
import sys

this_module = sys.modules[__name__]
for m in [
    m
    for m in inspect.getmembers(HP_2011_ds62, inspect.isclass)
    if m[1].__module__ == "burnman.minerals.HP_2011_ds62"
]:
    setattr(this_module, m[0], m[1])


"""
SOLID SOLUTIONS

The parameters in Jennings and Holland (2015) are given in the following units:
[kJ/mol], [kJ/K/mol], [kJ/kbar/mol]

N.B. The excess entropy terms in these solution models have the opposite sign
to the thermal parameters in Jennings and Holland, 2015.
This is consistent with its treatment as an excess entropy term
(W=W_H-T*W_S+P*W_V), rather than a thermal correction to the
interaction parameter (W=W_0+T*W_T+P*W_P).
"""


class ferropericlase(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "ferropericlase (FM)"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.per(), "[Mg]O"],
                [HP_2011_ds62.fper(), "[Fe]O"],
            ],
            energy_interaction=[[18.0e3]],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class plagioclase(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "plagioclase (NCAS)"
        self.solution_model = AsymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.an(), "[Ca][Al]2Si2O8"],
                [HP_2011_ds62.abh(), "[Na][Al1/2Si1/2]2Si2O8"],
            ],  # Al-avoidance model
            alphas=[0.39, 1.0],
            energy_interaction=[[22.4e3]],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class clinopyroxene(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "clinopyroxene (NCFMASCrO)"

        self.solution_model = AsymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.di(), "[Mg][Ca][Si]1/2O6"],
                [cfs(), "[Fe][Fe][Si]1/2O6"],
                [HP_2011_ds62.cats(), "[Al][Ca][Si1/2Al1/2]1/2O6"],
                [crdi(), "[Cr][Ca][Si1/2Al1/2]1/2O6"],
                [cess(), "[Fef][Ca][Si1/2Al1/2]1/2O6"],
                [HP_2011_ds62.jd(), "[Al][Na][Si]1/2O6"],
                [cen(), "[Mg][Mg][Si]1/2O6"],
                [cfm(), "[Mg][Fe][Si]1/2O6"],
            ],  # note cfm ordered endmember
            alphas=[1.2, 1.0, 1.9, 1.9, 1.9, 1.2, 1.0, 1.0],
            energy_interaction=[
                [20.0e3, 12.3e3, 8.0e3, 8.0e3, 26.0e3, 29.8e3, 18.0e3],
                [25.0e3, 34.0e3, 34.0e3, 24.0e3, 7.0e3, 4.0e3],
                [2.0e3, 2.0e3, 6.0e3, 45.7e3, 27.0e3],
                [2.0e3, 3.0e3, 48.0e3, 36.0e3],
                [3.0e3, 58.0e3, 36.0e3],
                [40.0e3, 40.0e3],
                [4.0e3],
            ],
            volume_interaction=[
                [0.0, -0.1e-5, 0.0, 0.0, 0.0, -0.03e-5, 0.0],
                [-0.1e-5, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, -0.29e-5, -0.1e-5],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0],
                [0.0],
            ],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class cfs(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="clinoferrosilite",
            mineral_list=[HP_2011_ds62.fs()],
            molar_amounts=[1.0],
            free_energy_adjustment=[3.8e3, 3.0, 0.03e-5],
        )


class crdi(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="chromium diopside",
            mineral_list=[HP_2011_ds62.cats(), HP_2011_ds62.kos(), HP_2011_ds62.jd()],
            molar_amounts=[1.0, 1.0, -1.0],
            free_energy_adjustment=[-3.0e3, 0.0, 0.0],
        )


class cess(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="ferric diopside",
            mineral_list=[HP_2011_ds62.cats(), HP_2011_ds62.acm(), HP_2011_ds62.jd()],
            molar_amounts=[1.0, 1.0, -1.0],
            free_energy_adjustment=[-6.0e3, 0.0, 0.0],
        )


class cen(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="clinoenstatite",
            mineral_list=[HP_2011_ds62.en()],
            molar_amounts=[1.0],
            free_energy_adjustment=[3.5e3, 2.0, 0.048e-5],
        )


class cfm(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="ordered clinoferroenstatite",
            mineral_list=[HP_2011_ds62.en(), HP_2011_ds62.fs()],
            molar_amounts=[0.5, 0.5],
            free_energy_adjustment=[-3.0e3, 0.0, 0.0],
        )


class olivine(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "olivine (FMS)"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.fo(), "[Mg]2SiO4"],
                [HP_2011_ds62.fa(), "[Fe]2SiO4"],
            ],
            energy_interaction=[[9.0e3]],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class spinel(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "disordered spinel (CFMASO)"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.sp(), "[Al2/3Mg1/3]3O4"],
                [HP_2011_ds62.herc(), "[Al2/3Fe1/3]3O4"],
                [HP_2011_ds62.mt(), "[Fef2/3Fe1/3]3O4"],
                [HP_2011_ds62.picr(), "[Cr2/3Mg1/3]3O4"],
            ],
            energy_interaction=[[4.0e3, 56.0e3, 39.0e3], [32.0e3, 27.0e3], [36.0e3]],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class garnet(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "garnet (CFMASCrO, low pressure)"

        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
                [HP_2011_ds62.gr(), "[Ca]3[Al]2Si3O12"],
                [HP_2011_ds62.andr(), "[Ca]3[Fef]2Si3O12"],
                [HP_2011_ds62.knor(), "[Mg]3[Cr]2Si3O12"],
            ],
            energy_interaction=[
                [4.0e3, 35.0e3, 91.0e3, 2.0e3],
                [4.0e3, 60.0e3, 6.0e3],
                [2.0e3, 47.0e3],
                [101.0e3],
            ],
            entropy_interaction=[
                [0.0, 0.0, -1.7, 0.0],
                [0.0, -1.7, 0.0],
                [0.0, 33.8],
                [32.1],
            ],  # note huge entropy additions! (and sign change from a + bT + cP format)
            volume_interaction=[
                [0.1e-5, 0.1e-5, 0.032e-5, 0.0],
                [0.1e-5, 0.032e-5, 0.01e-5],
                [0.0, 0.221e-5],
                [0.153e-5],
            ],
        )
        Solution.__init__(self, molar_fractions=molar_fractions)


class orthopyroxene(Solution):
    def __init__(self, molar_fractions=None):
        self.name = "orthopyroxene (CFMASCrO)"
        self.solution_model = AsymmetricRegularSolution(
            endmembers=[
                [HP_2011_ds62.en(), "[Mg][Mg][Si]0.5Si1.5O6"],
                [HP_2011_ds62.fs(), "[Fe][Fe][Si]0.5Si1.5O6"],
                [fm(), "[Fe][Mg][Si]0.5Si1.5O6"],
                [odi(), "[Mg][Ca][Si]0.5Si1.5O6"],
                [HP_2011_ds62.mgts(), "[Al][Mg][Si1/2Al1/2]0.5Si1.5O6"],
                [cren(), "[Cr][Mg][Si1/2Al1/2]0.5Si1.5O6"],
                [mess(), "[Fef][Mg][Si1/2Al1/2]0.5Si1.5O6"],
            ],  # fm ordered phase, fake T-site multiplicity
            alphas=[1.0, 1.0, 1.0, 1.2, 1.0, 1.0, 1.0],
            energy_interaction=[
                [5.2e3, 4.0e3, 32.2e3, 13.0e3, 8.0e3, 8.0e3],
                [4.0e3, 24.0e3, 7.0e3, 10.0e3, 10.0e3],
                [18.0e3, 2.0e3, 12.0e3, 12.0e3],
                [75.4e3, 30.0e3, 30.0e3],
                [2.0e3, 2.0e3],
                [2.0e3],
            ],
            volume_interaction=[
                [0.0, 0.0, 0.12e-5, -0.15e-5, 0.0, 0.0],
                [0.0, 0.0, -0.15e-5, 0.0, 0.0],
                [0.0, -0.15e-5, 0.0, 0.0],
                [-0.94e-5, 0.0, 0.0],
                [0.0, 0.0],
                [0.0],
            ],
        )

        Solution.__init__(self, molar_fractions=molar_fractions)


class fm(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="ordered ferroenstatite",
            mineral_list=[HP_2011_ds62.en(), HP_2011_ds62.fs()],
            molar_amounts=[0.5, 0.5],
            free_energy_adjustment=[-6.0e3, 0.0, 0.0],
        )


class odi(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="orthodiopside",
            mineral_list=[HP_2011_ds62.di()],
            molar_amounts=[1.0],
            free_energy_adjustment=[-0.1e3, -0.211, 0.005e-5],
        )  # note sign of *entropy* change.


class cren(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="chromium enstatite",
            mineral_list=[HP_2011_ds62.mgts(), HP_2011_ds62.kos(), HP_2011_ds62.jd()],
            molar_amounts=[1.0, 1.0, -1.0],
            free_energy_adjustment=[3.0e3, 0.0, 0.0],
        )


class mess(CombinedMineral):
    def __init__(self):
        CombinedMineral.__init__(
            self,
            name="ferrienstatite",
            mineral_list=[HP_2011_ds62.mgts(), HP_2011_ds62.acm(), HP_2011_ds62.jd()],
            molar_amounts=[1.0, 1.0, -1.0],
            free_energy_adjustment=[-15.0e3, 0.0, 0.15e-5],
        )


def construct_combined_covariance(
    original_covariance_dictionary, combined_mineral_list
):
    """
    This function takes a dictionary containing a list of endmember_names
    and a covariance_matrix, and a list of CombinedMineral instances,
    and creates an updated covariance dictionary containing those
    CombinedMinerals

    Parameters
    ----------
    original_covariance_dictionary : dictionary
        Contains a list of strings of endmember_names of length n
        and a 2D numpy array covariance_matrix of shape n x n

    combined_mineral_list : list of instances of :class:`burnman.CombinedMineral`
        List of minerals to be added to the covariance matrix

    Returns
    -------
    cov : dictionary
        Updated covariance dictionary, with the same keys as the original

    """
    cov_orig = original_covariance_dictionary

    # Update names
    cov = {"endmember_names": copy(cov_orig["endmember_names"])}
    for c in combined_mineral_list:
        cov["endmember_names"].append(c.name)

    # Update covariance matrix
    A = np.identity(len(cov_orig["endmember_names"]))
    for i, indices in enumerate(
        [
            [
                cov_orig["endmember_names"].index(name)
                for name in [mbr[0].params["name"] for mbr in c.mixture.endmembers]
            ]
            for c in combined_mineral_list
        ]
    ):
        B = np.zeros(len(cov_orig["endmember_names"]))
        B[indices] = combined_mineral_list[i].mixture.molar_fractions
        A = np.vstack((A, B))

    cov["covariance_matrix"] = A.dot(cov_orig["covariance_matrix"]).dot(A.T)

    return cov


def cov():
    """
    A function which returns the variance-covariance matrix of the
    zero-point energies of all the endmembers in the dataset.
    Derived from HP_2011_ds62, modified to include all
    the new CombinedMinerals.

    Returns
    -------
    cov : dictionary
        Dictionary keys are:
        - endmember_names: a list of endmember names, and
        - covariance_matrix: a 2D variance-covariance array for the endmember enthalpies of formation
    """
    return construct_combined_covariance(
        original_covariance_dictionary=HP_2011_ds62.cov(),
        combined_mineral_list=[
            cfs(),
            crdi(),
            cess(),
            cen(),
            cfm(),
            fm(),
            odi(),
            cren(),
            mess(),
        ],
    )
