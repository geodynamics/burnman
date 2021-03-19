# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_composition_fitting
----------------------

This example shows how to fit compositional data to a solution model,
how to partition a bulk composition between phases of known composition,
and how to assess goodness of fit.

*Uses:*

* :doc:`composition_fitting`
* :class:`burnman.Composition`
* :class:`burnman.solidsolution.SolidSolution`
* :class:`burnman.solutionmodel.SolutionModel`


*Demonstrates:*

* Fitting compositional data to a solution model
* Partitioning of a bulk composition between phases of known composition
* Assessing goodness of fit.

"""
from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

from burnman import minerals
from burnman.composition_fitting import fit_composition_to_solution
from burnman.composition_fitting import fit_phase_proportions_to_bulk_composition

import pandas as pd

if __name__ == "__main__":

    """
    CFMASCrO
    """
    gt = minerals.JH_2015.garnet()

    fitted_variables = ['Fe', 'Ca', 'Mg', 'Cr', 'Al', 'Si', 'Fe3+']
    variable_values = [1.1, 2., 0., 0, 1.9, 3., 0.1]
    variable_covariances = np.eye(7)
    variable_conversions = {'Fe3+': {'Fef_B': 1.}}
    popt, pcov, res = fit_composition_to_solution(gt,
                                                  fitted_variables,
                                                  variable_values,
                                                  variable_covariances,
                                                  variable_conversions)

    gt.set_composition(popt)
    print(gt.formula)
    print(popt, res)

    mars = pd.read_table('../burnman/data/input_fitting/Bertka_Fei_1997_mars_mantle.dat', sep=" ")

    samples = list(mars['sample'].unique())
    samples.remove("bulk_composition")

    pressures = [mars.loc[mars['sample'] == s]['pressure'].to_numpy()[0] for s in samples]

    phases = list(mars['phase'].unique())
    phases.remove("bulk")

    weight_proportions = np.zeros((len(samples), len(phases)))
    weight_proportions[:] = np.NaN

    bulk = mars.loc[mars['sample'] == "bulk_composition"]

    b = np.array([bulk['Na2O'], bulk['CaO'], bulk['FeO'],
                  bulk['MgO'], bulk['Al2O3'], bulk['SiO2']])[:,0]

    for i, sample in enumerate(samples):
        c = mars.loc[mars['sample'] == sample]

        A = np.array([c['Na2O'], c['CaO'], c['FeO'],
                      c['MgO'], c['Al2O3'], c['SiO2']])

        popt, pcov, res = fit_phase_proportions_to_bulk_composition(A, b)

        for j, phase in enumerate(c['phase']):
            weight_proportions[i, phases.index(phase)] = popt[j]

    for i, phase in enumerate(phases):
        plt.plot(pressures, weight_proportions[:, i], label=phase)
        plt.scatter(pressures, weight_proportions[:, i])

    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Weight fraction')
    plt.legend()
    plt.show()
