# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
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
import itertools

from collections import OrderedDict
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

from burnman import minerals
from burnman.composition_fitting import fit_composition_to_solution
from burnman.composition_fitting import fit_phase_proportions_to_bulk_composition

if __name__ == "__main__":

    """
    Example 1
    ---------
    Fitting single-phase EPMA data to a solid solution model

    The following example takes a set of variables
    (either elements such as "Fe", site_species such as "Fef_B",
    which would correspond to a species labelled Fef on the second site,
    or user-defined variables which are arithmetic sums of
    elements and/or site_species) and their values and covariances,
    and fits these values to a set of molar fractions of a solution model.
    """

    # Here, we use the Jennings and Holland (2015) garnet,
    # which has five endmembers: pyrope, almandine, grossular,
    # andradite and knorringite
    gt = minerals.JH_2015.garnet()

    print(f'Endmembers: {gt.endmember_names}')
    print(f'Elements: {gt.elements}')
    print('Stoichiometric matrix:')
    print(gt.stoichiometric_matrix)
    print()

    # The variables we are trying to fit are Fe (total), Ca, Mg
    # Cr, Al, Si and Fe3+, all given in mole amounts.
    fitted_variables = ['Fe', 'Ca', 'Mg', 'Cr', 'Al', 'Si', 'Fe3+']
    variable_values = np.array([1.1, 2., 0., 0, 1.9, 3., 0.1])
    variable_covariances = np.eye(7)*0.01*0.01

    # Add some noise.
    v_err = np.random.rand(7)
    np.random.seed(100)
    variable_values = np.random.multivariate_normal(variable_values,
                                                    variable_covariances)

    # Fe3+ isn't an element or a site-species of the solution model,
    # so we need to provide the linear conversion from Fe3+ to
    # elements and/or site species. In this case, Fe3+ resides only on
    # the second site, and the JH_2015.gt model has labelled Fe3+ on that
    # site as Fef. Therefore, the conversion is simply Fe3+ = Fef_B.
    variable_conversions = {'Fe3+': {'Fef_B': 1.}}

    # The next line does the heavy lifting
    popt, pcov, res = fit_composition_to_solution(gt,
                                                  fitted_variables,
                                                  variable_values,
                                                  variable_covariances,
                                                  variable_conversions)

    # We can set the composition of gt using the optimized parameters
    gt.set_composition(popt)

    # Print the optimized parameters and principal uncertainties
    print('Molar fractions:')
    for i in range(len(popt)):
        print(f'{gt.endmember_names[i]}: '
              f'{gt.molar_fractions[i]:.3f} +/- '
              f'{np.sqrt(pcov[i][i]):.3f}')
    print()
    print(f'Weighted residual: {res:.3f}')


    """
    Example 2
    ---------
    Fitting multiphase EPMA data to a bulk composition

    The following example takes a set of compositions (in wt %) of phases
    determined by EPMA experiments from high pressure experiments
    by Bertka and Fei (1997). It then calculates the mass fractions
    of the phases which best-fit the bulk composition reported in that paper,
    which corresponds to an estimate of the composition of Mars' mantle.
    """

    # Load and transpose input data
    filename = '../burnman/data/input_fitting/Bertka_Fei_1997_mars_mantle.dat'
    with open(filename) as f:
        column_names = f.readline().strip().split()[1:]
    data = np.genfromtxt(filename, dtype=None, encoding='utf8')
    data = list(map(list, itertools.zip_longest(*data, fillvalue=None)))

    # The first six columns are compositions given in weight % oxides
    compositions = np.array(data[:6])

    # The first row is the bulk composition
    bulk_composition = compositions[:,0]

    # Load all the data into a dictionary
    data = {column_names[i]: np.array(data[i])
            for i in range(len(column_names))}

    # Make ordered lists of samples (i.e. experiment ID) and phases
    samples = []
    phases = []
    for i in range(len(data['sample'])):
        if data['sample'][i] not in samples:
            samples.append(data['sample'][i])
        if data['phase'][i] not in phases:
            phases.append(data['phase'][i])

    samples.remove("bulk_composition")
    phases.remove("bulk")

    # Get the indices of all the phases present in each sample
    sample_indices = [[i for i in range(len(data['sample']))
                       if data['sample'][i] == sample]
                      for sample in samples]

    # Get the run pressures of each experiment
    pressures = [data['pressure'][indices[0]] for indices in sample_indices]

    # Create empty arrays to store the weight proportions of each phase,
    # and the principal uncertainties (we do not use the covariances here,
    # although they are calculated)
    weight_proportions = np.zeros((len(samples), len(phases)))*np.NaN
    weight_proportion_uncertainties = np.zeros((len(samples),
                                                len(phases)))*np.NaN

    # Loop over the samples, fitting phase proportions
    # to the provided bulk composition
    for i, sample in enumerate(samples):
        # This line does the heavy lifting
        popt, pcov, res = fit_phase_proportions_to_bulk_composition(compositions[:,sample_indices[i]],
                                                                    bulk_composition)

        # Fill the correct elements of the weight_proportions
        # and weight_proportion_uncertainties arrays
        sample_phases = [data['phase'][i] for i in sample_indices[i]]
        for j, phase in enumerate(sample_phases):
            weight_proportions[i, phases.index(phase)] = popt[j]
            weight_proportion_uncertainties[i, phases.index(phase)] = np.sqrt(pcov[j][j])

    # Plot the data
    for i, phase in enumerate(phases):
        ebar = plt.errorbar(pressures, weight_proportions[:, i],
                            yerr=weight_proportion_uncertainties[:, i],
                            fmt="none", zorder=2)
        plt.scatter(pressures, weight_proportions[:, i], label=phase, zorder=3)

    plt.xlim(0., 40.)
    plt.ylim(0., 1.)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Phase fraction (wt %)')
    plt.legend()
    plt.show()
