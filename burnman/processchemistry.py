# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This module provides the functions required to process the standard burnman formula compositions
# ProcessChemistry returns the number of atoms and molar mass of a compound given its unit formula as an argument.
# process_solution_chemistry returns information required to calculate
# solid solution properties from a set of endmember formulae

from __future__ import absolute_import
import re
import numpy as np
from fractions import Fraction
from collections import Counter
import pkgutil
from string import ascii_uppercase as ucase
from scipy.optimize import nnls

def simplify_matrix(arr):
    def f(i,j):
        return nsimplify(arr[i][j])
    return Matrix( len(arr), len(arr[0]), f )

def read_masses():
    """
    A simple function to read a file with a two column list of
    elements and their masses into a dictionary
    """
    datastream = pkgutil.get_data(
        'burnman', 'data/input_masses/atomic_masses.dat')
    datalines = [line.strip()
                 for line in datastream.decode('ascii').split('\n') if line.strip()]
    lookup = dict()
    for line in datalines:
        data = "%".join(line.split("%")[:1]).split()
        if data != []:
            lookup[data[0]] = float(data[1])
    return lookup

atomic_masses = read_masses()

def dictionarize_formula(formula):
    """
    A function to read a chemical formula string and
    convert it into a dictionary
    """
    f = dict()
    elements = re.findall('[A-Z][^A-Z]*', formula)
    for element in elements:
        element_name = re.split('[0-9][^A-Z]*', element)[0]
        element_atoms = re.findall('[0-9][^A-Z]*', element)
        if len(element_atoms) == 0:
            element_atoms = Fraction(1.0)
        else:
            element_atoms = Fraction(element_atoms[0])
        f[element_name] = f.get(element_name, 0.0) + element_atoms

    return f

def sum_formulae(formulae, amounts=None):
    if amounts == None:
        amounts = [1. for formula in formulae]
    else:
        assert (len(formulae) == len(amounts))

    summed_formula = Counter()
    for i, formula in enumerate(formulae):
        summed_formula.update(Counter({element: amounts[i] * n_atoms for (element, n_atoms) in formula.items()}))
    return summed_formula

def formula_mass(formula):
    """
    A function to take chemical formula and atomic mass
    dictionaries and compute the formula mass.
    """
    mass = sum(
        formula[element] * atomic_masses[element] for element in formula)
    return mass

def dictionarize_site_formula(formula):
    """
    A function to take a chemical formula with sites specified
    by square brackets and return a standard dictionary with
    element keys and atoms of each element per formula unit as items.
    """
    s = re.split(r'\[', formula)[1:]
    list_multiplicity = np.empty(shape=(len(s)))
    f = dict()

    for site in range(len(s)):
        site_occupancy = re.split(r'\]', s[site])[0]
        mult = re.split('[A-Z][^A-Z]*', re.split(r'\]', s[site])[1])[0]
        not_in_site = str(filter(None, re.split(r'\]', s[site])))[1]
        not_in_site = not_in_site.replace(mult, '', 1)
        if mult == '':
            list_multiplicity[site] = Fraction(1.0)
        else:
            list_multiplicity[site] = Fraction(mult)

        # Loop over elements on a site
        elements = re.findall('[A-Z][^A-Z]*', site_occupancy)
        for i in range(len(elements)):
            element_on_site = re.split('[0-9][^A-Z]*', elements[i])[0]
            proportion_element_on_site = re.findall(
                '[0-9][^A-Z]*', elements[i])
            if len(proportion_element_on_site) == 0:
                proportion_element_on_site = Fraction(1.0)
            else:
                proportion_element_on_site = Fraction(
                    proportion_element_on_site[0])
            n_element = float(mult) * proportion_element_on_site
            f[element_on_site] = f.get(element_on_site, 0.0) + n_element

        # Loop over elements not on a site
        for enamenumber in re.findall('[A-Z][^A-Z]*', not_in_site):
            element = str(filter(None, re.split(r'(\d+)', enamenumber)))
            f[element[0]] = f.get(element[0], 0.0) + float(element[1])

    return f

def solution_bounds(endmember_occupancies):
    """
    Parameters
    ----------
    endmember_occupancies : 2d array of floats
        A 1D array for each endmember in the solid solution,
        containing the number of atoms of each element on each site.

    Returns
    -------
    solution_bounds : 2d array of floats
        An abbreviated version of endmember_occupancies, 
        where the columns represent the independent compositional 
        bounds on the solution
    """
    # Find bounds for the solution
    i_sorted =zip(*sorted([(i,
                            sum([1 for val in endmember_occupancies.T[i]
                                 if val>1.e-10]))
                           for i in range(len(endmember_occupancies.T))
                                          if np.any(endmember_occupancies.T[i] > 1.e-10)],
                          key=lambda x: x[1]))[0]

    solution_bounds = endmember_occupancies[:,i_sorted[0],np.newaxis]
    for i in i_sorted[1:]:
        if np.abs(nnls(solution_bounds, endmember_occupancies.T[i])[1]) > 1.e-10:
            solution_bounds = np.concatenate((solution_bounds,
                                              endmember_occupancies[:,i,np.newaxis]),
                                             axis=1)
    return solution_bounds

def process_solution_chemistry(solution_model):
    """
    This function parses a class instance with a "formulas"
    attribute containing site information, e.g.

        [ '[Mg]3[Al]2Si3O12', '[Mg]3[Mg1/2Si1/2]2Si3O12' ]

    It outputs the bulk composition of each endmember
    (removing the site information), and also a set of
    variables and arrays which contain the site information.
    These are output in a format that can easily be used to
    calculate activities and gibbs free energies, given
    molar fractions of the phases and pressure
    and temperature where necessary.

    Parameters
    ----------
    solution_model : instance of class
        Class must have a "formulas" attribute, containing a
        list of chemical formulae with site information

    Creates attributes
    ------------------
    solution_formulae : list of dictionaries
        List of endmember formulae is output from site formula strings

    n_sites : integer
        Number of sites in the solid solution.
        Should be the same for all endmembers.

    sites : list of lists of strings
        A list of elements for each site in the solid solution

    site_names : list of strings
        A list of elements_site pairs in the solid solution, where
        each distinct site is given by a unique uppercase letter
        e.g. ['Mg_A', 'Fe_A', 'Al_A', 'Al_B', 'Si_B']

    n_occupancies : integer
        Sum of the number of possible elements on each of the sites
        in the solid solution.
        Example: A binary solution [[A][B],[B][C1/2D1/2]] would have
        n_occupancies = 5, with two possible elements on
        Site 1 and three on Site 2

    site_multiplicities : array of floats
        The number of each site per formula unit
        To simplify computations later, the multiplicities
        are repeated for each element on each site

    endmember_occupancies : 2d array of floats
        A 1D array for each endmember in the solid solution,
        containing the fraction of atoms of each element on each site.

    endmember_noccupancies : 2d array of floats
        A 1D array for each endmember in the solid solution,
        containing the number of atoms of each element on each site
        per mole of endmember.

    """
    formulae = solution_model.formulas
    n_sites = formulae[0].count('[')
    n_endmembers = len(formulae)

    # Check the number of sites is the same for all endmembers
    for i in range(n_endmembers):
        assert(formulae[i].count('[') == n_sites)

    solution_formulae = []
    sites = [[] for i in range(n_sites)]
    list_occupancies = []
    list_multiplicity = np.empty(shape=(n_sites))
    n_occupancies = 0

    # Number of unique site occupancies (e.g.. Mg on X etc.)
    for i_mbr in range(n_endmembers):
        solution_formula = dict()
        list_occupancies.append([[0] * len(sites[site])
                                for site in range(n_sites)])
        s = re.split(r'\[', formulae[i_mbr])[1:]

        for i_site, site_string in enumerate(s):
            site_split = re.split(r'\]', site_string)
            site_occupancy = site_split[0]

            mult = re.split('[A-Z][^A-Z]*', site_split[1])[0]
            if mult == '':
                list_multiplicity[i_site] = Fraction(1.0)
            else:
                list_multiplicity[i_site] = Fraction(mult)

            # Loop over elements on a site
            elements = re.findall('[A-Z][^A-Z]*', site_occupancy)

            for element in elements:

                # Find the element and proportion on the site
                element_split = re.split('([0-9][^A-Z]*)', element)
                element_on_site = element_split[0]
                if len(element_split) == 1:
                    proportion_element_on_site = Fraction(1.0)
                else:
                    proportion_element_on_site = Fraction(element_split[1])

                solution_formula[element_on_site] = solution_formula.get(
                    element_on_site, 0.0) + list_multiplicity[i_site] * proportion_element_on_site

                if element_on_site not in sites[i_site]:
                    n_occupancies += 1
                    sites[i_site].append(element_on_site)
                    i_el = sites[i_site].index(element_on_site)
                    for parsed_mbr in range(len(list_occupancies)):
                        list_occupancies[parsed_mbr][i_site].append(0)
                else:
                    i_el = sites[i_site].index(element_on_site)
                list_occupancies[i_mbr][i_site][i_el] = proportion_element_on_site

            # Loop over elements after site
            if len(site_split) != 1:
                not_in_site = str(filter(None, site_split[1]))
                not_in_site = not_in_site.replace(mult, '', 1)
                for enamenumber in re.findall('[A-Z][^A-Z]*', not_in_site):
                    element = list(
                        filter(None, re.split(r'(\d+)', enamenumber)))
                    # Look up number of atoms of element
                    if len(element) == 1:
                        nel = 1.
                    else:
                        nel = float(float(element[1]))
                    solution_formula[element[0]] = solution_formula.get(
                        element[0], 0.0) + nel

        solution_formulae.append(solution_formula)

    # Site occupancies and multiplicities
    endmember_occupancies = np.empty(shape=(n_endmembers, n_occupancies))
    site_multiplicities = np.empty(shape=(n_occupancies))
    for i_mbr in range(n_endmembers):
        n_element = 0
        for i_site in range(n_sites):
            for i_el in range(len(list_occupancies[i_mbr][i_site])):
                endmember_occupancies[i_mbr][
                    n_element] = list_occupancies[i_mbr][i_site][i_el]
                site_multiplicities[n_element] = list_multiplicity[i_site]
                n_element += 1

    # Site names
    solution_model.site_names = []
    for i, elements in enumerate(sites):
        for element in elements:
            solution_model.site_names.append('{0}_{1}'.format(element, ucase[i]))

    # Finally, make attributes for solution model instance:
    solution_model.solution_formulae = solution_formulae
    solution_model.n_sites = n_sites
    solution_model.sites = sites
    solution_model.n_occupancies = n_occupancies
    solution_model.site_multiplicities = site_multiplicities
    solution_model.endmember_occupancies = endmember_occupancies
    solution_model.endmember_noccupancies = np.einsum('ij, j->ij', endmember_occupancies, site_multiplicities)

def compositional_array(formulae):
    """
    Parameters
    ----------
    formulae : list of dictionaries
        List of chemical formulae

    Returns
    -------
    formula_array : 2D array of floats
        Array of endmember formulae

    elements : List of strings
        List of elements
    """
    elements = []
    for formula in formulae:
        for element in formula:
            if element not in elements:
                elements.append(element)

    formula_array = ordered_compositional_array(formulae, elements)

    return formula_array, elements


def ordered_compositional_array(formulae, elements):
    """
    Parameters
    ----------
    formulae : list of dictionaries
        List of chemical formulae

    elements : List of strings
        List of elements

    Returns
    -------
    formula_array : 2D array of floats
        Array of endmember formulae
    """
    formula_array = np.zeros(shape=(len(formulae), len(elements)))
    for idx, formula in enumerate(formulae):
        for element in formula:
            assert(element in elements)
            formula_array[idx][elements.index(element)] = formula[element]

    return formula_array
