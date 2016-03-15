# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This module provides the functions required to process the standard burnman formula compositions
# ProcessChemistry returns the number of atoms and molar mass of a compound given its unit formula as an argument.
# process_solution_chemistry returns information required to calculate
# solid solution properties from a set of endmember formulae

from __future__ import absolute_import
import re
import numpy as np
from fractions import Fraction
import pkgutil


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


def formula_mass(formula, atomic_masses):
    """
    A function to take chemical formula and atomic mass
    dictionaries and
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
    solution_formulae = dict()
    s = re.split(r'\[', formula)[1:]
    sites = [[] for i in range(len(s))]
    list_occupancies = []
    list_multiplicity = np.empty(shape=(len(s)))
    n_occupancies = 0
    f = dict()

    for site in range(len(s)):
        site_occupancy = re.split(r'\]', s[site])[0]
        mult = re.split('[A-Z][^A-Z]*', re.split(r'\]', s[site])[1])[0]
        not_in_site = str(filter(None, re.split(r'\]', s[site])))[1]
        not_in_site = not_in_site.replace(mult, '', 1)
        if mult == '':
            list_multiplicity[site] = 1.0
        else:
            list_multiplicity[site] = mult

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


def process_solution_chemistry(formulae):
    """
    This function parses a set of endmember formulae
    containing site information, e.g.

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
    formulae : list of strings
        List of chemical formulae with site information

    Returns
    -------
    solution_formulae : list of dictionaries
        List of endmember formulae is output from site formula strings

    n_sites : integer
        Number of sites in the solid solution.
        Should be the same for all endmembers.

    sites : list of lists of strings
        A list of elements for each site in the solid solution

    n_occupancies : integer
        Sum of the number of possible elements on each of the sites
        in the solid solution.
        Example: A binary solution [[A][B],[B][C1/2D1/2]] would have
        n_occupancies = 5, with two possible elements on
        Site 1 and three on Site 2

    endmember_occupancies : 2d array of floats
        A 1D array for each endmember in the solid solution,
        containing the number of atoms of each element on each site.

    site_multiplicities : array of floats
        The number of each site per formula unit
        To simplify computations later, the multiplicities
        are repeated for each element on each site

    """
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
    for endmember in range(n_endmembers):
        solution_formula = dict()
        list_occupancies.append([[0] * len(sites[site])
                                for site in range(n_sites)])
        s = re.split(r'\[', formulae[endmember])[1:]

        for site in range(n_sites):
            site_split = re.split(r'\]', s[site])
            site_occupancy = site_split[0]

            mult = re.split('[A-Z][^A-Z]*', site_split[1])[0]
            if mult == '':
                list_multiplicity[site] = 1.0
            else:
                list_multiplicity[site] = float(mult)

            # Loop over elements on a site
            elements = re.findall('[A-Z][^A-Z]*', site_occupancy)

            for i in range(len(elements)):

                # Find the element and proportion on the site
                element_split = re.split('([0-9][^A-Z]*)', elements[i])
                element_on_site = element_split[0]
                if len(element_split) == 1:
                    proportion_element_on_site = Fraction(1.0)
                else:
                    proportion_element_on_site = Fraction(element_split[1])

                solution_formula[element_on_site] = solution_formula.get(
                    element_on_site, 0.0) + list_multiplicity[site] * proportion_element_on_site

                if element_on_site not in sites[site]:
                    n_occupancies += 1
                    sites[site].append(element_on_site)
                    element_index = sites[site].index(element_on_site)
                    for parsed_mbr in range(len(list_occupancies)):
                        list_occupancies[parsed_mbr][site].append(0)
                else:
                    element_index = sites[site].index(element_on_site)
                list_occupancies[endmember][site][
                    element_index] = proportion_element_on_site

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
    for endmember in range(n_endmembers):
        n_element = 0
        for site in range(n_sites):
            for element in range(len(list_occupancies[endmember][site])):
                endmember_occupancies[endmember][
                    n_element] = list_occupancies[endmember][site][element]
                site_multiplicities[n_element] = list_multiplicity[site]
                n_element += 1

    return solution_formulae, n_sites, sites, n_occupancies, endmember_occupancies, site_multiplicities


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
