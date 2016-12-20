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
from collections import Counter

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
    s = re.split(r'\[', formula)[1:]
    list_multiplicity = np.empty(shape=(len(s)))
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


def component_to_atom_fractions(composition_dictionary, fraction_type):
    component_formulae = [dictionarize_formula(c) for c in composition_dictionary]

    if fraction_type == 'weight':
        component_masses = np.array([formula_mass(f, read_masses()) for f in component_formulae])
        n_moles = composition_dictionary.values() /component_masses
    elif fraction_type == 'molar':
        n_moles = composition_dictionary.values()
    else:
        raise Exception('Fraction type not recognised. Should be weight or molar')

    molar_composition = [{key: f[key]*m for key in f} for f, m in zip(*[component_formulae, n_moles])]
    atomic_fractions = sum((Counter(c) for c in molar_composition), Counter())
    
    n_atoms = sum(atomic_fractions.values())
    
    return {key: atomic_fractions[key]/n_atoms for key in atomic_fractions}


def binary_composition(composition1, composition2, x):
    """
    Returns the composition within a binary system that
    is defined by a fraction of the second composition.
    
    Parameters
    ----------
    composition1 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    composition2 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    x : float
        Composition as a fraction of composition2.

    Returns
    -------
    composition : dictionary of floats
        (1-x)*composition1 + x*composition2.
    """
    composition = {}
    for key, value in composition1.items():
        composition[key] = value*(1. - x)
    for key, value in composition2.items():
        if key in composition:
            composition[key] += value*x
        else:
            composition[key] = value*x

    return composition


def _cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    Shamelessly taken from stackoverflow (question 1208118)

    Parameters
    ----------
    arrays : list of 1-D arrays 
        Arrays from which to form the cartesian product.
    out : ndarray
        Array in which to place the cartesian product.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing 
        cartesian products formed of input arrays.
    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        _cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def _indicize(inp_arr):
    arr=np.zeros([len(inp_arr), 1+max(inp_arr.flat)])
    for i, site_indices in enumerate(inp_arr):
        for j in site_indices:
            arr[i][j] = 1
    return arr

def dependent_endmembers(formulae):
    """
    Find the set of possible dependent endmembers from a set of site formulae.
    Dependent endmembers are those which can be made from a set of independent 
    endmembers. For example, the simple alkali salts [Na][Cl], [K][Cl] 
    and [Na][I] form an independent set, with a single dependent endmember 
    [K][I] (equal to [K][Cl] + [Na][I] - [Na][Cl]). 
    If we also have the independent endmember [Cs][F], no more dependent 
    endmembers can be constructed; for example, [Cs][I] cannot be made because
    this would also require an additional independent endmember 
    ([Na][F] or [K][F]).

    Parameters
    ----------
    formulae : list of dictionaries
        Site formulae (see description in processchemistry)

    Returns
    -------
    dependent_endmembers : 2d numpy array
        2-D array containing the reaction coefficients of the independent 
        endmembers required to create the dependent endmembers.
    """

    solution_formulae, n_sites, sites, n_occupancies, endmember_occupancies, site_multiplicities = process_solution_chemistry(formulae)

    # First, we redefine our endmembers using unique site occupancy indices
    # We can then find the set of unique occupancies for each site
    # site indices contains the unique integer indices for the occupancies on
    # each site, indices contains the sets of indices for each site.
    i=0
    n_sites = 0
    indices = []
    site_indices = []
    
    for site in sites:
        site_indices.append([])
        site_occupancies = map(tuple, endmember_occupancies[:, i:i+len(site)])
        set_site_occupancies = map(np.array, set(site_occupancies))
        
        for site_occupancy in site_occupancies:
            site_indices[-1].append([i+n_sites for i, set_occupancy in enumerate(set_site_occupancies) if np.array_equal(set_occupancy, site_occupancy)][0])
        indices.append(np.arange(n_sites,n_sites+len(set_site_occupancies)))
        
        i += len(site)
        n_sites += len(set_site_occupancies)

    given_members = _indicize(np.array(zip(*site_indices)))

    # All the potential endmembers can be created from permutations of the unique site occupancies
    # Note that the solution may not allow all of these endmembers; for example, if there
    # is an endmember with unique occupancies on two sites,
    # (e.g. Ca2+ and Fe3+ in the X, Y sites in garnet)
    # then it will not be possible to exchange Ca2+ for Mg2+ on the X site
    # and retain Fe3+ on the Y site.
    all_members=_indicize(_cartesian(indices))

    # Now we return only the endmembers which are not contained within the user-provided
    # independent set
    a=np.concatenate((given_members, all_members), axis=0)
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx, cnts = np.unique(b, return_index=True, return_counts=True)
    potential_dependents = a[[i for i, cnt in zip(*[idx, cnts]) if cnt==1]]
    
    
    # We only want the dependent endmembers which can be described by a
    # linear combination of the independent set
    tol = 1e-12
    dependent_endmembers = np.array([]).reshape(0, len(formulae))
    for mbr in potential_dependents:
        a,resid,rank,s = np.linalg.lstsq(given_members.T, mbr)
        if resid[0] < tol:
            a.real[abs(a.real) < tol] = 0.0
            dependent_endmembers = np.concatenate((dependent_endmembers, np.array([a])), axis=0)
    return dependent_endmembers

def compositional_variables(assemblage, indices):
    """
    Takes an assemblage and outputs names for the
    compositional variables which describe the bulk composition
    and compositions of all the phases.

    Parameters
    ----------
    assemblage : composite
        The assemblage of minerals for which we want to find
        amounts and compositions.
    
    Returns
    -------
    var_names : list of strings
        Strings are provided in the same order as the X part
        of the PTX variable input to :func:`_set_eqns`. Phase
        amount names are given as 'x(phase.name)', while
        the molar fractions of endmembers in each phase are
        given as 'p(endmember.name)'
    """
    old_i = -1
    var_names=[]
    for (i, j) in indices:
        if i != old_i:
            var_names.append('x('+assemblage.phases[i].name+')')
            old_i = i
        else:
            var_names.append('p('+assemblage.phases[i].endmembers[j][0].name+')')

    return var_names
    
