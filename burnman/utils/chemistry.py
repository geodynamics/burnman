# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

# This module provides the functions required to process the
# standard burnman formula formats.
# tools.chemistry returns the number of atoms and molar mass of a compound
# given its unit formula as an argument.
# process_solution_chemistry returns information required to calculate
# solution properties from a set of endmember formulae

from __future__ import absolute_import
import re
import numpy as np
from fractions import Fraction
from collections import Counter
import pkgutil
from string import ascii_uppercase as ucase
from sympy import nsimplify


def read_masses():
    """
    A simple function to read a file with a two column list of
    elements and their masses into a dictionary
    """
    datastream = pkgutil.get_data("burnman", "data/input_masses/atomic_masses.dat")
    datalines = [
        line.strip() for line in datastream.decode("ascii").split("\n") if line.strip()
    ]
    lookup = dict()
    for line in datalines:
        data = "%".join(line.split("%")[:1]).split()
        if data != []:
            lookup[data[0]] = float(data[1])
    return lookup


"""
atomic_masses is a dictionary of atomic masses
"""
atomic_masses = read_masses()

"""
IUPAC_element_order provides a list of all the elements.
Element order is based loosely on electronegativity,
following the scheme suggested by IUPAC, except that H
comes after the Group 16 elements, not before them.
"""
IUPAC_element_order = [
    "v",
    "Og",
    "Rn",
    "Xe",
    "Kr",
    "Ar",
    "Ne",
    "He",  # Group 18
    "Fr",
    "Cs",
    "Rb",
    "K",
    "Na",
    "Li",  # Group 1 (not H)
    "Ra",
    "Ba",
    "Sr",
    "Ca",
    "Mg",
    "Be",  # Group 2
    "Lr",
    "No",
    "Md",
    "Fm",
    "Es",
    "Cf",
    "Bk",
    "Cm",
    "Am",
    "Pu",
    "Np",
    "U",
    "Pa",
    "Th",
    "Ac",  # Actinides
    "Lu",
    "Yb",
    "Tm",
    "Er",
    "Ho",
    "Dy",
    "Tb",
    "Gd",
    "Eu",
    "Sm",
    "Pm",
    "Nd",
    "Pr",
    "Ce",
    "La",  # Lanthanides
    "Y",
    "Sc",  # Group 3
    "Rf",
    "Hf",
    "Zr",
    "Ti",  # Group 4
    "Db",
    "Ta",
    "Nb",
    "V",  # Group 5
    "Sg",
    "W",
    "Mo",
    "Cr",  # Group 6
    "Bh",
    "Re",
    "Tc",
    "Mn",  # Group 7
    "Hs",
    "Os",
    "Ru",
    "Fe",  # Group 8
    "Mt",
    "Ir",
    "Rh",
    "Co",  # Group 9
    "Ds",
    "Pt",
    "Pd",
    "Ni",  # Group 10
    "Rg",
    "Au",
    "Ag",
    "Cu",  # Group 11
    "Cn",
    "Hg",
    "Cd",
    "Zn",  # Group 12
    "Nh",
    "Tl",
    "In",
    "Ga",
    "Al",
    "B",  # Group 13
    "Fl",
    "Pb",
    "Sn",
    "Ge",
    "Si",
    "C",  # Group 14
    "Mc",
    "Bi",
    "Sb",
    "As",
    "P",
    "N",  # Group 15
    "Lv",
    "Po",
    "Te",
    "Se",
    "S",
    "O",  # Group 16
    "H",  # hydrogen
    "Ts",
    "At",
    "I",
    "Br",
    "Cl",
    "F",
]  # Group 17


def dictionarize_formula(formula):
    """
    A function to read a chemical formula string and
    convert it into a dictionary

    :param formula: Chemical formula, written in the XnYm format, where
        the formula has n atoms of element X and m atoms of element Y
    :type formula: str

    :returns: The same chemical formula, but expressed as a dictionary.
    :rtype: dict
    """
    f = dict()
    elements = re.findall("[A-Z][^A-Z]*", formula)
    for element in elements:
        element_name = re.split("[0-9][^A-Z]*", element)[0]
        element_atoms = re.findall("[0-9][^A-Z]*", element)
        if len(element_atoms) == 0:
            element_atoms = Fraction(1.0)
        else:
            element_atoms = Fraction(element_atoms[0])
        f[element_name] = f.get(element_name, 0.0) + element_atoms

    return f


def sum_formulae(formulae, amounts=None):
    """
    Adds together a set of formulae.

    :param formulae: List of chemical formulae.
    :type formulae: list of dictionary or counter objects
    :param amounts: List of amounts of each formula.
    :type amounts: list of floats

    :returns: The sum of the user-provided formulae
    :rtype: Counter object
    """
    if amounts is None:
        amounts = [1.0 for formula in formulae]
    else:
        assert len(formulae) == len(amounts)

    summed_formula = Counter()
    for i, formula in enumerate(formulae):
        summed_formula.update(
            Counter(
                {
                    element: amounts[i] * n_atoms
                    for (element, n_atoms) in formula.items()
                }
            )
        )
    return summed_formula


def formula_mass(formula):
    """
    A function to take a chemical formula and compute the formula mass.

    :param formula: A chemical formula
    :type formula: dict or Counter object

    :returns: The mass per mole of formula [kg]
    :rtype: float
    """
    mass = sum(formula[element] * atomic_masses[element] for element in formula)
    return mass


def convert_formula(formula, to_type="mass", normalize=False):
    """
    Converts a chemical formula from one type (mass or molar)
    into the other. Renormalises amounts if normalize=True.

    :param formula: A chemical formula.
    :type formula: dict or Counter object

    :param to_type: Conversion type, one of 'mass' or 'molar'.
    :type to_type: str

    :param normalize: Whether or not to normalize the converted formula to 1.
    :type normalize: bool

    :returns: The converted formula.
    :rtype: dict
    """

    if to_type == "mass":
        f = {
            element: n_atoms * atomic_masses[element]
            for (element, n_atoms) in formula.items()
        }
    elif to_type == "molar":
        f = {
            element: n_atoms / atomic_masses[element]
            for (element, n_atoms) in formula.items()
        }
    else:
        raise Exception(
            "Value of parameter to_type not recognised. "
            'Should be either "mass" or "molar".'
        )

    if normalize:
        s = np.sum([n for (element, n) in f.items()])
        f = {element: n / s for (element, n) in f.items()}

    return f


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

    :param solution_model: Class must have a "formulas" attribute,
        containing a list of chemical formulae with site information
    :type solution model: instance of class

    :rtype: None

    .. note:: Nothing is returned from this function, but the solution_model
        object gains the following attributes:

        * solution_formulae [list of dictionaries]
            List of endmember formulae in dictionary form.
        * empty_formula [string]
            Abbreviated chemical formula with sites denoted by empty
            square brackets.
        * general_formula [string]
            General chemical formula with sites denoted by
            square brackets filled with a comma-separated list of species
        * n_sites [integer]
            Number of sites in the solution.
            Should be the same for all endmembers.
        * sites [list of lists of strings]
            A list of species for each site in the solution.
        * site_names [list of strings]
            A list of species_site pairs in the solution, where
            each distinct site is given by a unique uppercase letter
            e.g. ['Mg_A', 'Fe_A', 'Al_A', 'Al_B', 'Si_B'].
        * n_occupancies [integer]
            Sum of the number of possible species on each of the sites
            in the solution.
            Example: A binary solution [[A][B],[B][C1/2D1/2]] would have
            n_occupancies = 5, with two possible species on
            Site 1 and three on Site 2.
        * site_multiplicities [2D array of floats]
            A 1D array for each endmember in the solution,
            containing the multiplicities of each site per formula unit.
            To simplify computations later, the multiplicities
            are repeated for each species on each site, so the shape of
            this attribute is (n_endmembers, n_site_species).
        * endmember_occupancies [2d array of floats]
            A 1D array for each endmember in the solution,
            containing the fraction of atoms of each species on each site.
        * endmember_noccupancies [2d array of floats]
            A 1D array for each endmember in the solution,
            containing the number of atoms of each species on each site
            per mole of endmember.
    """
    formulae = solution_model.formulas
    n_sites = formulae[0].count("[")
    n_endmembers = len(formulae)

    # Check the number of sites is the same for all endmembers
    if not np.all(np.array([f.count("[") for f in formulae]) == n_sites):
        raise Exception("All formulae must have the same " "number of distinct sites.")

    solution_formulae = [{} for i in range(n_endmembers)]
    sites = [[] for i in range(n_sites)]
    list_occupancies = []
    list_multiplicities = np.empty(shape=(n_endmembers, n_sites))
    n_occupancies = 0

    # Number of unique site occupancies (e.g.. Mg on X etc.)
    for i_mbr in range(n_endmembers):
        list_occupancies.append([[0] * len(sites[site]) for site in range(n_sites)])
        s = re.split(r"\[", formulae[i_mbr])[1:]

        for i_site, site_string in enumerate(s):
            site_split = re.split(r"\]", site_string)
            site_occupancy = site_split[0]

            mult = re.split("[A-Z][^A-Z]*", site_split[1])[0]
            if mult == "":
                list_multiplicities[i_mbr][i_site] = Fraction(1.0)
            else:
                list_multiplicities[i_mbr][i_site] = Fraction(mult)

            # Loop over species on a site
            species = re.findall("[A-Z][^A-Z]*", site_occupancy)

            for sp in species:
                # Find the species and its proportion on the site
                species_split = re.split("([0-9][^A-Z]*)", sp)
                name_of_species = species_split[0]
                if len(species_split) == 1:
                    proportion_species_on_site = Fraction(1.0)
                else:
                    proportion_species_on_site = Fraction(species_split[1])

                solution_formulae[i_mbr][name_of_species] = solution_formulae[
                    i_mbr
                ].get(name_of_species, 0.0) + (
                    list_multiplicities[i_mbr][i_site] * proportion_species_on_site
                )

                if name_of_species not in sites[i_site]:
                    n_occupancies += 1
                    sites[i_site].append(name_of_species)
                    i_el = sites[i_site].index(name_of_species)
                    for parsed_mbr in range(len(list_occupancies)):
                        list_occupancies[parsed_mbr][i_site].append(0)
                else:
                    i_el = sites[i_site].index(name_of_species)
                list_occupancies[i_mbr][i_site][i_el] = proportion_species_on_site

            # Loop over species after site
            if len(site_split) != 1:
                not_in_site = str(filter(None, site_split[1]))
                not_in_site = not_in_site.replace(mult, "", 1)
                for enamenumber in re.findall("[A-Z][^A-Z]*", not_in_site):
                    sp = list(filter(None, re.split(r"(\d+)", enamenumber)))
                    # Look up number of atoms of element
                    if len(sp) == 1:
                        nel = 1.0
                    else:
                        nel = float(float(sp[1]))
                    solution_formulae[i_mbr][sp[0]] = (
                        solution_formulae[i_mbr].get(sp[0], 0.0) + nel
                    )

    # Site occupancies and multiplicities
    endmember_occupancies = np.empty(shape=(n_endmembers, n_occupancies))
    site_multiplicities = np.empty(shape=(n_endmembers, n_occupancies))

    for i_mbr in range(n_endmembers):
        n_species = 0
        for i_site in range(n_sites):
            for i_el in range(len(list_occupancies[i_mbr][i_site])):
                endmember_occupancies[i_mbr][n_species] = list_occupancies[i_mbr][
                    i_site
                ][i_el]
                site_multiplicities[i_mbr][n_species] = list_multiplicities[i_mbr][
                    i_site
                ]
                n_species += 1

    # Site names
    solution_model.site_names = []
    for i, species in enumerate(sites):
        for sp in species:
            solution_model.site_names.append("{0}_{1}".format(sp, ucase[i]))

    # Finally, make attributes for solution model instance:
    solution_model.solution_formulae = solution_formulae
    solution_model.n_sites = n_sites
    solution_model.sites = sites
    solution_model.site_multiplicities = site_multiplicities
    solution_model.n_occupancies = n_occupancies
    solution_model.endmember_occupancies = endmember_occupancies
    solution_model.endmember_noccupancies = np.einsum(
        "ij, ij->ij", endmember_occupancies, site_multiplicities
    )

    solution_model.empty_formula = re.sub(
        "([\\[]).*?([\\]])", "\\g<1>\\g<2>", solution_model.formulas[0]
    )
    split_empty = solution_model.empty_formula.split("[")
    solution_model.general_formula = split_empty[0]
    for i in range(n_sites):
        solution_model.general_formula += f"[{','.join(sites[i])}{split_empty[i+1]}"


def site_occupancies_to_strings(
    site_species_names, site_multiplicities, endmember_occupancies
):
    """
    Converts a list of endmember site occupancies into a list
    of string representations of those occupancies.

    :param site_species_names: A list of list of strings,
        giving the names of the species which reside on each site.
        List of sites, each of which contains a list of the species
        occupying each site.
    :type site_species_names: 2D list of strings

    :param site_multiplicities: List of floats giving the multiplicity
        of each site. If 2D, must have the same shape as endmember_occupancies.
        If 1D, must be either the same length as the number of sites, or
        the same length as site_species_names
        (with an implied repetition of the same
        number for each species on a given site).
    :type site_multiplicities: 1D or 2D numpy array of floats

    :param endmember_occupancies: A list of site-species occupancies
        for each endmember. The first dimension loops over the endmembers, and the
        second dimension loops over the site-species occupancies for that endmember.
        The total number and order of occupancies must
        be the same as the strings in site_species_names.
    :type endmember_occupancies: 2D numpy array of floats

    :returns: A list of strings in standard burnman format.
        For example, [Mg]3[Al]2 would correspond to the
        classic two-site pyrope garnet.
    :rtype: list of strings
    """

    site_multiplicities = np.array(site_multiplicities)
    endmember_occupancies = np.array(endmember_occupancies)
    n_endmembers = endmember_occupancies.shape[0]

    if len(site_multiplicities.shape) == 1:
        # Site multiplicities should either be given on a per-site basis,
        # or a per-species basis
        if len(site_species_names) == len(site_multiplicities):
            site_mults = []

            for i, site in enumerate(site_species_names):
                for species in site:
                    site_mults.append(site_multiplicities[i])

            site_multiplicities = np.array(site_mults)

        elif len(endmember_occupancies[0]) != len(site_multiplicities):
            raise Exception(
                "Site multiplicities should either be given "
                "on a per-site basis or a per-species basis"
            )

        site_multiplicities = np.einsum(
            "i, j->ij", np.ones(n_endmembers), site_multiplicities
        )
    elif len(site_multiplicities.shape) == 2:
        if site_multiplicities.shape != endmember_occupancies.shape:
            raise Exception(
                "If site_multiplicities is 2D, it should have "
                "the same shape as endmember_occupancies. "
                "They currently have shapes "
                f"{site_multiplicities.shape} and "
                f"{endmember_occupancies.shape}."
            )
    else:
        raise Exception("Site multiplicities should either be 1D or 2D.")

    site_formulae = []
    for i_mbr, mbr_occupancies in enumerate(endmember_occupancies):
        i = 0
        site_formulae.append("")
        for site in site_species_names:
            amounts = mbr_occupancies[i : i + len(site)]
            mult = site_multiplicities[i_mbr, i]
            if np.abs(mult - 1.0) < 1.0e-12:
                mult = ""
            else:
                mult = str(nsimplify(mult))
            amounts /= sum(amounts)
            site_occupancy = formula_to_string(dict(zip(site, amounts)))
            site_formulae[-1] += "[{0}]{1}".format(site_occupancy, mult)
            i += len(site)

    return site_formulae


def compositional_array(formulae):
    """
    :param formulae: List of chemical formulae
    :type formulae: list of dicts

    :returns: Array of endmember formulae and a list of elements.
    :rtype: 2D numpy.array of floats and a list of strs
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
    :param formulae: List of chemical formulae
    :type formulae: list of dicts

    :param elements : List of elements
    :type elements: list of strings

    :returns: Array of endmember formulae
    :rtype: 2D array of floats
    """
    formula_array = np.zeros(shape=(len(formulae), len(elements)))
    for idx, formula in enumerate(formulae):
        for element in formula:
            assert element in elements
            formula_array[idx][elements.index(element)] = formula[element]

    return formula_array


def formula_to_string(formula):
    """
    :param formula: Chemical formula
    :type formula: dict or Counter

    :returns: A formula string, with element order as given in the list
        IUPAC_element_order.
        If one or more keys in the dictionary are not one of the elements
        in the periodic table, then they are added at the end of the string.
    :rtype: str
    """

    formula_string = ""
    for e in IUPAC_element_order:
        if e in formula and np.abs(formula[e]) > 1.0e-12:
            if np.abs(formula[e] - 1.0) < 1.0e-12:
                formula_string += e
            else:
                formula_string += e + str(nsimplify(formula[e]))

    for e in formula:
        if e not in IUPAC_element_order:
            if e in formula and np.abs(formula[e]) > 1.0e-12:
                if np.abs(formula[e] - 1.0) < 1.0e-12:
                    formula_string += e
                else:
                    formula_string += e + str(nsimplify(formula[e]))

    return formula_string


def sort_element_list_to_IUPAC_order(element_list):
    """
    :param element_list : List of elements.
    :type element_list: list

    :returns: List of elements sorted into IUPAC order
    :rtype: list
    """
    sorted_list = [e for e in IUPAC_element_order if e in element_list]
    assert len(sorted_list) == len(element_list)
    return sorted_list


def convert_fractions(composite, phase_fractions, input_type, output_type):
    """
    Takes a composite with a set of user defined molar, volume
    or mass fractions (which do not have to be the fractions
    currently associated with the composite) and
    converts the fractions to molar, mass or volume.

    Conversions to and from mass require a molar mass to be
    defined for all phases. Conversions to and from volume
    require set_state to have been called for the composite.

    :param composite: Composite for which fractions are to be defined.
    :type composite: :class:`~burnman.Composite`

    :param phase_fractions: List of input phase fractions
        (of type input_type).
    :type phase_fractions: list of floats

    :param input_type: Input fraction type. One of 'molar', 'mass' or 'volume'.
    :type input_type: str

    :param output_type: Output fraction type. One of 'molar', 'mass' or 'volume'.
    :type output_type: str

    :returns: List of output phase fractions (of type output_type)
    :rtype: list of floats
    """
    if input_type == "volume" or output_type == "volume":
        if composite.temperature is None:
            raise Exception(
                composite.to_string()
                + ".set_state(P, T) has not been called, so volume fractions are currently undefined. Exiting."
            )

    if input_type == "molar":
        molar_fractions = phase_fractions
    if input_type == "volume":
        total_moles = sum(
            volume_fraction / phase.molar_volume
            for volume_fraction, phase in zip(phase_fractions, composite.phases)
        )
        molar_fractions = [
            volume_fraction / (phase.molar_volume * total_moles)
            for volume_fraction, phase in zip(phase_fractions, composite.phases)
        ]
    if input_type == "mass":
        total_moles = sum(
            mass_fraction / phase.molar_mass
            for mass_fraction, phase in zip(phase_fractions, composite.phases)
        )
        molar_fractions = [
            mass_fraction / (phase.molar_mass * total_moles)
            for mass_fraction, phase in zip(phase_fractions, composite.phases)
        ]

    if output_type == "volume":
        total_volume = sum(
            molar_fraction * phase.molar_volume
            for molar_fraction, phase in zip(molar_fractions, composite.phases)
        )
        output_fractions = [
            molar_fraction * phase.molar_volume / total_volume
            for molar_fraction, phase in zip(molar_fractions, composite.phases)
        ]
    elif output_type == "mass":
        total_mass = sum(
            molar_fraction * phase.molar_mass
            for molar_fraction, phase in zip(molar_fractions, composite.phases)
        )
        output_fractions = [
            molar_fraction * phase.molar_mass / total_mass
            for molar_fraction, phase in zip(molar_fractions, composite.phases)
        ]
    elif output_type == "molar":
        output_fractions = molar_fractions

    return output_fractions


def reaction_matrix_as_strings(reaction_matrix, compound_names):
    """
    Returns a list of string representations of all the reactions in
    reaction_matrix.

    :param reaction_matrix: Matrix of stoichiometric amounts
        of each compound j in reaction i.
    :type reaction_matrix: 2D numpy array

    :param compound_names: List of compound names.
    :type compound_names: list of strings

    :returns: List of strings corresponding to each reaction.
    :rtype: list of strings
    """
    reaction_strings = []
    for reaction in reaction_matrix:
        lhs, rhs = ("", "")
        for i, coefficient in enumerate(reaction):
            if coefficient < -1.0e-10:
                if len(lhs) > 0:
                    lhs += " + "
                lhs += f"{-coefficient} {compound_names[i]}"
            if coefficient > 1.0e-10:
                if len(rhs) > 0:
                    rhs += " + "
                rhs += f"{coefficient} {compound_names[i]}"
        reaction_strings.append(f"{lhs} = {rhs}")

    return reaction_strings
