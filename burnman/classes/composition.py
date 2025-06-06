# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2018 by the BurnMan team, released under the GNU
# GPL v2 or later.

from copy import deepcopy
import numpy as np
from scipy.optimize import nnls

from ..utils.chemistry import dictionarize_formula, formula_mass
from ..utils.misc import OrderedCounter


def file_to_composition_list(fname, unit_type, normalize):
    """
    Takes an input file with a specific format and returns a list of
    compositions (and associated comments) contained in that file.

    :param fname: Path to ascii file containing composition data.
        Lines beginning with a hash are not read.
        The first read-line of the datafile contains a list of tab or
        space-separated components (e.g. FeO or SiO2), followed by the
        word Comment.
        Following lines are lists of floats with the amounts of each component.
        After the component amounts, the user can write anything they like
        in the Comment section.
    :type fname: str
    :param unit_type: 'mass', 'weight' or 'molar'
        Specify whether the compositions in the file are given as
        mass (weight) or molar amounts.
    :type unit_type: str
    :param normalize: If False, absolute numbers of moles/grams of component are stored,
        otherwise the component amounts of returned compositions will
        sum to one (until Composition.renormalize() is used).
    :type normalize : bool

    """
    lines = list(
        filter(
            None, [line.rstrip("\n").split() for line in open(fname) if line[0] != "#"]
        )
    )
    n_components = lines[0].index("Comment")
    components = lines[0][:n_components]
    comments = [line[n_components:] for line in lines[1:]]
    compositions = np.array(
        [map(float, ln) for ln in list(zip(*(list(zip(*lines[1:]))[:n_components])))]
    )
    return [
        Composition(OrderedCounter(dict(zip(components, c))), unit_type, normalize)
        for c in compositions
    ], comments


class Composition(object):
    """
    Class for a composition object, which can be used
    to store, modify and renormalize compositions,
    and also convert between mass, molar
    and atomic amounts. Weight is provided as an alias
    for mass, as we assume that only Earthlings
    will use this software.

    This class is available as ``burnman.Composition``.
    """

    def __init__(self, composition_dictionary, unit_type="mass", normalize=False):
        """
        Create a composition using a dictionary and unit type.

        :param composition_dictionary: Dictionary of components
            (given as a string) and their amounts.
        :type composition_dictionary: dictionary
        :param unit_type: 'mass', 'weight' or 'molar' (optional, 'mass' as default)
            Specify whether the input composition is given as mass or
            molar amounts.
        :type unit_type: str
        :param normalize: If False, absolute numbers of kilograms/moles of component are
            stored, otherwise the component amounts of returned compositions
            will sum to one (until Composition.renormalize() is used).
        :type normalize: bool
        """

        self._cached = {}

        n_total = float(sum(composition_dictionary.values()))

        # Create the input dictionary, normalize if requested
        input_dictionary = OrderedCounter(deepcopy(composition_dictionary))
        if normalize:
            for k in composition_dictionary.keys():
                input_dictionary[k] = composition_dictionary[k] / n_total

        # Break component formulae into atomic dictionaries
        fl = self._component_formulae_and_element_lists(composition_dictionary.keys())
        self.component_formulae, self.element_list = fl

        if unit_type == "mass" or unit_type == "weight":
            self.mass_composition = input_dictionary
        elif unit_type == "molar":
            self.mass_composition = self._mole_to_mass_composition(input_dictionary)
        else:
            raise Exception(
                "Unit type not yet implemented. "
                "Should be either mass, weight or molar."
            )

    def _component_formulae_and_element_lists(self, component_strings):
        """
        Converts a list of components in string form into a list of components in dictionary form
        and a list of elements.

        :param component_strings: A list of strings containing a formula for each component
        :type component_strings: list of strings
        :return: a list of components in dictionary form and a list of elements.
        :rtype: _type_
        """

        # Break component formulae into atomic dictionaries
        component_formulae = {c: dictionarize_formula(c) for c in component_strings}

        # Create lists of elemental compositions of components
        element_list = OrderedCounter()
        for component in component_formulae.values():
            element_list += OrderedCounter(
                {element: n_atoms for (element, n_atoms) in component.items()}
            )
        element_list = list(element_list.keys())

        return component_formulae, element_list

    def renormalize(self, unit_type, normalization_component, normalization_amount):
        """
        Change the total amount of material in the composition
        to satisfy a given normalization condition
        (mass, weight, molar, or atomic)

        :param unit_type: 'mass', 'weight', 'molar' or 'atomic'
            Unit type with which to normalize the composition
        :type unit_type: str
        :param normalization_component: Component/element on which to renormalize.
            String must either be one of the components/elements
            already in the composition, or have the value 'total'.
        :type normalization_component: str
        :param normalization_amount: Amount of component in the
            renormalised composition.
        :type normalization_amount: float
        """

        if unit_type not in ["mass", "weight", "molar", "atomic"]:
            raise Exception(
                "unit_type not yet implemented."
                "Should be either mass, weight, molar or atomic."
            )

        c = self.composition(unit_type)

        if normalization_component == "total":
            f = normalization_amount / float(sum(c.values()))
        else:
            f = normalization_amount / c[normalization_component]

        new_mass_composition = OrderedCounter()
        for k in self.mass_composition.keys():
            new_mass_composition[k] = self.mass_composition[k] * f

        self.mass_composition = new_mass_composition

    def add_components(self, composition_dictionary, unit_type):
        """
        Add (or remove) components from the composition.
        The components are added to the current state of the
        (mass, weight or molar) composition; if the composition has
        been renormalised, then this should be taken into account.

        :param composition_dictionary: Components to add, and their amounts.
        :type composition_dictionary: dictionary
        :param unit_type: 'mass', 'weight' or 'molar'.
            Unit type of the components to be added.
        :type unit_type: str
        """
        if unit_type == "mass" or unit_type == "weight":
            composition = self.mass_composition
        elif unit_type == "molar":
            composition = self.molar_composition
        else:
            raise Exception(
                "Unit type not recognised. " "Should be either mass, weight or molar."
            )

        composition += OrderedCounter(composition_dictionary)

        # Reinitialize composition object
        self.__init__(composition, unit_type)

    def change_component_set(self, new_component_list):
        """
        Change the set of basis components without
        changing the bulk composition.

        Will raise an exception if the new component set is
        invalid for the given composition.

        :param new_component_list: New set of basis components.
        :type new_component_list: list of strings
        """

        old_element_list = self.element_list
        _, new_element_list = self._component_formulae_and_element_lists(
            new_component_list
        )

        self.element_list = list(set(old_element_list).union(set(new_element_list)))

        composition = np.array(
            [self.atomic_composition[element] for element in self.element_list]
        )
        component_matrix = np.zeros((len(new_component_list), len(self.element_list)))

        for i, component in enumerate(new_component_list):
            formula = dictionarize_formula(component)
            for element, n_atoms in formula.items():
                component_matrix[i][self.element_list.index(element)] = n_atoms

        sol = nnls(component_matrix.T, composition)
        if sol[1] < 1.0e-10:
            component_amounts = sol[0]
        else:
            raise Exception(
                "Failed to change component set. "
                "Could not find a non-negative "
                f"least squares solution (residual={sol[1]}). "
                "Can the bulk composition be described "
                "with this set of components?"
            )

        composition = OrderedCounter(dict(zip(new_component_list, component_amounts)))

        # Reinitialize the object
        self.__init__(composition, "molar")

    def remove_null_components(self, tol=1.0e-12):
        """
        Removes components with absolute concentrations less than a given tolerance.

        :param tol: zero tolerance, defaults to 1.e-12
        :type tol: float, optional
        """
        c = deepcopy(self.molar_composition)
        del_keys = []
        for key, value in c.items():
            if np.abs(value) < tol:
                del_keys.append(key)
        for key in del_keys:
            c.pop(key)

        # Reinitialize the object
        self.__init__(c, "molar")

    def _mole_to_mass_composition(self, molar_comp):
        """
        Hidden function to returns the mass composition as a counter [kg]
        """
        cf = self.component_formulae
        mass_composition = OrderedCounter(
            {c: molar_comp[c] * formula_mass(cf[c]) for c in molar_comp.keys()}
        )

        return mass_composition

    @property
    def weight_composition(self):
        """
        An alias for mass composition [kg].
        """
        return self.mass_composition

    @property
    def molar_composition(self):
        """
        Returns the molar composition as a counter [moles]
        """
        mass_comp = self.mass_composition
        cf = self.component_formulae

        return OrderedCounter(
            {c: mass_comp[c] / formula_mass(cf[c]) for c in mass_comp.keys()}
        )

    @property
    def atomic_composition(self):
        """
        Returns the atomic composition as a counter [moles]
        """

        return self._moles_to_atoms(self.molar_composition)

    def composition(self, unit_type):
        """
        Helper function to return the composition in the
        desired type.

        :param unit_type: One of 'mass', 'weight', 'molar' and 'atomic'.
        :type unit_type: str

        :returns: Mass (weight), molar or atomic composition.
        :rtype: OrderedCounter
        """
        return getattr(self, f"{unit_type}_composition")

    def _moles_to_atoms(self, molar_comp_dictionary):
        """
        Hidden function that converts a molar component
        dictionary into an atomic (elemental) dictionary
        """
        component_matrix = np.zeros(
            (len(self.component_formulae), len(self.element_list))
        )
        cf = self.component_formulae
        molar_composition_vector = np.zeros(len(cf))
        for i, (component, formula) in enumerate(cf.items()):
            molar_composition_vector[i] = molar_comp_dictionary[component]

            for element, n_atoms in formula.items():
                component_matrix[i][self.element_list.index(element)] = n_atoms

        atom_compositions = np.dot(molar_composition_vector, component_matrix)
        return OrderedCounter(dict(zip(self.element_list, atom_compositions)))

    def print(
        self,
        unit_type,
        significant_figures=1,
        normalization_component="total",
        normalization_amount=None,
    ):
        """
        Pretty-print function for the composition
        This does not renormalize the Composition object itself,
        only the printed values.

        :param unit_type: 'mass', 'weight', 'molar' or 'atomic'
            Unit type in which to print the composition.
        :type unit_type: str
        :param significant_figures: Number of significant figures
            for each amount.
        :type significant_figures: int
        :param normalization_component: Component/element on which to renormalize.
            String must either be one of the components/elements
            already in composite, or have the value 'total'.
            (default = 'total')
        :type normalization_component: str
        :param normalization_amount: Amount of component in the
            renormalised composition. If not explicitly set,
            no renormalization will be applied.
        :type normalization_amount: float
        """

        if unit_type not in ["mass", "weight", "molar", "atomic"]:
            raise Exception(
                "unit_type not yet implemented."
                "Should be either mass, weight, molar or atomic."
            )

        c = self.composition(unit_type)
        print(f"{unit_type.capitalize()} composition")

        if normalization_amount is None:
            f = 1
        elif normalization_component == "total":
            f = normalization_amount / float(sum(c.values()))
        else:
            f = normalization_amount / c[normalization_component]

        for key, value in sorted(c.items()):
            print(f"{key}: {value*f:0.{significant_figures}f}")
