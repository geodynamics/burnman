from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2018 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from collections import Counter, OrderedDict
from scipy.optimize import nnls
    
from .processchemistry import dictionarize_formula, formula_mass

class OrderedCounter(Counter, OrderedDict):
    """
    Counter that remembers the order elements are first encountered
    """

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)
    
def composition_property(func):
    """
    Decorator @composition_property to be used for cached properties of compositions.

    To be used on function in Composition or derived classes that should be exposed
    as read-only properties that are cached. 
    """
    class mat_obj():

        def __init__(self, func):
            self.func = func
            self.varname = self.func.__name__

        def get(self, obj):
            if not hasattr(obj, "_cached"):
                raise Exception("The composition_property decorator could not find class member _cached. "
                                "Did you forget to call Composition.__init__(self) in __init___?")
            cache_array = getattr(obj, "_cached")
            if self.varname not in cache_array:
                cache_array[self.varname] = self.func(obj)
            return cache_array[self.varname]

    return property(mat_obj(func).get, doc=func.__doc__)


def file_to_composition_list(fname, unit_type, normalize):
    """
    Takes an input file with a specific format and returns a list of compositions
    (and associated comments) contained in that file.

    Parameters
    ----------
    fname : string
        Path to ascii file containing composition data.
        Lines beginning with a hash are not read.
        The first read-line of the datafile contains a list of tab or
        space-separated components (e.g. FeO or SiO2), followed by the word Comment.
        Following lines are lists of floats with the amounts of each component.
        After the component amounts, the user can write anything they like
        in the Comment section.
    unit_type : 'weight' or 'molar'
        Specify whether the compositions in the file are given as weight or
        molar amounts.
    normalize : boolean
        If False, absolute numbers of moles/grams of component are stored,
        otherwise the component amounts of returned compositions will
        sum to one (until Composition.renormalize() is used).

    """
    lines = list(filter(None, [line.rstrip('\n').split() for line in open(fname) if line[0] != '#']))
    n_components = lines[0].index("Comment")
    components = lines[0][:n_components]
    comments = [line[n_components:] for line in lines[1:]]
    compositions = np.array([map(float, l) for l in list(zip(*(list(zip(*lines[1:]))[:n_components])))])
    return [Composition(OrderedCounter(dict(zip(components, c))), unit_type, normalize) for c in compositions], comments

    
class Composition(object):
    """
    Class for a composition object, which can be used
    to store, modify and renormalize compositions, 
    and also convert between molar, weight, 
    and atomic amounts.

    This class is available as ``burnman.Composition``.
    """

    def __init__(self, composition_dictionary, unit_type='weight', normalize=False):
        """
        Create a composition using a dictionary and unit type.

        Parameters
        ----------
        composition_dictionary : dictionary
            Dictionary of components (given as a string) and their amounts.
        unit_type : 'weight' or 'molar' (optional, 'weight' as standard)
            Specify whether the input composition is given as weight or
            molar amounts.
        normalize : boolean
            If False, absolute numbers of moles/grams of component are stored,
            otherwise the component amounts of returned compositions will 
            sum to one (until Composition.renormalize() is used).
        """

        self._cached = {}

        n_total = float(sum(composition_dictionary.values()))

        normalized_dictionary = {}
        for k in composition_dictionary.keys():
            normalized_dictionary[k] = composition_dictionary[k]/n_total
        
        self.normalization_component = {'weight': 'total',
                                        'molar': 'total',
                                        'atomic': 'total'}

        self.normalization_amount = {'weight': 1.,
                                     'molar': 1.,
                                     'atomic': 1.}


        # component formulae
        self.component_formulae = {c: dictionarize_formula(c)
                                   for c in composition_dictionary.keys()}
        
        # elemental compositions of components
        self.element_list = OrderedCounter()
        for component in self.component_formulae.values():
            self.element_list += OrderedCounter({element: n_atoms
                                          for (element, n_atoms) in component.items()})
        self.element_list = list(self.element_list.keys())
        
        if unit_type == 'weight':
            if normalize:
                self._cached['weight_composition'] = OrderedCounter(normalized_dictionary)
            else:
                self._cached['weight_composition'] = OrderedCounter(composition_dictionary)
                
                mole_total = sum([composition_dictionary[c] /
                                  formula_mass(self.component_formulae[c])
                                  for c in composition_dictionary.keys()])
                self.normalization_amount['weight'] = n_total
                self.normalization_amount['molar'] = mole_total
                self.normalization_amount['atomic'] = sum(self._moles_component_to_atoms(self.molar_composition).values())
                
        elif unit_type == 'molar':
            if normalize:
                self._cached['molar_composition'] = OrderedCounter(normalized_dictionary)
            else:
                self._cached['molar_composition'] = OrderedCounter(composition_dictionary)
                
                weight_total = sum([composition_dictionary[c] *
                                    formula_mass(self.component_formulae[c])
                                    for c in composition_dictionary.keys()])
                
                self.normalization_amount['weight'] = weight_total
                self.normalization_amount['molar'] = n_total
                self.normalization_amount['atomic'] = sum(self._moles_component_to_atoms(self.molar_composition).values())
            
        else:
            raise Exception('Unit type not yet implemented. '
                            'Should be either weight or molar.')

        
        
            

    def renormalize(self, unit_type, normalization_component, normalization_amount):
        """
        Change the normalization for a given unit type 
        (weight, molar, or atomic)
        Resets cached composition only for that unit type

        Parameters
        ----------
        unit_type : 'weight', 'molar' or 'atomic'
            Unit type composition to be renormalised
        normalization_component: string 
            Component/element on which to renormalize.
            String must either be one of the components/elements
            already in composite, or have the value 'total'
        normalization_amount: float
            Amount of component in the renormalised composition 
        """
        if unit_type == 'weight':
            s = self.weight_composition
        elif unit_type == 'molar':
            s = self.molar_composition
        elif unit_type == 'atomic':
            s = self.atomic_composition
        else:
            raise Exception('Unit type not recognised. '
                            'Should be one of weight, molar and atomic')

        self.normalization_component[unit_type] = normalization_component
        self.normalization_amount[unit_type] = normalization_amount
        self._cached[unit_type+'_composition'] = self._normalize_to_basis(s, unit_type)
            
    def add_components(self, composition_dictionary, unit_type):
        """
        Add (or remove) components from the composition.
        The components are added to the current state of the 
        (weight or molar) composition; if the composition has 
        been renormalised, then this should be taken into account.

        Parameters
        ----------
        composition_dictionary : dictionary
            Components to add, and their amounts, in dictionary form
        unit_type : 'weight' or 'molar'
            Unit type of the components to be added
        """
        if unit_type == 'weight':
            composition = self.weight_composition
        elif unit_type == 'molar':
            composition = self.molar_composition
        else:
            raise Exception('Unit type not recognised. '
                            'Should be either weight or molar.')

        composition += OrderedCounter(composition_dictionary)

        self.__init__(composition, unit_type)

    def change_component_set(self, new_component_list):
        """
        Change the set of basis components without 
        changing the bulk composition. 

        Will raise an exception if the new component set is 
        invalid for the given composition.

        Parameters
        ----------
        new_component_list : list of strings
            New set of basis components.
        """
        composition = np.array([self.atomic_composition[element]
                                for element in self.element_list])
        component_matrix = np.zeros((len(new_component_list), len(self.element_list)))
        
        for i, component in enumerate(new_component_list):
            formula = dictionarize_formula(component)
            for element, n_atoms in formula.items():
                component_matrix[i][self.element_list.index(element)] = n_atoms

        sol = nnls(component_matrix.T, composition)
        if sol[1] < 1.e-12:
            component_amounts = sol[0]
        else:
            raise Exception('Failed to change component set. '
                            'Could not find a non-negative least squares solution. '
                            'Can the bulk composition be described with this set of components?')

        composition = OrderedCounter(dict(zip(new_component_list, component_amounts)))
        self.__init__(composition, 'molar')
    
    def _normalize_to_basis(self, composition, unit_type):
        if self.normalization_component[unit_type] == 'total':
            n_orig = float(sum(composition.values()))
        else:
            n_orig = composition[self.normalization_component[unit_type]]
            
        for k in composition.keys():
            composition[k] *= self.normalization_amount[unit_type]/n_orig
    
        return composition

    @composition_property
    def molar_composition(self):
        """
        Returns the molar composition as a counter [moles]
        """
        mole_compositions = OrderedCounter({c: self.weight_composition[c] /
                                            formula_mass(self.component_formulae[c])
                                            for c in self.weight_composition.keys()})

        return self._normalize_to_basis(mole_compositions, 'molar')

    @composition_property
    def weight_composition(self):
        """
        Returns the weight composition as a counter [g]
        """
        weight_compositions = OrderedCounter({c: self.molar_composition[c] *
                                              formula_mass(self.component_formulae[c])
                                              for c in self.molar_composition.keys()})

        return self._normalize_to_basis(weight_compositions, 'weight')
    
    @composition_property
    def atomic_composition(self):
        """
        Returns the atomic composition as a counter [moles]
        """
        atom_compositions = self._moles_component_to_atoms(self.molar_composition)
        
        return self._normalize_to_basis(atom_compositions, 'atomic')

    def _moles_component_to_atoms(self, molar_composition_dictionary):
        component_matrix = np.zeros((len(self.component_formulae), len(self.element_list)))
        molar_composition_vector = np.zeros(len(self.component_formulae))
        for i, (component, formula) in enumerate(self.component_formulae.items()):
            molar_composition_vector[i] = molar_composition_dictionary[component]
            
            for element, n_atoms in formula.items():
                component_matrix[i][self.element_list.index(element)] = n_atoms
        
        atom_compositions = np.dot(molar_composition_vector, component_matrix)
        return OrderedCounter(dict(zip(self.element_list, atom_compositions)))
        

    def print(self, unit_type, significant_figures=1,
              normalization_component='total', normalization_amount=100.):
        """
        Pretty-print function for the composition
        This does not renormalize the Composition internally

        Parameters
        ----------
        unit_type : 'weight', 'molar' or 'atomic'
            Unit type in which to print the composition
        significant_figures : integer
            Number of significant figures for each amount
        normalization_component: string 
            Component/element on which to renormalize.
            String must either be one of the components/elements
            already in composite, or have the value 'total'.
            (default = 'total')
        normalization_amount: float
            Amount of component in the renormalised composition.
            (default = '100.')
        """
        if unit_type == 'weight':
            print('Weight composition')
            
            if normalization_component == 'total':
                total_stored = float(sum(self.weight_composition.values()))
            else:
                total_stored = self.weight_composition[normalization_component]
            f = normalization_amount/total_stored
            
            for (key, value) in sorted(self.weight_composition.items()):
                print('{0}: {1:0.{sf}f}'.format(key, value*f, sf=significant_figures))
        elif unit_type == 'molar':
            print('Molar composition')

            if normalization_component == 'total':
                total_stored = float(sum(self.molar_composition.values()))
            else:
                total_stored = self.molar_composition[normalization_component]
            f = normalization_amount/total_stored
            
            for (key, value) in sorted(self.molar_composition.items()):
                print('{0}: {1:0.{sf}f}'.format(key, value*f, sf=significant_figures))
        elif unit_type == 'atomic':
            print('Atomic composition')

            if normalization_component == 'total':
                total_stored = float(sum(self.atomic_composition.values()))
            else:
                total_stored = self.atomic_composition[normalization_component]
            f = normalization_amount/total_stored
            
            for (key, value) in sorted(self.atomic_composition.items()):
                print('{0}: {1:0.{sf}f}'.format(key, value*f, sf=significant_figures))
        else:
            raise Exception('unit_type not yet implemented. Should be either weight,  molar or atomic.')
        
