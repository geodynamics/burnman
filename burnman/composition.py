from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from collections import Counter
from scipy.optimize import nnls

from .processchemistry import dictionarize_formula, formula_mass

def composition_property(func):
    """
    Decorator @material_property to be used for cached properties of materials.

    To be used on function in Material or derived classes that should be exposed
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


def file_to_composition_list(fname, unit_type):
    lines = filter(None, [line.rstrip('\n').split() for line in open(fname)])
    n_components = lines[0].index("Comment")
    components = lines[0][:n_components]
    comments = [line[n_components:] for line in lines[1:]]
    compositions = np.array([map(float, l) for l in zip(*(zip(*lines[1:])[:n_components]))])
    return [Composition(Counter(dict(zip(components, c))), unit_type) for c in compositions], comments

    
class Composition(object):
    def __init__(self, composition_dictionary, unit_type='weight'):
        self._cached = {}

        n_total = np.sum(composition_dictionary.values())
        for k in composition_dictionary.keys():
            composition_dictionary[k] /= n_total
        
        if unit_type == 'weight':
            self._cached['weight_composition'] = Counter(composition_dictionary)
            
        elif unit_type == 'molar':
            self._cached['molar_composition'] = Counter(composition_dictionary)
            
        else:
            raise Exception('unit_type not yet implemented. Should be either weight or molar.')


        self.component_formulae = {c: dictionarize_formula(c) for c in composition_dictionary.keys()}
        
        # atom compositions
        self.element_list = Counter()
        for component in self.component_formulae.values():
            self.element_list += Counter({element: n_atoms for (element, n_atoms) in component.items()})
        self.element_list = list(self.element_list.keys())
            
        self.basis = {'weight': 'total', 'molar': 'total', 'atomic': 'total'}
        self.n_basis = {'weight': 1., 'molar': 1., 'atomic': 1.}
        

    def change_basis(self, unit_type, basis, n_basis):

        if unit_type == 'weight':
            s = self.weight_composition
        elif unit_type == 'molar':
            s = self.molar_composition
        elif unit_type == 'atomic':
            s = self.atomic_composition
        else:
            raise Exception('unit type not recognised. should be one of weight, molar and atomic')

        self.basis[unit_type] = basis
        self.n_basis[unit_type] = n_basis
        self._cached[unit_type+'_composition'] = self._normalize_to_basis(s, unit_type)
            
    def add_components(self, composition_dictionary, unit_type):
        if unit_type == 'weight':
            composition = self.weight_composition
        elif unit_type == 'molar':
            composition = self.molar_composition
        else:
            raise Exception('unit type not recognised. should be either weight or molar.')

        composition += Counter(composition_dictionary)

        self.__init__(composition, unit_type)

    def change_component_set(self, new_component_list):
        composition = np.array([self.atomic_composition[element] for element in self.element_list])
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

        composition = Counter(dict(zip(new_component_list, component_amounts)))
        self.__init__(composition, 'molar')
    
    def _normalize_to_basis(self, composition, unit_type):
        if self.basis[unit_type] == 'total':
            n_orig = np.sum(composition.values())
        else:
            n_orig = composition[self.basis[unit_type]]
            
        for k in composition.keys():
            composition[k] *= self.n_basis[unit_type]/n_orig
    
        return composition

    @composition_property
    def molar_composition(self):
        mole_compositions = Counter({c: self.weight_composition[c] /
                                     formula_mass(self.component_formulae[c])
                                     for c in self.weight_composition.keys()})

        return self._normalize_to_basis(mole_compositions, 'molar')

    @composition_property
    def weight_composition(self):
        weight_compositions = Counter({c: self.molar_composition[c] *
                                     formula_mass(self.component_formulae[c])
                                     for c in self.molar_composition.keys()})

        return self._normalize_to_basis(weight_compositions, 'weight')
    
    @composition_property
    def atomic_composition(self):

        component_matrix = np.zeros((len(self.component_formulae), len(self.element_list)))
        molar_composition_vector = np.zeros(len(self.component_formulae))
        for i, (component, formula) in enumerate(self.component_formulae.items()):
            molar_composition_vector[i] = self.molar_composition[component]
            
            for element, n_atoms in formula.items():
                component_matrix[i][self.element_list.index(element)] = n_atoms
        
        atom_compositions = np.dot(molar_composition_vector, component_matrix)
        atom_compositions = (atom_compositions.T / np.sum(atom_compositions)).T
        atom_compositions = Counter(dict(zip(self.element_list, atom_compositions)))

        return self._normalize_to_basis(atom_compositions, 'atomic')
        

    def print(self, unit_type, significant_figures=1, total=100.):
        if unit_type == 'weight':
            print('\nWeight composition')
            for (key, value) in self.weight_composition.items():
                print('{0}: {1:0.{sf}f}'.format(key, value*total, sf=significant_figures))
        elif unit_type == 'molar':
            print('\nMolar composition')
            for (key, value) in self.molar_composition.items():
                print('{0}: {1:0.{sf}f}'.format(key, value*total, sf=significant_figures))
        elif unit_type == 'atomic':
            print('\nAtomic composition')
            for (key, value) in self.atomic_composition.items():
                print('{0}: {1:0.{sf}f}'.format(key, value*total, sf=significant_figures))
        else:
            raise Exception('unit_type not yet implemented. Should be either weight,  molar or atomic.')
        
