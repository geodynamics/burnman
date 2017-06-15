from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


import burnman
from burnman.composition import file_to_composition_list, Composition
from burnman.processchemistry import dictionarize_formula

print('1) Converting weight percents to molar and atomic fractions')
compositions, comments = file_to_composition_list('../burnman/data/input_compositions/dhz_mineral_compositions.dat', 'weight')
for i, composition in enumerate(compositions):
    n_O, name, reference = comments[i]
    formatted_weight_composition = [float('{0:.5f}'.format(c)) for c in composition.weight_composition.values()]
    formatted_molar_composition = [float('{0:.5f}'.format(c)) for c in composition.molar_composition.values()]
    print('{0}:\n {1}\n {2} (wt)\n {3} (molar)\n'.format(name, composition.weight_composition.keys(),
                                            formatted_weight_composition, formatted_molar_composition))
    
    compositions[i].change_basis('atomic', 'O', float(n_O))
    formatted_composition = [float('{0:.5f}'.format(c)) for c in composition.atomic_composition.values()]
    print(' {0}\n {1} \n (atoms, {2} oxygen basis)\n\n'.format(composition.atomic_composition.keys(), formatted_composition, float(n_O)))


print('2) Fayalite starting mix calculations:')
composition = Composition(dictionarize_formula('Fe2SiO4'), 'molar')

# The first step is to split the desired starting mix into a set of starting oxides
# (alternatively, we could have initialised the composition with a dictionary of these oxides)
composition.change_component_set(['FeO', 'SiO2'])

# Let's check the atomic, molar and weight composition of this compound
composition.change_basis('atomic', 'Si', 1.)
print('Atomic composition: {0} {1}'.format(composition.atomic_composition.keys(), composition.atomic_composition.values()))
composition.change_basis('molar', 'SiO2', 1.)
print('Molar composition (1 mole SiO2): {0} {1}'.format(composition.molar_composition.keys(), composition.molar_composition.values()))
composition.change_basis('weight', 'SiO2', 1.)
print('Weight composition (1g SiO2): {0} {1}'.format(composition.weight_composition.keys(), composition.weight_composition.values()))

print('')
print('FeO doesn\'t exist as a stoichiometric compound, so we commonly create \na starting mix from magnetite powder and then reduce the mix.')
print('Here, we add one third of an oxygen for every FeO in the required bulk composition.')
composition.add_components({'O1/3': composition.molar_composition['FeO']},
                           'molar')

composition.change_component_set(['Fe3O4', 'SiO2'])
composition.change_basis('atomic', 'Si', 1.)
print('Atomic composition: {0} {1}\n'.format(composition.atomic_composition.keys(), composition.atomic_composition.values()))
composition.change_basis('weight', 'SiO2', 1.)

print('We can now print the modified starting mix:')
print('Weight composition (1g SiO2): {0} {1}'.format(composition.weight_composition.keys(), composition.weight_composition.values()))
