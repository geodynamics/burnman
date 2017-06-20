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
from burnman.processchemistry import dictionarize_formula

# Let's import a function to read compositions from file, and the Composition class
# The function assumes that the first line is a header of components, and then the word "Comment"
# The following lines are either molar or weight proportions of the components given in the header,
# followed by any string
from burnman.composition import file_to_composition_list, Composition

# Let's read in composition data from an example file
print('1) Converting weight percents to molar and atomic fractions')
compositions, comments = file_to_composition_list('../burnman/data/input_compositions/dhz_mineral_compositions.dat', 'weight')

# The hard work has already been done in that one line: each of the compositions is stored as an instance of
# the Composition class in a list "compositions", and the postceding comments in the list "comments".

# In the following lines we parse the comments (which in this case contain useful information about the number of
# oxygens per formula unit, and the names of the minerals to which each component corresponds
for i, composition in enumerate(compositions):
    # Parsing the comments list
    n_O, name, reference = comments[i]

    # Each composition has the attributes "weight_composition", "molar_composition" and "atomic_composition",
    # which are dictionaries. Here we just format the weight and molar dictionaries for printing.
    formatted_weight_composition = [float('{0:.5f}'.format(c)) for c in composition.weight_composition.values()]
    formatted_molar_composition = [float('{0:.5f}'.format(c)) for c in composition.molar_composition.values()]
    print('{0}:\n {1}\n {2} (wt)\n {3} (molar)\n'.format(name, composition.weight_composition.keys(),
                                            formatted_weight_composition, formatted_molar_composition))

    # By default, the compositions are normalised to 1 g (weight) or 1 mole (of components or atoms)
    # We can change this with the class function "change_basis". Here we change the number of atoms in the
    # atomic composition so that the number of oxygens is equal to n_O.
    # This does not change the molar or weight bases.
    compositions[i].change_basis('atomic', 'O', float(n_O))
    formatted_composition = [float('{0:.5f}'.format(c)) for c in composition.atomic_composition.values()]
    print(' {0}\n {1} \n (atoms, {2} oxygen basis)\n\n'.format(composition.atomic_composition.keys(), formatted_composition, float(n_O)))

    
# Let's do something a little more complicated.
# When we're making a starting mix for petrological experiments, we often have to add additional components.
# For example, sodium is typically added as NaCO3 even if we don't want carbon in our starting mix. This is because
# Na2O is reactive and hygroscopic. Similarly, we add iron as Fe3O4 even if we want a reduced starting mix, because
# FeO is not a stable stoichiometric compound. Let's see how easy it is to make a general structure for creating such mixes.

# We start with a fayalite starting composition
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

# Here we modify the bulk composition by adding oxygen to the starting mix equal to one third the total FeO (on a molar basis)
# This is equivalent to adding FeO as Fe3O4 (1/3 Fe3O4 = FeO + 1/3 O)
print('')
print('FeO doesn\'t exist as a stoichiometric compound, so we commonly create \na starting mix from magnetite powder and then reduce the mix.')
print('Here, we add one third of an oxygen for every FeO in the required bulk composition.')
composition.add_components({'O1/3': composition.molar_composition['FeO']},
                           'molar')

# Now we can change the component set again, this time into the set of compounds that we'll use
# to make the starting mix
composition.change_component_set(['Fe3O4', 'SiO2'])

# Let's print out the new atomic composition
composition.change_basis('atomic', 'Si', 1.)
print('Atomic composition: {0} {1}\n'.format(composition.atomic_composition.keys(), composition.atomic_composition.values()))

# Finally, let's print out the starting composition that we'll use,
# assuming that we want to start with 2 g of Fe3O4 and SiO2 powder
print('We can now print the modified starting mix:')
composition.change_basis('weight', 'total', 2.)
print('Weight composition (1g SiO2): {0} {1}'.format(composition.weight_composition.keys(), composition.weight_composition.values()))
