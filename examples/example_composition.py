# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_composition
-------------------

This example script demonstrates the use of BurnMan's Composition class.

*Uses:*

* :class:`burnman.Composition`

*Demonstrates:*

* Creating an instance of the Composition class with a molar
  or weight composition
* Printing weight, molar, atomic compositions
* Renormalizing compositions
* Modifying the independent set of components
* Modifying compositions by adding and removing components
"""
from __future__ import absolute_import
from __future__ import print_function

from burnman.tools.chemistry import dictionarize_formula

# Import the Composition class and a function to read compositions from file
# The function assumes that the first line is a header of components,
# followed by the word "Comment"
# All lines after the first are either molar or weight proportions
# of the components given in the header, followed by the comment string
from burnman import Composition
from burnman.classes.composition import file_to_composition_list


if __name__ == "__main__":
    print('1) Creating a Composite instance, printing '
          'molar, weight and weight compositions:\n')

    print('Olivine (fo90):')
    forsterite_composition = Composition({'MgO': 1.8,
                                          'FeO': 0.2,
                                          'SiO2': 1.}, 'molar')

    # The composition class has a pretty-print method called print.
    # It can output the currently stored composition in
    # absolute moles, mass (weight) or atoms:
    forsterite_composition.print('molar', significant_figures=4)

    # Alternatively, the amount can be normalized to a desired total
    # or amount of a given component / element.
    # The print method does not renormalize the composition itself,
    # only the printed values:
    forsterite_composition.print('weight', significant_figures=4,
                                 normalization_component='total',
                                 normalization_amount=1.)
    forsterite_composition.print('atomic', significant_figures=4,
                                 normalization_component='total',
                                 normalization_amount=7.)

    # Let's read in composition data from an example file
    print('\n\n2) Reading compositions from file and renormalizing them:\n')
    d = file_to_composition_list('../burnman/data/input_compositions/'
                                 'dhz_mineral_compositions.dat',
                                 unit_type='weight',
                                 normalize=True)
    compositions, comments = d

    # The hard work has already been done in that one line:
    # each of the compositions is stored as an instance of
    # the Composition class in a list "compositions", and the
    # postceding comments in the list "comments".

    # Now, let's print out the information in a useful form
    for i, composition in enumerate(compositions):
        # Parsing the comments list
        # (which in this case contain useful information about the number of
        # oxygens per formula unit, and the names of the minerals
        # to which each component corresponds)
        n_O, name, reference = comments[i]

        # Each composition has the dictionary attributes
        # "weight_composition", "molar_composition" and "atomic_composition".
        # Here we just format the weight and molar dictionaries for printing.
        components = sorted(composition.weight_composition.keys())
        composition.renormalize('weight', 'total', 100.)
        wf = [float('{0:.3f}'.format(composition.weight_composition[c]))
              for c in components]
        composition.renormalize('molar', 'total', 100.)
        mf = [float('{0:.3f}'.format(composition.molar_composition[c]))
              for c in components]
        print('{0}:\n {1}\n {2} (wt %)\n {3} (mol %)\n'.format(name,
                                                               components,
                                                               wf, mf))

        # In the above, we used the renormalize method to renormalize
        # the composition object to 100 kg or 100 moles total.

        # We can also renormalize to a given amount of
        # a component or element. We now renormalize
        # so that the number of oxygens is equal to n_O.
        compositions[i].renormalize('atomic', 'O', float(n_O))
        components = sorted(composition.atomic_composition.keys())
        af = [float('{0:.3f}'.format(composition.atomic_composition[c]))
              for c in components]
        print(' {0}\n {1} (atoms, {2} oxygen basis)\n'.format(components, af,
                                                              float(n_O)))

    # Let's do something a little more complicated.
    # When we're making a starting mix for petrological experiments,
    # we often have to add additional components.
    # For example, we add iron as Fe2O3 even if we want a reduced
    # oxide starting mix, because FeO is not a stable stoichiometric compound.

    # Here we show how to use BurnMan to create such mixes.

    # We start with a fayalite starting composition
    print('\n3) Fayalite starting mix calculations:\n')
    composition = Composition(dictionarize_formula('Fe2SiO4'), 'molar')

    # The composition currently has the elements Fe, Si and O as
    # distinct components. First, we change this component set into a
    # set of starting oxides
    # (alternatively, we could have directly initialized the
    # composition with a dictionary of these oxides)
    composition.change_component_set(['FeO', 'SiO2'])

    # Let's check the molar composition of this composition
    composition.print('molar', significant_figures=4, normalization_amount=1.)

    # Now we modify the bulk composition by adding oxygen to the
    # starting mix equal to one third the total FeO (on a molar basis)
    # This is equivalent to converting all the FeO into Fe2O3
    # (1/2 Fe2O3 = FeO + 1/2 O)
    print('')
    print('FeO doesn\'t exist as a stoichiometric compound, '
          'but we can create\n'
          'a starting mix from hematite powder and then reduce the mix.\n'
          'Here, we add one half of an oxygen for every FeO '
          'in the required bulk composition.\n\n'
          'The modified starting mix')
    composition.add_components({'O1/2': composition.molar_composition['FeO']},
                               'molar')

    # Now we can change the component set again, this time into
    # the set of compounds that we'll use to make the starting mix
    composition.change_component_set(['Fe2O3', 'SiO2'])

    # Let's print out the new atomic composition
    composition.renormalize('atomic', 'total', 8.)
    elements = sorted(list(composition.atomic_composition.keys()))
    v = ['{0:.3f}'.format(composition.atomic_composition[e])
         for e in elements]
    print('Atomic composition\n{0}\n{1}\n'.format(elements, v))

    # Finally, let's print out the starting composition that we'll use,
    # assuming that we want to start with 2 g of Fe2O3 and SiO2 powder
    composition.print('weight', significant_figures=4, normalization_amount=2.)
    composition.renormalize('weight', 'total', 2.)

    # Now let's do something similar, but for a more complicated
    # starting mix involving calcium and sodium.
    # Calcium and sodium are typically added as CaCO3 and Na2CO3
    # even if we don't want carbonate in our starting mix.
    # This is because CaO and Na2O are reactive and hygroscopic
    # (we don't want an unknown amount of water screwing up our weights).

    print('\n\n4) KLB-1 starting mix calculations '
          '(carbonated, Fe57 as metallic):\n')

    KLB1 = Composition({'SiO2': 39.4,
                        'Al2O3': 2.0,
                        'CaO': 3.3,
                        'MgO': 49.5,
                        'FeO': 5.2,
                        'Na2O': 0.26}, 'molar')  # from Holland et al., 2013

    # Now we need to modify the composition so that we can make it from
    # standard starting materials:

    # 1) We want to add calcium and sodium as CaCO3 and Na2CO3.

    # 2) We want to make a fraction f of total iron Fe^{57},
    # where the Fe57 is added as metallic Fe and the Fe(natural) is Fe2O3
    # We want to retain the correct total amount of iron, but we are free
    # to vary the amount of oxygen, so
    # FeO + m O -> f Fe(57) + (1-f)/2 Fe2O3

    # Therefore m = 3*(1-f)/2 - 1 = 0.5 - 1.5f
    f = 0.5
    m = 0.5 - 1.5*f

    # Here's where we add the carbonate and oxygen
    CO2_molar = KLB1.molar_composition['CaO'] + KLB1.molar_composition['Na2O']
    O_molar = KLB1.molar_composition['FeO']*m
    KLB1.add_components(composition_dictionary={'CO2': CO2_molar,
                                                'O': O_molar},
                        unit_type='molar')

    # Here's where we change the components:
    KLB1.change_component_set(['Na2CO3', 'CaCO3', 'Fe',
                               'Fe2O3', 'MgO', 'Al2O3', 'SiO2'])
    KLB1.print('weight', significant_figures=4, normalization_amount=1.)
