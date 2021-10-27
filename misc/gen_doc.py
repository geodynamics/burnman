
""" generates a list with the examples """
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import re

sys.path.insert(1, os.path.abspath('../examples'))

# order the examples sensibly:
ordered_examples = ['example_mineral.py',
                    'example_gibbs_modifiers.py',
                    'example_solid_solution.py',
                    'example_composite.py',
                    'example_geotherms.py',
                    'example_composition.py',
                    'example_fit_composition.py',
                    'example_beginner.py',
                    'example_seismic.py',
                    'example_anisotropy.py',
                    'example_anisotropic_mineral.py',
                    'example_spintransition.py',
                    'example_spintransition_thermal.py',
                    'example_user_input_material.py',
                    'example_averaging.py',
                    'example_woutput.py',
                    'example_compare_all_methods.py',
                    'example_optimize_pv.py',
                    'example_fit_data.py',
                    'example_fit_eos.py',
                    'example_grid.py',
                    'example_chemical_potentials.py',
                    'example_perplex.py',
                    'example_polytopetools.py',
                    'example_equilibrate.py',
                    'example_olivine_binary.py',
                    'example_geodynamic_adiabat.py',
                    'example_layer.py',
                    'example_build_planet.py'
                    ]

for ex in ordered_examples:
    # print "*",ex
    print(__import__(ex.replace(".py", "")).__doc__)


# check we do not forget an example
all = os.listdir('../examples/')
examples = list(filter(lambda k: re.match("^example(.*)\.py$", k), all))
not_listed = list(filter(lambda x: (x not in ordered_examples), examples))

sys.stderr.write('The following examples are not included '
                 'in the automatically generated documentation:\n')
for example in not_listed:
    sys.stderr.write(f'   {example}\n')
