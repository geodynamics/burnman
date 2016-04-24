
""" generates a list with the examples """
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import os
import re

sys.path.insert(1, os.path.abspath('../examples'))

# order the examples sensibly:

ordered_examples = ['example_beginner.py',
                    'example_geotherms.py',
                    'example_seismic.py',
                    'example_composition.py',
                    'example_user_input_material.py',
                    'example_averaging.py',
                    'example_woutput.py',
                    'example_compare_all_methods.py',
                    'example_optimize_pv.py',
                    'example_fit_data.py',
                    'example_grid.py',
                    'example_chemical_potentials.py',
                    'example_solid_solution.py',
                    'example_build_planet.py'
                    ]

for ex in ordered_examples:
    # print "*",ex
    print(__import__(ex.replace(".py", "")).__doc__)


# check we do not forget an example
all = os.listdir('../examples/')
examples = list(filter(lambda k: re.match("^example(.*)\.py$", k), all))
not_listed = list(filter(lambda x: (not x in ordered_examples), examples))

for l in not_listed:
    sys.stderr.write("WARNING EXAMPLE NOT LISTED: " + l + "\n")
