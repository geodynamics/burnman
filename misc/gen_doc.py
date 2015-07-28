
""" generates a list with the examples """
from __future__ import print_function

import os, sys
import os
import re

sys.path.insert(1,os.path.abspath('../examples'))

# order the examples sensibly:

ordered_examples = ['example_beginner.py', \
                    'example_geotherms.py', \
                    'example_seismic.py', \
                    'example_composition.py', \
                    'example_user_input_material.py', \
                    'example_averaging.py', \
                    'example_woutput.py', \
                    'example_compare_all_methods.py', \
                    'example_spintransition.py', \
                    'example_partition_coef.py', \
                    'example_compare_enstpyro.py', \
                    'example_optimize_pv.py', \
                        ]

for ex in ordered_examples:
    #print "*",ex
    print(__import__(ex.replace(".py","")).__doc__)



#check we do not forget an example

all = os.listdir('../examples/')
examples = [k for k in all if re.match("^example(.*)\.py$", k)]
not_listed = [x for x in examples if (not x in ordered_examples)]

for l in not_listed:
    sys.stderr.write("WARNING EXAMPLE NOT LISTED: "+l+"\n")
