# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


""" Generates a text table with mineral properties. Run 'python table.py latex' to write a tex version of the table to mytable.tex """
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
import inspect

from burnman import minerals
from burnman import tools

if __name__ == "__main__":

    def create_list(name, mineral):
        ownname = mineral.to_string().replace(
            "'", "").replace("burnman.minerals.", "")
        if name != ownname:
            name = name + " (" + ownname + ")"
        row = [name]
        for param in params:

            row.append(str(p.params[param] if param in p.params else ""))
        return row

    libs = dir(minerals)
    for l in libs:

        names = set()
        phasenames = []

        mineralgroup = getattr(minerals, l)

        if mineralgroup.__class__.__name__ == "module":
            for m in dir(mineralgroup):
                mineral = getattr(mineralgroup, m)
                if inspect.isclass(mineral) and mineral != burnman.Mineral and issubclass(mineral, burnman.Mineral) \
                        and issubclass(mineral, burnman.mineral_helpers.HelperSpinTransition) == False:
                    if issubclass(mineral, burnman.SolidSolution):
                        continue
                    # print mineral.__module__ + mineral.__name__
                    name1 = str(mineralgroup.__name__) + "." + str(m)
                    name = name1.replace("burnman.minerals.", "")
                    # name = mineral.__module__.replace("minlib_","").replace("burnman.","").replace("minerals.","") + "." + mineral.__name__
                    # print mineral, name, name1
                    if not name in names:
                        names.add(name)
                        try:
                            x = mineral()
                            phasenames.append((name, x))
                        except:
                            print("Could not create '%s'" % name)

            eos = phasenames[0][1].params['equation_of_state']
            if eos == 'hp_tmt':
                params = [
                    'V_0', 'K_0', 'Kprime_0', 'Kdprime_0', 'molar_mass', 'n', 'Cp']
            elif eos == 'slb2' or eos == 'slb3' or eos == 'mgd2' or eos == 'mgd3':
                params = ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0',
                          'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0', 'eta_s_0']
            elif eos == 'cork':
                params = ['cork_params', 'cork_T', 'cork_P', 'Cp']

            table = [['Name'] + params]
            tablel = []

            sortedlist = sorted(phasenames, key=lambda x: x[0])

            for (name, p) in sortedlist:
                p.set_state(1e9, 300)
                row = create_list(name, p)
                table.append(row)
                tablel.append(row)

            if (len(sys.argv) == 1):
                tools.pretty_print_table(table, False)
            elif sys.argv[1] == "tab":
                tools.pretty_print_table(table, True)
