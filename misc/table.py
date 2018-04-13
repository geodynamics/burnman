# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


""" Generates a text table with mineral properties. Run 'python table.py latex' to write a tex version of the table to mytable.tex """
from __future__ import absolute_import
from __future__ import print_function


import os
import sys
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

    # Create a dictionary of the equation of state parameters to insert into the table
    # Equations of state not included in this dictionary will be ignored
    eos_params = {'hp_tmt': ['V_0', 'K_0', 'Kprime_0', 'Kdprime_0', 'molar_mass', 'n', 'Cp'],
                  'slb2': ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0',
                           'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0', 'eta_s_0'],
                  'slb3': ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0',
                           'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0', 'eta_s_0'],
                  'mgd2': ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0',
                           'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0', 'eta_s_0'],
                  'mgd3': ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0',
                           'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0', 'eta_s_0'],
                  'cork': ['cork_params', 'cork_T', 'cork_P', 'Cp'],
                  'dks_s': ['V_0', 'T_0', 'E_0', 'S_0', 'K_0',
                            'Kprime_0', 'Kdprime_0', 'n', 'Cv', 'grueneisen_0', 'q_0'],
                  'dks_l': ['V_0', 'T_0', 'O_theta', 'O_f', 'm', 'zeta_0', 'xi', 'Tel_0', 'eta'],
                  'aa': ['T_0', 'S_0', 'V_0', 'K_S', 'Kprime_S', 'Kprime_prime_S', 'grueneisen_0', 'grueneisen_prime', 'grueneisen_n']}

    libs = dir(minerals)
    for l in libs:

        names = set()
        phasenames = []

        mineralgroup = getattr(minerals, l)

        if mineralgroup.__class__.__name__ == "module":
            for m in dir(mineralgroup):
                mineral = getattr(mineralgroup, m)
                if inspect.isclass(mineral) and mineral != burnman.Mineral and mineral != burnman.combinedmineral.CombinedMineral and issubclass(mineral, burnman.Mineral) \
                        and not issubclass(mineral, burnman.mineral_helpers.HelperSpinTransition):
                    if issubclass(mineral, burnman.SolidSolution):
                        continue
                    # print mineral.__module__ + mineral.__name__
                    name1 = str(mineralgroup.__name__) + "." + str(m)
                    name = name1.replace("burnman.minerals.", "")
                    # name = mineral.__module__.replace("minlib_","").replace("burnman.","").replace("minerals.","") + "." + mineral.__name__
                    # print mineral, name, name1
                    if name not in names:
                        names.add(name)
                        try:
                            x = mineral()
                            phasenames.append((name, x))
                        except:
                            print("Could not create '%s'" % name)

            # The following groups minerals from each module into groups based on eos and
            # prints a separate table for each eos group.
            list_eoses = [phasename[1].params['equation_of_state'] for phasename in phasenames]
            for (eos, params) in sorted(eos_params.items(), key=lambda x: x[0]):
                eos_phasenames = [phasenames[i] for i, e in enumerate(list_eoses) if e == eos]
                if len(eos_phasenames) > 0:
                    table = []
                    tablel = []
                    
                    table.append(['Name ({0} equation of state)'.format(eos)] + params)
                    tablel.append([])
                    
                    sortedlist = sorted(eos_phasenames, key=lambda x: x[0])

                    for (name, p) in sortedlist:
                        p.set_state(1e9, 300)
                        row = create_list(name, p)
                        table.append(row)
                        tablel.append(row)

                    if (len(sys.argv) == 1):
                        tools.pretty_print_table(table, False)
                    elif sys.argv[1] == "tab":
                        tools.pretty_print_table(table, True)
                        
