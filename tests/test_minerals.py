from __future__ import absolute_import
import unittest
import inspect
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals
from util import BurnManTest

# Instantiate every mineral in the mineral library


class instantiate_minerals(BurnManTest):

    def test_instantiate_minerals(self):

        P = 25.e9
        T = 1500.
        mineral_libraries = dir(
            minerals)  # Get a list of mineral libraries available

        for minlib_name in mineral_libraries:

            minlib_ = getattr(minerals, minlib_name)

            # Pull all the members of the given mineral library
            member_list = [getattr(minlib_, m) for m in dir(minlib_)]

            # Now restrict the members list only to minerals that can be
            # instantiated
            mineral_list = [m for m in member_list if inspect.isclass(m)
                            and issubclass(m, burnman.Mineral)
                            and m is not burnman.mineral.Mineral
                            and m is not burnman.solidsolution.SolidSolution
                            and m is not burnman.combinedmineral.CombinedMineral]

            for mineral_ in mineral_list:
                m = mineral_()  # instantiate

                # Call set_composition if necessary
                if isinstance(m, burnman.solidsolution.SolidSolution):
                    m.set_composition([1. / m.n_endmembers] * m.n_endmembers)

                # test that it works
                m.set_state(P, T)
                V = m.molar_volume


if __name__ == '__main__':
    unittest.main()
