from __future__ import absolute_import
import unittest
from util import BurnManTest

import burnman_path
from burnman import minerals

from burnman.tools.chemistry import chemical_potentials
from burnman.tools.chemistry import relative_fugacity
from burnman.tools.chemistry import fugacity

assert burnman_path  # silence pyflakes warning


class ChemicalPotentials(BurnManTest):

    def test_chemical_potentials(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        per = minerals.SLB_2011.periclase()

        assemblage = [bdg, per]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        mu_SiO2 = chemical_potentials(assemblage, [{'Si': 1., 'O': 2.}])[0]

        self.assertFloatEqual(bdg.gibbs - per.gibbs, mu_SiO2)

    def test_fugacity(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        stv = minerals.SLB_2011.stishovite()

        assemblage = [bdg, stv]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        self.assertFloatEqual(1., fugacity(stv, assemblage))

    def test_relative_fugacity(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        stv = minerals.SLB_2011.stishovite()

        assemblage = [bdg, stv]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        self.assertFloatEqual(1., relative_fugacity(stv, assemblage,
                                                    assemblage))


if __name__ == '__main__':
    unittest.main()
