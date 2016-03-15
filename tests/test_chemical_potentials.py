from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


class ChemicalPotentials(BurnManTest):

    def test_chemical_potentials(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        per = minerals.SLB_2011.periclase()
        stv = minerals.SLB_2011.stishovite()

        assemblage = [bdg, per]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        mu_SiO2 = burnman.chemicalpotentials.chemical_potentials(
            assemblage, [{'Si': 1., 'O': 2.}])[0]

        self.assertFloatEqual(bdg.gibbs - per.gibbs, mu_SiO2)

    def test_fugacity(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        stv = minerals.SLB_2011.stishovite()

        assemblage = [bdg, stv]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        fugacity = burnman.chemicalpotentials.fugacity(stv, assemblage)
        self.assertFloatEqual(1., fugacity)

    def test_relative_fugacity(self):
        bdg = minerals.SLB_2011.mg_fe_perovskite()
        stv = minerals.SLB_2011.stishovite()

        assemblage = [bdg, stv]

        bdg.set_composition([1.0, 0.0, 0.0])

        for phase in assemblage:
            phase.set_state(25.e9, 2000.)

        relative_fugacity = burnman.chemicalpotentials.relative_fugacity(
            stv, assemblage, assemblage)
        self.assertFloatEqual(1., relative_fugacity)

if __name__ == '__main__':
    unittest.main()
