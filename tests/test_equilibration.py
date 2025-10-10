import unittest
from util import BurnManTest

import burnman
from burnman import equilibrate
from burnman.minerals import HP_2011_ds62, SLB_2011

import numpy as np


class equilibration(BurnManTest):
    def test_univariant_line(self):
        andalusite = HP_2011_ds62.andalusite()
        kyanite = HP_2011_ds62.ky()
        composition = kyanite.formula

        pressures = np.linspace(1.0e5, 1.0e9, 11)

        assemblage = burnman.Composite([andalusite, kyanite])
        equality_constraints = [
            ("P", pressures),
            ("phase_fraction", (andalusite, np.array([0.0]))),
        ]
        sols, prm = equilibrate(composition, assemblage, equality_constraints)

        Ts = [sol.assemblage.temperature for sol in sols]
        Ts_ref = [
            464.95783870074905,
            544.0828494722284,
            623.6421086627697,
            703.7607424193357,
            784.4975351197959,
            865.8789417005959,
            947.9152419887187,
            1030.608906228574,
            1113.95923718214,
            1197.9651099824446,
            1282.6266945380332,
        ]
        self.assertArraysAlmostEqual(Ts, Ts_ref)

    def test_invariant(self):
        sillimanite = HP_2011_ds62.sill()
        andalusite = HP_2011_ds62.andalusite()
        kyanite = HP_2011_ds62.ky()

        composition = sillimanite.formula
        assemblage = burnman.Composite([sillimanite, andalusite, kyanite])
        equality_constraints = [
            ("phase_fraction", (kyanite, np.array([0.0]))),
            ("phase_fraction", (sillimanite, np.array([0.0]))),
        ]

        sol, prm = equilibrate(composition, assemblage, equality_constraints)

        self.assertArraysAlmostEqual(
            sol.x, [4.30671426e08, 8.09342875e02, 0.0, 1.0, 0.0]
        )

    def test_univariant_line_fo(self):
        forsterite = SLB_2011.forsterite()
        periclase = SLB_2011.periclase()
        bridgmanite = SLB_2011.mg_bridgmanite()
        composition = forsterite.formula

        temperatures = np.linspace(300.0, 2000.0, 11)

        assemblage = burnman.Composite([forsterite, periclase, bridgmanite])
        equality_constraints = [
            ("T", temperatures),
            ("phase_fraction", (forsterite, np.array([0.0]))),
        ]
        sols, prm = equilibrate(composition, assemblage, equality_constraints)

        Ps = np.array([sol.assemblage.pressure for sol in sols])
        Ps_ref = [
            17405868549,
            17657896452,
            17941848050,
            18231743600,
            18517274303,
            18793435710,
            19057314295,
            19306923194,
            19540715743,
            19757356711,
            19955602937,
        ]
        self.assertArraysAlmostEqual(Ps, Ps_ref)

    def test_ol_wad_eqm(self):
        ol = SLB_2011.mg_fe_olivine()
        wad = SLB_2011.mg_fe_wadsleyite()

        assemblage = burnman.Composite([ol, wad], [0.7, 0.3])
        ol.set_composition([0.5, 0.5])
        wad.set_composition([0.6, 0.4])

        assemblage.set_state(10.0e9, 1200.0)
        equality_constraints = [
            ("P", 10.0e9),
            (
                "phase_composition",
                (ol, (["Mg_A", "Fe_A"], [0.0, 1.0], [1.0, 1.0], 0.45)),
            ),
        ]
        composition = {"Mg": 1.0, "Fe": 1.0, "Si": 1.0, "O": 4.0}

        sol, prm = equilibrate(composition, assemblage, equality_constraints)
        self.assertArraysAlmostEqual(
            [assemblage.temperature, ol.molar_fractions[1], wad.molar_fractions[1]],
            [1620.532183457096, 0.45, 0.6791743],
        )

    def test_binary_solution_convergence(self):
        mg_bdg = SLB_2011.mg_bridgmanite()
        fe_bdg = SLB_2011.fe_bridgmanite()
        mg_ppv = SLB_2011.mg_post_perovskite()
        fe_ppv = SLB_2011.fe_post_perovskite()
        bdg_endmembers = [[mg_bdg, "[Mg]"], [fe_bdg, "[Fe]"]]
        ppv_endmembers = [[mg_ppv, "[Mg]"], [fe_ppv, "[Fe]"]]

        bdg = burnman.Solution(
            "bdg",
            solution_model=burnman.classes.solutionmodel.IdealSolution(bdg_endmembers),
        )
        ppv = burnman.Solution(
            "ppv",
            solution_model=burnman.classes.solutionmodel.IdealSolution(ppv_endmembers),
        )

        composition = {"Mg": 0.9, "Fe": 0.1, "Si": 1.0, "O": 3.0}
        assemblage = burnman.Composite(
            phases=[bdg, ppv], fractions=[0.9, 0.1], name="MgSiO3-pv-ppv-assemblage"
        )
        bdg.set_composition([0.9, 0.1])
        ppv.set_composition([0.9, 0.1])

        temperatures = np.linspace(1000, 4000, 4)
        assemblage.set_state(120.0e9, temperatures[0])

        equality_constraints = [("T", temperatures), ("phase_fraction", (ppv, 0.0))]

        sol, prm = equilibrate(composition, assemblage, equality_constraints)
        self.assertFalse(any(not s.success for s in sol))


if __name__ == "__main__":
    unittest.main()
