import unittest
from util import BurnManTest
import numpy as np

import burnman
from burnman.tools import output_seismo


class test_seismic(BurnManTest):
    def test_internal_depth_list(self):
        models = [
            burnman.seismic.PREM(),
            burnman.seismic.STW105(),
            burnman.seismic.AK135(),
            burnman.seismic.IASP91(),
        ]

        ref_depth_lists = {
            "PREM": [0.0, 6371000.0, 94],
            "STW105": [0.0, 6371000.0, 750],
            "AK135": [0.0, 6371000.0, 145],
            "IASP91": [0.0, 6371000.0, 152],
        }

        for model in models:
            dl = model.internal_depth_list()
            name = model.__class__.__name__
            stats = [min(dl), max(dl), len(dl)]
            # print model.__class__.__name__, stats
            self.assertArraysAlmostEqual(stats, ref_depth_lists[name])

    def test_evaluate(self):
        models = [
            burnman.seismic.PREM(),
            burnman.seismic.Fast(),
            burnman.seismic.Slow(),
            burnman.seismic.STW105(),
            burnman.seismic.AK135(),
            burnman.seismic.IASP91(),
        ]

        ref = {
            "PREM": [
                12817.6924,
                6932.8549000000003,
                10010.358087588364,
                5120.6290999999992,
            ],
            "Fast": [
                12795.360742611414,
                6941.4225568201909,
                9973.8053779446909,
                5120.6290999999992,
            ],
            "Slow": [
                12795.360742611414,
                6904.3291880013148,
                10008.07546360332,
                5120.6290999999992,
            ],
            "STW105": [
                12817.85012987013,
                6932.9342065251822,
                10010.486818122406,
                5120.6749255622426,
            ],
            "AK135": [
                12798.468686868688,
                6920.5212121212116,
                9997.1320353488773,
                5103.8969696969698,
            ],
            "IASP91": [12794.4, 6921.0, 9991.4805389391604],
        }

        for model in models:
            dl = model.internal_depth_list()
            name = model.__class__.__name__
            depth = 2000e3
            vars = ["v_p", "v_s", "v_phi", "density"]  # skip gravity
            if name == "IASP91":
                vars = vars[0:-1]  # skip density
            result = model.evaluate(vars, [depth])
            result = list(result.T[0])
            # print "'%s': %s," % (name, result)
            self.assertArraysAlmostEqual(result, ref[name])

    def test_output(self):
        mg_fe_perovskite = burnman.minerals.SLB_2011.mg_fe_perovskite()
        mg_fe_perovskite.set_composition([0.9, 0.1, 0])  # frac_mg, frac_fe, frac_al
        rock = burnman.Composite([mg_fe_perovskite], [1.0])
        depths = np.linspace(2890e3, 670e3, 20)
        lower_mantle = burnman.Layer(name="Pv LM", radii=6371.0e3 - depths)
        lower_mantle.set_material(rock)
        lower_mantle.set_temperature_mode(
            temperature_mode="adiabatic", temperature_top=1900.0
        )
        lower_mantle.set_pressure_mode(
            pressure_mode="self-consistent", pressure_top=2.4e10, gravity_bottom=10.0
        )
        lower_mantle.make()

        # Make sure that all the writing functions return a value
        tvel = output_seismo.tvel_formatted_data_and_header(
            lower_mantle, background_model=burnman.seismic.PREM()
        )
        axisem = output_seismo.axisem_formatted_data_and_reference([lower_mantle])
        mineos = output_seismo.mineos_formatted_data_and_reference([lower_mantle])

        # Write to nothing
        output_seismo.write_tvel_file(lower_mantle, [], burnman.seismic.PREM())
        output_seismo.write_axisem_input([lower_mantle], [], verbose=False)
        output_seismo.write_mineos_input([lower_mantle], [], verbose=False)


if __name__ == "__main__":
    unittest.main()
