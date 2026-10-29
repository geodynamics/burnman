import unittest
from util import BurnManTest
import inspect

import burnman
from burnman import minerals


# Instantiate every mineral in the mineral library
class instantiate_minerals(BurnManTest):
    def test_instantiate_minerals(self):
        P = 25.0e9
        T = 1500.0
        mineral_libraries = dir(minerals)  # Get a list of mineral libraries available

        for minlib_name in mineral_libraries:
            minlib_ = getattr(minerals, minlib_name)

            # Pull all the members of the given mineral library
            member_list = [getattr(minlib_, m) for m in dir(minlib_)]

            # Now restrict the members list only to minerals that can be
            # instantiated
            mineral_list = [
                m
                for m in member_list
                if inspect.isclass(m)
                and issubclass(m, burnman.Mineral)
                and not issubclass(m, burnman.RelaxedSolution)
                and m is not burnman.Mineral
                and m is not burnman.SolidSolution
                and m is not burnman.CombinedMineral
            ]

            for mineral_ in mineral_list:
                m = mineral_()  # instantiate

                # Call set_composition if necessary
                if isinstance(m, burnman.SolidSolution):
                    m.set_composition([1.0 / m.n_endmembers] * m.n_endmembers)

                # Respect the common domain of the material and all of its
                # endmembers (including components of CombinedMinerals).
                parameter_sets = []
                pending = [m]
                while pending:
                    component = pending.pop()
                    parameter_sets.append(component.params)
                    if isinstance(component, burnman.CombinedMineral):
                        pending.append(component.mixture)
                    elif isinstance(component, burnman.SolidSolution):
                        pending.extend(
                            endmember for endmember, _ in component.endmembers
                        )
                pressure = min(
                    P, min(params.get("P_max", P) for params in parameter_sets)
                )
                temperature = min(
                    T, min(params.get("T_max", T) for params in parameter_sets)
                )
                temperature = max(
                    temperature,
                    max(params.get("T_min", 0.0) for params in parameter_sets),
                )
                with self.subTest(library=minlib_name, mineral=mineral_.__name__):
                    m.set_state(pressure, temperature)
                    V = m.molar_volume


class SLB2022LandauVolumes(BurnManTest):
    def test_critical_volumes_are_converted_from_cubic_centimetres(self):
        # HeFESTo parameter set 010121 gives critical volumes in cm^3/mol.
        # The molar quantities already use one chemical formula unit.
        pressure, temperature = 10.0e9, 800.0
        for name, critical_volume, entropy, critical_temperature in (
            ("neph", 0.8, 10.0, 467.0),
            ("qtz", 1.35936, 5.76, 847.0),
        ):
            with self.subTest(mineral=name):
                phase = getattr(minerals.SLB_2022, name)()
                params = phase.property_modifiers[0][1]
                self.assertAlmostEqual(
                    params["V_D"] * 1.0e6, critical_volume, places=12
                )
                # HeFESTo uses pressure in GPa and volume in cm^3/mol.
                tc = critical_temperature + critical_volume * 10.0 / (0.001 * entropy)
                self.assertAlmostEqual(
                    params["Tc_0"] + params["V_D"] * pressure / params["S_D"], tc
                )
                phase.set_state(pressure, temperature)
                self.assertGreater(phase.molar_volume, 0.0)
                self.assertGreater(phase.density, 0.0)
                self.assertGreater(phase.shear_wave_velocity, 0.0)


class VolumeStateModifiers(BurnManTest):
    def test_volume_state_does_not_add_modifier_volume_twice(self):
        mineral = minerals.SLB_2011.forsterite()
        mineral.property_modifiers = [
            ["linear", {"delta_V": 2.0e-7, "delta_S": 0.5, "delta_E": 100.0}]
        ]
        mineral.set_state(3.0e9, 700.0)
        volume = mineral.molar_volume
        gibbs, helmholtz = mineral.molar_gibbs, mineral.molar_helmholtz
        mineral.set_state_with_volume(volume, 700.0)
        self.assertAlmostEqual(mineral.molar_volume / volume, 1.0, places=12)
        self.assertAlmostEqual(mineral.pressure / 3.0e9, 1.0, places=10)
        self.assertAlmostEqual(mineral.molar_gibbs / gibbs, 1.0, places=12)
        self.assertAlmostEqual(mineral.molar_helmholtz / helmholtz, 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
