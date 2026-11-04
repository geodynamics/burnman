import unittest
import numpy as np

TOTALS = (
    "mass",
    "volume",
    "entropy",
    "enthalpy",
    "gibbs",
    "helmholtz",
    "internal_energy",
    "heat_capacity_p",
    "heat_capacity_v",
)
INTENSIVES = tuple("molar_" + name for name in TOTALS) + (
    "density",
    "thermal_expansivity",
    "isothermal_bulk_modulus_reuss",
    "isentropic_bulk_modulus_reuss",
    "isothermal_compressibility_reuss",
    "isentropic_compressibility_reuss",
    "grueneisen_parameter",
)
TENSORS = (
    "cell_vectors",
    "isothermal_stiffness_tensor",
    "isentropic_compliance_tensor",
    "thermal_expansivity_tensor",
    "grueneisen_tensor",
    "molar_isometric_heat_capacity",
)


class BurnManTest(unittest.TestCase):
    def assertFloatEqual(self, a, b, tol=1e-5, tol_zero=1e-16):
        self.assertAlmostEqual(a, b, delta=max(tol_zero, max(abs(a), abs(b)) * tol))

    def assertArraysAlmostEqual(self, a, b, tol=1e-5, tol_zero=1e-16):
        self.assertEqual(len(a), len(b))
        for i1, i2 in zip(a, b):
            self.assertFloatEqual(i1, i2, tol, tol_zero)

    def assertMaterialAmountScaling(
        self, material, state_setter=None, extra_properties=()
    ):
        """Check sample scaling, thermodynamic identities, copying and evaluation."""
        evaluate = state_setter is None
        if state_setter is None:
            state_setter = lambda: material.set_state(10e9, 1500.0)
        state_setter()
        original_amount = material.number_of_moles
        properties = INTENSIVES + tuple(extra_properties)
        expected = {
            name: np.array(getattr(material, name), copy=True)
            for name in properties + TOTALS
        }
        for factor in (0.1, 1.0, 2.5):
            with self.subTest(factor=factor):
                material.number_of_moles = factor * original_amount
                state_setter()
                for name in properties:
                    np.testing.assert_allclose(
                        getattr(material, name),
                        expected[name],
                        rtol=1e-7,
                        atol=1e-12,
                        err_msg=name,
                    )
                for name in TOTALS:
                    value = getattr(material, name)
                    self.assertTrue(np.isfinite(value), name)
                    self.assertFloatEqual(value, factor * expected[name], tol=1e-7)
                    self.assertFloatEqual(
                        value,
                        material.number_of_moles * getattr(material, "molar_" + name),
                        tol=1e-10,
                    )
                for alias, name in {
                    "V": "volume",
                    "S": "entropy",
                    "H": "enthalpy",
                    "C_p": "heat_capacity_p",
                    "C_v": "heat_capacity_v",
                }.items():
                    self.assertFloatEqual(
                        getattr(material, alias), getattr(material, name)
                    )
                self.assertFloatEqual(material.mass, material.density * material.volume)
                self.assertFloatEqual(
                    material.gibbs,
                    material.enthalpy - material.temperature * material.entropy,
                )
                self.assertFloatEqual(
                    material.helmholtz,
                    material.internal_energy - material.temperature * material.entropy,
                )
        copied = material.copy()
        self.assertFloatEqual(copied.number_of_moles, material.number_of_moles)
        for name in TOTALS:
            self.assertFloatEqual(getattr(copied, name), getattr(material, name))
        if evaluate:
            values = material.evaluate(TOTALS, [10e9], [1500.0])[:, 0]
            self.assertArraysAlmostEqual(
                values, [getattr(material, name) for name in TOTALS]
            )

    def assertReferenceAmountIndependence(
        self, material, references, state_setter=None
    ):
        """Internal endmembers and reference materials supply molar properties."""
        if state_setter is None:
            state_setter = lambda: material.set_state(10e9, 1500.0)
        state_setter()
        expected = [getattr(material, name) for name in INTENSIVES]
        amounts = [reference.number_of_moles for reference in references]
        try:
            for i, reference in enumerate(references):
                reference.number_of_moles = 0.3 + i
            state_setter()
            self.assertArraysAlmostEqual(
                [getattr(material, name) for name in INTENSIVES], expected
            )
        finally:
            for reference, amount in zip(references, amounts):
                reference.number_of_moles = amount
            state_setter()


class Huh(BurnManTest):
    def test(self):
        self.assertFloatEqual(5200.01, 5200.015)

    def test_to_zero(self):
        self.assertFloatEqual(0, 1e-20)


if __name__ == "__main__":
    unittest.main()
