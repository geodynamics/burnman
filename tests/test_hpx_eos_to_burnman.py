import ast
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from burnman import Solution
from burnman.classes.solutionmodel import IdealSolution
from burnman.data.input_raw_endmember_datasets.hpx_eos_to_burnman import convert_dataset
from burnman.minerals import HGP_2018_ds633, HP_2011_fluids


def solution_block(name, label="shared", recipe="fo 1", adjustment="0 0 0", n_make=1):
    """A two-endmember HPx-eos solution with one generated endmember."""
    return f"""{name} 2 1
q({name}) 0.5
p({label}) 1 1 0 1 1 q
p(fo) 1 1 1 1 -1 q
sf
W({label},fo) 1 0 0
2
xMg 1 1 0 1 1 q
xFe 1 1 1 1 -1 q
{label} 1 1 xMg 1
check 1
make {n_make} {recipe}
delG(make) {adjustment}
fo 1 1 xFe 1
check 0
end
"""


class hpx_eos_porter(unittest.TestCase):
    def convert(self, *blocks):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "tc-test.txt"
            path.write_text("".join(blocks))
            source = convert_dataset(path, "HGP_2018_ds633")

        # Every generated object and import must have its own module identifier.
        definitions = []
        for statement in ast.parse(source).body:
            if isinstance(statement, ast.ClassDef):
                definitions.append(statement.name)
            elif isinstance(statement, ast.Assign):
                definitions.extend(target.id for target in statement.targets)
            elif isinstance(statement, (ast.Import, ast.ImportFrom)):
                definitions.extend(
                    alias.asname or alias.name for alias in statement.names
                )
        self.assertEqual(len(definitions), len(set(definitions)))

        namespace = {"__package__": "burnman.minerals"}
        exec(compile(source, "ported_dataset.py", "exec"), namespace)
        return namespace

    def test_distinct_recipes_and_identical_definitions(self):
        namespace = self.convert(
            solution_block("first"),
            solution_block("second", recipe="en 1"),
            solution_block("third"),
            solution_block("fourth", recipe="fo 1/2"),
        )
        first = namespace["first"]().endmembers[0][0]
        second = namespace["second"]().endmembers[0][0]
        third = namespace["third"]().endmembers[0][0]
        fourth = namespace["fourth"]().endmembers[0][0]
        self.assertEqual(first.formula, HGP_2018_ds633.fo().formula)
        self.assertEqual(second.formula, HGP_2018_ds633.en().formula)
        self.assertIs(first, third)
        self.assertIsNot(first, second)
        self.assertIsNot(first, fourth)
        self.assertEqual(fourth.formula["Mg"], 1.0)
        self.assertEqual(first.name, "shared")
        self.assertEqual(second.name, "shared")

    def test_distinct_energy_adjustments(self):
        namespace = self.convert(
            solution_block("first"),
            solution_block("second", adjustment="1 0 0"),
            solution_block("third", adjustment="0 1 0"),
            solution_block("fourth", adjustment="0 0 1"),
        )
        endmembers = [
            namespace[name]().endmembers[0][0]
            for name in ("first", "second", "third", "fourth")
        ]
        self.assertEqual(len({id(m) for m in endmembers}), 4)
        for mineral in endmembers:
            mineral.set_state(1.0e5, 300.0)
        baseline = endmembers[0].molar_gibbs
        for mineral, difference in zip(endmembers[1:], (1000.0, 300000.0, 1.0)):
            self.assertAlmostEqual(mineral.molar_gibbs - baseline, difference)

    def test_identical_recipes_with_reordered_components(self):
        namespace = self.convert(
            solution_block("first", recipe="fo 1/2 en 1/2", n_make=2),
            solution_block("second", recipe="en 1/2 fo 1/2", n_make=2),
        )
        self.assertIs(
            namespace["first"]().endmembers[0][0],
            namespace["second"]().endmembers[0][0],
        )

    def test_suffixes_do_not_shadow_source_names(self):
        namespace = self.convert(
            solution_block("first"),
            solution_block("second", recipe="en 1"),
            solution_block("third", label="shared_second", recipe="fs 1"),
            solution_block("second", recipe="fa 1"),
        )
        for name, mineral in (
            ("first", HGP_2018_ds633.fo()),
            ("second", HGP_2018_ds633.en()),
            ("third", HGP_2018_ds633.fs()),
            ("second_solution", HGP_2018_ds633.fa()),
        ):
            self.assertEqual(
                namespace[name]().endmembers[0][0].formula, mineral.formula
            )
        self.assertEqual(
            namespace["shared_second"].formula, HGP_2018_ds633.fs().formula
        )

    def test_endmembers_do_not_shadow_solution_classes(self):
        namespace = self.convert(
            solution_block("first", label="second"),
            solution_block("second", recipe="en 1"),
        )
        self.assertTrue(issubclass(namespace["second"], Solution))
        self.assertEqual(
            namespace["first"]().endmembers[0][0].formula, HGP_2018_ds633.fo().formula
        )

    def test_import_names_and_keywords(self):
        for name in (
            "np",
            "array",
            "nan",
            "Mineral",
            "Solution",
            "CombinedMineral",
            "SymmetricRegularSolution",
            "AsymmetricRegularSolution",
            "HGP_2018_ds633",
            "HP_2011_fluids",
            "class",
        ):
            with self.subTest(name=name):
                namespace = self.convert(solution_block(name, label=name))
                mineral = namespace[f"{name}_solution"]().endmembers[0][0]
                self.assertEqual(mineral.name, name)
                self.assertEqual(mineral.formula, HGP_2018_ds633.fo().formula)

    def test_ordering_helpers_do_not_shadow_other_objects(self):
        namespace = self.convert(
            solution_block("first", label="san_nood"),
            solution_block("san_nood", recipe="ordered san 1"),
            solution_block("third", recipe="disordered san 1"),
        )
        for name in ("san_nood", "third"):
            mineral = namespace[name]().endmembers[0][0]
            self.assertEqual(mineral.formula, HGP_2018_ds633.san().formula)
            mineral.set_state(1.0e5, 300.0)
            self.assertTrue(mineral.molar_gibbs < 0.0)
        self.assertEqual(
            namespace["first"]().endmembers[0][0].formula, HGP_2018_ds633.fo().formula
        )

    def test_ambiguous_local_endmember_names_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate endmember names.*ambiguous"):
            self.convert(solution_block("ambiguous", label="fo"))

    def test_variable_site_multiplicities(self):
        namespace = self.convert("""melt 3 1
x(melt) 0.2
w(melt) 0.3
p(fo) 1 1 1 2 -1 x -1 w
p(fa) 1 1 0 1 1 x
p(qL) 1 1 0 1 1 w
sf
W(fo,fa) 0 0 0
W(fo,qL) 0 0 0
W(fa,qL) 0 0 0
4
fac 1 1 1 1 -1 w
pmg 1 1 1 2 -1 x -1 w
pfe 1 1 0 1 1 x
pw 1 1 0 1 1 w
fo 1 2 fac 1 pmg 1
fa 1 2 fac 1 pfe 1
qL 1 1 pw 2
end
""")
        model = namespace["melt"]().solution_model
        activities = IdealSolution.activities(model, 1.0e9, 1500.0, [0.5, 0.2, 0.3])
        for actual, expected in zip(activities, [0.35, 0.14, 0.09]):
            self.assertAlmostEqual(actual, expected)

    def test_modern_headers_and_explicit_activity_denominators(self):
        namespace = self.convert("""melt 3 2 silicate_melt 100 % TC 3.51
x(melt) -
w(melt) -
p(fo) 1 1 1 2 -1 x -1 w
p(fa) 1 1 0 1 1 x
p(qL) 1 1 0 1 1 w
sf
W(fo,fa) 0 0 0
W(fo,qL) 0 0 0
W(fa,qL) 0 0 0
8
pmg 1 1 1 2 -1 x -1 w
pfe 1 1 0 1 1 x
mgM 1 1 1 2 -1 x -1 w
feM 1 1 0 1 1 x
sumM 1 1 1 1 -1 w
sumT 1 1 1 1 -1 w
xv 1 1 1 1 -1 w
xh 1 1 0 1 1 w
fo 1 5 pmg 1 mgM 1 sumM -1 sumT -1 xv 2
fa 1 5 pfe 1 feM 1 sumM -1 sumT -1 xv 2
qL 1 1 xh 2
end
""")
        model = namespace["melt"]().solution_model
        activities = IdealSolution.activities(model, 1.0e9, 1400.0, [0.5, 0.2, 0.3])
        for actual, expected in zip(activities, [0.25, 0.04, 0.09]):
            self.assertAlmostEqual(actual, expected)

    def test_buffer_and_solution_name_collisions(self):
        namespace = self.convert(
            solution_block("shared"),
            """shared 1 buffer 1
make 1 equilibrium fo 2
DQF 1 0 0
""",
        )
        self.assertTrue(issubclass(namespace["shared"], Solution))
        buffer = namespace["shared_buffer"]
        self.assertEqual(buffer.formula["Mg"], 4.0)
        self.assertEqual(buffer.property_modifiers[0][1]["delta_E"], 1000.0)

    def test_temperature_dependent_interaction_units(self):
        namespace = self.convert(
            solution_block("first").replace(
                "W(shared,fo) 1 0 0", "W(shared,fo) 1 -0.01 0.1"
            )
        )
        model = namespace["first"]().solution_model
        pressure, temperature, p = 1.0e9, 1400.0, [0.5, 0.5]
        ideal = IdealSolution(model.endmembers).excess_gibbs_energy(
            pressure, temperature, p
        )
        excess = model.excess_gibbs_energy(pressure, temperature, p)
        self.assertAlmostEqual(
            excess - ideal, 0.25 * (1000.0 - 10.0 * temperature + 1.0e-6 * pressure)
        )

    def test_multiple_energy_corrections_are_added(self):
        block = solution_block("first", adjustment="1 0 0")
        block = block.replace("delG(make) 1 0 0", "delG(make) 1 0 0\nDQF 2 0 0")
        namespace = self.convert(block)
        correction = namespace["first"]().endmembers[0][0].property_modifiers[0][1]
        self.assertEqual(correction["delta_E"], 3000.0)

    def test_make_without_an_energy_correction(self):
        block = solution_block("first", recipe="fo 1/2")
        block = block.replace("check 1\n", "").replace("delG(make) 0 0 0\n", "")
        namespace = self.convert(block)
        self.assertEqual(namespace["first"]().endmembers[0][0].formula["Mg"], 1.0)

    def test_dqf_energy_correction(self):
        namespace = self.convert(
            solution_block("first", adjustment="-8 0 -0.25").replace(
                "delG(make)", "DQF"
            )
        )
        correction = namespace["first"]().endmembers[0][0].property_modifiers[0][1]
        self.assertEqual(correction["delta_E"], -8000.0)
        self.assertAlmostEqual(correction["delta_V"], -2.5e-6)

    def test_water_reference_in_direct_and_made_endmembers(self):
        direct = solution_block("direct", label="H2O")
        direct = direct.replace("check 1\nmake 1 fo 1\ndelG(make) 0 0 0\n", "")
        namespace = self.convert(direct, solution_block("made", recipe="H2O 1"))
        self.assertIsInstance(
            namespace["direct"]().endmembers[0][0], HP_2011_fluids.H2O
        )
        self.assertEqual(namespace["made"]().endmembers[0][0].formula, {"H": 2, "O": 1})

    def test_dqf_without_a_make_instruction(self):
        block = solution_block("first", label="foL", adjustment="-8 0 0")
        block = block.replace("check 1\nmake 1 fo 1\n", "").replace("delG(make)", "DQF")
        namespace = self.convert(block)
        mineral = namespace["first"]().endmembers[0][0]
        self.assertEqual(mineral.formula, HGP_2018_ds633.foL().formula)
        self.assertEqual(mineral.property_modifiers[0][1]["delta_E"], -8000.0)


if __name__ == "__main__":
    unittest.main()
