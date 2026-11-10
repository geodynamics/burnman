#!/usr/bin/env python3
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
This script converts the original HHPH2013 THERMOCALC files into a BurnMan
readable dataset.

Run from any directory; input defaults are relative to this script:

    python /path/to/HHPH2013data_to_burnman.py -o HHPH_2013.py

The endmember dataset file tc-ds62.txt supplies endmembers and their modifiers.
The shared HPx-eos converter reads the solution models from tc-NCFMAS.txt,
contained in
https://hpxeosandthermocalc.org/wp-content/uploads/2022/06/hollandetal2013.zip
Selected endmembers are those in Table 1 of the original paper, plus coesite.

With -o, HHPH_2013_cov.py is also written beside the output module. Without
-o, the endmember/solution source is printed to stdout. Nothing is downloaded
and importing this converter performs no conversion. Generated modules use SI.
"""

import argparse
import ast
from decimal import Decimal
import math
from pathlib import Path
import pprint
import sys
from types import ModuleType

import numpy as np

SOURCE_DIR = Path(__file__).resolve().parent
THERMOCALC_DIR = SOURCE_DIR / "tc-thermoinput-hhph2013"
ENDMEMBER_FILE = THERMOCALC_DIR / "tc-ds62.txt"
SOLUTION_FILE = THERMOCALC_DIR / "tc-NCFMAS.txt"
ARCHIVE_URL = (
    "https://hpxeosandthermocalc.org/wp-content/uploads/2022/06/hollandetal2013.zip"
)
COMPONENTS = (
    "Si",
    "Ti",
    "Al",
    "Fe",
    "Mg",
    "Mn",
    "Ca",
    "Na",
    "K",
    "O",
    "H",
    "C",
    "Cl",
    "e-",
    "Ni",
    "Zr",
    "S",
    "Cu",
    "Cr",
)
# Table 1 endmembers in published order, followed by coesite. This selects
# the subset; all thermodynamic parameters come from tc-ds62.txt.
ENDMEMBER_NAMES = (
    "fo",
    "fa",
    "mwd",
    "fwd",
    "mrw",
    "frw",
    "mpv",
    "fpv",
    "apv",
    "npv",
    "cpv",
    "mak",
    "fak",
    "maj",
    "nagt",
    "py",
    "alm",
    "gr",
    "en",
    "cen",
    "hen",
    "hfs",
    "fs",
    "mgts",
    "di",
    "hed",
    "jd",
    "cats",
    "stv",
    "macf",
    "mscf",
    "fscf",
    "nacf",
    "cacf",
    "manal",
    "nanal",
    "msnal",
    "fsnal",
    "canal",
    "per",
    "fper",
    "cor",
    "mcor",
    "coe",
)


def _number(value, scale="1"):
    """Scale decimal values once, without binary conversion artefacts."""
    result = float(Decimal(str(value)) * Decimal(scale))
    if not math.isfinite(result):
        raise ValueError(f"Non-finite dataset value: {value}")
    return result


def _read_thermocalc_file(path):
    """Return format-5 endmember records and the remaining covariance lines."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < 3 or lines[2].strip() != "5":
        raise ValueError(f"{path}: expected THERMOCALC ds62 format 5")
    count = int(lines[0].split()[0])
    if len(lines) < 3 + 4 * count:
        raise ValueError(f"{path}: truncated endmember records")
    records = {}
    for index in range(count):
        row = [line.split() for line in lines[3 + 4 * index : 7 + 4 * index]]
        if not row[0]:
            raise ValueError(f"{path}: missing endmember header at record {index}")
        name = row[0][0]
        if name in records:
            raise ValueError(f"{path}: duplicate endmember {name!r}")
        records[name] = row
    return records, lines[3 + 4 * count :]


def read_thermocalc_endmembers(path=ENDMEMBER_FILE, names=None):
    """Read four-line records in the authors' ds62 (format 5) source file.

    S and all Cp coefficients use kJ-based units: unlike printed Table 1,
    Cp_b has no additional column scaling. Expansion is already in 1/K.
    Landau and Bragg-Williams modifiers are retained. The complete covariance
    is read separately by read_thermocalc_covariance().
    """
    records, _ = _read_thermocalc_file(path)
    if names is None:
        names = ENDMEMBER_NAMES
    result = []
    for name in names:
        if name not in records:
            raise ValueError(f"{path}: missing endmember {name!r}")
        header, reference, cp, eos = records[name]
        if header[-1] != "0" or len(header[2:-1]) % 2:
            raise ValueError(f"{path}: malformed composition for {name!r}")
        formula = {}
        for component, amount in zip(header[2:-1:2], header[3:-1:2]):
            component_index = int(component)
            if not 1 <= component_index <= len(COMPONENTS):
                raise ValueError(f"{path}: unknown component {component}")
            element = COMPONENTS[component_index - 1]
            formula[element] = formula.get(element, 0.0) + _number(amount)
        if len(reference) != 3 or len(cp) != 4 or len(eos) < 5:
            raise ValueError(f"{path}: malformed thermodynamic record for {name!r}")
        flag, modifier = int(eos[4]), eos[5:]
        property_modifiers = []
        if flag == 1 and len(modifier) == 3:
            property_modifiers = [
                [
                    "landau_hp",
                    dict(
                        P_0=1e5,
                        T_0=298.15,
                        Tc_0=_number(modifier[0]),
                        S_D=_number(modifier[1], "1e3"),
                        V_D=_number(modifier[2], "1e-5"),
                    ),
                ]
            ]
        elif flag == 2 and len(modifier) == 6:
            property_modifiers = [
                [
                    "bragg_williams",
                    dict(
                        deltaH=_number(modifier[0], "1e3"),
                        deltaV=_number(modifier[1], "1e-5"),
                        Wh=_number(modifier[2], "1e3"),
                        Wv=_number(modifier[3], "1e-5"),
                        n=_number(modifier[4]),
                        factor=_number(modifier[5]),
                    ),
                ]
            ]
        elif flag != 0 or modifier:
            raise ValueError(f"{path}: unsupported modifier flag {flag} for {name!r}")
        result.append(
            dict(
                name=name,
                formula=formula,
                params=dict(
                    H_0=_number(reference[0], "1e3"),
                    S_0=_number(reference[1], "1e3"),
                    V_0=_number(reference[2], "1e-5"),
                    Cp=[_number(x, "1e3") for x in cp],
                    a_0=_number(eos[0]),
                    K_0=_number(eos[1], "1e8"),
                    Kprime_0=_number(eos[2]),
                    Kdprime_0=_number(eos[3], "1e-8"),
                ),
                property_modifiers=property_modifiers,
            )
        )
    return result


def read_thermocalc_covariance(path=ENDMEMBER_FILE):
    """Return all raw labels and their formation-energy covariance [J^2/mol^2].

    Format 5 stores a scalar followed by the packed upper triangle, in
    endmember-record order. No covariance for EOS or mixing parameters is
    supplied in this block. Preserve all source labels, including endmembers
    outside the HHPH2013 subset, and convert (kJ/mol)² to (J/mol)².
    """
    records, lines = _read_thermocalc_file(path)
    names = list(records)
    values = []
    for line in lines:
        if line.lstrip().startswith("tc-ds"):
            break
        values.extend(line.split())
    count = len(names)
    expected = 1 + count * (count + 1) // 2
    if len(values) != expected:
        raise ValueError(
            f"{path}: invalid triangular covariance size: "
            f"expected {expected} values, found {len(values)}"
        )
    upper = [_number(value, "1e6") for value in values[1:]]
    indices = np.triu_indices(count)
    matrix = np.zeros((count, count))
    matrix[indices] = upper
    matrix[(indices[1], indices[0])] = upper
    return {"endmember_names": names, "covariance_matrix": matrix}


def convert_covariance_dataset(endmember_file=ENDMEMBER_FILE):
    """Generate a separate, lazily imported covariance module."""
    covariance = read_thermocalc_covariance(endmember_file)
    names = covariance["endmember_names"]
    upper = covariance["covariance_matrix"][np.triu_indices(len(names))]
    return (
        "# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit\n"
        "# for the Earth and Planetary Sciences\n"
        "# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU\n"
        "# GPL v2 or later.\n\n"
        '"""HHPH2013 reference formation-energy covariance [J^2/mol^2].\n\n'
        "Contains all labels from the authors' archived tc-ds62.txt.\n"
        'Autogenerated using HHPH2013data_to_burnman.py.\n"""\n\n'
        "import numpy as np\n\n"
        f"endmember_names = {names!r}\n"
        f"upper = np.array({upper.tolist()!r})\n"
        f"matrix = np.zeros(({len(names)}, {len(names)}))\n"
        "indices = np.triu_indices(len(endmember_names))\n"
        "matrix[indices] = upper\n"
        "matrix[(indices[1], indices[0])] = upper\n"
        "cov = {'endmember_names': endmember_names, 'covariance_matrix': matrix}\n"
        "del upper, indices, matrix\n"
    )


def _endmember_source(records):
    source = ""
    for record in records:
        name = record["name"]
        source += f"\n\nclass {name}(Mineral):\n    def __init__(self):\n"
        if isinstance(record["formula"], str):
            source += f"        formula = dictionarize_formula({record['formula']!r})\n"
        else:
            source += f"        formula = {record['formula']!r}\n"
        source += "        self.params = {\n"
        source += f"            'name': {name!r},\n            'formula': formula,\n"
        source += "            'equation_of_state': 'hp_tmt',\n"
        for key, value in record["params"].items():
            source += f"            {key!r}: {value!r},\n"
        source += "            'n': sum(formula.values()),\n"
        source += "            'molar_mass': formula_mass(formula),\n        }\n"
        if record["property_modifiers"]:
            source += (
                "        self.property_modifiers = "
                + pprint.pformat(record["property_modifiers"], sort_dicts=False)
                + "\n"
            )
        source += "        Mineral.__init__(self)\n"
    return source


def _preamble(title, description):
    return (
        "# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit\n"
        "# for the Earth and Planetary Sciences\n"
        "# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU\n"
        "# GPL v2 or later.\n\n"
        f'"""\n{title}\n{"^" * len(title)}\n\n{description}\n\n'
        "All generated parameters use SI units.\n"
        'Autogenerated using HHPH2013data_to_burnman.py.\n"""\n\n'
        "import numpy as np\nfrom numpy import array, nan\n"
        "from ..classes.mineral import Mineral\n"
        "from ..classes.solution import Solution\n"
        "from ..classes.solutionmodel import SymmetricRegularSolution, AsymmetricRegularSolution\n"
        "from ..classes.combinedmineral import CombinedMineral\n"
        "from ..utils.chemistry import dictionarize_formula, formula_mass\n"
    )


def convert_dataset(endmember_file=ENDMEMBER_FILE, solution_file=SOLUTION_FILE):
    """Generate endmembers and solutions from both original THERMOCALC files."""
    from burnman.data.input_raw_endmember_datasets.hpx_eos_to_burnman import (
        convert_dataset as convert_solutions,
    )

    module_name = "HHPH_2013"
    source = _preamble(
        module_name,
        "Upper mantle, transition zone and shallow-lower mantle dataset from\n"
        "Holland, Powell & Harte (2013; doi:10.1093/petrology/egt035).\n\n"
        "The parameters in this dataset are derived from the authors' archived THERMOCALC\n"
        f"archive ({ARCHIVE_URL}),\n"
        "which supersede the printed Table 1 values.\n"
        "The classes mpv, cpv and cor are pure endmembers; their corresponding solutions\n"
        "are mpv_solution, cpv_solution and cor_solution.\n\n"
        "No composition limits are imposed by this dataset, but users should be aware that\n"
        "the hpx and opx solutions need to be limited to Ca-poor compositions\n"
        "to match the published papers (Ca(M2) < 0.2 is a reasonable bound).",
    ) + _endmember_source(read_thermocalc_endmembers(endmember_file))
    # Use freshly converted data; never import a stale generated dataset.
    module = ModuleType(f"burnman.minerals.{module_name}")
    module.__package__ = "burnman.minerals"
    exec(compile(source, str(endmember_file), "exec"), module.__dict__)
    solutions = convert_solutions(
        solution_file, module_name, endmember_module=module, inline=True
    )
    source += (
        "\n\n# Solid solutions, read directly from the THERMOCALC activity models.\n"
    )
    source += solutions
    classes = []
    for statement in ast.parse(solutions).body:
        if isinstance(statement, ast.ClassDef):
            assignment = next(
                node
                for node in ast.walk(statement)
                if isinstance(node, ast.Assign)
                and any(
                    isinstance(t, ast.Attribute) and t.attr == "name"
                    for t in node.targets
                )
            )
            classes.append((ast.literal_eval(assignment.value), statement.name))
    source += "\n\nSOLUTION_CLASSES = {\n"
    source += "".join(f"    {label!r}: {name},\n" for label, name in classes)
    return source + '''}


def cov():
    """Return all 273 raw endmember labels and formation-energy covariance.

    The dictionary keys are ``endmember_names`` and ``covariance_matrix``.
    Names follow the original tc-ds62.txt record order, including endmembers
    outside the HHPH2013 subset. Matrix entries use (J/mol)².

    Covariance for CombinedMineral recipes can be propagated with A C A.T,
    as in BurnMan's thermobarometry tools.
    """
    from .HHPH_2013_cov import cov

    return cov
'''


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endmembers", type=Path, help="Input THERMOCALC ds62 file.")
    parser.add_argument(
        "--solutions", type=Path, help="Input THERMOCALC activity-model file."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output module; also writes HHPH_2013_cov.py beside it. Defaults to stdout.",
    )
    args = parser.parse_args(argv)
    source = convert_dataset(
        args.endmembers or ENDMEMBER_FILE, args.solutions or SOLUTION_FILE
    )
    if args.output is None:
        sys.stdout.write(source)
    else:
        covariance_source = convert_covariance_dataset(
            args.endmembers or ENDMEMBER_FILE
        )
        args.output.write_text(source, encoding="utf-8")
        args.output.with_name("HHPH_2013_cov.py").write_text(
            covariance_source, encoding="utf-8"
        )


if __name__ == "__main__":
    sys.path.insert(0, str(SOURCE_DIR.parents[2]))
    main()
