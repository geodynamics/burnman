# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""Convert the THERMOCALC 3.51 dataset format supplied with W24 and G25.

Usage: python HP636data_to_burnman.py path/to/tc-ds636.txt
The output contains solids, liquids, PS1994 water and the complete covariance
matrix. Other gases and aqueous species require EOSs outside this conversion.
"""

import argparse
from pathlib import Path
from pprint import pformat
import subprocess
import sys

import numpy as np

from burnman.utils.chemistry import formula_mass

COMPONENTS = [
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
]


def read_dataset(path):
    """Return supported endmember records and the covariance in SI units."""
    lines = [
        line.split() for line in Path(path).read_text().splitlines() if line.strip()
    ]
    if lines[0][0] != "6":
        raise ValueError("Expected THERMOCALC dataset format 6.")
    count = int(lines[0][1])
    records, names = [], []
    for i in range(count):
        composition, reference, cp, elastic = lines[1 + 4 * i : 5 + 4 * i]
        name, long_name = composition[:2]
        names.append(name)
        flag = int(reference[0])
        if flag not in (0, 1, 2, 7) and name != "H2O":
            continue
        formula = {}
        j = 3
        while composition[j] != "0":
            element = COMPONENTS[int(composition[j]) - 1]
            formula[element] = float(composition[j + 1])
            j += 2
        params = {
            "name": name,
            "formula": formula,
            "equation_of_state": "hp_tmtL" if flag == 7 else "hp_tmt",
            "H_0": float(reference[1]) * 1.0e3,
            "S_0": float(reference[2]) * 1.0e3,
            "V_0": float(reference[3]) * 1.0e-5,
            "Cp": [float(c) * 1.0e3 for c in cp],
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }
        modifiers = []
        if name != "H2O":
            expansion, bulk, derivative = map(float, elastic[:3])
            params.update(
                {
                    "a_0": expansion,
                    "K_0": bulk * 1.0e8,
                    "Kprime_0": derivative,
                    # Format 6 omits K''. The HP approximation is K'' = -K'/K.
                    "Kdprime_0": -derivative / (bulk * 1.0e8),
                }
            )
            if flag == 7:
                params["dKdT_0"] = float(elastic[3]) * 1.0e8
            elif flag == 1:
                tc, entropy, volume = map(float, elastic[3:6])
                modifiers = [
                    [
                        "landau_hp",
                        {
                            "P_0": 1.0e5,
                            "T_0": 298.15,
                            "Tc_0": tc,
                            "S_D": entropy * 1.0e3,
                            "V_D": volume * 1.0e-5,
                        },
                    ]
                ]
            elif flag == 2:
                dh, dv, wh, wv, n, factor = map(float, elastic[3:9])
                modifiers = [
                    [
                        "bragg_williams",
                        {
                            "deltaH": dh * 1.0e3,
                            "deltaV": dv * 1.0e-5,
                            "Wh": wh * 1.0e3,
                            "Wv": wv * 1.0e-5,
                            "n": n,
                            "factor": factor,
                        },
                    ]
                ]
        records.append((name, params, modifiers))
    if len(names) != len(set(names)):
        raise ValueError("Duplicate endmember identifiers in dataset 6.36.")
    values = []
    for row in lines[1 + 4 * count :]:
        if row[0].startswith("tc-ds"):
            break
        values.extend(map(float, row))
    # THERMOCALC precedes the packed upper triangle with a scalar, as in
    # the older dataset formats handled by HPdata_to_burnman.py.
    if len(values) != 1 + count * (count + 1) // 2:
        raise ValueError("Invalid triangular covariance matrix size.")
    return records, names, np.asarray(values[1:]) * 1.0e6


def convert_dataset(path):
    """Return source for the endmember and covariance modules."""
    records, names, covariance = read_dataset(path)
    preamble = (
        "# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit\n"
        "# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU\n"
        "# GPL v2 or later.\n\n"
    )
    source = preamble + '''"""
HPx_ds636
^^^^^^^^^

THERMOCALC endmember dataset 6.36, distributed with Weller et al. (2024,
doi:10.1093/petrology/egae098) and Green et al. (2025,
doi:10.1093/petrology/egae079). Values are in SI units.
The solids and liquids use the Holland-Powell (2011) equations of state;
the parameter values come from dataset 6.36.
The supplied tc-ds636.txt records generation on 3 May 2023.
Contains solids, liquids and PS1994 water. Other gases and aqueous species
are omitted; the covariance retains all 291 original endmember labels.
Autogenerated by HP636data_to_burnman.py from tc-ds636.txt.
"""

from ..classes.mineral import Mineral
from . import HP_2011_fluids

'''
    for name, params, modifiers in records:
        if name == "H2O":
            reference = {key: params[key] for key in ("H_0", "S_0", "Cp")}
            source += (
                "class H2O(HP_2011_fluids.H2O):\n"
                "    def __init__(self):\n"
                "        super().__init__()\n"
                f"        self.params.update({reference!r})\n"
                "        self.set_method('pitzer-sterner')\n\n"
            )
            continue
        identifier = "andalusite" if name == "and" else name
        source += f"class {identifier}(Mineral):\n    def __init__(self):\n"
        source += "        self.params = " + pformat(params, sort_dicts=False) + "\n"
        if modifiers:
            source += (
                "        self.property_modifiers = "
                + pformat(modifiers, sort_dicts=False)
                + "\n"
            )
        source += "        Mineral.__init__(self)\n\n"
    source += '''def cov():
    """Return all source labels and their formation-energy covariance [J²/mol²]."""
    from .HPx_ds636_cov import cov
    return cov
'''
    covariance_source = preamble + (
        '"""Dataset 6.36 formation-energy covariance, generated by HP636data_to_burnman.py."""\n'
        "import numpy as np\n\n"
        f"endmember_names = {names!r}\n"
        f"upper = np.array({covariance.tolist()!r})\n"
        f"matrix = np.zeros(({len(names)}, {len(names)}))\n"
        "indices = np.triu_indices(len(endmember_names))\n"
        "matrix[indices] = upper\n"
        "matrix[(indices[1], indices[0])] = upper\n"
        "cov = {'endmember_names': endmember_names, 'covariance_matrix': matrix}\n"
        "del upper, indices, matrix\n"
    )
    return source, covariance_source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    args = parser.parse_args()
    output = Path(__file__).resolve().parents[2] / "minerals"
    files = [output / "HPx_ds636.py", output / "HPx_ds636_cov.py"]
    for path, source in zip(files, convert_dataset(args.dataset)):
        path.write_text(source)
    subprocess.run(
        [sys.executable, "-m", "black", *(str(path) for path in files)], check=True
    )


if __name__ == "__main__":
    main()
