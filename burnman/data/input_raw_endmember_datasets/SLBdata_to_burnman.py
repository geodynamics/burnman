# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This is a standalone program that converts a tabulated version of the
# Stixrude and Lithgow-Bertelloni data format into the standard burnman
# format (printed to stdout)


import os.path

if os.path.isfile("stx11ver.dat") == False:
    print("This code requires the data file stx11ver.dat.")
    print("This file is bundled with the Perple_X data files, which can be found here:")
    print("http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/")
    print("")
    print("Please download the file and place it in this directory.")
    exit()


def read_dataset(datafile):
    # Perple_X distributes stx11ver.dat with Windows-1252 punctuation.
    with open(datafile, "r", encoding="cp1252") as f:
        return [line.split() for line in f]


def read_dataset_utf8(datafile):
    with open(datafile, "r", encoding="utf-8") as f:
        return [line.split() for line in f]


def process_stixrude():
    perplex_ds = read_dataset("stx11ver.dat")
    alias_mapping = read_dataset("perplex_slb_2011_names.txt")

    alias = {}
    accurate_values = {}
    configurational_entropies = {}
    for line in alias_mapping:
        alias[line[0]] = line[2]

    process = 0  # flag to process data or not (to cut out first lines)
    n = 0  # line number for each record
    for line in perplex_ds:
        if len(line) != 0:
            if process == 1:
                if line[0] == "end":
                    n = 0
                else:
                    n = n + 1
                if n == 1:
                    name = alias[line[0]]
                    data = []
                    configurational_entropy = "None"
                elif n == 3:
                    # F0, (n, ignored), -V0
                    data.extend([float(line[2]), -1.0 * float(line[8])])
                elif n == 4:
                    # K_0, K_prime, Debye_0, grueneisen_0, q_0, eta_s0,
                    # [configurational_and_magnetic_entropy]
                    data.extend(
                        [
                            float(line[2]),
                            float(line[5]),
                            float(line[8]),
                            float(line[11]),
                            float(line[14]),
                            float(line[17]),
                        ]
                    )
                    if len(line) == 21:
                        configurational_entropy = float(line[20])
                elif n == 5:
                    # mu_S0, mu_S0_prime
                    data.extend([float(line[2]), float(line[5])])
                    # Sort into correct order for processing: F, V, K, K',
                    #                        Debye, gruen, q, G, Gprime, etaS0
                    accurate_values[name] = [
                        data[0] / 1.0e3,
                        data[1] * 10.0,
                        data[2] / 1.0e4,
                        data[3],
                        data[4],
                        data[5],
                        data[6],
                        data[8] / 1.0e4,
                        data[9],
                        data[7],
                    ]
                    configurational_entropies[name] = configurational_entropy

            if line[0] == "end_components":
                process = 1

    return accurate_values, configurational_entropies


accurate_values, configurational_entropies = process_stixrude()

# Write the fixed site entropies as occupancy expressions. Perple_X supplies
# rounded entropy offsets for spinel, hercynite and Ca-Tschermaks; the
# jadeite-majorite offset is absent there. The hercynite expression also
# retains the magnetic entropy included in the Perple_X offset.
configurational_entropies.update(
    {
        # [Mg3/4Al1/4]4[Al7/8Mg1/8]8O16
        "spinel": (
            "-_R * (3.0 * np.log(3.0 / 4.0) + 1.0 * np.log(1.0 / 4.0) "
            "+ 7.0 * np.log(7.0 / 8.0) + 1.0 * np.log(1.0 / 8.0))"
        ),
        # Same site occupancies as spinel ([Fe3/4Al1/4]4[Al7/8Fe1/8]8O16),
        # plus the Perple_X Fe2+ magnetic entropy: 4 high-spin Fe2+ ions
        # (S=2, 2S+1=5 spin states) per formula unit (Fe4Al8O16).
        "hercynite": (
            "-_R * (3.0 * np.log(3.0 / 4.0) + 1.0 * np.log(1.0 / 4.0) "
            "+ 7.0 * np.log(7.0 / 8.0) + 1.0 * np.log(1.0 / 8.0)) "
            "+ 4.0 * _R * np.log(5.0)"
        ),
        # [Ca][Al][Si1/2Al1/2]2O6
        "ca_tschermaks": "-_R * (1.0 * np.log(1.0 / 2.0) + 1.0 * np.log(1.0 / 2.0))",
        # [Na2/3Al1/3]3[Al][Si]Si3O12
        "jd_majorite": "-_R * (2.0 * np.log(2.0 / 3.0) + 1.0 * np.log(1.0 / 3.0))",
    }
)


ds = read_dataset_utf8("slb_2011.txt")
landau = read_dataset_utf8("slb_2011_landau.txt")
landau_params = {}
for line in landau:
    phase = str(line[0]).lower()
    if phase != "#":
        landau_params[phase] = [float(line[1]), float(line[2]) * 1.0e-6, float(line[3])]

print("# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit")
print("# for the Earth and Planetary Sciences")
print("# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU")
print("# GPL v2 or later.")
print("")
print("")
print('"""')
print("SLB_2011")
print("^^^^^^^^")
print("")
print("Minerals from Stixrude & Lithgow-Bertelloni 2011 and references therein.")
print("File autogenerated using SLBdata_to_burnman.py.")
print('"""')
print("")
print("")
print("")
print("import numpy as np")
print("")
print("from ..classes.mineral import Mineral")
print("from ..classes.solution import Solution")
print("from ..classes.solutionmodel import (")
print("    IdealSolution,")
print("    SymmetricRegularSolution,")
print("    AsymmetricRegularSolution,")
print(")")
print("from ..utils.chemistry import dictionarize_formula, formula_mass")
print("")
print("# Gas constant (J/mol/K) used by HeFESTo (CODATA 2006) to compute the fixed")
print("# mixed-site configurational entropies below; differs negligibly from the")
print("# current exact SI value in burnman.constants.gas_constant.")
print("_R = 8.314472")
print("")

param_scales = [
    -1.0,
    -1.0,  # not numbers, so we won't scale
    1.0e3,
    1.0e3,  # KJ -> J
    1.0e-6,
    1.0e-6,  # cm^3/mol -> m^3/mol
    1.0e9,
    1.0e9,  # GPa -> Pa
    1.0,
    1.0,  # no scale for K'
    1.0,
    1.0,  # no scale for Debye
    1.0,
    1.0,  # no scale for gruneisen
    1.0,
    1.0,  # no scale for q
    1.0e9,
    1.0e9,  # GPa -> Pa
    1.0,
    1.0,  # no scale for G'
    1.0,
    1.0,
]  # no scale for eta_s


solutionfile = "slb_2011_solutions.txt"
with open(solutionfile, "r") as fin:
    print(fin.read())
fin.close()

print('"""')
print("ENDMEMBERS")
print('"""')
print("")
formula = "0"
for idx, m in enumerate(ds):
    if idx == 0:
        param_names = m
    else:
        print("class", m[0].lower(), "(Mineral):")
        print("    def __init__(self):")
        print("".join(["        formula='", m[1], "'"]))
        print("        formula = dictionarize_formula(formula)")
        print("        self.params = {")
        print("".join(["            'name': '", m[0], "',"]))
        print("            'formula': formula,")
        print("            'equation_of_state': 'slb3',")
        for pid, param in enumerate(accurate_values[m[0].lower()]):
            print(
                "            '" + param_names[(pid + 1) * 2] + "':",
                float(param) * param_scales[(pid + 1) * 2],
                ",",
            )
        print("            'n': sum(formula.values()),")
        print("            'molar_mass': formula_mass(formula)}")
        print("")
        # The SLB2011 HeFESTo Landau term uses the disordered reference state.
        # The generic BurnMan modifier uses the ordered state; an extra linear
        # term cancels that reference offset for the published SLB2011 phases.
        # Stishovite's negative Tc_0 encodes the high-pressure CaCl2 transition.
        has_landau = m[0].lower() in landau_params
        if has_landau or configurational_entropies[m[0].lower()] != "None":
            print(
                "        self.property_modifiers = [",
            )
            if has_landau:
                print(
                    "['landau', {'Tc_0':",
                    landau_params[m[0].lower()][0],
                    ", 'S_D':",
                    landau_params[m[0].lower()][2],
                    ", 'V_D':",
                    landau_params[m[0].lower()][1],
                    "}],",
                )
                if m[0].lower() in ("quartz", "nepheline", "stishovite"):
                    tc0, vd, sd = landau_params[m[0].lower()]
                    print(
                        "['linear', {'delta_E':",
                        f"-2.0 * {sd} * {tc0} / 3.0",
                        ", 'delta_S':",
                        f"-{sd}",
                        ", 'delta_V':",
                        f"-{vd}",
                        "}],",
                    )
            if configurational_entropies[m[0].lower()] != "None":
                print(
                    "['linear', {'delta_E':",
                    0.0,
                    ", 'delta_S':",
                    configurational_entropies[m[0].lower()],
                    ", 'delta_V':",
                    0.0,
                    "}],",
                )
            print("]")
            print("")
        print("        self.uncertainties = {")
        for pid, param in enumerate(m):
            if pid > 1 and pid % 2 == 1 and pid < 21:
                print(
                    "            '" + param_names[pid] + "':",
                    float(param) * param_scales[pid],
                    ",",
                )
        pid = 21
        param = m[pid]
        print(
            "            '" + param_names[pid] + "':",
            float(param) * param_scales[pid],
            "}",
        )
        print("        Mineral.__init__(self)")
        print("")


aliasfile = "slb_2011_aliases.txt"
with open(aliasfile, "r") as fin:
    print(fin.read())
fin.close()
