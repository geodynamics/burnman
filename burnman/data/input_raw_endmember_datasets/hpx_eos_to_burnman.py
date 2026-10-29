# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


# This is a script that converts a HPx-eos data format into the standard burnman
# format (written to burnman/minerals/). For example, the dataset file
# tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASHTOCr.txt
# can be processed via the command
# python hpx_eos_to_burnman.py
# the text "from . import ig50NCKFMASHTOCr must then be added to the minerals __init__.py

import numpy as np
import logging
import sys
from sympy import Symbol, prod, sympify, cancel, linsolve, Poly
from sympy.parsing.sympy_parser import parse_expr
from burnman.constants import gas_constant
from datetime import date
import importlib
import keyword
import re
import subprocess
from pathlib import Path

ds = [
    ["tc-thermoinput-igW24-2025-06-28/tc-ig51W24.txt", "HPx_ds636"],
    ["tc-thermoinput-igG25-H18corr-2024-12-21/tc-ig51G25.txt", "HPx_ds636"],
    ["tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASHTOCr.txt", "HGP_2018_ds633"],
    ["tc-thermoinput-igneous-2022-01-23/tc-ig50NCKFMASTOCr.txt", "HGP_2018_ds633"],
    ["tc-thermoinput-metabasite-2022-01-30/tc-mb50NCKFMASHTO.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50KFMASH.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50MnNCKFMASHTO.txt", "HP_2011_ds62"],
    ["tc-thermoinput-metapelite-2022-01-23/tc-mp50NCKFMASHTO.txt", "HP_2011_ds62"],
]


def dataset_name(path):
    return Path(path).stem.removeprefix("tc-")


# The original H18 melt model was withdrawn following an algebraic error.
# Its corrected G25 replacement uses dataset 6.36, rather than 6.33.
# https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/hpx-eos-igneous-sets/
unsupported_solutions = {
    "ig50NCKFMASHTOCr": {"liq"},
    "ig50NCKFMASTOCr": {"liq"},
}

dataset_sources = {
    "ig51W24": (
        "Weller et al. (2024), doi:10.1093/petrology/egae098.\n"
        "Source: tc-ig51W24.txt, W24 release of 2025-06-28.\n"
        "https://hpxeosandthermocalc.org/wp-content/uploads/2025/06/"
        "tc-thermoinput-igw24-2025-06-28.zip\n"
    ),
    "ig51G25": (
        "Green et al. (2025), doi:10.1093/petrology/egae079.\n"
        "Source: tc-ig51G25.txt, G25 release uploaded 2025-05-07.\n"
        "https://hpxeosandthermocalc.org/wp-content/uploads/2025/05/"
        "tc-thermoinput-igg25-h18corr-2024-12-21.zip\n"
    ),
}

# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
current_year = str(date.today())[:4]


def reverse_polish(lines):
    # get symbols
    syms = []
    for line in lines:
        n_terms = int(line[0])
        i = 2
        for term in range(n_terms):
            n_add = int(line[i])
            syms.extend(line[i + 2 : i + 2 * n_add + 2][::2])
            i += 2 * n_add + 2
    unique_syms = list(np.unique(sorted(syms)))
    sympy_syms = {u: Symbol(u) for u in unique_syms}

    # get expressions
    expr = sympify(0)
    for line in lines:
        n_terms = int(line[0])
        i = 2
        terms = []
        for term in range(n_terms):
            n_add = int(line[i])
            en = parse_expr(line[i - 1])  # this is the constant
            for j in range(n_add):
                en += parse_expr(line[i + 1 + 2 * j]) * sympy_syms[line[i + 2 + 2 * j]]
            terms.append(en)
            # syms.extend(l[i+2:i+2*n_add + 2][::2])
            # print("a", n_add)
            i += 2 * n_add + 2
        product = prod(terms)
        expr += product
    return expr


def ordering_modifier(property_modifier, ordered):
    mtype, params = property_modifier
    if mtype == "bragg_williams":
        if ordered:
            return np.array([0.0, 0.0, 0.0])
        else:
            n = params["n"]
            if params["factor"] > 0.0:
                f = [params["factor"], params["factor"]]
            else:
                f = [1.0, -params["factor"]]

            S_disord = (
                -gas_constant
                * (
                    f[0] * (np.log(1.0 / (n + 1.0)) + n * np.log(n / (n + 1.0)))
                    + f[1]
                    * (n * np.log(1.0 / (n + 1.0)) + n * n * np.log(n / (n + 1.0)))
                )
                / (n + 1.0)
            )
            return np.array([params["deltaH"], S_disord, params["deltaV"]])
    else:
        if params["T_0"] < params["Tc_0"]:
            Q_0 = np.power((params["Tc_0"] - params["T_0"]) / params["Tc_0"], 0.25)
        else:
            Q_0 = 0.0

        if ordered:
            E_ord = (
                params["Tc_0"] * params["S_D"] * (Q_0 * Q_0 - np.power(Q_0, 6.0) / 3.0)
                - (
                    params["S_D"] * params["Tc_0"] * 2.0 / 3.0
                    - params["V_D"] * params["P_0"]
                )
                - params["P_0"] * params["V_D"] * Q_0 * Q_0
            )
            S_ord = params["S_D"] * (Q_0 * Q_0 - 1.0)
            V_ord = params["V_D"] * (Q_0 * Q_0 - 1.0)
            return np.array([E_ord, S_ord, V_ord])

        else:
            E_disord = (
                params["Tc_0"] * params["S_D"] * (Q_0 * Q_0 - np.power(Q_0, 6.0) / 3.0)
                - params["P_0"] * params["V_D"] * Q_0 * Q_0
            )
            S_disord = params["S_D"] * Q_0 * Q_0
            V_disord = params["V_D"] * Q_0 * Q_0
            return np.array([E_disord, S_disord, V_disord])


digit_map = {
    "0": "zero",
    "1": "one",
    "2": "two",
    "3": "three",
    "4": "four",
    "5": "five",
    "6": "six",
    "7": "seven",
    "8": "eight",
    "9": "nine",
}


def replace_numbers(match):
    num_str = match.group()
    return "".join(digit_map[d] for d in num_str)


def temkin_site_formulae(mbr_proportions, site_expressions, endmember_rows):
    """Recover variable site multiplicities from the ideal activity expressions.

    Each activity factor is a ratio of linear combinations of endmember
    amounts. Factors with proportional denominators occupy the same site.
    Their activity exponents give the species counts on that site, and their
    sum gives its multiplicity, which can differ between endmembers.
    """
    proportions = []
    for lines in mbr_proportions:
        expression = [line[:] for line in lines]
        expression[0] = expression[0][2:]
        proportions.append(reverse_polish(expression))

    variables = sorted(set().union(*(p.free_symbols for p in proportions)), key=str)
    coefficients = [Symbol(f"_coefficient_{i}") for i in range(len(proportions))]
    linear_combination = sum(c * p for c, p in zip(coefficients, proportions))
    exponents = {species: [] for species in site_expressions}
    for row in endmember_rows:
        powers = {
            row[3 + 2 * i]: parse_expr(row[4 + 2 * i]) for i in range(int(row[2]))
        }
        for species in exponents:
            exponents[species].append(powers.get(species, sympify(0)))

    # Modern THERMOCALC files give species amounts and explicit denominator
    # factors (e.g. mgM**4 / sumM**4). Turn those pairs into site fractions.
    # Start with the largest powers to distinguish microscopic counts from
    # molecular proportions, which may share the same expression.
    expressions = dict(site_expressions)
    denominators = [s for s, counts in exponents.items() if any(c < 0 for c in counts)]
    denominators.sort(key=lambda s: min(exponents[s]))
    paired_species = set()
    for denominator in denominators:
        remaining = [-c for c in exponents[denominator]]
        if any(c < 0 for c in remaining):
            raise ValueError(f"Mixed signs in activity factor '{denominator}'")
        for species, counts in exponents.items():
            if species in denominators or species in paired_species or not any(counts):
                continue
            if all(c == 0 or c == r for c, r in zip(counts, remaining)):
                expressions[species] = cancel(
                    expressions[species] / site_expressions[denominator]
                )
                remaining = [r - c for r, c in zip(remaining, counts)]
                paired_species.add(species)
        if any(remaining):
            raise ValueError(f"Cannot pair the activity denominator '{denominator}'")

    sites = {}
    for species, expression in expressions.items():
        counts = exponents[species]
        if species in denominators:
            continue
        if not any(counts):
            continue
        if any(count < 0 for count in counts):
            raise ValueError(f"Negative activity exponent for '{species}'")
        numerator = sum(c * p for c, p in zip(counts, proportions))
        denominator = cancel(numerator / expression)
        equations = Poly(linear_combination - denominator, *variables).coeffs()
        solutions = list(linsolve(equations, coefficients))
        if len(solutions) != 1 or any(c.free_symbols for c in solutions[0]):
            raise ValueError(f"Cannot recover the mixing site for '{species}'")
        denominator_counts = solutions[0]
        if any(c < 0 for c in denominator_counts):
            raise ValueError(f"Negative site multiplicity for '{species}'")
        scale = next(c for c in denominator_counts if c != 0)
        site = tuple(c / scale for c in denominator_counts)
        sites.setdefault(site, []).append((species, counts))

    formulae = ["" for _ in proportions]
    for denominator_counts, species in sites.items():
        multiplicities = [
            sum(counts[i] for _, counts in species) for i in range(len(proportions))
        ]
        scale = next(
            m / d for m, d in zip(multiplicities, denominator_counts) if d != 0
        )
        if any(m != scale * d for m, d in zip(multiplicities, denominator_counts)):
            raise ValueError("Activity factors have incompatible site multiplicities")
        for i, multiplicity in enumerate(multiplicities):
            if multiplicity == 0:
                formulae[i] += "[]0"
                continue
            occupancies = []
            for label, counts in species:
                if counts[i] == 0:
                    continue
                label = re.sub(r"\d+", replace_numbers, label)
                label = "".join(filter(str.isalpha, label)).title()
                fraction = counts[i] / multiplicity
                occupancies.append(label + (str(fraction) if fraction != 1 else ""))
            formulae[i] += "[" + "".join(occupancies) + "]"
            if multiplicity != 1:
                formulae[i] += str(multiplicity)
    return formulae


class UniqueNames:
    """Allocate module identifiers without overwriting imports or other objects."""

    def __init__(self, reserved_names, source_names):
        self.used_names = set(reserved_names) | set(keyword.kwlist)
        self.source_names = set(source_names)

    def allocate(self, name, context):
        identifier = name
        if identifier in self.used_names:
            stem = f"{name}_{context}"
            identifier = stem
            suffix = 2
            while identifier in self.used_names or identifier in self.source_names:
                identifier = f"{stem}_{suffix}"
                suffix += 1
        self.used_names.add(identifier)
        return identifier


def convert_dataset(solution_file, mbr_dataset):
    """Return the Python source for one HPx-eos solution dataset."""
    logging.debug(f"{solution_file}, {mbr_dataset}")

    dataset = importlib.import_module(f"burnman.minerals.{mbr_dataset}")
    extra_datasets = set()

    def endmember_reference(label):
        if label == "H2O" and not hasattr(dataset, label):
            extra_datasets.add("HP_2011_fluids")
            return "HP_2011_fluids.H2O()"
        if label == "and":
            label = "andalusite"
        return f"{mbr_dataset}.{label}()"

    out_ss = ""
    out_make = ""
    solution_dataset = dataset_name(solution_file)
    excluded = unsupported_solutions.get(solution_dataset, set())

    i = 0
    data = []
    with open(solution_file, "r") as file:
        # Read each line in the file one by one
        process = True
        for line in file:
            # Process the line (for example, print it)
            ln = line.split("%", 1)[0].strip().split()
            if len(ln) > 0:
                if ln[0] == "verbatim" or ln[0] == "header":
                    process = not process
                if (
                    ln[0][0] != "%"
                    and ln[0][0] != "*"
                    and ln[0] != "header"
                    and ln[0] != "verbatim"
                    and process
                ):
                    data.append(ln)

        start_indices = []
        for i_d, d in enumerate(data):
            legacy_header = len(d) == 3 and d[1].isdigit() and d[2] in ("1", "2")
            modern_header = (
                len(d) >= 5
                and d[1].isdigit()
                and d[2] == "2"
                and re.fullmatch(r"[A-Za-z][\w-]*", d[3]) is not None
                and re.fullmatch(r"\w+", d[0]) is not None
            )
            if (
                (legacy_header or modern_header)
                and int(d[1]) > 1
                and d[0] not in excluded
            ):
                start_indices.append(i_d)

    identifiers = UniqueNames(
        [
            "np",
            "array",
            "nan",
            mbr_dataset,
            "HP_2011_fluids",
            "Mineral",
            "Solution",
            "SymmetricRegularSolution",
            "AsymmetricRegularSolution",
            "CombinedMineral",
        ],
        [line[0] for line in data],
    )
    # Allocate classes first so that endmembers cannot shadow solution names.
    solution_identifiers = {
        ind: identifiers.allocate(data[ind][0], "solution") for ind in start_indices
    }
    noods = {}
    combined_endmembers = {}

    def resolve_endmember(rows, context):
        nonlocal out_make
        label = rows[0][0]
        commands = [row[0] for row in rows[1:]]
        if not any(c == "make" or c == "DQF" or "delG" in c for c in commands):
            return endmember_reference(label)
        if "make" in commands:
            make = next(row[1:] for row in rows if row[0] == "make")
        else:
            make = ["1", label, "1"]
        corrections = [
            np.array(row[1:4], dtype=float)
            for row in rows
            if "delG" in row[0] or row[0] == "DQF"
        ]
        delG = np.sum(corrections, axis=0) if corrections else np.zeros(3)
        delG *= [1.0e3, -1.0e3, 1.0e-5]
        combined_mbrs, combined_amounts = [], []
        el = 1
        for _ in range(int(make[0])):
            if make[el] in ("ordered", "disordered"):
                ordered = make[el] == "ordered"
                el += 1
                m_string = make[el]
                amount = float(parse_expr(make[el + 1]))
                mineral = getattr(dataset, m_string)()
                delG += (
                    ordering_modifier(mineral.property_modifiers[0], ordered) * amount
                )
                if m_string not in noods:
                    noods[m_string] = identifiers.allocate(
                        f"{m_string}_nood", "endmember"
                    )
                combined_mbrs.append(noods[m_string])
            else:
                if make[el] == "equilibrium":
                    el += 1
                combined_mbrs.append(endmember_reference(make[el]))
            combined_amounts.append(float(parse_expr(make[el + 1])))
            el += 2
        # Labels are local to each solution. Reuse identical definitions,
        # including reordered recipes, without overwriting distinct objects.
        recipe = tuple(sorted(zip(combined_mbrs, combined_amounts)))
        definition = (label, recipe, tuple(delG))
        if definition not in combined_endmembers:
            identifier = identifiers.allocate(label, context)
            combined_endmembers[definition] = identifier
            out_make += f'{identifier} = CombinedMineral([{", ".join(combined_mbrs)}], {combined_amounts}, {list(delG)}, "{label}")\n'
        return combined_endmembers[definition]

    names = []
    for ind in start_indices:
        i = ind
        name = data[i][0]
        names.append(name)
        n_mbrs = int(data[i][1])
        logging.debug(f"Solution: {name}")
        logging.debug(f"Number of endmembers: {n_mbrs}")

        vars = data[i + 1 : i + n_mbrs]
        logging.debug(vars)
        i += n_mbrs

        mbr_proportions = []
        mbr_names = []
        for j in range(n_mbrs):
            n_lines = int(data[i][1])
            mbr_names.append(data[i][0].split("(")[1][:-1])
            mbr_proportions.append(data[i : i + n_lines])
            i += n_lines
        logging.debug(mbr_proportions)

        if len(set(mbr_names)) != len(mbr_names):
            raise ValueError(f"Duplicate endmember names in solution '{name}'")

        formulation = data[i][0]
        logging.debug(f"formulation: {formulation}")
        i += 1
        n_ints = int((n_mbrs * (n_mbrs - 1.0)) / 2.0)
        interactions = data[i : i + n_ints]
        logging.debug(interactions)

        nonideal_entropies = False
        nonideal_volumes = False

        We = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        Ws = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        Wv = [[0.0 for j in mbr_names[:-i]] for i in range(1, n_mbrs)]
        for ints in interactions:
            m = ints[0].split("(")[1].replace(")", "").split(",")
            m0 = mbr_names.index(m[0])
            m1 = mbr_names.index(m[1])
            We[m0][m1 - m0 - 1] = float(ints[1]) * 1.0e3
            Ws[m0][m1 - m0 - 1] = -float(ints[2]) * 1.0e3
            Wv[m0][m1 - m0 - 1] = float(ints[3]) * 1.0e-5

            if np.abs(Ws[m0][m1 - m0 - 1]) > 1.0e-10:
                nonideal_entropies = True
            if np.abs(Wv[m0][m1 - m0 - 1]) > 1.0e-10:
                nonideal_volumes = True

        i += n_ints

        if formulation == "asf":
            text_alphas = data[i : i + n_mbrs]
            logging.debug(text_alphas)
            i += n_mbrs

            Ae = [0.0 for i in range(n_mbrs)]
            As = [0.0 for i in range(n_mbrs)]
            Av = [0.0 for i in range(n_mbrs)]
            for alpha in text_alphas:
                try:
                    m = alpha[0].split("(")[1].replace(")", "").split(",")[0]
                except IndexError:
                    m = alpha[0]
                m0 = mbr_names.index(m)
                Ae[m0] = float(alpha[1])
                As[m0] = -float(alpha[2])
                Av[m0] = float(alpha[3]) * 1.0e6

        n_site_species = int(data[i][0])

        logging.debug(f"Number of site species: {n_site_species}")
        i += 1

        site_fractions = []
        for j in range(n_site_species):
            n_lines = int(data[i][1])
            site_fractions.append(data[i : i + n_lines])
            i += n_lines

        summed = 0.0
        sites = [[]]
        site_expressions = {}
        for site_species_and_expression in site_fractions:
            site_species = site_species_and_expression[0][0]
            expression = [line[:] for line in site_species_and_expression]
            expression[0] = expression[0][2:]
            expression = reverse_polish(expression)
            site_expressions[site_species] = expression
            summed += expression
            summed = sympify(summed)
            sites[-1].append(site_species)
            if summed == 1:
                summed = 0
                sites.append([])
        sites = sites[:-1]
        site_species = {}
        for n_site, site in enumerate(sites):
            for species in site:
                site_species[species] = n_site
        n_sites = len(sites)

        endmember_site_fractions_and_checks = []
        for j in range(n_mbrs):
            n_lines = 1
            while i + n_lines < len(data):
                if (
                    data[i + n_lines][0] == "make"
                    or data[i + n_lines][0] == "delG(tran)"
                    or data[i + n_lines][0] == "check"
                    or data[i + n_lines][0] == "delG(od)"
                    or data[i + n_lines][0] == "delG(make)"
                    or data[i + n_lines][0] == "delG(mod)"
                    or data[i + n_lines][0] == "delG(rcal)"
                    or data[i + n_lines][0] == "DQF"
                ):
                    n_lines += 1
                else:
                    break
            endmember_site_fractions_and_checks.append(data[i : i + n_lines])
            i += n_lines

        temkin_formulae = None
        if any(
            f[0][3 + 2 * j] not in site_species or parse_expr(f[0][4 + 2 * j]) < 0
            for f in endmember_site_fractions_and_checks
            for j in range(int(f[0][2]))
        ):
            temkin_formulae = temkin_site_formulae(
                mbr_proportions,
                site_expressions,
                [f[0] for f in endmember_site_fractions_and_checks],
            )

        mbr_initializations = []
        for i_mbr, f in enumerate(endmember_site_fractions_and_checks):
            mbr = resolve_endmember(f, name)
            if temkin_formulae is not None:
                mbr_initializations.append(f'[{mbr}, "{temkin_formulae[i_mbr]}"]')
                continue
            n_occs = int(f[0][2])
            site_occ = [[] for ln in range(n_sites)]
            site_occ_list = f[0][3:]
            ss = [[] for idx in range(n_sites)]
            site_amounts = [[] for idx in range(n_sites)]
            for j in range(n_occs):
                site_sp = site_occ_list[2 * j]
                i_site = site_species[site_sp]
                # rename for insertion into string
                if site_sp[0] == "x":
                    site_sp = site_sp[1:]
                    site_sp = re.sub(r"\d+(\.\d+)?", replace_numbers, site_sp)
                    site_sp = "".join(filter(str.isalpha, site_sp)).title()

                ss[i_site].append(site_sp)
                site_amounts[i_site].append(parse_expr(site_occ_list[2 * j + 1]))
            total_site_amounts = [sum(amounts) for amounts in site_amounts]
            site_fractions = [
                [amount / total_site_amounts[j] for amount in amounts]
                for j, amounts in enumerate(site_amounts)
            ]

            site_formula = ""
            for m, site in enumerate(ss):
                site_formula += "["
                for n, species in enumerate(site):
                    site_formula += species
                    if site_fractions[m][n] != 1:
                        site_formula += str(site_fractions[m][n])
                site_formula += "]"
                if str(total_site_amounts[m]) != "1":
                    site_formula += str(total_site_amounts[m])

            mbr_initializations.append(f'[{mbr}, "{site_formula}"]')

        # Create docstring
        doc = '        """\n'
        doc += f"        Initialisation for a {name} solution object.\n"

        doc += "        Contains the following endmembers with associated site occupancies:\n"

        for mbr_initialization in mbr_initializations:
            doc += f"        * {mbr_initialization}\n"

        if formulation == "asf":
            doc += "\n        This is implemented as an asymmetric solution.\n"
        else:
            doc += "\n        This is implemented as a symmetric solution.\n"
        doc += '        """\n'

        out_ss += f"class {solution_identifiers[ind]}(Solution):\n"
        out_ss += "    def __init__(self, molar_fractions=None):\n"
        out_ss += doc
        out_ss += f'        self.name = "{name}"\n'
        if formulation == "asf":
            out_ss += "        self.solution_model = AsymmetricRegularSolution(\n"
        else:
            out_ss += "        self.solution_model = SymmetricRegularSolution(\n"
        out_ss += "            endmembers=[\n"
        for mbr_initialization in mbr_initializations:
            out_ss += f"                {mbr_initialization},\n"
        out_ss += "            ],\n"
        if formulation == "asf":
            out_ss += f"            alphas={Ae},\n"
        out_ss += f"            energy_interaction={We},\n"
        if nonideal_entropies:
            out_ss += f"            entropy_interaction={Ws},\n"
        if nonideal_volumes:
            out_ss += f"            volume_interaction={Wv},\n"

        out_ss += "        )\n"
        out_ss += "        Solution.__init__(self, molar_fractions=molar_fractions)\n\n"

        logging.debug(site_fractions)
        logging.debug(endmember_site_fractions_and_checks)
        logging.debug("")

    buffers = []
    for i, row in enumerate(data):
        if len(row) >= 3 and row[1:3] == ["1", "buffer"]:
            rows = [row]
            for command in data[i + 1 :]:
                if command[0] not in ("make", "DQF") and "delG" not in command[0]:
                    break
                rows.append(command)
            buffers.append(resolve_endmember(rows, "buffer"))

    out_noods = ""
    for nood in sorted(noods):
        mineral = getattr(dataset, nood)()
        params = mineral.params
        out_noods += f"{noods[nood]} = Mineral({params})\n\n"

    # Preamble
    underline = "^" * len(solution_dataset)
    preamble = (
        "# This file is part of BurnMan - a thermoelastic\n"
        "# and thermodynamic toolkit for the Earth and Planetary Sciences\n"
        f"# Copyright (C) 2012 - {current_year} by the BurnMan team, released under the GNU\n"
        "# GPL v2 or later.\n"
        "\n"
        '"""\n'
        f"{solution_dataset}\n"
        f"{underline}\n"
        "\n"
        "HPx-eos solutions using endmembers from\n"
        f"dataset {mbr_dataset}.\n\n"
        "Contains the following solutions:\n"
    )
    for name in names:
        preamble += f"* {name}\n"

    if buffers:
        preamble += "\nOxygen buffers: " + ", ".join(buffers) + ".\n"

    if solution_dataset in dataset_sources:
        preamble += "\n" + dataset_sources[solution_dataset]

    if excluded:
        preamble += (
            "\nThe original H18 liq model is omitted because its authors withdrew it\n"
            "after discovering an algebraic error. The corrected G25 replacement\n"
            "requires endmember dataset 6.36. See:\n"
            "https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/hpx-eos-igneous-sets/\n"
        )

    preamble += (
        "\nThe values in this document are all in S.I. units,\n"
        "unlike those in the original THERMOCALC file.\n\n"
        "This file is autogenerated using hpx_eos_to_burnman.py\n"
        '"""\n\n'
        f"import numpy as np\n"
        f"from numpy import array, nan\n"
        f"from . import {mbr_dataset}\n"
        + "".join(f"from . import {d}\n" for d in sorted(extra_datasets))
        + "from ..classes.mineral import Mineral\n"
        "from ..classes.solution import Solution\n"
        "from ..classes.solutionmodel import SymmetricRegularSolution\n"
        "from ..classes.solutionmodel import AsymmetricRegularSolution\n"
        "from ..classes.combinedmineral import CombinedMineral\n\n"
    )

    output_string = f"{preamble}\n{out_noods}\n{out_make}\n{out_ss}"

    return output_string


def main():
    source_directory = Path(__file__).resolve().parent
    output_directory = source_directory.parent.parent / "minerals"
    output_files = []
    for solution_file, mbr_dataset in ds:
        output_string = convert_dataset(source_directory / solution_file, mbr_dataset)
        solution_dataset = dataset_name(solution_file)
        output_file = output_directory / f"{solution_dataset}.py"
        output_file.write_text(output_string)
        output_files.append(str(output_file))
        print(f"from . import {solution_dataset}")

    print("Copy the above to ../../minerals/__init__.py.\n")
    subprocess.run([sys.executable, "-m", "black", *output_files], check=True)


if __name__ == "__main__":
    main()
