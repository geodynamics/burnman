# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2023 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import os
import subprocess
import shutil
from copy import deepcopy

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata

from .material import Material, material_property

from ..utils.misc import copy_documentation
from ..utils.misc import extract_lines_between_markers
from ..utils.misc import run_cli_program_with_input

# Some hard-coded data
databases = {}

databases["hp-igneous"] = {
    "data_file": "hp634ver.dat",
    "solution_model_file": "solution_model.dat",
    "solutions_list": [],
    "excluded_endmember_list": [],
    "transformations": [
        ["Fe", "FeO", ["O2"], "1 -0.5"],
        ["Ni", "NiO", ["02"], "1 -0.5"],
        ["Cu", "CuO", ["O2"], "1 -0.5"],
        ["H2", "H2O", ["O2"], "1 -0.5"],
        ["C", "CO2", ["O2"], "1 -1"],
    ],
    "component_set": [
        "Na2O",
        "MgO",
        "Al2O3",
        "SiO2",
        "K2O",
        "CaO",
        "TiO2",
        "MnO",
        "Fe",
        "Ni",
        "ZrO2",
        "Cl2",
        "O2",
        "H2",
        "C",
        "Cu",
        "Cr2O3",
        "S2",
        "F2",
        "N2",
    ],
}

databases["stx21"] = {
    "data_file": "stx21ver.dat",
    "solution_model_file": "stx21_solution_model.dat",
    "solutions_list": [
        "C2/c",
        "Wus",
        "Pv",
        "Pl",
        "Sp",
        "O",
        "Wad",
        "Ring",
        "Opx",
        "Cpx",
        "Aki",
        "Gt",
        "Ppv",
        "CF",
        "NaAl",
    ],
    "excluded_endmember_list": [],
    "transformations": [
        ["FeO", "FEO", [], "1"],
        ["Na2O", "NA2O", [], "1"],
        ["MgO", "MGO", [], "1"],
        ["Al2O3", "AL2O3", [], "1"],
        ["SiO2", "SIO2", [], "1"],
        ["CaO", "CAO", [], "1"],
    ],
    "component_set": ["Na2O", "MgO", "Al2O3", "SiO2", "CaO", "FeO"],
}

databases["stx24"] = {
    "data_file": "stx24ver.dat",
    "solution_model_file": "stx24_solution_model.dat",
    "solutions_list": [
        "C2/c",
        "Wus",
        "Pv",
        "Pl",
        "Sp",
        "O",
        "Wad",
        "Ring",
        "Opx",
        "Cpx",
        "Aki",
        "Gt",
        "Ppv",
        "CF",
        "NaAl",
    ],
    "excluded_endmember_list": [],
    "transformations": [],
    "component_set": ["Na2O", "MgO", "Al2O3", "SiO2", "CaO", "Fe", "O", "Cr2O3"],
}


def ctransf_stdin(database):
    """
    Generate the standard input string for the Perple_X 'ctransf' utility.

    :param database: A dictionary containing database-related keys, including 'data_file' and 'transformations',
                     where 'transformations' is a list of component transformations to apply to the data file.
    :type database: dict
    :return: A string formatted for use as stdin to the 'ctransf' program.
    :rtype: str
    """

    # check if there are saturated components
    potentially_saturated_components = extract_lines_between_markers(
        database["data_file"],
        "begin_special_components",
        "end_special_components",
        inclusive=False,
    )

    stdin = f"{database['data_file']}\n"
    for t in database["transformations"]:
        if t[1] in potentially_saturated_components:
            sat = "N\n"
        else:
            sat = ""
        comps = "\n".join(t[2])

        stdin += f"Y\n{t[0]}\n{t[1]}\n{sat}{comps}\n\n{t[3]}\nY\n"
    stdin += "\n"
    return stdin


def build_stdin(
    project_name,
    database,
    perplex_option_file,
    composition,
    pressure_range,
    temperature_range,
):
    """
    Generate the standard input string required for the Perple_X 'build' utility.

    :param project_name: Name of the project; used to name output files.
    :type project_name: str
    :param database: Dictionary containing keys like 'solution_model_file', 'dataset_file', and possibly others.
    :type database: dict
    :param perplex_option_file: Path to the Perple_X option file containing user-selected build parameters.
    :type perplex_option_file: str
    :param composition: burnman.Composition object containing the desired composition.
    :type composition: dict
    :param pressure_range: Tuple specifying the pressure range as (min_pressure, max_pressure).
    :type pressure_range: tuple[float, float]
    :param temperature_range: Tuple specifying the temperature range as (min_temperature, max_temperature).
    :type temperature_range: tuple[float, float]
    :return: A formatted string to be passed to the Perple_X 'build' executable as stdin.
    :rtype: str
    """

    minimum_pressure_bar = pressure_range[0] / 1.0e5
    maximum_pressure_bar = pressure_range[1] / 1.0e5

    # check if there are saturated components
    if (
        len(
            extract_lines_between_markers(
                "ctransf.dat",
                "begin_special_components",
                "end_special_components",
                inclusive=False,
            )
        )
        > 0
    ):
        saturated_string = "N\n"
    else:
        saturated_string = ""

    c = deepcopy(composition)
    c.change_component_set(database["component_set"])
    c.remove_null_components()
    c.renormalize("mass", "total", 100.0)
    components = "\n".join(c.mass_composition.keys())
    component_amounts = "\n".join(
        format(x, "10.5f") for x in c.mass_composition.values()
    )

    if len(database["excluded_endmember_list"]) > 0:
        exclude_phases = "\n".join(database["excluded_endmember_list"]) + "\n"
    else:
        exclude_phases = ""
    if len(database["solutions_list"]) > 0:
        solutions = "\n".join(database["solutions_list"]) + "\n"
    else:
        solutions = ""

    stdin = (
        f"{project_name}\n"
        "ctransf.dat\n"
        f"{perplex_option_file}\n"
        "N\n"  # transform components
        "2\n"  # minimization on a 2D grid
        "N\n"  # no saturated fluids
        f"{saturated_string}"  # no saturated components
        "N\n"  # Use chemical potentials, ... as independent variables
        f"{components}\n\n"  # last \n is because a blank is required to end the list
        "n\n"  # P is not dependent on T
        "1\n"  # P represents the x-axis
        f"{minimum_pressure_bar} {maximum_pressure_bar}\n"
        f"{temperature_range[0]} {temperature_range[1]}\n"
        "Y\n"  # compositions as mass amounts
        f"{component_amounts}\n"
        "N\n"  # no print file
        "Y\n"  # exclude endmember phases
        "N\n"  # don't prompt for endmember phases
        f"{exclude_phases}\n"
        "Y\n"  # include solution models
        f"{database["solution_model_file"]}\n"
        f"{solutions}\n"
        f"{project_name}\n"  # title, program ends without needing another carriage return
    )
    return stdin


def make_build_file(
    perplex_bindir,
    project_name,
    database,
    perplex_option_file,
    composition,
    pressure_range,
    temperature_range,
    verbose=False,
):
    """
    Generates a Perple_X build file for a specified project.

    :param project_name: Name of the Perple_X project.
    :type project_name: str
    :param database: A dictionary containing database information, including 'data_file' and 'transformations'.
    :type database: dict
    :param perplex_bindir: Path to the Perple_X binary directory.
    :type perplex_bindir: str
    :param perplex_option_file: Path to the Perple_X option file.
    :type perplex_option_file: str
    :param composition: Dictionary of component compositions for the system.
    :type composition: dict
    :param pressure_range: Tuple indicating pressure range (min, max).
    :type pressure_range: tuple[float, float]
    :param temperature_range: Tuple indicating temperature range (min, max).
    :type temperature_range: tuple[float, float]
    :param verbose: Whether to print verbose output, defaults to False.
    :type verbose: bool, optional
    """

    ctransf_path = os.path.join(perplex_bindir, "ctransf")
    build_path = os.path.join(perplex_bindir, "build")

    if verbose:
        print("Transforming data file...")

    str_ctransf = ctransf_stdin(database)
    if not database.get("transformations"):
        shutil.copyfile(database["data_file"], "ctransf.dat")
    else:
        run_cli_program_with_input(ctransf_path, str_ctransf, verbose)

    if verbose:
        print("Generating build input...")

    str_build = build_stdin(
        project_name,
        database,
        perplex_option_file,
        composition,
        pressure_range,
        temperature_range,
    )

    if verbose:
        print(f"Removing old build file if exists: {project_name}.dat")
    try:
        os.remove(f"{project_name}.dat")
    except FileNotFoundError:
        pass

    if verbose:
        print("Generating build file...")

    run_cli_program_with_input(build_path, str_build, verbose)


def run_vertex(perplex_bindir, project_name, verbose=False):
    """
    Run Perple_X vertex on a project.

    :param vertex_path: Path to Perple_X vertex
    :type vertex_path: string
    :param project_name: Name of project defined during build
    :type project_name: string
    """
    vertex_path = os.path.join(perplex_bindir, "vertex")
    run_cli_program_with_input(vertex_path, project_name, verbose)


def run_pssect(perplex_bindir, project_name, convert_to_pdf=True, verbose=False):
    """
    Runs the Perple_X 'pssect' program to generate a PostScript section plot,
    and optionally converts the resulting .ps file to a .pdf file.

    :param perplex_bindir: Path to the directory containing Perple_X executables.
    :type perplex_bindir: str
    :param project_name: Name of the Perple_X project (used to locate input/output files).
    :type project_name: str
    :param convert_to_pdf: Whether to convert the output .ps file to .pdf using Ghostscript, defaults to True.
    :type convert_to_pdf: bool, optional
    :param verbose: Whether to print detailed status messages, defaults to False.
    :type verbose: bool, optional
    """
    pssect_path = os.path.join(perplex_bindir, "pssect")
    input_str = f"{project_name}\nn\n"

    if verbose:
        print("Running pssect...")
    run_cli_program_with_input(pssect_path, input_str, verbose)

    ps_file = f"{project_name}.ps"
    pdf_file = f"{project_name}.pdf"

    if convert_to_pdf:
        if verbose:
            print(f"Converting {ps_file} to {pdf_file}...")

        try:
            subprocess.run(["ps2pdf", ps_file, pdf_file], check=True)
            if verbose:
                print(f"Conversion successful: {pdf_file}")
        except FileNotFoundError:
            print("Error: 'ps2pdf' not found. Please install Ghostscript.")
        except subprocess.CalledProcessError as e:
            print(f"Conversion failed: {e}")


def create_perplex_table(
    perplex_bindir,
    project_name,
    outfile,
    n_pressures,
    n_temperatures,
    pressure_range=None,
    temperature_range=None,
    verbose=False,
):
    """
    Uses Perple_X's 'werami' program to output a thermodynamic property table
    for a previously computed model.

    Properties output:
        2 - Density (kg/m³)
        4 - Expansivity (1/K)
        5 - Compressibility (1/bar)
        10 - Adiabatic bulk modulus (bar)
        11 - Adiabatic shear modulus (bar)
        12 - Sound velocity (km/s)
        13 - P-wave velocity (Vp, km/s)
        14 - S-wave velocity (Vs, km/s)
        17 - Entropy (J/K/kg)
        18 - Enthalpy (J/kg)
        19 - Heat Capacity (J/K/kg)
        22 - Molar Volume (J/bar)

    :param perplex_bindir: Path to directory containing Perple_X executables.
    :param project_name: Name of the Perple_X project (input for werami).
    :param outfile: Desired output filename for the table.
    :param n_pressures: Number of pressure steps.
    :param n_temperatures: Number of temperature steps.
    :param pressure_range: Optional. Tuple of (Pmin, Pmax) in Pascals. Defaults to the range given in the build file.
    :param temperature_range: Optional. Tuple of (Tmin, Tmax) in Kelvin. Defaults to the range given in the build file.
    :param verbose: Print detailed process info (default: False).
    """
    if verbose:
        print(f"Creating a {n_pressures}x{n_temperatures} P-T table using werami...")

    werami_path = os.path.join(perplex_bindir, "werami")

    if pressure_range and temperature_range:
        pt_input = (
            "y\n"
            f"{pressure_range[0] / 1e5:.6f} {pressure_range[1] / 1e5:.6f}\n"
            f"{temperature_range[0]} {temperature_range[1]}\n"
        )
    else:
        if verbose:
            print("Using original P-T grid from vertex run.")
        pt_input = "n\n"

    properties = [2, 4, 5, 10, 11, 12, 13, 14, 17, 18, 19, 22]
    stdin_lines = [project_name, "2"]  # table mode
    for prop in properties:
        stdin_lines.append(str(prop))
        stdin_lines.append("n")  # default output format

    stdin_lines += ["0", pt_input.strip(), f"{n_pressures} {n_temperatures}", "0"]
    stdin_str = "\n".join(stdin_lines) + "\n"

    try:
        process = subprocess.Popen(
            werami_path,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8",
        )

        process.stdin.write(stdin_str)
        process.stdin.close()

        output_filename = None
        stdout_data = []

        for line in process.stdout:
            stdout_data.append(line)
            if "missing *.tof file" in line:
                raise RuntimeError(
                    f"You must run vertex before using werami on the project '{project_name}'."
                )
            if "set sample_on_grid to F" in line:
                raise RuntimeError(
                    "Set 'sample_on_grid = F' in the perplex option file (e.g., perplex_option.dat)."
                )
            if "Output has been written to the" in line:
                output_filename = line.strip().split()[-1]

        process.wait()
        output_text = "".join(stdout_data)

        if verbose:
            print(output_text)

        if not output_filename or not os.path.isfile(output_filename):
            raise FileNotFoundError(
                "werami did not produce a recognizable output file."
            )

        os.rename(output_filename, outfile)
        if verbose:
            print(f"Output file renamed to '{outfile}'")
            print("Processing complete.")

    except FileNotFoundError:
        raise RuntimeError(f"'werami' not found at path: {werami_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to run werami: {e}")


class PerplexMaterial(Material):
    """
    This is the base class for a PerpleX material. States of the material
    can only be queried after setting the pressure and temperature
    using set_state().

    Instances of this class are initialised with
    a 2D PerpleX tab file. This file should be in the standard format
    (as output by werami), and should have columns with the following names:
    'rho,kg/m3', 'alpha,1/K', 'beta,1/bar', 'Ks,bar', 'Gs,bar', 'v0,km/s',
    'vp,km/s', 'vs,km/s', 's,J/K/kg', 'h,J/kg', 'cp,J/K/kg', 'V,J/bar/mol'.
    The order of these names is not important.

    Properties of the material are determined by linear interpolation from
    the PerpleX grid. They are all returned in SI units on a molar basis,
    even though the PerpleX tab file is not in these units.

    This class is available as ``burnman.PerplexMaterial``.
    """

    def __init__(self, tab_file, name="Perple_X material"):
        self.name = name
        self.params = {"name": name}
        (
            self._property_interpolators,
            self.params["molar_mass"],
            self.bounds,
        ) = self._read_2D_perplex_file(tab_file)
        Material.__init__(self)

    def _read_2D_perplex_file(self, filename):
        with open(filename, "r") as f:
            datastream = f.read()

        lines = [
            line.strip().split() for line in datastream.split("\n") if line.strip()
        ]

        if lines[2][0] != "2":
            raise Exception("This is not a 2D PerpleX table")

        bounds = [
            (float(lines[4][0]), float(lines[5][0]), int(lines[6][0])),
            (float(lines[8][0]), float(lines[9][0]), int(lines[10][0])),
        ]

        if lines[3][0] == "P(bar)" and lines[7][0] == "T(K)":
            Pmin, Pint, nP = bounds[0]
            Tmin, Tint, nT = bounds[1]
        elif lines[3][0] == "T(K)" and lines[7][0] == "P(bar)":
            Pmin, Pint, nP = bounds[1]
            Tmin, Tint, nT = bounds[0]
        else:
            raise Exception(
                "This file does not have a recognised PerpleX structure.\n"
                "Are the independent variables P(bar) and T(K)?"
            )

        Pmin = Pmin * 1.0e5  # bar to Pa
        Pint = Pint * 1.0e5  # bar to Pa
        Pmax = Pmin + Pint * (nP - 1.0)
        Tmax = Tmin + Tint * (nT - 1.0)
        pressures = np.linspace(Pmin, Pmax, nP)
        temperatures = np.linspace(Tmin, Tmax, nT)
        n_properties = int(lines[11][0])
        property_list = lines[12]

        # property_table[i][j][k] returns the kth property at the ith pressure and jth temperature
        table = np.array(
            [[float(string) for string in line] for line in lines[13 : 13 + nP * nT]]
        )

        if lines[3][0] == "P(bar)":
            property_table = np.swapaxes(table.reshape(nT, nP, n_properties), 0, 1)
        else:
            property_table = table.reshape(nP, nT, n_properties)

        ordered_property_list = [
            "rho,kg/m3",
            "alpha,1/K",
            "beta,1/bar",
            "Ks,bar",
            "Gs,bar",
            "v0,km/s",
            "vp,km/s",
            "vs,km/s",
            "s,J/K/kg",
            "h,J/kg",
            "cp,J/K/kg",
            "V,J/bar/mol",
        ]
        p_indices = [
            i
            for i, p in enumerate(property_list)
            for ordered_p in ordered_property_list
            if p == ordered_p
        ]

        properties = {}
        for i, p_idx in enumerate(p_indices):
            # Fill in NaNs as long as they aren't in the corners of the P-T grid
            a = np.array(property_table[:, :, [p_idx]][:, :, 0])
            x, y = np.indices(a.shape)
            a[np.isnan(a)] = griddata(
                (x[~np.isnan(a)], y[~np.isnan(a)]),  # points we know
                a[~np.isnan(a)],  # values we know
                (x[np.isnan(a)], y[np.isnan(a)]),
            )

            # Fill any remaining NaNs with zeros
            properties[ordered_property_list[i]] = np.nan_to_num(a, 0.0)

        densities = properties["rho,kg/m3"]
        volumes = 1.0e-5 * properties["V,J/bar/mol"]
        molar_masses = densities * volumes
        molar_mass = np.mean(molar_masses)

        property_interpolators = {
            "rho": RegularGridInterpolator(
                (pressures, temperatures), densities, bounds_error=True
            ),
            "alpha": RegularGridInterpolator(
                (pressures, temperatures), properties["alpha,1/K"], bounds_error=True
            ),
            "K_T": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e5 / properties["beta,1/bar"],
                bounds_error=True,
            ),
            "K_S": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e5 * properties["Ks,bar"],
                bounds_error=True,
            ),
            "G_S": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e5 * properties["Gs,bar"],
                bounds_error=True,
            ),
            "bulk_sound_velocity": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e3 * properties["v0,km/s"],
                bounds_error=True,
            ),
            "p_wave_velocity": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e3 * properties["vp,km/s"],
                bounds_error=True,
            ),
            "s_wave_velocity": RegularGridInterpolator(
                (pressures, temperatures),
                1.0e3 * properties["vs,km/s"],
                bounds_error=True,
            ),
            "S": RegularGridInterpolator(
                (pressures, temperatures),
                properties["s,J/K/kg"] * molar_masses,
                bounds_error=True,
            ),
            "H": RegularGridInterpolator(
                (pressures, temperatures),
                properties["h,J/kg"] * molar_masses,
                bounds_error=True,
            ),
            "C_p": RegularGridInterpolator(
                (pressures, temperatures),
                properties["cp,J/K/kg"] * molar_masses,
                bounds_error=True,
            ),
            "V": RegularGridInterpolator(
                (pressures, temperatures), volumes, bounds_error=True
            ),
        }

        bounds = [[Pmin, Pmax], [Tmin, Tmax]]
        return property_interpolators, molar_mass, bounds

    @copy_documentation(Material.set_state)
    def set_state(self, pressure, temperature):
        if not np.logical_and(
            np.all(self.bounds[0][0] <= pressure), np.all(pressure <= self.bounds[0][1])
        ):
            raise ValueError(
                "The set_state pressure ({0:.4f}) is outside the bounds of this rock ({1:.4f}-{2:.4f} GPa)".format(
                    pressure, self.bounds[0][0] / 1.0e9, self.bounds[0][1] / 1.0e9
                )
            )
        if not np.logical_and(
            np.all(self.bounds[1][0] <= temperature),
            np.all(temperature <= self.bounds[1][1]),
        ):
            raise ValueError(
                "The set_state temperature ({0:.1f}) is outside the bounds of this rock ({1:.1f}-{2:.1f} K)".format(
                    temperature, self.bounds[1][0], self.bounds[1][1]
                )
            )
        Material.set_state(self, pressure, temperature)

    """
    Properties by linear interpolation of Perple_X output
    """

    @material_property
    @copy_documentation(Material.molar_volume)
    def molar_volume(self):
        return self._property_interpolators["V"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.molar_enthalpy)
    def molar_enthalpy(self):
        return self._property_interpolators["H"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.molar_entropy)
    def molar_entropy(self):
        return self._property_interpolators["S"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.isothermal_bulk_modulus_reuss)
    def isothermal_bulk_modulus_reuss(self):
        return self._property_interpolators["K_T"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.isentropic_bulk_modulus_reuss)
    def isentropic_bulk_modulus_reuss(self):
        return self._property_interpolators["K_S"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.molar_heat_capacity_p)
    def molar_heat_capacity_p(self):
        return self._property_interpolators["C_p"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.thermal_expansivity)
    def thermal_expansivity(self):
        return self._property_interpolators["alpha"]([self.pressure, self.temperature])[
            0
        ]

    @material_property
    @copy_documentation(Material.shear_modulus)
    def shear_modulus(self):
        return self._property_interpolators["G_S"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.p_wave_velocity)
    def p_wave_velocity(self):
        return self._property_interpolators["p_wave_velocity"](
            [self.pressure, self.temperature]
        )[0]

    @material_property
    @copy_documentation(Material.bulk_sound_velocity)
    def bulk_sound_velocity(self):
        return self._property_interpolators["bulk_sound_velocity"](
            [self.pressure, self.temperature]
        )[0]

    @material_property
    @copy_documentation(Material.shear_wave_velocity)
    def shear_wave_velocity(self):
        return self._property_interpolators["s_wave_velocity"](
            [self.pressure, self.temperature]
        )[0]

    """
    Properties from mineral parameters,
    Legendre transformations
    or Maxwell relations
    """

    @material_property
    @copy_documentation(Material.molar_gibbs)
    def molar_gibbs(self):
        return self.molar_enthalpy - self.temperature * self.molar_entropy

    @material_property
    @copy_documentation(Material.molar_mass)
    def molar_mass(self):
        if "molar_mass" in self.params:
            return self.params["molar_mass"]
        else:
            raise ValueError(
                "No molar_mass parameter for mineral " + self.to_string + "."
            )

    @material_property
    @copy_documentation(Material.density)
    def density(self):
        return self._property_interpolators["rho"]([self.pressure, self.temperature])[0]

    @material_property
    @copy_documentation(Material.molar_internal_energy)
    def molar_internal_energy(self):
        return (
            self.molar_gibbs
            - self.pressure * self.molar_volume
            + self.temperature * self.molar_entropy
        )

    @material_property
    @copy_documentation(Material.molar_helmholtz)
    def molar_helmholtz(self):
        return self.molar_gibbs - self.pressure * self.molar_volume

    @material_property
    @copy_documentation(Material.isothermal_compressibility_reuss)
    def isothermal_compressibility_reuss(self):
        return 1.0 / self.isothermal_bulk_modulus_reuss

    @material_property
    @copy_documentation(Material.isentropic_compressibility_reuss)
    def isentropic_compressibility_reuss(self):
        return 1.0 / self.isentropic_bulk_modulus_reuss

    @material_property
    @copy_documentation(Material.molar_heat_capacity_v)
    def molar_heat_capacity_v(self):
        return (
            self.molar_heat_capacity_p
            - self.molar_volume
            * self.temperature
            * self.thermal_expansivity
            * self.thermal_expansivity
            * self.isothermal_bulk_modulus_reuss
        )

    @material_property
    @copy_documentation(Material.grueneisen_parameter)
    def grueneisen_parameter(self):
        return (
            self.thermal_expansivity
            * self.molar_volume
            * self.isentropic_bulk_modulus_reuss
            / self.molar_heat_capacity_p
        )
