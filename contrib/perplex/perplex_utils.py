# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.
import os
import subprocess
import shutil
from copy import deepcopy
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from scipy.signal import savgol_filter
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from shapely.ops import polygonize, unary_union
from shapely.prepared import prep
from shapely.strtree import STRtree
from shapely.errors import GEOSException
import itertools

from burnman.utils.misc import extract_lines_between_markers
from burnman.utils.misc import run_cli_program_with_input
from burnman.utils.plotting import visual_center_of_polygon

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
        "Y\n"  # no print file
        "Y\n"  # exclude endmember phases
        "N\n"  # don't prompt for endmember phases
        f"{exclude_phases}\n"
        "Y\n"  # include solution models
        f"{database['solution_model_file']}\n"
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

    return run_cli_program_with_input(build_path, str_build, verbose)


def run_vertex(perplex_bindir, project_name, verbose=False):
    """
    Run Perple_X vertex on a project.

    :param vertex_path: Path to Perple_X vertex
    :type vertex_path: string
    :param project_name: Name of project defined during build
    :type project_name: string
    """
    vertex_path = os.path.join(perplex_bindir, "vertex")
    return run_cli_program_with_input(vertex_path, project_name, verbose)


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


def create_werami_table(
    perplex_bindir,
    project_name,
    outfile,
    n_pressures,
    n_temperatures,
    property_ids,
    pressure_range=None,
    temperature_range=None,
    verbose=False,
):
    """
    Uses Perple_X's 'werami' program to output a thermodynamic property table
    for a previously computed model.

    :param perplex_bindir: Path to directory containing Perple_X executables.
    :param project_name: Name of the Perple_X project (input for werami).
    :param outfile: Desired output filename for the table.
    :param n_pressures: Number of pressure steps.
    :param n_temperatures: Number of temperature steps.
    :param property_ids: List of property IDs to include in the output table.
    :type property_ids: list[int]
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

    stdin_lines = [project_name, "2"]  # table mode
    for prop in property_ids:
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


def create_mode_table(
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
    Uses Perple_X's 'werami' program to output a table
    containing all the phase modes from a vertex calculation.

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

    stdin_lines = [project_name, "2"]  # table mode
    stdin_lines.append("25")
    stdin_lines.append("n")  # default output format

    stdin_lines += ["n", f"{n_pressures} {n_temperatures}", "0"]
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


def create_perplex_class_table(
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

    property_ids = [2, 4, 5, 10, 11, 12, 13, 14, 17, 18, 19, 22]
    create_werami_table(
        perplex_bindir,
        project_name,
        outfile,
        n_pressures,
        n_temperatures,
        property_ids,
        pressure_range,
        temperature_range,
        verbose,
    )


def load_phase_mode_tab_file(filename):
    """
    Load phase mode data from a tab file.

    :param filename: Path to the tab-delimited file containing phase mode data.
    :type filename: str
    :return: Tuple of (data array reshaped to (n_P, n_T, -1).T, column names list)
    :rtype: (numpy.ndarray, list[str])
    """
    with open(filename, "r") as file:
        lines = file.readlines()

    n_P = int(lines[6].strip().split()[0])
    n_T = int(lines[10].strip().split()[0])
    header_index = 12
    column_names = lines[header_index].split()

    data = [
        [float(val) if val.lower() != "nan" else 0.0 for val in line.split()]
        for line in lines[header_index + 1 :]
        if line.strip()
    ]

    return np.array(data).reshape(n_T, n_P, -1).T, column_names


def _split_paths_by_gap(paths, dx, dy, max_gap=2):
    """
    Split paths into segments where gaps between points exceed a threshold.

    :param paths: List of numpy arrays representing point sequences.
    :type paths: list[numpy.ndarray]
    :param dx: Scale factor for x-coordinate differences.
    :type dx: float
    :param dy: Scale factor for y-coordinate differences.
    :type dy: float
    :param max_gap: Maximum allowed gap between consecutive points
    (in scaled units).
    :type max_gap: float
    :return: List of numpy arrays with split path segments.
    :rtype: list[numpy.ndarray]
    """
    split_paths = []
    for path in paths:
        if len(path) == 0:
            continue
        dists = np.sqrt(
            (np.diff(path[:, 0]) / dx) ** 2 + (np.diff(path[:, 1]) / dy) ** 2
        )
        split_indices = np.where(dists > max_gap)[0]

        start = 0
        for idx in split_indices:
            if idx + 1 > start + 1:
                split_paths.append(path[start : idx + 1])
            start = idx + 1
        if len(path[start:]) > 1:
            split_paths.append(path[start:])
    return split_paths


def _find_contours(x, y, z, levels=[1.0e-10]):
    """
    Find contour lines at specified levels from gridded data.

    :param x: 2D array of x coordinates.
    :type x: numpy.ndarray
    :param y: 2D array of y coordinates.
    :type y: numpy.ndarray
    :param z: 2D array of data values.
    :type z: numpy.ndarray
    :param levels: Contour levels to find.
    :type levels: list[float]
    :return: List of contour segments split by gaps.
    :rtype: list[numpy.ndarray]
    """
    tmpfig = plt.figure()
    tmpax = tmpfig.add_subplot(111)
    CS = tmpax.contour(x, y, z, levels=levels, zorder=0)
    intersection_curves = [path.vertices for path in CS._paths]
    plt.close(tmpfig)

    dx, dy = x[1, 0] - x[0, 0], y[0, 1] - y[0, 0]
    return _split_paths_by_gap(intersection_curves, dx, dy, max_gap=2)


def _savgol_filter_points(points, window_length, polyorder, mode="nearest"):
    """
    Smooth points using a Savitzky-Golay filter.

    :param points: Array of 2D points.
    :type points: numpy.ndarray
    :param window_length: Length of the filter window (must be odd).
    :type window_length: int
    :param polyorder: Polynomial order to fit.
    :type polyorder: int
    :param mode: Mode by which to smooth the data.
    :type mode: string
    :return: Smoothed points array.
    :rtype: numpy.ndarray
    """
    if len(points) < 2:
        return points.copy()

    # Make filter independent of direction traversed along line
    x, y = points[:, 0], points[:, 1]
    xhat = (
        savgol_filter(x, window_length, polyorder, mode=mode)
        + savgol_filter(x[::-1], window_length, polyorder, mode=mode)[::-1]
    ) / 2.0
    yhat = (
        savgol_filter(y, window_length, polyorder, mode=mode)
        + savgol_filter(y[::-1], window_length, polyorder, mode=mode)[::-1]
    ) / 2.0
    xhat[0], xhat[-1] = x[0], x[-1]
    yhat[0], yhat[-1] = y[0], y[-1]
    return np.vstack((xhat, yhat)).T


def _strip_colinear_points(points, eps=1e-6):
    """
    Remove colinear and nearly colinear points from a sequence.

    :param points: Array of 2D points.
    :type points: numpy.ndarray
    :param eps: Tolerance for colinearity.
    :type eps: float
    :return: Filtered points array with colinear points removed.
    :rtype: numpy.ndarray
    """
    if len(points) <= 2:
        return points.copy()
    keep = [0]
    for i in range(1, len(points) - 1):
        a, b, c = points[i - 1], points[i], points[i + 1]
        ab = np.append(b - a, 0)
        bc = np.append(c - b, 0)
        cross = np.linalg.norm(np.cross(ab, bc))
        if cross > eps:
            keep.append(i)
    keep.append(len(points) - 1)
    return points[keep]


def _smooth_polygon_edges_between_intersections(
    polygon, intersection_tree, window_length=6, polyorder=2, eps=1e-4
):
    """
    Smooth polygon edges between intersection points using
    a Savitzky-Golay filter.

    :param polygon: Shapely Polygon to smooth.
    :type polygon: shapely.geometry.Polygon
    :param intersection_tree: STRtree spatial index of intersection points.
    :type intersection_tree: shapely.strtree.STRtree
    :param window_length: Savitzky-Golay filter window length.
    :type window_length: int
    :param polyorder: Savitzky-Golay polynomial order.
    :type polyorder: int
    :param eps: Buffer radius to detect intersections.
    :type eps: float
    :return: Tuple of smoothed Polygon and list of smoothed edge point arrays.
    :rtype: (shapely.geometry.Polygon, list[numpy.ndarray])
    """
    if not polygon.is_valid or polygon.is_empty:
        return polygon

    coords = np.array(polygon.exterior.coords)
    smoothed_coords = []
    segment = [coords[0]]

    intersecting_nodes = [
        i
        for i in range(len(coords))
        if len(intersection_tree.query(Point(coords[i]).buffer(eps))) > 0
    ]

    # assert intersecting_nodes[0] == 0
    # assert intersecting_nodes[-1] == len(coords) - 1

    segments = [
        coords[start : end + 1]
        for start, end in zip(intersecting_nodes[:-1], intersecting_nodes[1:])
    ]
    smoothed_segments = [
        _savgol_filter_points(segment, window_length, polyorder) for segment in segments
    ]

    if len(smoothed_segments) > 0:
        smoothed_coords = [smoothed_segments[0][0]]
        for segment in smoothed_segments:
            smoothed_coords.extend(segment[1:])

        try:
            return Polygon(smoothed_coords), smoothed_segments
        except ValueError:
            return polygon, coords
    else:
        return polygon, coords


def get_fields_assemblages_and_bounds(
    filename, phase_name_replacements, smoothing_window, smoothing_order
):
    """
    Extract phase assemblage fields and plot bounds from a mode tab file.

    :param filename: Path to the mode tab file.
    :type filename: str
    :param phase_name_replacements: Dictionary mapping original phase names
    to replacement names.
    :type phase_name_replacements: dict[str, str]
    :param smoothing_window: Window length for smoothing polygon edges.
    :type smoothing_window: int
    :param smoothing_order: Polynomial order for smoothing polygon edges.
    :type smoothing_order: int
    :return: Tuple containing list of polygon data dicts and bounds
    [[P_min, P_max], [T_min, T_max]].
    :rtype: (list[dict], list[list[float]])
    """
    data, column_names = load_phase_mode_tab_file(filename)

    P, T = data[:2]
    pressures = P[:, 0]
    temperatures = T[0]

    phase_modes = data[2:]
    num_phases = phase_modes.shape[0]
    phase_names = [phase_name_replacements.get(name, name) for name in column_names[2:]]

    eps_edge_T = 1.0e-5
    eps_edge_P = 1.0e-5

    lines = []
    for i in range(num_phases):
        Z = phase_modes[i]
        contours = _find_contours(P, T, Z, levels=[1e-10])
        for contour in contours:
            contour = _strip_colinear_points(contour)

            # Shift all contours at the edge of the domain
            # to be a tiny bit wider so that the intersection works
            mask = contour[:, 0] < pressures[0] + eps_edge_P
            contour[mask, 0] = contour[mask, 0] - eps_edge_P
            mask = contour[:, 0] > pressures[-1] - eps_edge_P
            contour[mask, 0] = contour[mask, 0] + eps_edge_P

            mask = contour[:, 1] < temperatures[0] + eps_edge_T
            contour[mask, 1] = contour[mask, 1] - eps_edge_T
            mask = contour[:, 1] > temperatures[-1] - eps_edge_T
            contour[mask, 1] = contour[mask, 1] + eps_edge_T

            lines.append(contour)

    # Shift all contours at the edge of the domain
    # to be a tiny bit wider but less so than
    # the other contour lines so that the intersection works
    P_min = pressures[0] - eps_edge_P / 2.0
    P_max = pressures[-1] + eps_edge_P / 2.0
    T_min = temperatures[0] - eps_edge_T / 2.0
    T_max = temperatures[-1] + eps_edge_T / 2.0
    lines.extend(
        [
            LineString([[P_min, T_min], [P_max, T_min]]),  # bottom
            LineString([[P_min, T_max], [P_max, T_max]]),  # top
            LineString([[P_min, T_min], [P_min, T_max]]),  # left
            LineString([[P_max, T_min], [P_max, T_max]]),  # right
        ]
    )

    lines = [LineString(arr) for arr in lines]
    multi_line = unary_union(lines)
    polygons = list(polygonize(multi_line))

    # Now get all the intersections
    unique_intersections = [
        Point(x, y)
        for x, y in {
            (round(p.x, 10), round(p.y, 10))
            for line1, line2 in itertools.combinations(lines, 2)
            for inter in [line1.intersection(line2)]
            if not inter.is_empty
            for p in (
                [inter]
                if inter.geom_type == "Point"
                else inter.geoms if inter.geom_type == "MultiPoint" else []
            )
        }
    ]
    intersection_tree = STRtree(unique_intersections)

    grid_points = np.array([P.ravel(), T.ravel()]).T

    phase_modes_flat = phase_modes.reshape(phase_modes.shape[0], -1)

    polygon_data = []

    for poly in polygons:
        smoothed_poly, smoothed_edges = _smooth_polygon_edges_between_intersections(
            poly,
            intersection_tree,
            window_length=smoothing_window,
            polyorder=smoothing_order,
            eps=1e-8,
        )

        select_smoothed_edges = []
        for edge in smoothed_edges:
            # Don't plot if 1D or at the edges of the domain
            if len(edge.shape) == 2:
                if (
                    np.max(edge[:, 0]) > pressures[0]
                    and np.min(edge[:, 0]) < pressures[-1]
                    and np.max(edge[:, 1]) > temperatures[0]
                    and np.min(edge[:, 1]) < temperatures[-1]
                ):

                    select_smoothed_edges.append(edge)

        x, y = smoothed_poly.exterior.xy
        contour = np.vstack((x, y)).T

        # Build assemblage string and count phases
        prepared_poly = prep(poly)
        inside_mask = np.array([prepared_poly.contains(Point(p)) for p in grid_points])
        if not inside_mask.any():
            continue

        eps_mode = 1.0e-6
        sum_phase_modes = np.sum(phase_modes_flat[:, inside_mask], axis=1)
        if np.max(sum_phase_modes) <= eps_mode:
            continue

        phases_present = [
            phase_names[i] for i, v in enumerate(sum_phase_modes) if v > eps_mode
        ]
        assemblage_string = " ".join(sorted(phases_present))
        n_phases = len(phases_present)

        # Scale pressures to GPa
        # We do this here rather than above to avoid misidentification of
        # empty fields
        contour[:, 0] = contour[:, 0] / 1.0e4
        for i in range(len(select_smoothed_edges)):
            select_smoothed_edges[i][:, 0] = select_smoothed_edges[i][:, 0] / 1.0e4

        polygon_data.append(
            {
                "contour": contour,
                "edges": select_smoothed_edges,
                "label": assemblage_string,
                "n_phases": n_phases,
            }
        )

    bounds = [
        [pressures[0] / 1.0e4, pressures[-1] / 1.0e4],
        [temperatures[0], temperatures[-1]],
    ]
    return polygon_data, bounds


def get_label_position_from_polygon_contour(contour, x_range, y_range, label_scaling):
    """
    Compute label position for polygon contour based on the visual center
    of the polygon after scaling by given ranges.

    :param contour: Nx2 array of polygon contour points.
    :type contour: numpy.ndarray
    :param x_range: Range of x values.
    :type x_range: float
    :param y_range: Range of y values.
    :type y_range: float
    :param label_scaling: Scaling factor corresponding to an approximate value
    of the width:height ratio of the label.
    :type label_scaling: float
    :return: Position (x, y) to place label, float indicating scaled distance to the nearest edge.
    :rtype: numpy.ndarray, float
    """

    scaled = contour.copy()
    scaled[:, 0] /= x_range * label_scaling
    scaled[:, 1] /= y_range

    pos, dist = visual_center_of_polygon([scaled], with_distance=True)
    pos[0] *= x_range * label_scaling
    pos[1] *= y_range

    return pos, dist


def merge_polygons_by_label(polygon_data):

    grouped = defaultdict(list)

    # Group dicts by label
    for polygon in polygon_data:
        grouped[polygon["label"]].append(polygon)

    merged_polygons = []

    for label, group in grouped.items():
        # Merge all polygons with the same label
        polygons = [
            Polygon(c) for c in [item["contour"] for item in group] if len(c) >= 3
        ]

        try:
            merged_polygon = unary_union(polygons).buffer(0)
            if type(merged_polygon) is MultiPolygon:
                merged_polygons.extend(group)
            else:
                x, y = merged_polygon.exterior.xy
                merged_contour = np.vstack((x, y)).T
                merged_entry = {
                    "contour": merged_contour,
                    "label": label,
                    "n_phases": group[0]["n_phases"],
                    "edges": list(
                        itertools.chain.from_iterable(
                            item.get("edges", []) for item in group
                        )
                    ),
                }

                merged_polygons.append(merged_entry)
        except GEOSException:
            merged_polygons.extend(group)

    return merged_polygons


def merge_bounds(bounds_list):
    """
    Given a list of bounds (each as [[P_min, P_max], [T_min, T_max]]),
    return an overarching bound in the same structure.
    """
    p_mins = [b[0][0] for b in bounds_list]
    p_maxs = [b[0][1] for b in bounds_list]
    t_mins = [b[1][0] for b in bounds_list]
    t_maxs = [b[1][1] for b in bounds_list]

    return [
        [min(p_mins), max(p_maxs)],
        [min(t_mins), max(t_maxs)],
    ]


def pretty_plot_phase_diagram(
    ax,
    werami_mode_tab_filenames,
    phase_name_replacements,
    bounding_colors=["#44015a", "#ffffff"],
    smoothing_window=4,
    smoothing_order=1,
    linewidth=0.5,
    label_scaling=3.0,
    label_clearance=0.01,
    n_phases_bounds=[2, 10],
):
    """
    Plot the Perple_X calculated phase diagram on a matplotlib axis.

    :param ax: Matplotlib Axes object to plot on.
    :type ax: matplotlib.axes.Axes
    :param werami_mode_tab_filenames: Paths to mode tab files.
    :type werami_mode_tab_filenames: List[str]
    :param phase_name_replacements: Dict mapping phase names to replacement names.
    :type phase_name_replacements: dict[str, str]
    :param bounding_colors: List of two colors for bounding colormap.
    :type bounding_colors: [str, str]
    :param smoothing_window: Savitzky-Golay smoothing window length.
    :type smoothing_window: int
    :param smoothing_order: Savitzky-Golay smoothing polynomial order.
    :type smoothing_order: int
    :param label_scaling: Scaling factor corresponding to an
    approximate value of the width:height ratio of the label.
    :type label_scaling: float
    :param label_clearance: Only plot the label if there is a certain
    clearance to the edge of the field. Scaled to the height of the domain.
    :type label_clearance: float
    :return: None
    """

    polygon_data = []
    bounds = []

    for werami_mode_tab_filename in werami_mode_tab_filenames:
        polygon_data_i, bounds_i = get_fields_assemblages_and_bounds(
            werami_mode_tab_filename,
            phase_name_replacements,
            smoothing_window,
            smoothing_order,
        )
        polygon_data.extend(polygon_data_i)
        bounds.append(bounds_i)

    polygon_data = merge_polygons_by_label(polygon_data)
    bounds = merge_bounds(bounds)

    P_range = bounds[0][1] - bounds[0][0]
    T_range = bounds[1][1] - bounds[1][0]

    ax.set_facecolor("white")

    n_phases_min, n_phases_max = n_phases_bounds
    n_phases_range = np.arange(n_phases_min, n_phases_max + 1)

    # Generate N colors from the original colormap
    base_cmap = LinearSegmentedColormap.from_list("base", bounding_colors)
    discrete_colors = [
        base_cmap(i / (len(n_phases_range) - 1)) for i in range(len(n_phases_range))
    ]

    # Create a discrete ListedColormap
    cmap = ListedColormap(discrete_colors)
    boundaries = np.arange(
        n_phases_min - 0.5, n_phases_max + 1.5, 1
    )  # one bin per phase count
    norm = BoundaryNorm(boundaries, cmap.N)

    # ScalarMappable for colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)

    for polygon in polygon_data:

        color = cmap(norm(polygon["n_phases"]))

        contour = polygon["contour"]
        ax.fill(
            contour[:, 0],
            contour[:, 1],
            color=color,
            linewidth=0,
            edgecolor="none",
            zorder=2,
        )

        for edge in polygon["edges"]:
            ax.plot(edge[:, 0], edge[:, 1], c="black", linewidth=linewidth)

        label_pos, label_dist = get_label_position_from_polygon_contour(
            contour, P_range, T_range, label_scaling
        )
        if label_dist > label_clearance:
            ax.text(*label_pos, polygon["label"], fontsize=6, ha="center", va="center")

    cbar = plt.colorbar(sm, ax=ax, boundaries=boundaries, ticks=n_phases_range)
    cbar.set_label("n phases")

    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Temperature (K)")

    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.grid(True, linestyle="--", color="gray", alpha=0.7)
    ax.set_axisbelow(False)
