import os
import shutil

import burnman
from perplex_utils import databases, make_build_file
from perplex_utils import run_vertex, run_pssect
from perplex_utils import create_mode_table
from perplex_utils import pretty_plot_phase_diagram

import matplotlib.pyplot as plt


def run_perplex(
    perplex_bindir,
    project_name,
    database,
    perplex_option_file,
    composition,
    pressure_range,
    temperature_range,
):
    # Make the build file
    # This file contains the information about the project, such as the
    # composition, the database, and the options for Perple_X.
    print("* Making build file...")
    make_build_file(
        perplex_bindir,
        project_name,
        database,
        perplex_option_file,
        composition,
        pressure_range,
        temperature_range,
        verbose=False,
    )

    # Run the Perple_X vertex command to do the Gibbs minimization
    # and create output files containing the equilibrium assemblages.
    print("* Running vertex...")
    run_vertex(perplex_bindir, project_name, verbose=True)

    # Run the Perple_X pssect command to create a postscript
    # file containing the phase diagram.
    print("* Running pssect...")
    run_pssect(perplex_bindir, project_name, convert_to_pdf=False, verbose=False)


# The following compositions are from Xu et al. (2008)
# (https://doi.org/10.1016/j.epsl.2008.08.012)
compositions = {
    "pyrolite": burnman.Composition(
        {
            "SiO2": 38.71,
            "MgO": 49.85,
            "FeO": 6.17,
            "CaO": 2.94,
            "Al2O3": 2.22,
            "Na2O": 0.11,
        },
        "molar",
    ),
    "basalt": burnman.Composition(
        {
            "SiO2": 51.75,
            "MgO": 14.94,
            "FeO": 7.06,
            "CaO": 13.88,
            "Al2O3": 10.19,
            "Na2O": 2.18,
        },
        "molar",
    ),
    "harzburgite": burnman.Composition(
        {
            "SiO2": 36.07,
            "MgO": 56.51,
            "FeO": 6.07,
            "CaO": 0.81,
            "Al2O3": 0.53,
            "Na2O": 0.001,
        },
        "molar",
    ),
    "modified_harzburgite": burnman.Composition(
        {
            "SiO2": 36.04,
            "MgO": 56.54,
            "FeO": 5.97,
            "CaO": 0.79,
            "Al2O3": 0.65,
            "Na2O": 0.001,
        },
        "molar",
    ),
}


if __name__ == "__main__":
    # Define the project parameters
    project_name = "basalt"
    database = databases["stx21"]
    composition = compositions[project_name]
    pressure_range_total = [1.0e5, 140.0e9]
    temperature_range_total = [200.0, 3000.0]

    # Split pressure and temperature so that PerpleX
    # no temperature splits seems to make diagram with
    # less prominent discontinuities
    n_pressures_per_split = 121
    n_temperatures_per_split = 601
    n_splits_pressure = 14
    n_splits_temperature = 1

    # If this script has already been run, and you just want to
    # tweak the figure, perplex should not be run again.
    perplex_should_be_run = True

    # End of project definitions

    outfile_base = f"{project_name}_table"
    perplex_dir = os.path.join(os.getcwd(), "perplex-installer/Perple_X")
    perplex_bindir = os.path.join(os.getcwd(), "perplex-installer/Perple_X/bin")

    delta_P = (pressure_range_total[1] - pressure_range_total[0]) / n_splits_pressure
    delta_T = (
        temperature_range_total[1] - temperature_range_total[0]
    ) / n_splits_temperature

    if perplex_should_be_run:
        # Check if the project directory already exists
        # If it does, delete it to start fresh.
        try:
            shutil.rmtree(project_name)
            print(f"Deleted old project directory ({project_name}).")
        except FileNotFoundError:
            pass

        os.mkdir(project_name)
        os.chdir(project_name)
        print(f"Now working in newly created directory ({project_name}).")

        print("Creating local copies of requested thermodynamic data files...")
        endmember_file = os.path.join(perplex_dir, "datafiles", database["data_file"])
        solution_file = os.path.join(
            perplex_dir, "datafiles", database["solution_model_file"]
        )
        shutil.copyfile(
            endmember_file, os.path.join(os.getcwd(), database["data_file"])
        )
        shutil.copyfile(
            solution_file, os.path.join(os.getcwd(), database["solution_model_file"])
        )
        perplex_option_file = "burnman_perplex_options.dat"

        # Define the number of exploratory nodes.
        # These are the nodes over which the full Gibbs minimization is performed.
        # In the refinement stage, the number of nodes is increased.
        # Any solutions not present on any of the bounding nodes will be excluded
        # from the calculations.
        n_P_exploratory = int((n_pressures_per_split + 1) / 2)
        n_T_exploratory = int((n_temperatures_per_split + 1) / 2)

        # Create the Perple_X options file
        with open(perplex_option_file, "w") as f:
            f.write(
                "sample_on_grid            F\n"
            )  # Do not force sampling on grid when running werami
            f.write("auto_refine  auto\n")  # Automatically refine the grid
            f.write("grid_levels  1 1\n")  # Do not use adaptive grid refinement
            f.write(
                f"x_nodes {n_P_exploratory} {n_pressures_per_split}\n"
            )  # Exploratory and final number of pressure nodes
            f.write(
                f"y_nodes {n_T_exploratory} {n_temperatures_per_split}\n"
            )  # Exploratory and final number of temperature nodes

        # Create the table(s) using the defined parameters
        for i_P in range(n_splits_pressure):
            for i_T in range(n_splits_temperature):
                pressure_range = [
                    pressure_range_total[0] + i_P * delta_P,
                    pressure_range_total[0] + (i_P + 1) * delta_P,
                ]
                temperature_range = [
                    temperature_range_total[0] + i_T * delta_T,
                    temperature_range_total[0] + (i_T + 1) * delta_T,
                ]
                if n_splits_pressure > 1 or n_splits_temperature > 1:
                    split_project_name = f"{project_name}_P{i_P + 1:02d}_T{i_T + 1:02d}"
                else:
                    split_project_name = project_name

                print(
                    f"Creating table [{i_P + 1}/{n_splits_pressure}, {i_T + 1}/{n_splits_temperature}]"
                )
                # Create the table for the current split
                run_perplex(
                    perplex_bindir,
                    split_project_name,
                    database,
                    perplex_option_file,
                    composition,
                    pressure_range,
                    temperature_range,
                )

                modes_filename = split_project_name + "_modes.dat"

                print("* Creating phase mode table...")
                create_mode_table(
                    perplex_bindir,
                    split_project_name,
                    modes_filename,
                    n_pressures_per_split,
                    n_temperatures_per_split,
                    pressure_range,
                    temperature_range,
                )
    else:
        os.chdir(project_name)

    print("Plotting phase diagram")
    phase_name_replacements = {
        "O": "ol",
        "Wad": "wad",
        "Ring": "rw",
        "Cpx": "cpx",
        "Opx": "opx",
        "Gt": "gt",
        "Aki": "aki",
        "Sp": "sp",
        "C2/c": "hpx",
        "Pl": "pl",
        "Pv": "bdg",
        "Ppv": "ppv",
        "NaAl": "nal",
        "CF": "cfer",
        "Wus": "fper",
        "st": "stv",
    }
    bounding_colors = ["#44015a", "#ffffff"]

    modes_filenames = []
    for i_P in range(n_splits_pressure):
        for i_T in range(n_splits_temperature):
            pressure_range = [
                pressure_range_total[0] + i_P * delta_P,
                pressure_range_total[0] + (i_P + 1) * delta_P,
            ]
            temperature_range = [
                temperature_range_total[0] + i_T * delta_T,
                temperature_range_total[0] + (i_T + 1) * delta_T,
            ]
            if n_splits_pressure > 1 or n_splits_temperature > 1:
                split_project_name = f"{project_name}_P{i_P + 1:02d}_T{i_T + 1:02d}"
            else:
                split_project_name = project_name

            modes_filenames.append(split_project_name + "_modes.dat")

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    small_fields = pretty_plot_phase_diagram(
        ax,
        modes_filenames,
        phase_name_replacements,
        bounding_colors,
        n_phases_bounds=[2.0, 7.0],
        smoothing_window=8,
        smoothing_order=1,
        linewidth=0.5,
        label_scaling=3.0,
        label_clearance=0.01,
        number_small_fields=True,
    )

    with open(f"{project_name}_field_ids.txt", "w") as f:
        for i, small_field in enumerate(small_fields):
            line = f"{i+1}: {small_field}"
            print(line)
            f.write(f"{line}\n")

    fig.savefig(f"{project_name}.pdf")
    plt.show()

    print("Processing complete.")
