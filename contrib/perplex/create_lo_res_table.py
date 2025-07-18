import os
import shutil

import burnman
from perplex_utils import databases, make_build_file
from perplex_utils import run_vertex, run_pssect
from perplex_utils import create_perplex_table

if __name__ == "__main__":
    # Define the project parameters
    project_name = "iron_olivine_lo_res"
    database = databases["stx24"]
    composition = burnman.Composition(
        {"MgO": 1.8, "FeO": 0.2, "SiO2": 1.0, "Fe": 0.1}, "molar"
    )
    pressure_range = [1.0e5, 10.0e9]
    temperature_range = [200.0, 3000.0]
    n_pressures = 11
    n_temperatures = 15
    outfile = "iron_olivine_lo_res_table.dat"

    perplex_dir = os.path.join(os.getcwd(), "perplex-installer/Perple_X")
    perplex_bindir = os.path.join(os.getcwd(), "perplex-installer/Perple_X/bin")

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
    shutil.copyfile(endmember_file, os.path.join(os.getcwd(), database["data_file"]))
    shutil.copyfile(
        solution_file, os.path.join(os.getcwd(), database["solution_model_file"])
    )
    perplex_option_file = "burnman_perplex_options.dat"

    # Define the number of exploratory nodes.
    # These are the nodes over which the full Gibbs minimization is performed.
    # In the refinement stage, the number of nodes is increased.
    # Any solutions not present on any of the bounding nodes will be excluded
    # from the calculations.
    n_P_exploratory = int((n_pressures + 1) / 2)
    n_T_exploratory = int((n_temperatures + 1) / 2)

    # Create the Perple_X options file
    with open(perplex_option_file, "w") as f:
        f.write(
            "sample_on_grid            F\n"
        )  # Do not force sampling on grid when running werami
        f.write("auto_refine  auto\n")  # Automatically refine the grid
        f.write("grid_levels  1 1\n")  # Do not use adaptive grid refinement
        f.write(
            f"x_nodes {n_P_exploratory} {n_pressures}\n"
        )  # Exploratory and final number of pressure nodes
        f.write(
            f"y_nodes {n_T_exploratory} {n_temperatures}\n"
        )  # Exploratory and final number of temperature nodes

    # Make the build file
    # This file contains the information about the project, such as the
    # composition, the database, and the options for Perple_X.
    print("Making build file...")
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
    print("Running vertex...")
    run_vertex(perplex_bindir, project_name, verbose=False)

    # Run the Perple_X pssect command to create a postscript
    # file containing the phase diagram.
    print("Running pssect...")
    run_pssect(perplex_bindir, project_name, convert_to_pdf=False, verbose=False)

    # Create the BurnMan-readable table from the vertex output
    print("Creating BurnMan-readable table...")
    create_perplex_table(
        perplex_bindir,
        project_name,
        outfile,
        n_pressures,
        n_temperatures,
        pressure_range,
        temperature_range,
    )
    print("Processing complete.")
