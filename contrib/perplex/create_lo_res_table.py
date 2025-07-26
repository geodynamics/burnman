import os
import shutil

import burnman
from perplex_utils import databases, make_build_file
from perplex_utils import run_vertex, run_pssect
from perplex_utils import create_perplex_class_table


def run_perplex(
    perplex_bindir,
    project_name,
    database,
    perplex_option_file,
    composition,
    pressure_range,
    temperature_range,
    n_pressures,
    n_temperatures,
    outfile,
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
    run_vertex(perplex_bindir, project_name, verbose=False)

    # Run the Perple_X pssect command to create a postscript
    # file containing the phase diagram.
    print("* Running pssect...")
    run_pssect(perplex_bindir, project_name, convert_to_pdf=False, verbose=False)

    # Create the BurnMan-readable table from the vertex output
    print("* Creating BurnMan-readable table...")
    create_perplex_class_table(
        perplex_bindir,
        project_name,
        outfile,
        n_pressures,
        n_temperatures,
        pressure_range,
        temperature_range,
    )


if __name__ == "__main__":
    # Define the project parameters
    project_name = "iron_olivine_lo_res"
    database = databases["stx24"]
    composition = burnman.Composition(
        {"MgO": 1.8, "FeO": 0.2, "SiO2": 1.0, "Fe": 0.1}, "molar"
    )
    pressure_range_total = [1.0e5, 10.0e9]
    temperature_range_total = [200.0, 3000.0]
    n_pressures_per_split = 6
    n_temperatures_per_split = 8
    n_splits_pressure = 2
    n_splits_temperature = 2

    outfile_base = "iron_olivine_lo_res_table"

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
    delta_P = (pressure_range_total[1] - pressure_range_total[0]) / n_splits_pressure
    delta_T = (
        temperature_range_total[1] - temperature_range_total[0]
    ) / n_splits_temperature

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
                outfile = f"{outfile_base}_P{i_P + 1:02d}_T{i_T + 1:02d}.dat"
            else:
                split_project_name = project_name
                outfile = f"{outfile_base}.dat"

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
                n_pressures_per_split,
                n_temperatures_per_split,
                outfile,
            )

    # If there are splits, merge the output files into one
    if n_splits_pressure > 1 or n_splits_temperature > 1:
        print("Merging output files...")

        # Read the header lines from the first split file
        first_split_outfile = f"{outfile_base}_P01_T01.dat"
        with open(first_split_outfile, "r") as infile:
            header_lines = [next(infile) for _ in range(13)]

        total_n_pressures = (n_pressures_per_split - 1) * n_splits_pressure + 1
        total_n_temperatures = (n_temperatures_per_split - 1) * n_splits_temperature + 1
        header_lines[6] = f"          {total_n_pressures}\n"
        header_lines[10] = f"          {total_n_temperatures}\n"

        # Load all other lines from the split files into a single list of lists
        all_lines = []
        for i_P in range(n_splits_pressure):
            for i_T in range(n_splits_temperature):
                split_outfile = f"{outfile_base}_P{i_P + 1:02d}_T{i_T + 1:02d}.dat"
                with open(split_outfile, "r") as infile:
                    lines = infile.readlines()[13:]
                    all_lines.extend(lines)

        # Sort the lines by pressure and temperature, looping over pressures first
        all_lines.sort(key=lambda x: (float(x.split()[1]), float(x.split()[0])))

        # Remove duplicate pressures and temperatures from the sorted list
        unique_lines = []
        seen = set()
        for line in all_lines:
            pressure_temp = (float(line.split()[0]), float(line.split()[1]))
            if pressure_temp not in seen:
                seen.add(pressure_temp)
                unique_lines.append(line)

        # Check the total number of unique lines is as expected
        if len(unique_lines) != total_n_pressures * total_n_temperatures:
            raise ValueError(
                f"Expected {total_n_pressures * total_n_temperatures} unique lines, "
                f"but found {len(unique_lines)}."
            )
        # Write the merged output file
        merged_outfile = f"{outfile_base}.dat"
        with open(merged_outfile, "w") as outfile:
            outfile.writelines(header_lines)
            outfile.writelines(unique_lines)
        print(f"Merged output file created: {merged_outfile}")

    print("Processing complete.")
