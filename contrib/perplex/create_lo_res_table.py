from __future__ import absolute_import
from __future__ import print_function

import os
import shutil

import burnman
from perplex_utils import databases, make_build_file
from perplex_utils import run_vertex, run_pssect
from perplex_utils import create_perplex_table

if __name__ == "__main__":
    project_name = "iron_olivine_lo_res"
    database = databases["stx24"]
    composition = burnman.Composition(
        {"MgO": 1.8, "FeO": 0.2, "SiO2": 1.0, "Fe": 0.1}, "molar"
    )
    pressure_range = [1.0e5, 10.0e9]
    temperature_range = [200.0, 3000.0]
    n_pressures = 11
    n_temperatures = 11
    outfile = "iron_olivine_lo_res_table.dat"

    perplex_dir = os.path.join(os.getcwd(), "perplex-installer/Perple_X")
    perplex_bindir = os.path.join(os.getcwd(), "perplex-installer/Perple_X/bin")

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

    with open(perplex_option_file, "w") as f:
        f.write("sample_on_grid            F")

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

    print("Running vertex...")
    run_vertex(perplex_bindir, project_name, verbose=False)

    print("Running pssect...")
    run_pssect(perplex_bindir, project_name, convert_to_pdf=False, verbose=False)

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
