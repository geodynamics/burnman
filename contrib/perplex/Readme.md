# Perple_X P-T Table Generation

This repository contains scripts and resources to download, install, and use **Perple_X** to generate and read a low-resolution thermodynamic table for an iron-olivine system. It also contains a script to
create a smoothed 1D adiabatic profile for use in ASPECT (https://github.com/geodynamics/aspect/).

These scripts can be used as a template to create, read and process your own tables.

## Contents

- **download_and_install_perplex.sh**  
  A shell script to automatically download and install Perple_X in the local directory.

- **create_lo_res_table.py**  
  A Python script to generate a low-resolution thermodynamic table using Perple_X. This script calls Perple_X programs and configures the run for a simplified system.

- **read_lo_res_table.py**  
  A Python script to read and parse the generated low-resolution table for further
  analysis and plotting.

- **generate_aspect_compatible_1D_adiabat_table.py**
  A Python script to generate a smoothed table of properties along a 1D isentrope using
  a generated table.

- **perplex_utils.py**  
  A python file containing useful Perple_X-related functions. The start of the file defines
  some database dictionaries for the Holland and Powell (2018),
  and Stixrude and Lithgow-Bertelloni (2021, 2024) datasets which can be passed to the other
  functions or used as a basis for your own dictionaries.

## Getting Started

1. **Install Perple_X:**

   ```bash
   ./download_and_install_perplex.sh
   ```

   This will download the latest Perple_X release and compile it locally.

2. **Create the Thermodynamic Table:**

   ```bash
   python create_lo_res_table.py
   ```

   This will use Perple_X to generate the low-resolution table in a newly created project directory.

3. **Read and Analyze the Table:**

   ```bash
   python read_lo_res_table.py
   ```

   This script parses the output table for further use (e.g., evaluating properties or plotting).

4. **Create ASPECT compatible Table:**

   ```bash
   python generate_aspect_compatible_1D_adiabat_table.py
   ```

   This script creates an ASPECT-compatible table.

## Notes

- The resolution and thermodynamic scope in this example are intentionally simplified. You will
  need to change the scripts if you want to increase the resolution. 
- Make sure to review Perple_X license and citation guidelines when using in published work.
