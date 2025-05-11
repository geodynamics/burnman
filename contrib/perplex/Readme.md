# Perple_X P-T Table Generation

This repository contains scripts and resources to download, install, and use **Perple_X** to generate and read a low-resolution thermodynamic table for an iron-olivine system. It can be used as a template to
create and read your own tables.

## Contents

- **download_and_install_perplex.sh**  
  A shell script to automatically download and install Perple_X in the local directory.

- **create_lo_res_table.py**  
  A Python script to generate a low-resolution thermodynamic table using Perple_X. This script calls Perple_X programs and configures the run for a simplified system.

- **read_lo_res_table.py**  
  A Python script to read and parse the generated low-resolution table (e.g., from `table.txt`) for further analysis or plotting.


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

   This script parses the output table for further use (e.g., plotting or machine learning workflows).

## Notes

- This setup is designed for educational or rapid prototyping use; the resolution and thermodynamic scope are intentionally simplified.
- Make sure to review Perple_X license and citation guidelines when using in published work.
