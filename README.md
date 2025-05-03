# Project-IV
# Modelling Traumatic Brain Injury

This project simulates filament mechanics relevant to traumatic brain injury using a Fortran-based LBFGS optimisation algorithm. The simulation models bundles of elastic rods and minimises a custom energy functional to understand filament behaviour under mechanical stress.

Originally developed by Chris Prior and Jack Panter, this version has been adapted and extended for dissertation research.

## Overview

The simulation minimises energy to determine the equilibrium configuration of rod-like filaments. This process involves:

- Computing bending energy and gradients (in `potential.f90`)
- Managing rod and segment data (in `data.f90`)
- Running an LBFGS optimisation loop
- Outputting the minimum-energy configuration (as a `lowest` file)
- Visualising filament geometry (using `viewcoords.py`)

## File Structure

Top-level files and folders:

- `potential.f90` – Custom energy and gradient definitions (main modifications here)
- `data.f90` – Input and structural data definitions
- `make_runs.sh` – Bash script to compile and run the simulation
- `viewcoords.py` – Python script to visualise output from the `lowest` file
- Other `.f90` files – Supporting modules for geometry, memory, etc.
- `lowest` – Output file containing optimised coordinates
- `README.md` – This file

## Running the Simulation

To compile and run the code, execute the following in your terminal:

bash make_runs.sh

This script handles compilation and launches the simulation.

## Visualising Results

Once the simulation finishes, the results will be saved in a file called `lowest`. To visualise the output, run:

python3 viewcoords.py

Make sure you have Python 3 and the `matplotlib` library installed.

## Editing the Model

You can customise the simulation by editing the following:

- `potential.f90` – Modify energy terms and gradient computations
- `data.f90` – Change number of rods, segment lengths, initial angles, or parameters

Other files may require edits if you introduce new features or constraints.

## Acknowledgements

- Original code developed by **Chris Prior** and **Jack Panter**
- This version includes modifications made for a **dissertation on traumatic brain injury modelling**

## License

This project is intended for academic use only. Redistribution is not permitted without permission from the original authors.

