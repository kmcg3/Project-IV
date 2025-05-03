# Modelling Traumatic Brain Injury

This project simulates the behaviour of filament bundles — relevant to understanding traumatic brain injury — using a Fortran-based LBFGS algorithm. The model minimises energy functions to determine the equilibrium configuration of elastic rods under loads.

Originally developed by **Chris Prior** and **Jack Panter**, this version includes modifications for extensibility of the model and introduces some new energies and gradients.

> **Note:** Some file paths are hardcoded and may need to be edited.

---

## Overview

This simulation:
- Minimises a filament bundle's energy function using the LBFGS algorithm;
- Uses angles and lengths as parameters;
- Outputs the lowest energy configuration;
- Includes a Python script to visualise results.

---
**Some of the files involved in the simulation:**
| File | Description |
|------|-------------|
| [`potential.f90`](temp/NewGmin/source/override/lamina_new/potential.f90) | Energy and gradient definitions |
| [`data.f90`](temp/NewGmin/user/lamina_phase_diagram/data.f90) | inputs data for [`data.in`](temp/NewGmin/user/lamina_phase_diagram/template/data.in) |
| [`make_runs.sh`](temp/NewGmin/user/lamina_phase_diagram/make_runs.sh) | Shell script to compile and run simulation, can edit rod number, lengths, etc. here |
| [`view_coords.py`](temp/NewGmin/user/lamina_phase_diagram/view_coords.py) | Python script to visualise the bundle |
| [`lowest`](temp/NewGmin/lowest) | Output file containing lowest energy coordinates |

---

## Modifications

This version includes the following changes as part of a dissertation:

- Changes to the files mentioned above and multiple other `f.90` files;
- Introduced the `LENGTHS()` variable;
- Added new energy and gradient terms to `potential.f90`;
- Developing (in progress) a new curvature for the bending energy and gradient.

---

## How to Run

Open a terminal, navigate to the `lamina_phase_diagram` directory, and run:

```bash
./make_runs.sh
```

This will compile the Fortran code and begin the simulation process.

---

## Visualising Results

After the simulation, the lowest energy configuration is saved in the `lowest` file. You can visualise it with:

```bash
python3 view_coords.py
```

Make sure:
- You are in `lamina_phase_diagram`;
- File paths in `view_coords.py` are correct;
- You are using the correct `lowest` file generated through running `make_runs.sh`.

---

## Acknowledgements

- Original code by **Chris Prior** and **Jack Panter**;
- Modified and extended for undergraduate work on  **modelling traumatic brain injury**.

---

## License

This project is intended for academic use only. Redistribution is not permitted without permission from the original authors.

