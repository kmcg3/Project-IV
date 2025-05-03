# Modelling Traumatic Brain Injury

This project simulates the mechanical behaviour of bundled filaments — relevant to understanding traumatic brain injury — using a Fortran-based LBFGS optimisation algorithm. The model minimises a custom energy functional to determine the equilibrium configuration of elastic rods under stress.

Originally developed by **Chris Prior** and **Jack Panter**, this version includes major modifications for academic research.

> **Note:** Some file paths are hardcoded and may need to be edited to suit your system. Check:
> - [`make_runs.sh`](temp/NewGmin/user/lamina_phase_diagram/make_runs.sh)
> - [`view_coords.py`](temp/NewGmin/user/lamina_phase_diagram/view_coords.py)
> - Any I/O in Fortran source files (especially in `potential.f90` and `data.f90`)

---

## Overview

This simulation:
- Minimises a filament energy functional using LBFGS
- Uses segment angles and lengths as parameters
- Outputs the lowest-energy configuration to a file
- Includes a Python script to visualise results

---

## Key Components

| File | Description |
|------|-------------|
| [`potential.f90`](temp/NewGmin/source/override/lamina_new/potential.f90) | Energy and gradient definitions (core logic and edits) |
| [`data.f90`](temp/NewGmin/user/lamina_phase_diagram/data.f90) | Defines rod structure and input parameters |
| [`make_runs.sh`](temp/NewGmin/user/lamina_phase_diagram/make_runs.sh) | Shell script to compile and run simulation |
| [`view_coords.py`](temp/NewGmin/user/lamina_phase_diagram/view_coords.py) | Python script to plot final filament configuration |
| [`lowest`](temp/NewGmin/lowest) | Output file containing optimised coordinates |

---

## Custom Modifications

This version includes the following changes as part of a dissertation:

- Created a `LENGTHS()` variable to model segment lengths dynamically
- Added new gradient terms in `potential.f90` to support extended energy models
- Introduced a `test_gradient()` function to compare analytical and numerical gradients
- Developed (in progress) a new curvature-based energy and gradient formulation

---

## How to Run

Open a terminal, navigate to the `lamina_phase_diagram` directory, and run:

```bash
./make_runs.sh
```

This will compile the Fortran code and begin the simulation process.

---

## Visualising Results

After the simulation, the optimised configuration is saved in the `lowest` file. You can visualise it with:

```bash
python3 view_coords.py
```

Make sure:
- You are in `lamina_phase_diagram`
- Python 3 and `matplotlib` are installed

---

## Editing the Model

To modify the simulation:

- Use [`potential.f90`](temp/NewGmin/source/override/lamina_new/potential.f90) to change energy/gradient calculations or add new terms
- Use [`data.f90`](temp/NewGmin/user/lamina_phase_diagram/data.f90) to change the number of rods, initial conditions, segment lengths, etc.

Other `.f90` files may also need edits depending on your extensions.

---

## Acknowledgements

- Original code by **Chris Prior** and **Jack Panter**
- Modified and extended for a **dissertation on traumatic brain injury modelling**

---

## License

This project is intended for academic use only. Redistribution is not permitted without permission from the original authors.

