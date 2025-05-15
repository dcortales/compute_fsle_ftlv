# Computation of 2D Lagrangian Variables from Velocity Fields

This repository provides a Fortran module (`ftlv_fsle_pdua_iv.f90`) for computing Lagrangian trajectories, finite-time Lagrangian vorticity (FTLV), and the finite-size Lyapunov exponent (FSLE) from velocity fields.
The scripts `compute_fsle.py` and `compute_ftlv.py` read and select the velocity fields before passing them to the Fortran module for processing and save snapshots of the newly created lagrangian variables.

## Installation and Usage

The Fortran module used independently or compiled into a Python-compatible library using `f2py`. To compile with `f2py`, run:

```bash
f2py -c -m ftlv_fsle_computation ftlv_fsle_computation.f90
```
or can be simply downloaded within this package: `ftlv_fsle_ekst_computation.cpython-312-x86_64-linux-gnu.so`

Once compiled, the module can be imported as a library and used in `compute_fsle.py` and `compute_ftlv.py` scripts.

## Requirements

- Linux OP
- Fortran compiler
- numpy (for f2py functionality)
- Python (compatible with f2py)

## Authors

Diego Cortés Morales

Ismael Hérnandez Carrasco
