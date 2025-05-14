# Computation of 2D Lagrangian Variables from Velocity Fields

This repository provides a Fortran module, ftlv_fsle_pdua_iv.f90, for computing Lagrangian trajectories, finite-time Lagrangian vorticity (FTLV), and the finite-size Lyapunov exponent (FSLE) from velocity fields.
The script compute_ftlv.py reads and selects the velocity fields before passing them to the Fortran module for processing.

## Installation and Usage

The Fortran module used independently or compiled into a Python-compatible library using f2py with the following command:

f2py -c -m ftlv_fsle_pdua_iv ftlv_fsle_pdua_iv.f90

Once compiled, the module can be imported as a library and used in the compute_fsle.py and compute_ftlv.py scripts.

## Requirements

- Linux OP
- Fortran compiler
- numpy (for f2py functionality)
- Python (compatible with f2py)

## License

## Author

Diego Cortés Morales
