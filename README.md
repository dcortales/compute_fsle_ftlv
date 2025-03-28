Computation of Lagrangian Variables from Velocity Fields

This repository provides a Fortran module, ./ftlv_fsle_pdua_iv.f90, for computing Lagrangian trajectories, finite-time Lagrangian vorticity (FTLV), and the finite-size Lyapunov exponent (FSLE) from velocity fields.
Installation and Usage

The Fortran module can be compiled into a Python-compatible library using f2py with the following command:

f2py -c -m ftlv_fsle_pdua_iv ftlv_fsle_pdua_iv.f90

Once compiled, the module can be imported as a library and used in the compute_ftlv.py script.
Requirements

- Linux OP
- Fortran compiler
- numpy (for f2py functionality)
- Python (compatible with f2py)

License

[Include license information if applicable.]
Author

[Your name or institution if needed.]
