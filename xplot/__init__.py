# Created on Tue 28 Nov 2023

"""
--------------------------------- atomsplot ------------------------------------

Small program to plot images of atomic structures.

Key features:
- Read atomic structures from various formats (ASE, VASP, etc.)
- Generate high-quality images using PovRay
- Colorcode atoms by forces, coordination numbers, magnetic moments
- Draw arrows to visualize forces and magnetic moments
- Depth cueing to better see the depth using fog
- Generate movies from trajectory files
- Plot charge density isosurfaces

@author: Enrico Pedretti
--------------------------------------------------------------------------------

"""

# ------------------------------ Useful links ------------------------------------
# Repository: https://github.com/Enry99/atomsplot

# Realized as a spin-off of Xsorb, a tool to automate molecular adsorption:
# E. Pedretti, P. Restuccia, M.C. Righi, Comput. Phys. Commun. 291 (2023), 108827
# https://doi.org/10.1016/j.cpc.2023.108827 (https://github.com/Enry99/xsorb)

__version__ = '1.0'
