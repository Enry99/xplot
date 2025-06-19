'''
Custom I/O functions for ASE Atoms objects
'''

from ase.io import read, write

import atomsplot.ase_custom.espresso
import atomsplot.ase_custom.extxyz
import atomsplot.ase_custom.vasp
from atomsplot.ase_custom import AtomsCustom


def ase_custom_read(filename, **kwargs):
    """
    Modified version of ase.io.read that returns an ase_custom.AtomsCustom object
    """
    atoms_or_atoms_list = read(filename=filename, **kwargs)
    if isinstance(atoms_or_atoms_list, list):
        return [AtomsCustom(at) for at in atoms_or_atoms_list]
    else:
        return AtomsCustom(atoms_or_atoms_list)
