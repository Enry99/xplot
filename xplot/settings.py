'''
General settings and rendering parameters class for the xplot package.
'''

import os
import json
import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
from ase.data.colors import jmol_colors, cpk_colors

logger = logging.getLogger(__name__)


ATOMIC_RADIUS_DEFAULT = 0.7
BOND_RADIUS_DEFAULT = 0.8
BOND_LINE_WIDTH_DEFAULT = 0.15
CELL_LINE_WIDTH_DEFAULT = 0.1

# VESTA colors for elements, normalized to [0, 1] range
vesta_colors = np.array([
    [255, 0, 0],    # None - Placeholder for no element
    [112, 171, 250],  # Ac
    [192, 192, 192],  # Ag
    [129, 178, 214],  # Al
    [84, 92, 242],    # Am
    [207, 254, 196],  # Ar
    [116, 208, 87],   # As
    [117, 79, 69],    # At
    [255, 209, 35],   # Au
    [31, 162, 15],    # B
    [0, 201, 0],      # Ba
    [94, 215, 123],   # Be
    [224, 0, 56],     # Bh
    [158, 79, 181],   # Bi
    [138, 79, 227],   # Bk
    [126, 49, 2],     # Br
    [76, 76, 76],     # C
    [90, 150, 189],   # Ca
    [255, 217, 143],  # Cd
    [255, 255, 199],  # Ce
    [161, 54, 212],   # Cf
    [49, 252, 2],     # Cl
    [120, 92, 227],   # Cm
    [0, 0, 175],      # Co
    [0, 0, 158],      # Cr
    [87, 23, 143],    # Cs
    [34, 71, 220],    # Cu
    [209, 0, 79],     # Db
    [31, 255, 199],   # Dy
    [0, 230, 117],    # Er
    [179, 31, 212],   # Es
    [97, 255, 199],   # Eu
    [176, 185, 230],  # F
    [181, 113, 0],    # Fe
    [179, 31, 186],   # Fm
    [66, 0, 102],     # Fr
    [158, 227, 115],  # Ga
    [69, 255, 199],   # Gd
    [126, 110, 166],  # Ge
    [255, 204, 204],  # H
    [252, 232, 206],  # He
    [77, 194, 255],   # Hf
    [184, 184, 208],  # Hg
    [0, 255, 156],    # Ho
    [230, 0, 46],     # Hs
    [148, 0, 148],    # I
    [166, 117, 115],  # In
    [23, 84, 135],    # Ir
    [161, 33, 246],   # K
    [250, 193, 243],  # Kr
    [90, 196, 73],    # La
    [134, 223, 115],  # Li
    [199, 0, 102],    # Lr
    [0, 171, 36],     # Lu
    [179, 13, 166],   # Md
    [251, 123, 21],   # Mg
    [167, 8, 157],    # Mn
    [84, 181, 181],   # Mo
    [235, 0, 38],     # Mt
    [176, 185, 230],  # N
    [249, 220, 60],   # Na
    [115, 194, 201],  # Nb
    [199, 255, 199],  # Nd
    [254, 55, 181],   # Ne
    [183, 187, 189],  # Ni
    [189, 13, 135],   # No
    [0, 128, 255],    # Np
    [254, 3, 0],      # O
    [38, 102, 150],   # Os
    [192, 156, 194],  # P
    [0, 161, 255],    # Pa
    [87, 89, 97],     # Pb
    [0, 105, 133],    # Pd
    [163, 255, 199],  # Pm
    [171, 92, 0],     # Po
    [217, 255, 199],  # Pr
    [208, 208, 224],  # Pt
    [0, 107, 255],    # Pu
    [0, 125, 0],      # Ra
    [112, 46, 176],   # Rb
    [38, 125, 171],   # Re
    [204, 0, 89],     # Rf
    [10, 125, 140],   # Rh
    [66, 130, 150],   # Rn
    [36, 143, 143],   # Ru
    [255, 250, 0],    # S
    [158, 99, 181],   # Sb
    [181, 99, 171],   # Sc
    [154, 239, 15],   # Se
    [217, 0, 69],     # Sg
    [27, 59, 250],    # Si
    [143, 255, 199],  # Sm
    [154, 142, 185],  # Sn
    [0, 255, 0],      # Sr
    [77, 166, 255],   # Ta
    [48, 255, 199],   # Tb
    [59, 158, 158],   # Tc
    [212, 122, 0],    # Te
    [0, 186, 255],    # Th
    [120, 202, 255],  # Ti
    [166, 84, 77],    # Tl
    [0, 212, 82],     # Tm
    [0, 143, 255],    # U
    [229, 25, 0],     # V
    [33, 148, 214],   # W
    [66, 158, 176],   # Xe
    [148, 255, 255],  # Y
    [0, 191, 56],     # Yb
    [143, 143, 129],  # Zn
    [0, 255, 0]       # Zr
]) / 255.0


@dataclass
class CustomSettings:
    """Configuration settings for rendering atoms
    """

    atomic_colors : dict = {}
    molecule_colors : dict = {}
    color_scheme : np.ndarray = jmol_colors.copy()  # Use ASE's Jmol colors by default

    mol_indices : Optional[list[int]] = None  # Indices of molecules to highlight
    nontransparent_atoms : list[int] = []

    atomic_radius : float = ATOMIC_RADIUS_DEFAULT
    bond_radius : float = BOND_RADIUS_DEFAULT
    bond_line_width : float = BOND_LINE_WIDTH_DEFAULT
    cell_line_width : float = CELL_LINE_WIDTH_DEFAULT


    def __post_init__(self):
        """Post-initialization to set custom settings."""

        custom_settings_path = "image_settings.json"
        if os.path.isfile("image_settings.json"):
            with open(custom_settings_path, "r") as f:
                custom_settings : dict = json.load(f)
                logger.info("Custom colors read from %s file.", custom_settings_path)

            if custom_settings.pop("povray_old_style", None):
                os.environ['POVRAY_OLD_STYLE'] = '1'

            color_scheme = custom_settings.pop("color_scheme", "jmol").lower()
            if color_scheme == "vesta":
                self.color_scheme = vesta_colors.copy()
                logger.info("Using VESTA colors for elements.")
            elif color_scheme == "cpk":
                self.color_scheme = cpk_colors.copy()
                logger.info("Using CPK colors for elements.")
            elif color_scheme == "jmol":
                self.color_scheme = jmol_colors.copy()
                logger.info("Using Jmol colors for elements.")
            else:
                logger.warning("Unknown color scheme '%s'. Using Jmol colors by default.", color_scheme)

            for key, value in custom_settings.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    logger.warning("Custom setting '%s' not recognized.", key)
