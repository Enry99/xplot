'''
Mods to ase library to handle custom labels, and for rendering with POVRAY.
Need to be monkey patched at runtime by importin this module.
'''

#CHANGELOG
#-17 Oct 2024: taken from atomsplot, commit 743d1a9. Update to ase 3.23.0
#-23 Oct 2024: small fix for removing custom labels number in AtomsCustom __init__,
#              added ground fog height and changed pale texture definition,
#              small fix for povray, that was not monkey-patching correctly the
#              write_pov method of POVRAY class.
#-24 Oct 2024: added fix for read_espresso_out to skip initial positions in restarts
#              (taken from script_dinamiche)
#-17 Jun 2025: align with 11 Jun commit to xsorb:
#              -update to ase 3.25.0
#              -explictly marked the custom parts with #### CUSTOM ...
# -19 Jun 2025: improved arrows in povray, removed monkey patching of
#               PlottingVariables, now using directly POVRAY class.

# Runtime patch for read/write
import atomsplot.ase_custom.espresso
import atomsplot.ase_custom.vasp
import atomsplot.ase_custom.extxyz
from atomsplot.ase_custom.atoms import AtomsCustom
from atomsplot.ase_custom.extxyz import write_xyz_custom

# Runtime patch for rendering is not called here but only explicitly
# with rom xsorb.ase_custom import povray when necessary
