'''
Mods to ase library to handle custom labels, and for rendering with POVRAY.
Need to be monkey patched at runtime by importin this module.
'''

#CHANGELOG
#-17 Oct 2024: taken from xplot, commit 743d1a9. Update to ase 3.23.0
#-23 Oct 2024: small fix for removing custom labels number in AtomsCustom __init__,
#              added ground fog height and changed pale texture definition,
#              small fix for povray, that was not monkey-patching correctly the
#              write_pov method of POVRAY class.

# Runtime patch for read/write
import xplot.ase_custom.espresso
import xplot.ase_custom.vasp
import xplot.ase_custom.xyz
from xplot.ase_custom.atoms import AtomsCustom
from xplot.ase_custom.xyz import write_xyz_custom

# Runtime patch for rendering is not called here but only explicitly
# with from xsorb.ase_custom import povray when necessary
