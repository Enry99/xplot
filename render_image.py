#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Nov 2023

@author: Enrico Pedretti

Function to render an image using ASE and povray
"""

import os, sys
from ase.io.pov import get_bondpairs
from ase.io import read, write
import json
import numpy as np

ATOMIC_RADIUS_DEFAULT = 0.65
BOND_RADIUS_DEFAULT = 0.8


def render_image(atoms, label = str, povray : bool = True, witdth_res : int = 700, rotations : str = '-90z, -90x'):

    image_out_filename = label + '.png'
     
    
    if witdth_res is None and povray: 
        witdth_res = 500  # I used 3000 for the Xsorb paper. From 1500 is still quite good. 2000 maybe best compromise (still very high res)
    

    #cut the vacuum above the top slab
    #cell = atoms.cell.lengths()
    #cell[2] = np.max([atom.z for atom in atoms]) + 2
    #atoms.set_cell(cell, scale_atoms = False)


    from ase.data.colors import jmol_colors
    ATOM_COLORS = jmol_colors.copy()

    if os.path.isfile("custom_colors.json"):
        
        with open("custom_colors.json", "r") as f:
            custom_colors = json.load(f)

        #print("Custom colors read from file.")

        USER_COLORS       = custom_colors["atomic_colors"] if "atomic_colors" in custom_colors else []
        ATOMIC_RADIUS     = custom_colors["bond_radius"] if "bond_radius" in custom_colors else ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = custom_colors["bond_radius"] if "bond_radius" in custom_colors else BOND_RADIUS_DEFAULT
        CELLLINEWIDTH     = custom_colors["cell_line_width"] if "cell_line_width" in custom_colors else 0
    else:
        USER_COLORS  = []
        ATOMIC_RADIUS     = ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = BOND_RADIUS_DEFAULT
        CELLLINEWIDTH     = 0

    for color in USER_COLORS:

        ATOM_COLORS[color[0]] = color[1]

    colors = [ ATOM_COLORS[atom.number] for atom in atoms]

    #fading color for lower layers in top view
    #zmax = max([atom.z for atom in atoms])
    #zmin = min([atom.z for atom in atoms])
    #delta = zmax - zmin
    #colors_top = [ (colors[atom.index] + (np.array([1,1,1]) - colors[atom.index])*(zmax - atom.z)/delta).round(4) for atom in atoms ]


    if(povray): #use POVray renderer (high quality, CPU intensive)
        config_copy = atoms.copy()
        config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

        write('{0}.pov'.format(label), 
            atoms, 
            format='pov',
            radii = ATOMIC_RADIUS, 
            rotation=rotations,
            colors=colors,
            povray_settings=dict(canvas_width=witdth_res, celllinewidth=CELLLINEWIDTH, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS))
            #camera_type='perspective'
        ).render()
        os.remove('{0}.pov'.format(label))
        os.remove('{0}.ini'.format(label))

    else: # use ASE renderer (low quality, does not draw bonds)
        write(image_out_filename, atoms, rotation=rotations, scale = 100, colors=colors)



def read_file(filename : str, index : str = '-1', movie = False, framerate = 10):
    
    atoms = read(filename, index=index) #read xyz file
    label = filename.split('.xyz')[0]


    if type(atoms) is list:
        for i, atom in enumerate(atoms):
            render_image(atom, label + '_{:04d}'.format(i))
        if movie:
            os.system(f'ffmpeg -framerate {framerate} -i {label}_%04d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {label}.mp4')
    else:
        render_image(atoms, label)


if (len(sys.argv) > 2):
    read_file(sys.argv[1], sys.argv[2], sys.argv[3] == 'movie' if len(sys.argv) > 3 else False, framerate=5)
else:
    read_file(sys.argv[1])