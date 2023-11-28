#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Nov 2023

@author: Enrico Pedretti

Functions to render an image and to setup the rendering
"""


from ase.io.pov import get_bondpairs
from ase.io import read, write
import json
import os

#scaling factors for drawing atoms and bonds
ATOMIC_RADIUS_DEFAULT = 0.65
BOND_RADIUS_DEFAULT = 0.8


def render_image(atoms, label = str, 
                 povray : bool = True, 
                 width_res : int = 700, 
                 rotations : str = '-90x', 
                 custom_colors = None):
     
    #cut the vacuum above the top slab
    #cell = atoms.cell.lengths()
    #cell[2] = np.max([atom.z for atom in atoms]) + 2
    #atoms.set_cell(cell, scale_atoms = False)


    #set custom colors if present ################################################
    from ase.data.colors import jmol_colors
    ATOM_COLORS = jmol_colors.copy()

    if custom_colors is not None:
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
    ############################################################################


    if(povray): #use POVray renderer (high quality, CPU intensive)
        config_copy = atoms.copy()
        #config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

        write('{0}.pov'.format(label), 
            atoms, 
            format='pov',
            radii = ATOMIC_RADIUS, 
            rotation=rotations,
            colors=colors,
            povray_settings=dict(canvas_width=width_res, 
                                 celllinewidth=CELLLINEWIDTH, 
                                 transparent=False, 
                                 camera_type='orthographic', 
                                 camera_dist=50., 
                                 bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS)
                                 #camera_type='perspective'
                                )                                
        ).render()
        os.remove('{0}.pov'.format(label))
        os.remove('{0}.ini'.format(label))

    else: # use ASE renderer (low quality, does not draw bonds)
        write(label + '.png', atoms, rotation=rotations, scale = 100, colors=colors)


def start_rendering(filename : str, 
                    index : str = '-1', 
                    movie : bool = False, 
                    framerate : int = 10, 
                    povray : bool = True, 
                    width_res : int = 700, 
                    rotations : str = '90x'
                    ):
    
    atoms = read(filename, index=index) #read file (any format supported by ASE)
    label = os.path.splitext(filename)[0]
    print('File was read successfully.')

    if os.path.isfile("custom_colors.json"):
        with open("custom_colors.json", "r") as f:
            custom_colors = json.load(f)
            print("Custom colors read from file.")
    else:
        custom_colors = None

    if type(atoms) is list:

        print(f'Rendering {len(atoms)} images...')
        for i, atom in enumerate(atoms):
            render_image(atom, label + '_{:05d}'.format(i), povray=povray, width_res=width_res, rotations=rotations, custom_colors=custom_colors)
        print('Rendering complete.')

        if movie:
            print('Generating movie...')
            os.system(f'ffmpeg -framerate {framerate} -i {label}_%05d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {label}.mp4')
            print('Movie generated.')
    else:
        print('Rendering image...')
        render_image(atoms, label, povray=povray, width_res=width_res, rotations=rotations, custom_colors=custom_colors)
        print('Rendering complete.')

    print('Job done.')