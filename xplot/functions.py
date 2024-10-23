#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Nov 2023

@author: Enrico Pedretti

Functions to render an image and to setup the rendering
"""

import os

from tqdm import tqdm
from ase.io import read

from xplot import ase_custom
from xplot.render import render_image



def start_rendering(filename : str,
                    outfile : str = None,
                    index : str = '-1',
                    movie : bool = False,
                    framerate : int = 10,
                    custom_settings : dict = None,
                    **kwargs
                    ):


    if index == '-1' and movie: #if we want to render a movie, we need to read the whole trajectory
        index = ':'
    atoms = read(filename, index=index) #read file (any format supported by ASE)
    label = os.path.splitext(outfile if outfile is not None else os.path.basename(filename))[0]
    print('File was read successfully.')

    if type(atoms) is list:

        print(f'Rendering {len(atoms)} images...')

        os.makedirs('rendered_frames', exist_ok=True)
        main_dir = os.getcwd()
        os.chdir('rendered_frames')

        for i, atoms_frame in enumerate(tqdm(atoms)):
            render_image(atoms=atoms_frame,
                         outfile=f'{label}_{i:05d}.png',
                         custom_settings=custom_settings,
                         **kwargs)
        print('Rendering complete.')

        os.chdir(main_dir)

        if movie:
            print('Generating movie...')
            success = os.system(f'ffmpeg -framerate {framerate} '\
                        f'-i rendered_frames/{label}_%05d.png  '\
                        '-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" '\
                        f'-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '\
                        f'{label}.mp4')

            if success:
                print('Movie generated.')
            else:
                print('Error generating movie, '\
                        'however the frames are still present in the  rendered_frames folder.')
    else:
        print('Rendering image...')
        render_image(atoms=atoms,
                     outfile=f'{label}.png',
                     custom_settings=custom_settings,
                     **kwargs)
        print('Rendering complete.')

    print('Job done.')
