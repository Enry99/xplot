#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Nov 2023

@author: Enrico Pedretti

Function to setup and start the rendering, and generate movies
"""

from __future__ import annotations

import os
import shutil
import subprocess
import logging

from tqdm import tqdm
from ase.io import read

from xplot import ase_custom # monkey patch. pylint: disable=unused-import
from xplot.render import render_image

logger = logging.getLogger(__name__)


def setup_rendering(filename : str,
                    outfile : str | None = None,
                    index : str = '-1',
                    movie : bool = False,
                    framerate : int = 10,
                    custom_settings : dict | None = None,
                    **kwargs):
    """
    Setup the rendering of an atomic structure or a trajectory.

    Parameters
    ----------
    filename : str
        Path to the file containing the atomic structure or trajectory.
    outfile : str, optional
        Name of the output file. If not provided, it will be derived from the filename.
    index : str, optional
        Index of the frame to be rendered. Default is '-1' (last frame) for single image,
        and ':' (all frames) for movie.
    movie : bool, optional
        If True, generate a movie from the frames. Default is False.
    framerate : int, optional
        Framerate of the movie (frames per second). Default is 10.
    custom_settings : dict, optional
        Custom settings for rendering, such as colors and styles.
        If not provided, default settings will be used.
    **kwargs : dict
        Additional keyword arguments for rendering, such as
        rotations, supercell, wrapping, depth cueing, range cut,
        color coding, arrows, and PovRay settings.
    """

    if index == '-1' and movie: #if we want to render a movie, we need to read the whole trajectory
        index = ':'

    atoms = read(filename, index=index)
    label = os.path.splitext(outfile if outfile is not None else os.path.basename(filename))[0]
    logger.info('File was read successfully.')

    if isinstance(atoms, list): # multiple frames

        if os.path.exists('rendered_frames'):
            logger.info('Removing old rendered_frames folder...')
            shutil.rmtree('rendered_frames')
        os.makedirs('rendered_frames')
        main_dir = os.getcwd()
        os.chdir('rendered_frames')

        for i, atoms_frame in enumerate(tqdm(atoms, desc='Rendering frames:',)):
            render_image(atoms=atoms_frame,
                         outfile=f'{label}_{i:05d}.png',
                         custom_settings=custom_settings,
                         **kwargs)
        logger.info('Rendering complete.')

        os.chdir(main_dir)

        if movie:
            logger.info('Generating movie...')

            # first, try to use ffmpeg:
            ffmpeg_cmd = ['ffmpeg', '-y', '-framerate', str(framerate), '-i',
                    f'rendered_frames/{label}_%05d.png', '-c:v', 'libx264',
                    '-pix_fmt', 'yuv420p', f'{label}.mp4']
            ret = subprocess.run(ffmpeg_cmd, check=True, capture_output=True)
            success = ret.returncode == 0

            if not success:
                logger.error('ffmpeg failed, trying to use imagemagick...')
                # if ffmpeg fails, try imagemagick:
                ret = subprocess.run(['magick', 'convert', '-delay', str(1000 // framerate),
                                     '-loop', '0', f'rendered_frames/{label}_*.png',
                                     f'{label}.gif'], check=True, capture_output=True)
                success = ret.returncode == 0

                if not success:
                    logger.error('Imagemagick also failed. '\
                          'Please check your installation of ffmpeg or imagemagick.')

            if success:
                logger.info('Movie generated.')
            else:
                logger.error('Error generating movie, '\
                        'however the frames are still present in the rendered_frames folder.')

    else: # single frame
        logger.info('Rendering image...')
        render_image(atoms=atoms,
                     outfile=f'{label}.png',
                     custom_settings=custom_settings,
                     **kwargs)
        logger.info('Rendering complete.')

    logger.info('Job done.')
