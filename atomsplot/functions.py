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
import numpy as np
from ase.io import read
from ase.units import Bohr

from atomsplot import ase_custom # monkey patch. pylint: disable=unused-import
from atomsplot.render import render_image
from atomsplot.settings import CustomSettings

logger = logging.getLogger(__name__)


def _deduce_chg_format(filename: str) -> str | None:

    if filename.endswith('.cube'):
        return 'cube'
    elif filename.startswith('CHG'):
        return 'vasp'
    else:
        return None

def _read_charge_file(filename, fmt='cube', upscale : int | None =None):
    """
    Read charge density file (cube of VASP CHGCAR/CHG format).

    Parameters
    ----------
    filename : str
        Path to the file containing the isosurfaces.
    fmt : str, optional
        'cube' or VASP CHGCAR/CHG.
    upscale : int, optional
        Upscale factor for the density grid.

    Returns
    -------
    atoms : ase.Atoms
        Atoms object containing the atomic structure.
    density_grid : np.ndarray
        3D numpy array representing the charge density grid.
    """

    if fmt == 'cube':
        data_dict = read(filename, read_data=True, full_output=True)
        atoms = data_dict["atoms"]
        density_grid = data_dict["data"]

    elif fmt == 'vasp':
        from ase.calculators.vasp import VaspChargeDensity # pylint: disable=import-outside-toplevel
        vcd = VaspChargeDensity(filename)

        atoms = vcd.atoms[0]
        density_grid = np.array(vcd.chg[0]) * (Bohr ** 3) # convert volume in Angstrom^3 to bohr^3

        logging.debug('Charge density grid shape: %s', density_grid.shape)
    else:
        raise ValueError("Unsupported format. Use 'cube' or 'vasp'.")

    if upscale is not None and upscale > 1:
        logging.info('Upscaling charge density grid...')
        try:
            from scipy.ndimage import zoom #pylint: disable=import-outside-toplevel
            density_grid = zoom(density_grid, upscale, order=3)
        except ImportError:
            logging.error('scipy is not installed. Cannot upscale the charge density grid.')

    return atoms, density_grid


def setup_rendering(filename : str,
                    outfile : str | None = None,
                    index : str = '-1',
                    movie : bool = False,
                    framerate : int = 10,
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
    **kwargs : dict
        Additional keyword arguments for rendering
    """

    custom_settings = CustomSettings()

    chg_format =kwargs.pop('chg_format', None)
    if chg_format is None:
        chg_format = _deduce_chg_format(filename)


    if chg_format is not None:
        # read charge density file
        atoms, chg_grid = _read_charge_file(filename=filename,
                                    fmt=chg_format,
                                    upscale=kwargs.pop('chg_upscale', 1))
        kwargs['chg_grid'] = chg_grid
    else:
        if index == '-1' and movie: #if we want to render a movie, we need the whole trajectory
            index = ':'

        atoms = read(filename, index=index)


    label = os.path.splitext(outfile if outfile is not None else os.path.basename(filename))[0]
    logger.info('File was read successfully.')

    # remove None from kwargs
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    if isinstance(atoms, list): # multiple frames

        kwargs['fixed_bounds'] = True

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
            success = False

            # first, try to use ffmpeg:
            ffmpeg_cmd = f'ffmpeg -y -framerate {framerate} '\
                f'-i rendered_frames/{label}_%05d.png '\
                '-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" '\
                f'-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '\
                f'{label}.mp4'
            try:
                ret = subprocess.run([ffmpeg_cmd], check=True, capture_output=True, shell=True)
                success = ret.returncode == 0
            except (subprocess.CalledProcessError, FileNotFoundError) as ffmpeg_error:
                logger.error('ffmpeg failed: %s. trying to use imagemagick...', ffmpeg_error)
                # if ffmpeg fails, try imagemagick:
                try:
                    convert_cmd = f'convert -delay {1000 // framerate} '\
                        f'-loop 0 rendered_frames/{label}_*.png {label}.gif'
                    ret = subprocess.run([convert_cmd], check=True, capture_output=True, shell=True)
                    success = ret.returncode == 0
                except (subprocess.CalledProcessError, FileNotFoundError) as imagemagick_error:
                    logger.error('Imagemagick also failed: %s', imagemagick_error)

            if success:
                logger.info('Movie generated.')
            else:
                logger.error('Error generating movie, '
                        'however the frames are still present in the rendered_frames folder.')

    else: # single frame
        logger.info('Rendering image...')
        render_image(atoms=atoms,
                     outfile=f'{label}.png',
                     custom_settings=custom_settings,
                     **kwargs)
        logger.info('Rendering complete.')

    logger.info('Job done.')
