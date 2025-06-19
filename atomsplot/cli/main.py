#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging

from atomsplot import __version__
from atomsplot.cli.parser import cli_parse

logger = logging.getLogger(__name__)

def main():
    '''
    Main entry point for the atomsplot command line interface.
    '''

    # set up logging to print only to console
    logging.basicConfig(level=logging.INFO,
                        format='%(message)s')

    logger.info('atomsplot v%s', __version__)

    if len(sys.argv) == 1:
        logger.info("No command provided. The program will terminate.")
        return

    # parse command line arguments
    args = cli_parse()

    # import here to not impact time to display cli help
    from atomsplot.functions import setup_rendering #pylint: disable=import-outside-toplevel
    setup_rendering(filename=args.filename,
                    index=args.index,
                    outfile=args.output,
                    movie=args.movie,
                    framerate=args.framerate,
                    rotations=args.rotations,
                    supercell=args.supercell,
                    # repeat_slab=args.repeat_slab,
                    # center_molecule=args.center_molecule,
                    wrap=args.wrap,
                    hide_cell=args.hide_cell,
                    depth_cueing=args.depth_cueing,
                    range_cut=args.range_cut,
                    cut_vacuum=args.cut_vacuum,
                    colorcode=args.colorcode,
                    ccrange=args.ccrange,
                    arrows=args.arrows,
                    arrows_scale=args.arrows_scale,
                    bonds=args.bonds,
                    #highlight_mol=args.highlight_mol,
                    chg_format=args.chg_format,
                    chg_iso_threshold=args.chg_iso_threshold,
                    chg_upscale=args.chg_upscale,
                    povray=not args.no_povray,
                    width_res=args.width_res,
                    fixed_bounds=args.fixed_bounds
                )
