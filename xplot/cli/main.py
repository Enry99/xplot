#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import logging

from xplot import __version__
from xplot.cli.parser import cli_parse

logger = logging.getLogger(__name__)

def main():
    '''
    Main entry point for the xplot command line interface.
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
    from xplot.functions import setup_rendering #pylint: disable=import-outside-toplevel
    setup_rendering(filename=args.filename,
                    index=args.index,
                    outfile=args.output,
                    movie=args.movie,
                    framerate=args.framerate,
                    rotations=args.rotations,
                    supercell=args.supercell,
                    wrap=args.wrap,
                    depth_cueing=args.depth_cueing,
                    range_cut=args.range_cut,
                    cut_vacuum=args.cut_vacuum,
                    colorcode=args.colorcode,
                    ccrange=args.ccrange,
                    arrows=args.arrows,
                    bonds=args.bonds,
                    highlight_mol=args.highlight_mol,
                    chg_format=args.chg_format,
                    chg_file=args.chg_file,
                    chg_iso_threshold=args.chg_iso_threshold,
                    povray=args.povray,
                    width_res=args.width_res,
                )
