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
    from xplot.settings import read_custom_settings #pylint: disable=import-outside-toplevel

    #try to read custom settings
    custom_settings = read_custom_settings()
    if custom_settings.get("povray_old_style", None):
        os.environ['POVRAY_OLD_STYLE'] = '1'


    setup_rendering(filename = args.filename,
                    index = args.index,
                    movie = args.movie,
                    framerate = args.framerate,
                    custom_settings = custom_settings,
                    outfile = args.output,
                    rotations = args.rotations,
                    supercell = args.supercell,
                    wrap=args.wrap,
                    depth_cueing=args.depth_cueing,
                    range_cut=args.range_cut,
                    cut_vacuum=args.cut_vacuum,
                    colorcode=args.colorcode,
                    ccrange=args.ccrange,
                    arrows=args.arrows,
                    nobonds = args.nobonds,
                    povray = not args.nopov,
                    width_res = args.width_resolution,
                )
