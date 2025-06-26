'''
Argument parser for atomsplot command line interface.
'''

import argparse

import atomsplot

def _positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not positive")
    return ivalue

def _positive_float(value):
    fvalue = float(value)
    if fvalue < 1e-10:
        raise argparse.ArgumentTypeError(f"{value} is not positive")
    return fvalue

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    '''
    Combine the formatters to have a more informative help message.
    '''

def cli_parse():
    """
    Parse command line arguments for atomsplot.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """

    parser = argparse.ArgumentParser(
        prog='atomsplot',
        description=atomsplot.__doc__,
        formatter_class=CustomFormatter,
        allow_abbrev=False)

    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+atomsplot.__version__)

    # input/output options
    parser.add_argument('filename',
                        help='Can be in any format readable by ASE')
    parser.add_argument('-i','--index',
                        type=str,
                        default='-1',
                        help="Index of the frame to be rendered in a trajectory, "\
                        "e.g. '0' (first), '-1' (last), ':' (all), '::10' (one every 10 frames).")
    parser.add_argument('-o','--output',
                        type=str,
                        help='Output file name. Default = [filename].png')

    # structure options
    parser.add_argument('-r','--rotations',
                        type=str,
                        default='0x,0y,0z',
                        help="List of rotations for the visualization, where 0x,0y,0z = top view "\
                            "Presets for front view: front (= -90x), front2 (90z,-90x). "\
                            "If the first rotation has a negative angle, preceed it "
                            "with a dummy rotation, e.g. 0z,-90x ")
    parser.add_argument('-s','--supercell',
                    nargs = 3,
                    type=_positive_int,
                    metavar=('nx', 'ny', 'nz'),
                    help="Replicate the cell nx ny nz times along the three cell vectors.")
    # parser.add_argument('rs', '--repeat-slab',
    #                     action='store_true',
    #                     default=False,
    #                     help='Replicate only the slab, not the molecule.')
    # parser.add_argument('-cmol', '--center-molecule',
    #                     action='store_true',
    #                     default=False,
    #                     help='Center the molecule in the middle of the slab.')
    parser.add_argument('-wr', '--wrap',
                        action='store_true',
                        default=False,
                        help='Wrap atoms according to pbc.')
    parser.add_argument('-rc','--range-cut',
                    nargs=2,
                    type=float,
                    metavar=('zmin', 'zmax'),
                    help='range to be displayed in the z direction.')
    parser.add_argument('-cv','--cut-vacuum',
                    action='store_true',
                    help='Cut vacuum above and below the slab (avoid white empty region).')
    parser.add_argument('-b','--bonds',
                        type=str,
                        choices=['none', 'single', 'multiple'],
                        default='single',
                        help='Draw bonds between atoms.')

    # color and style options
    parser.add_argument('-hc', '--hide-cell',
                    action='store_true',
                    default=False,
                    help='Hide the cell box.')
    parser.add_argument('-dc','--depth-cueing',
                    nargs='?',
                    type=_positive_float,
                    const=1.0,
                    help='Enable depth cueing. Optional parameter: intensity (>0).')
    parser.add_argument('-fo', '--fog-offset',
                    type=float,
                    nargs='?',
                    const=2.0,
                    help='''Offset in Angstrom for where depth cueing starts,
                    starting from the camera position (level of the topmost atom).
                    If the argument is used with no value, the default is 2.''')
    # parser.add_argument('-hm', '--highlight-mol',
    #                 action='store_true',
    #                 default=False,
    #                 help='Highlight the molecule (draw atoms in a different color).')
    parser.add_argument('-cc', '--colorcode',
                    type=str,
                    choices=['forces', 'magmoms', 'coordnum'],
                    help='''Color atoms according to a property
                    (e.g. forces, magnetic moments, coordination number).''')
    parser.add_argument('--ccrange',
                    nargs=2,
                    type=float,
                    metavar=('min', 'max'),
                    help='range for the colorcode.')
    parser.add_argument('-arr', '--arrows',
                    type=str,
                    choices=['forces', 'magmoms'],
                    help='''Draw arrows representing the vectors,
                    with lenghth proportional to the magnitude.''')
    parser.add_argument('-as', '--arrows-scale',
                    type=_positive_float,
                    default=1.0,
                    help='''Scale factor for the arrows.
                    Default = 1.0 (no scaling).''')
    parser.add_argument('-chgfm', '--chg-format',
                    type=str,
                    choices=['cube', 'vasp'],
                    help='''Format of the charge density file.
                    Options: 'cube' or 'vasp' (CHGCAR/CHG).''')
    parser.add_argument('-chgu', '--chg-upscale',
                    type=_positive_int,
                    help='''Upscale the charge density grid by this factor.
                    (1 = no upscaling).''')
    parser.add_argument('-iso', '--chg-iso-threshold',
                    type=_positive_float,
                    help='''Iso-surface threshold for the charge density.
                    If not specified, VESTA default is used''')

    # rendering options
    parser.add_argument('-w', '--width-res',
                    type=_positive_int,
                    default=700,
                    help='Horizontal resolution in pixels.')
    parser.add_argument('-fb', '--fixed-bounds',
                    action='store_true',
                    default=False,
                    help='Fix the canvas white space to the unit cell. Always true for trajs.')
    parser.add_argument('-nopov','--no-povray',
                    action='store_true',
                    default=False,
                    help='Use povray for rendering (much better quality).')

    # movie options
    parser.add_argument('-m','--movie',
                    action='store_true',
                    default=False,
                    help='Create movie from the frames.')
    parser.add_argument('-f', '--framerate',
                    type=float,
                    default=10.0,
                    help='Framerate of the movie (frames per second).')



    args = parser.parse_args()

    if args.rotations == 'front': args.rotations = '-90x' # pylint: disable=multiple-statements
    if args.rotations == 'front2': args.rotations = '90z,-90x' # pylint: disable=multiple-statements

    return args
