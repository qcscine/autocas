"""Main module."""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import argparse
import sys

from scine_autocas.main_functions import MainFunctions


def main():
    parser = argparse.ArgumentParser(description="SCINE autoCAS")
    parser.add_argument(
        "-y", "--yaml_input", default=None, type=str, help="The yaml input file."
    )
    parser.add_argument(
        "-x",
        "--xyz_file",
        metavar="xyz_file",
        default=None,
        type=str,
        help="The molecule to process.",
    )
    parser.add_argument(
        "-b",
        "--basis_set",
        type=str,
        default="cc-pvtz",
        help="Basis set for the calculation. Please use the format required by your interface.",
    )
    parser.add_argument(
        "-p", "--plot", action="store_true", help="Creates an entanglement diagram."
    )
    parser.add_argument(
        "-i",
        "--interface",
        default="molcas",
        type=str,
        help="The interface to use for the calculations. Available: Molcas,  Default: Molcas",
    )
    args = vars(parser.parse_args())
    if args["xyz_file"] is None and args["yaml_input"] is None:
        parser.print_help()
        sys.exit(1)
    main_class = MainFunctions()
    main_class.main(args)


if __name__ == "__main__":
    main()
