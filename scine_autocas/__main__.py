# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import sys

# Remove Pyscf warning on libxc
try:
    from pyscf import __config__
    __config__.B3LYP_WITH_VWN5 = False
except ModuleNotFoundError:
    pass

from typing import Any, Dict

from .io import FileHandler, logger, make_parser
from .io.actions import analyze_action, create_action, initcas_action, run_action
from .io.input_handler import InputHandler
from .utils.exceptions import UnkownActionError


def _print_settings(settings: Dict[str, Any], key: str):
    """Print specific settings from user input.

    Parameters
    ----------
    settings : Dict[str, Any]
        settings dict
    key : str
        name of the setting, can be 'Interface', 'AutoCAS', 'Molecule'
    """
    if key in settings:
        print(f"{key}:")
        for setting in settings[key].items():
            print(f"\t{setting[0]}: {setting[1]}")


def print_all_settings(settings: Dict[str, Any]):
    """Print all settings.

    Parameters
    ----------
    settings : Dict[str, Any]
        settings dict
    """
    logger.frame("Settings")
    _print_settings(settings, "Interface")
    logger.empty_line()
    _print_settings(settings, "AutoCAS")
    logger.empty_line()
    _print_settings(settings, "Molecule")
    logger.empty_line()


def main():
    """Main function.

    Setup the logger, cli parser and select the correct actions to run correct workflows.
    """

    parser = make_parser()
    args = vars(parser.parse_args())
    try:
        user_input = InputHandler(args)
    except UnkownActionError:
        print()
        print("*" * 80)
        print("Please specify one of the actions described below.")
        print("*" * 80)
        print()
        parser.print_help()
        sys.exit(1)

    # Do nothing after create
    if user_input.action == "create":
        print("Starting <create> action")
        create_action(user_input)
        return

    # enable logging
    FileHandler.setup_project()
    logger.banner()

    print_all_settings(user_input.get_settings())
    if user_input.action == "run":
        print("Starting <run> action")
        run_action(user_input)
    elif user_input.action == "analyze":
        print("Starting <analyze> action")
        analyze_action(user_input)
    elif user_input.action == "initcas":
        print("Starting <initcas> action")
        initcas_action(user_input)
    else:
        raise UnkownActionError()


if __name__ == "__main__":
    main()
