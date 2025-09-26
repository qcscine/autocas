# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
import os
from typing import Optional

import yaml

from scine_autocas.io.input_handler import InputHandler
from scine_autocas.utils.defaults import Defaults
from scine_autocas.utils.helper_functions import get_all_params


def create_action(input_handler: Optional[InputHandler] = None):
    """Create the base yaml file.

    Parameters
    ----------
    input_handler : Optional[InputHandler]
        class to handle the input from the cli.

    Notes
    -----
    if no arguement, default yaml name is 'autocas.yaml'
    """
    if input_handler:
        yaml_file_name = input_handler.get_settings()["output_name"]
    else:
        yaml_file_name = "autocas.yaml"
    defaults_dict = get_all_params(Defaults)
    del defaults_dict["Defaults"]["DirName"]
    del defaults_dict["Defaults"]["QcMaquisNames"]
    del defaults_dict["Defaults"]["PlotNames"]
    if os.path.isfile(yaml_file_name):
        print(f"Autocas config {yaml_file_name} exists.")
        print("File is NOT overwritten, please select a new name to create a new config file")
        return

    print(f"Storing autocas config at {yaml_file_name}")
    with open(yaml_file_name, "w") as file:
        yaml.dump(defaults_dict, file)


if __name__ == "__main__":
    create_action()
