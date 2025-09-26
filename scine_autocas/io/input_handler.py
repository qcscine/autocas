# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
import os
from typing import Any, Dict, Union

import yaml

from scine_autocas.interfaces.interface import Interface
from scine_autocas.utils.defaults import AvailableActions, AvailableInteraces, Defaults
from scine_autocas.utils.exceptions import InputError, UnkownActionError
from scine_autocas.utils.helper_functions import get_all_params


class InputHandler:
    """Class to convert cli input to dicts.

    Parses the arguments from Argument parser and create options dicts to initialize
    required objects with.

    Attributes
    ----------
    args : Dict[str, Any]
        stores all arguments from Argument parser
    """
    __slots__ = ("args", "action", "_settings", "_yaml_path", "interface")

    def __init__(self, args: Dict[str, Any]):
        """Initialize InputHandler.

        Parameters
        ----------
        args : Dict[str, Any]
            arguments from Argument parser
        """
        self.args = args
        """Arguments from argparser."""
        self.action: str = ""
        """The action to run. See ActionAliases for defined actions."""
        self._settings: Dict[str, Any] = {}
        """Handles all parameters."""
        self._yaml_path: str = ""
        """Path to the yaml file"""
        self.interface: AvailableInteraces
        """name of the interface"""

        # parse the arguments
        self._parse_input()

    def _read_yaml_input(self):
        """Read yaml input and fill settings."""
        input_file = self._yaml_path
        # logging.debug(f"Reading YAML: {self._yaml_path}")
        with open(input_file, encoding="utf-8") as file:
            dict_of_endless_dicts = yaml.load(file, Loader=yaml.FullLoader)

        # first fill defaults
        for key in dict_of_endless_dicts:
            if key == "Defaults":
                for other_key in dict_of_endless_dicts[key]:
                    for bla in dict_of_endless_dicts[key][other_key]:
                        self._settings[other_key][bla] = dict_of_endless_dicts[key][other_key][bla]

        # overwrite defaults
        for key in dict_of_endless_dicts:
            if key != "Defaults":
                for other_key in dict_of_endless_dicts[key]:
                    self._settings[key][other_key] = dict_of_endless_dicts[key][other_key]

    def _get_action(self):
        """Get the action from input."""
        if "action" not in self.args:
            raise UnkownActionError()
        if self.args["action"] is None:
            raise UnkownActionError()
        if not AvailableActions.has(self.args["action"].lower()):
            raise UnkownActionError()
        self.action = self.args["action"].lower()

    def _get_interface(self):
        """Get the interface form input."""
        self.interface = AvailableInteraces.get(self.get_interface_options()["Interface"]["interface"])
        # logging.info(f"Interface: {self.interface}")

    def _check_minimal_input(self):
        """Check that minimum requirements for calculation exist."""
        if self.args["yaml"] is None and self.args["xyz_file"] is None:
            raise InputError("The input requires at least a yaml input or a molecule (xyz-file)")

    def _parse_input(self):
        """Parse input.

        Notes
        -----
        Function remove parsed elements from self.args. At the end every unparsed element is printed,
        to debug unparsed elements.
        """
        # logging.debug(f"All arguments: \n{self.args}")
        # get action, e.g. run, analyze, ...
        self._get_action()

        # with create only a yaml with autocas defaults is created
        if self.action == "create":
            self._settings["output_name"] = self.args["output_name"]
            del self.args["output_name"]
            return
        del self.args["action"]

        self._correct_path_to_file("xyz_file")
        if "molden_file" in self.args:
            self._correct_path_to_file("molden_file")

        # self._strip_defaults_from_cli()

        if self.action in ("run", "initcas"):
            self._check_minimal_input()

        # set defaults
        self._set_defaults()

        # if yaml exists parse yaml and overwrite defaults
        if self.args["yaml"] is not None:
            self._yaml_path = self.args["yaml"]
            self._read_yaml_input()
        try:
            del self.args["yaml"]
        except KeyError:
            pass

        # update settings from CLI and overwrite if they exist
        self._update_settings()

        # for item in self.args.items():
        # logging.warning(f"Unused argument: {item}")

        self._get_interface()

    def _correct_path_to_file(self, param_name: str):
        """Set correct path to xyz file."""
        if self.args[param_name].startswith("/") or self.args[param_name].startswith("~"):
            return
        self.args[param_name] = os.getcwd() + "/" + self.args[param_name]

    def _get_value_from_dict(self, key_str: str) -> Union[bool, Any]:
        """Get value from args.

        Parameters
        ----------
        key_str : str
            name of the argument

        Returns
        -------
        value : Any
            if key exists
        False : bool
            else
        """
        value = False
        try:
            value = self.args[key_str]
        except KeyError:
            return value
        return value

    def get_workflow(self) -> str:
        """Get correct workflow from _settings"""
        if self._settings["AutoCAS"]["large_cas"] and self._settings["Interface"]["n_excited_states"] > 0:
            return "large_cas_excited_states"
        if self._settings["AutoCAS"]["large_cas"]:
            return "large_cas"
        if self._settings["Interface"]["n_excited_states"]:
            return "excited_states"
        return "conventional"

    def _set_defaults(self):
        """Fill _settings with default values for Molecule, AutoCAS, Interface.

        See also
        --------
        scine_autocas.utils.defaults
        """
        default_params = get_all_params(Defaults)["Defaults"]
        self._settings["Molecule"] = default_params["Molecule"]
        self._settings["AutoCAS"] = default_params["AutoCAS"]
        self._settings["Interface"] = default_params["Interface"]
        self._settings["DirName"] = default_params["DirName"]

    def _update_settings(self):
        """Update _settings for Molecule, AutoCAS, Interface."""
        self._overwrite_settings("Molecule")
        self._overwrite_settings("AutoCAS")
        self._overwrite_settings("Interface")

    def _overwrite_settings(self, key: str):
        """Overwrite _settings for a give key, with corresponding values from
        self.args and remove key, value from self.args.

        Parameters
        ----------
        key : str
            Either Molecule, AutoCAS or Interface
        """
        keys = list(self.args)
        interface_settings = Interface.Settings.__slots__
        for val in keys:
            if val in self._settings[key]:
                self._settings[key][val] = self.args[val]
                del self.args[val]
            # xyz file has no default
            elif val in ("molden_file", "xyz_file") and key == "Molecule":
                self._settings[key][val] = self.args[val]
                del self.args[val]

            elif val in interface_settings and key == "Interface":
                self._settings[key][val] = self.args[val]
                del self.args[val]

    def get_molecule_options(self) -> Dict[str, Any]:
        """Get options to create Molecule object.

        Returns
        -------
        molecule_dict: Dict[str, Any]
            Only key is Molecule, with dict as value, which stores all options.
        """
        # logging.info(f"Molecule options: {self._settings['Molecule']}")
        return {"Molecule": self._settings["Molecule"]}

    def get_autocas_options(self) -> Dict[str, Any]:
        """Get options to create AutoCAS object.

        Returns
        -------
        autocas_dict: Dict[str, Any]
            Only key is AutoCAS, with dict as value, which stores all options.
        """
        # logging.info(f"AutoCAS options: {self._settings['AutoCAS']}")
        return {"AutoCAS": self._settings["AutoCAS"]}

    def get_interface_options(self) -> Dict[str, Any]:
        """Get options to create Interface object.

        Returns
        -------
        interface_dict: Dict[str, Any]
            Only key is Interface, with dict as value, which stores all options.
        """
        # logging.info(f"Interface options: {self._settings['Interface']}")
        return {"Interface": self._settings["Interface"]}

    def get_settings(self) -> Dict[str, Any]:
        """Return the full settings dict.

        Returns
        -------
        _settings : Dict[str, Any]
            the settings from yaml, cli and defaults
        """
        return self._settings
