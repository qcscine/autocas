"""Handle input files for autoCAS.

This module provides the InputHandler class, which parses and modifies
yaml input files.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from typing import Any, Dict

import yaml

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces import Interface
from scine_autocas.interfaces.molcas import Molcas


class InputHandler:
    """A class to handle the input for autoCAS.

    AutoCASs input is formated in the yaml format, so that the provided input
    is modified into a dictionary by the yaml library, which can be used easiely to
    store the current state of an autoCAS calculation or set up a calculation from an
    input.

    Attributes
    ----------
    input_file : str
        location of the yaml input file
    settings_dir : Dict[str, Any]
        the dict created from the yaml input
    interface : Interface
        the interface providing electronic structure caluclations
    molecule : Molecule
        the molecule object, created from the yaml input
    autocas : Autocas
        the Autocas object, created from the yaml input
    large_cas : bool
        flag to enable the large active space protocol
    """

    __slots__ = (
        "input_file",
        "settings_dir",
        "interface",
        "molecule",
        "autocas",
        "large_cas",
    )

    def __init__(self, yaml_input: str):
        """Construct the InputHandler.

        The constructor directly read the provided yaml input and sets
        all class attributes.

        Parameters
        ----------
        yaml_input : str
            path to the yaml input file
        """
        self.input_file: str
        """path to the yaml input"""
        self.settings_dir: Dict[str, Any]
        """all settings provided by the yaml input are stored in a dict"""
        self.interface: Interface
        """providing the electronic structure software"""
        self.molecule: Molecule
        """stores molecular system information"""
        self.autocas: Autocas
        """handles the active space search"""
        self.large_cas: bool = False
        """enables the large active space protocol"""
        if yaml_input != "":
            self.read_input(yaml_input)

    def print_settings(self):
        """Print all settings from yaml file"""
        for main_settings in self.settings_dir:
            if not isinstance(self.settings_dir[main_settings], dict):
                print(f"{main_settings}: {self.settings_dir[main_settings]}")
                print(f"{'-' * (len(main_settings) + 1)}")
                print("")
            else:
                print(f"{main_settings}:")
                print(f"{'-' * (len(main_settings) + 1)}")
                settings = yaml.dump(self.settings_dir[main_settings])
                print(settings)

    def read_input(self, input_file: str):
        """Read yaml input and set up class attributes.

        Parameters
        ----------
        input_file : str
            path to the yaml input file
        """
        self.input_file = input_file
        with open(self.input_file, encoding="utf-8") as file:
            self.settings_dir = yaml.load(file, Loader=yaml.FullLoader)

    def get_molecule(self) -> Molecule:
        """Create molecule object from input.

        Returns
        -------
        molecule : Molecule
            the molecule object, based on provided input file.
        """
        molecule_settings = self.settings_dir["molecule"]
        self.molecule = Molecule(settings_dict=molecule_settings)
        return self.molecule

    def get_autocas(self) -> Autocas:
        """Set up Autocas object from input.

        Returns
        -------
        autocas : Autocas
            the Autocas object, based on provided input file.

        Notes
        -----
        The autocas object relies on information of the molecule.
        """
        if self.molecule is None:
            self.get_molecule()
        try:
            autocas_settings = self.settings_dir["autocas"]
        except KeyError:
            autocas_settings = None
        self.autocas = Autocas(molecule=self.molecule, settings_dict=autocas_settings)
        if autocas_settings is not None:
            if "large_cas" in autocas_settings:
                self.large_cas = autocas_settings["large_cas"]
        return self.autocas

    def get_interface(self) -> Interface:
        """Set up interface from input.

        Returns
        -------
        interface : Interface
            an object, which is based on a inherited class of Interface, based on the provided input.

        Notes
        -----
        The interface object relies on information of the molecule.
        """
        if self.molecule is None:
            self.get_molecule()
        interface_settings = self.settings_dir["interface"]
        self.interface = Molcas(
            molecules=[self.molecule], settings_dict=interface_settings
        )
        return self.interface
