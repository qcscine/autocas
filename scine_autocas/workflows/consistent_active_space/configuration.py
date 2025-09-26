# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import os.path
from typing import List, Dict
import yaml
from scine_autocas.utils.defaults import Defaults
from scine_autocas.io import FileHandler


class ConsistentActiveSpaceConfiguration:
    """
    Configuration class for the consistent active space workflow.
    """
    __slots__ = (
        "load_orbitals",
        "cas_method",
        "autocas_indices",
        "basis_set",
        "large_active_space",
        "exclude_core",
        "unmappable",
        "xyz_files",
        "system_names",
        "base_load_path",
        "project_name"
    )

    def __init__(self):
        self.load_orbitals: bool = False
        """
        bool
            If true, the orbitals are loaded from the Serenity HDF5 files.
        """
        self.cas_method: str = "CASPT2"
        """
        str
            The post CAS method to be used. Default is "CASPT2".
        """
        self.autocas_indices: List[int] = []
        """
        List[int]
            The system indicies to run the autoCAS workflow on. If None, all systems are used.
        """
        self.basis_set: str = Defaults.Interface.basis_set
        """
        str
            The AO basis set set.
        """
        self.large_active_space: bool = False
        """
        bool
            If true, the large active space workflow is used.
        """
        self.exclude_core: bool = True
        """
        bool
            If true, core orbitals are always excluded from the active space.
        """
        self.unmappable: bool = False
        """
        bool
            If true, unmappable orbitals are always included in the active space.
        """
        self.xyz_files: List[str] = []
        """
        List[str]
            The structure files in XYZ format to be used for the workflow.
        """
        self.system_names: List[str] = []
        """
        List[str]
            The system names to be used for the workflow.
        """
        self.base_load_path: str = os.getcwd()
        """
        str
            The base path where the Serenity orbitals are loaded from. This is used when load_orbitals is True.
            By default, it is set to the current directory.
        """
        self.project_name: str = Defaults.DirName.project_name
        """
        str
            The project name. Serenity and Molcas files will be written to this directory.
        """

    def write_yaml_file(self, file_name: str = "consistent_cas.configuration.yaml") -> str:
        """
        Write the configuration to a YAML file.

        Returns
        -------
        str
            The configuration is written to a file in the current directory.
            The filename is returned.
        """
        config_dict = {attr: getattr(self, attr) for attr in self.__slots__}
        with open(file_name, 'w') as file:
            yaml.dump(config_dict, file, default_flow_style=False)
        return file_name

    @staticmethod
    def input_sanity_checks(config) -> None:
        """
        Perform sanity checks on the input configuration and complete missing/incomplete fields.

        Parameters
        ----------
        config : ConsistentActiveSpaceConfiguration
            The configuration object to be checked and completed.
        """
        n_systems = max(len(config.xyz_files), len(config.system_names))
        if n_systems < 2:
            raise ValueError("At least two XYZ files or system names are required for the consistent active"
                             " space workflow.")
        if not config.autocas_indices:
            config.autocas_indices = list(range(n_systems))

        if max(config.autocas_indices) > n_systems - 1:
            raise ValueError("The autocas_indices must not exceed the number of XYZ files/system names minus"
                             " one.")
        if not config.system_names:
            config.system_names = [f"system_{i}" for i in range(n_systems)]
        if config.load_orbitals:
            config.xyz_files = [os.path.join(config.base_load_path, name, name + ".xyz")
                                for name in config.system_names]
        if len(config.xyz_files) != len(config.system_names):
            raise ValueError("The number of XYZ files must match the number of system names.")
        for p in config.xyz_files:
            if not os.path.isfile(p):
                raise FileNotFoundError(f"The file {p} does not exist.")

    @classmethod
    def from_options(cls, options, xyz_files: List[str]):
        """
        Create a configuration from command line options.

        Parameters
        ----------
        options
            The command line options.
        xyz_files : List[str]
            The list of XYZ files to be used in the workflow.

        Returns
        -------
        ConsistentActiveSpaceConfiguration
            The configuration created from the options.
        """
        config = cls()
        config.load_orbitals = options.load_orbitals
        config.cas_method = options.cas_method
        config.autocas_indices = [int(i) for i in options.autocas_molecule_indices.split()]
        config.basis_set = options.basis_set
        config.large_active_space = options.large_active_space
        config.exclude_core = options.exclude_core
        config.unmappable = options.always_include_unmapables
        config.xyz_files = [os.path.abspath(p) if p[0] != "/" else p for p in xyz_files]
        ConsistentActiveSpaceConfiguration.input_sanity_checks(config)
        config.base_load_path = os.path.join(*config.xyz_files[0].split("/")[:-1])  # type: ignore
        return config

    @classmethod
    def from_file(cls, filename: str):
        """
        Load the configuration from a file.

        Parameters
        ----------
        filename : str
            The path to the configuration file.

        Returns
        -------
        ConsistentActiveSpaceConfiguration
            The loaded configuration.
        """
        with open(filename, 'r') as file:
            configuration_dict = yaml.safe_load(file)
        config = cls()
        for key, value in configuration_dict.items():
            setattr(config, key, value)
        config.xyz_files = [os.path.abspath(p) if p[0] != "/" else p for p in config.xyz_files]
        ConsistentActiveSpaceConfiguration.input_sanity_checks(config)
        return config

    def get_serenity_interface_settings(self) -> Dict:
        """
        Get the Serenity interface settings for the consistent active space workflow.

        Returns
        -------
        Dict
            The settings for the Serenity interface.
        """
        settings = {
            "Interface": {
                "uhf": False,
                "localisation_method": "IBO",
                "alignment": True,
                "localize_virtuals": True,
                "optimized_mapping": True,
                "work_dir": "serenity/",
                "basis_set_set": self.basis_set,
                "score_start": 1.0,
                "skip_localization": False,
                "system_names": self.system_names,
                "molcas_orbital_files": [FileHandler.get_project_path() + "/initial/" for _ in self.system_names]
            }
        }
        if self.xyz_files:
            settings["Interface"]["xyz_files"] = self.xyz_files
        return settings

    def get_serenity_load_paths(self) -> List[str]:
        """
        Get the paths to load Serenity orbitals from.

        Returns
        -------
        List[str]
            The list of paths to load Serenity orbitals from.
        """
        return [self.base_load_path for _ in range(len(self.system_names))]
