"""The Molcas environment.

This module handles the molcas environment, while respecting already set
environment variable, without overwriting them.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
from typing import Dict


class Environment:
    """Class to control and set up the environment for Molcas.

    The environment respects already set system variables. So if you have
    already set variables in the current session, these are not overwritten
    by the Environment class.

    Attributes
    ----------
    environment : Dict[str, str]
        the system environment
    project_name : str
        the molcas project name
    molcas_binary : str
        path to the pymolcas binary
    molcas_flags : str, default = "-b 1 -f"
        all pymolcas flags
    molcas_memory : str, default = "12000"
        the amount of memory accessible by molcas in [mb]
    molcas_scratch_dir : str
        the molcas scratch dir, NOT the work dir
    molcas_nprocs : str
        the number of molcas mpi nodes and qcmaquis omp_threads

    Notes
    -----
    if no molcas_binary is set and the "MOLCAS" variable is found in the environment, e.g. through a .bashrc,
    the environemnt sets the binary to the path $MOLCAS/pymolcas
    """

    __slots__ = (
        "environment",
        "project_name",
        "molcas_binary",
        "molcas_flags",
        "molcas_memory",
        "molcas_scratch_dir",
        "molcas_nprocs",
    )

    def __init__(self, settings_dict=None):
        """Init function."""
        self.environment: Dict[str, str] = {}
        """The environment"""
        self.project_name: str = ""
        """name of the project. Will be used by Molcas to name some files"""
        self.molcas_binary: str = ""
        """name of the Molcas binar"""
        self.molcas_flags: str = "-b 1 -f"
        """flags for the Molcas binar"""
        self.molcas_memory: str = "12000"
        """amount available memory for Molca"""
        self.molcas_scratch_dir: str = ""
        """dir for molcas to store internal stuff"""
        self.molcas_nprocs: str = "1"
        """number of mpi nodes"""
        if settings_dict is not None:
            for key in settings_dict:
                if hasattr(self, key):
                    setattr(self, key, str(settings_dict[key]))

    def __str__(self) -> str:
        """Print environment.

        Returns
        -------
        output : str
            the output if class is printed
        """
        output = ""
        output += f"MOLCAS_PROJECT:  {self.project_name}\n"
        output += f"MOLCAS_MEM:      {self.molcas_memory}\n"
        output += f"WorkDir:         {self.molcas_scratch_dir}\n"
        output += f"MOLCAS_NPROCS:   {self.molcas_nprocs}\n"
        output += f"OMP_NUM_THREADS: {self.molcas_nprocs}\n"
        output += f"QCMaquis_CPU:    {self.molcas_nprocs}"
        return output

    def _set_environment(self, envir_string: str, environment_variable: str):
        """Set the environment variable.

        Parameters
        ----------
        envir_string : str
            Name of the environment to set
        environment_variable : str
            if environment not set, the value to set
        """
        if envir_string not in self.environment:
            self.environment[envir_string] = environment_variable

    def make_environment(self, calc_dir: str) -> Dict[str, str]:
        """Set up the environment for Molcas.

        Parameters
        ----------
        calc_dir : str
            Path to the current work directory

        Returns
        -------
        environment : Dict[str, str]
            environment object

        Raises
        ------
        EnvironmentError
            if neither molcas_binary, nor the system variable MOLCAS is defined.
        """
        self.environment = os.environ.copy()
        if self.molcas_binary != "":
            pass
        elif "MOLCAS" in self.environment:
            self.molcas_binary = self.environment["MOLCAS"] + "/pymolcas"
        else:
            raise EnvironmentError(
                """
                no molcas_binary is set or Molcas is not installed or 'MOLCAS' system variable is not set.
                Please use the following command: \n
                export MOLCAS=/path/to/molcas/build \n
                and set the correct path to the molcas build directory.
                """
            )
        # if not self.molcas_binary:
        self._set_environment("MOLCAS_PROJECT", self.project_name)
        self._set_environment("MOLCAS_MEM", self.molcas_memory)
        # self._set_environment("MOLCAS_WORKDIR", "/tmp")
        self._set_environment("WorkDir", self.molcas_scratch_dir)
        self._set_environment("MOLCAS_OUTPUT", calc_dir)
        if self.molcas_nprocs != "1":
            self._set_environment("MOLCAS_NPROCS", self.molcas_nprocs)
            self._set_environment("OMP_NUM_THREADS", self.molcas_nprocs)
            self._set_environment("QCMaquis_CPUS", self.molcas_nprocs)
        if self.environment["WorkDir"] == "":
            self.environment["WorkDir"] = "/tmp"
        return self.environment
