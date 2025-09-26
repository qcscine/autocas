# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
from typing import Dict, Optional

from scine_autocas.io import FileHandler


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

    def __init__(self) -> None:  # , settings_dict: Dict[Any] = None):
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
        self.molcas_scratch_dir: str = f"{FileHandler.get_project_path()}/scratch"
        """dir for molcas to store internal stuff"""
        self.molcas_nprocs: str = "1"
        """number of mpi nodes"""

    def get_nthreads(self) -> str:
        """Get largest number of defined threads.

        Returns
        -------
        str
            number of threads
        """
        # default is single core
        defaults = int(self.molcas_nprocs)
        # environment variabels are optional
        try:
            omp_threads = int(self._get_value("OMP_NUM_THREADS"))
        except ValueError:
            omp_threads = 1
        try:
            molcas_nprocs = int(self._get_value('MOLCAS_NPROCS'))
        except ValueError:
            molcas_nprocs = 1
        try:
            qcmaquis_cpus = int(self._get_value('QCMaquis_CPU'))
        except ValueError:
            qcmaquis_cpus = 1
        n_threads = max(defaults, omp_threads, molcas_nprocs, qcmaquis_cpus)
        return str(n_threads)

    def _get_value(self, variable_name: str) -> str:
        """Get value from local environment.

        Parameters
        ----------
        environment : Dict[str, str]
            local environemnt
        variabel_name: str
            name of the variable

        Returns
        -------
        str
            variable value
        """
        try:
            return self.environment[variable_name]
        except KeyError:
            return ""

    def __str__(self) -> str:
        """Print environment.

        Returns
        -------
        output : str
            the output if class is printed
        """

        output = ""
        output += f"MOLCAS_PROJECT:  {self._get_value('MOLCAS_PROJECT')}\n"
        output += f"MOLCAS_MEM:      {self._get_value('MOLCAS_MEM')}\n"
        output += f"WorkDir:         {self._get_value('WorkDir')}\n"
        output += f"MOLCAS_NPROCS:   {self._get_value('MOLCAS_NPROCS')}\n"
        output += f"OMP_NUM_THREADS: {self._get_value('OMP_NUM_THREADS')}\n"
        output += f"QCMaquis_CPU:    {self._get_value('QCMaquis_CPU')}"
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

    def check_molcas_exists(self) -> bool:
        """Validate that pymolcas is in the correct path.

        Returns
        -------
        bool
            True if molcas is found

        Raises
        ------
        OSError
            if pymolcas is not found
        """
        if os.path.isfile(self.molcas_binary) and os.access(self.molcas_binary, os.X_OK):
            return True
        print("Cannot find molcas binary")
        print(f"molcas_binary variable is set to {self.molcas_binary}")
        return False

    def make_environment(self, calc_dir: Optional[str] = None) -> Dict[str, str]:
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
        OSError
            if neither molcas_binary, nor the system variable MOLCAS is defined.
        """
        self.environment = os.environ.copy()

        if self.molcas_binary != "":
            pass
        elif "MOLCAS" in self.environment:
            self.molcas_binary = self.environment["MOLCAS"] + "/pymolcas"
        else:
            raise OSError(
                """
                no molcas_binary is set or Molcas is not installed or 'MOLCAS' system variable is not set.
                Please use the following command: \n
                export MOLCAS=/path/to/molcas/build \n
                and set the correct path to the molcas build directory.
                """
            )
        if not self.check_molcas_exists():
            raise OSError(
                """
                no molcas_binary is set or Molcas is not installed or 'MOLCAS' system variable is not set.
                Please use the following command: \n
                export MOLCAS=/path/to/molcas/build \n
                and set the correct path to the molcas build directory.
                Ensure that a <pymolcas> binary is inside the build folder.
                """
            )
        print(f"""Found <pymolcas> binary in {self.molcas_binary}""")
        # if not self.molcas_binary:
        self._set_environment("MOLCAS_PROJECT", self.project_name)
        self._set_environment("MOLCAS_MEM", self.molcas_memory)
        # self._set_environment("MOLCAS_WORKDIR", "/tmp")
        self._set_environment("WorkDir", self.molcas_scratch_dir)
        if calc_dir:
            self._set_environment("MOLCAS_OUTPUT", calc_dir)
        else:
            self._set_environment("MOLCAS_OUTPUT", os.getcwd())

        if self.molcas_nprocs != "1":
            self._set_environment("MOLCAS_NPROCS", self.molcas_nprocs)
            self._set_environment("OMP_NUM_THREADS", self.molcas_nprocs)
            self._set_environment("QCMaquis_CPUS", self.molcas_nprocs)
        if self.environment["WorkDir"] == "":
            self.environment["WorkDir"] = "/tmp"
        return self.environment
