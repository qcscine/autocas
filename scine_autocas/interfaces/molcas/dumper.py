"""Takes care of dumping molcas related data."""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import shutil

# import subprocess


class Dumper:
    """Class to prevent calculations to overide previous ones.

    Attributes
    ----------
    initial_orbital_dir : str, default = "initial"
        name of the dir to store the initial hf calculation
    initial_dmrg_dir : str, default = "dmrg"
        name of the dir to store the initial dmrg calculation
    final_calculation_dir : str, default = "final"
        name of the dir to store the final calculation
    project_name : str
        name of the project
    project_dir : str
        path to the project
    orbital_file : str
        name of the current orbital file
    large_cas_counter : int
        counter of dmrg calculations for the large active space protocol
    current_dir : str
        path of the current dir
    large_cas : bool, default = False
        flag to enable dumping for the large active space protocol
    """

    __slots__ = (
        "initial_orbital_dir",
        "initial_dmrg_dir",
        "final_calculation_dir",
        "project_name",
        "project_dir",
        "orbital_file",
        "large_cas_counter",
        "current_dir",
        "large_cas",
    )

    def __init__(self):
        """Construct dumper."""
        self.initial_orbital_dir: str = "initial"
        """directory to store the calculation of the initial orbitals"""
        self.initial_dmrg_dir: str = "dmrg"
        """directory to store the initial DMRG calculation which provides the entropies"""
        self.final_calculation_dir: str = "final"
        """dir to store the final calculation"""
        self.project_name: str = ""
        """name of the project"""
        self.project_dir: str = ""
        """root directory of the project"""
        self.orbital_file: str = ""
        """name of the Molcas orbital file"""
        self.large_cas_counter: int = 0
        """counts the number of DMRG calculations for the large CAS protocol to name directories"""
        self.current_dir: str = ""
        """path of the current dir"""
        self.large_cas: bool = False
        """flag to enable dumping for the large active space protocol"""

    def copy_orbital_file(self) -> str:
        """Copy the orbital file from the initial SCF calculation to a new
        location.

        Raises
        ------
        FileNotFoundError
            if no orbital file can be found
        """
        new_orbital_file = self.current_dir + "/" + self.project_name + ".scf.h5_sel"
        try:
            shutil.copyfile(self.orbital_file, new_orbital_file)
            # subprocess.call(
            #     f"cp {self.orbital_file} {new_orbital_file}", shell=True
            # )
        except EnvironmentError as exc:
            raise FileNotFoundError(
                f"could not find {self.orbital_file}. Did Molcas run?"
            ) from exc
        # self.orbital_file = new_orbital_file
        return new_orbital_file

    def create_project_dir(self):
        """Create project directory."""
        if self.project_dir.split("/")[-1] != self.project_name:
            self.project_dir = self.project_dir + "/" + self.project_name
            # print("Create project directory: %s" % self.project_dir)
            os.makedirs(self.project_dir, exist_ok=True)

    def setup_sub_dir(self):
        """Create all subdirectories."""
        copy_orbitals = True
        if self.current_dir == "":
            self.current_dir = self.project_dir + "/" + self.initial_orbital_dir
            self.orbital_file = self.current_dir + "/" + self.project_name + ".scf.h5"
            copy_orbitals = False
        elif self.current_dir.split("/")[-1] == self.initial_orbital_dir or (
            self.current_dir.split("/")[-1][:4] == self.initial_dmrg_dir
            and self.large_cas
        ):
            self.current_dir = self.project_dir + "/" + self.initial_dmrg_dir
            if self.large_cas:
                self.large_cas_counter += 1
                self.current_dir += "_" + str(self.large_cas_counter)
        elif (
            self.current_dir.split("/")[-1][:4] == self.initial_dmrg_dir
            and not self.large_cas
        ):
            self.current_dir = self.project_dir + "/" + self.final_calculation_dir
        else:
            print("ups")
        # print("Creating sub dir:", self.current_dir)
        os.makedirs(self.current_dir, exist_ok=True)
        if copy_orbitals:
            self.copy_orbital_file()
