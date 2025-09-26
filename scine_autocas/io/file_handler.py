# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
from typing import List, Optional

from scine_autocas.utils.defaults import Defaults


class FileHandler:
    """Handle all io files and path"""

    class DirectoryNames:
        """Handle all directory names."""
        initial_orbs: str = Defaults.DirName.initial_orbs
        """Initial orbital dir name"""
        initial_dmrg: str = Defaults.DirName.initial_dmrg
        """Initial dmrg dir name"""
        final_calc: str = Defaults.DirName.final_calc
        """Final calc dir name"""
        project_name: str = Defaults.DirName.project_name
        """Autocas project dir name"""
        dumper_name: str = Defaults.DirName.dumper_name
        """Autocas dump dir name"""
        scratch: str = Defaults.DirName.scratch
        """scratch dir for all calculations"""

    class QcMaquisNames:
        """Handle all filenames"""
        qcmaquis_result_file: str = Defaults.QcMaquisNames.qcmaquis_result_file
        """Name of qcmaquis result file"""
        qcmaquis_checkpoint_dir: str = Defaults.QcMaquisNames.qcmaquis_checkpoint_dir
        """Name of qcmaquis checkpoint dir"""

    class PlotNames:
        """Handle all plot names"""
        entanglement_file: str = Defaults.PlotNames.entanglement_file
        """Name of the entanglement diagram"""
        threshold_file: str = Defaults.PlotNames.threshold_file
        """Name of the threshold diagram"""

    current_dir: str = ""
    """Path to current dir"""
    project_dir: str = ""
    """Path to base autocas project dir"""

    @staticmethod
    def setup_project():
        """Convenience function to setup a standard autocas project."""
        if FileHandler.project_dir != "":
            FileHandler.set_project_dir(FileHandler.project_dir)
        else:
            FileHandler.set_project_dir()
        FileHandler.make_project_dir()

    # -------------------------------------- Setters --------------------------------------
    # Set the name or path of a dir without changing the path to it.

    @staticmethod
    def set_current_dir(current_dir: Optional[str] = None):
        """Set current directory for file handler.

        Parameters
        ----------
        current_dir : str, default=os.getcwd()
            path to current directory
        """
        if current_dir:
            if not current_dir.startswith("/") and not current_dir.startswith("~"):
                current_dir = os.getcwd() + f"/{current_dir}"
            FileHandler.current_dir = current_dir if not current_dir.endswith("/") else current_dir[:-1]
        else:
            FileHandler.current_dir = os.getcwd()
        print(f"Current directory: {FileHandler.current_dir}")

    @staticmethod
    def set_project_dir(project_dir: Optional[str] = None):
        """Set project directory.

        Parameters
        ----------
        project_dir : str, optional
            path to project dir, default is a dir in os.getcwd()
        """

        if project_dir:
            if not project_dir.startswith("/") and not project_dir.startswith("~"):
                project_dir = os.getcwd() + f"/{project_dir}"
            FileHandler.project_dir = project_dir if not project_dir.endswith("/") else project_dir[:-1]
        else:
            FileHandler.project_dir = os.getcwd()

    # -------------------------------------- Makers --------------------------------------
    # Create a directory in project path and change into it

    @staticmethod
    def make_scratch_dir():
        """Create the project directory."""
        if FileHandler.check_dir_exists(FileHandler.DirectoryNames.scratch):
            return

        os.mkdir(FileHandler.get_scratch_path())

    @staticmethod
    def make_project_dir():
        """Create the project directory."""
        if FileHandler.check_project_dir_exists():
            os.chdir(FileHandler.get_project_path())
            FileHandler.set_current_dir()
            return

        os.mkdir(FileHandler.get_project_path())
        os.chdir(FileHandler.get_project_path())
        FileHandler.set_current_dir()

    @staticmethod
    def make_initial_orbital_dir():
        """Create initial orbital dir in the project root dir.

        Function changes to that dir and set current dir to this location.
        The name is set by 'FileHandler.DirectoryNames.initial_orbs' and
        always relative to the project dir path.
        """
        if FileHandler.check_project_dir_exists():
            FileHandler.make_scratch_dir()
            os.chdir(FileHandler.get_project_path())
            FileHandler._make_calc_dir(FileHandler.DirectoryNames.initial_orbs)
            os.chdir(FileHandler.DirectoryNames.initial_orbs)
            FileHandler.set_current_dir()

    @staticmethod
    def make_initial_dmrg_dir():
        """Create initial dmrg dir in the project root dir.

        Function changes to that dir and set current dir to this location.
        The name is set by 'FileHandler.DirectoryNames.initial_orbs' and
        always relative to the project dir path.
        """
        if FileHandler.check_project_dir_exists():
            FileHandler.make_scratch_dir()
            os.chdir(FileHandler.get_project_path())
            FileHandler._make_calc_dir(FileHandler.DirectoryNames.initial_dmrg)
            os.chdir(FileHandler.DirectoryNames.initial_dmrg)
            FileHandler.set_current_dir()

    @staticmethod
    def make_final_calc_dir():
        """Create final calculation dir in project root dir.

        Function changes to that dir and set current dir to this location.
        The name is set by 'FileHandler.DirectoryNames.initial_orbs' and
        always relative to the project dir path.
        """
        if FileHandler.check_project_dir_exists():
            FileHandler.make_scratch_dir()
            os.chdir(FileHandler.get_project_path())
            FileHandler._make_calc_dir(FileHandler.DirectoryNames.final_calc)
            os.chdir(FileHandler.DirectoryNames.final_calc)
            FileHandler.set_current_dir()

    @staticmethod
    def _make_calc_dir(calc_dir: str):
        """Create a calculation directory.

        Parameters
        ----------
        calc_dir: str
            name of the directory
        """
        try:
            os.mkdir(calc_dir)
        except FileExistsError:
            pass

    @staticmethod
    def make_directory_from_project_root(new_dir: str):
        """Create any directory in the project root dir.

        Function changes to that dir and set current dir to this location.
        The name is set by 'FileHandler.DirectoryNames.initial_orbs' and
        always relative to the project dir path.

        Parameters
        ----------
        new_dir : str
            Path to directory.
        """
        if new_dir.startswith("/") or new_dir.startswith("~"):
            raise ValueError("You can only create realtive directories to the project path")
        if FileHandler.check_project_dir_exists():
            FileHandler.make_scratch_dir()
            os.chdir(FileHandler.get_project_path())
            FileHandler._make_calc_dir(new_dir)
            os.chdir(new_dir)
            FileHandler.set_current_dir()

    # -------------------------------------- Checkers --------------------------------------
    # just check if a directory or so exists
    @staticmethod
    def check_dir_exists(dir_name: str) -> bool:
        """Check if dir in autocas project exists.

        Parameters
        ----------
        dir_name : str
            name of the directory in project path.

        Returns
        -------
        bool
            True, if dir exists
            False, else
        """
        if dir_name.startswith("/") or dir_name.startswith("~"):
            raise ValueError("You can only check directories inside project path")
        if os.path.exists(FileHandler.get_project_path() + f"/{dir_name}"):
            return True
        return False

    @staticmethod
    def check_project_dir_exists() -> bool:
        """Check if autocas project exists.

        Returns
        -------
        bool
            True, if dir exists
            False, else
        """
        if os.path.exists(FileHandler.get_project_path()):
            return True
        return False

    # -------------------------------------- Changers --------------------------------------
    # just change to a directory
    @staticmethod
    def ch_to_initial_orbital_dir():
        """Change directory to initial orbital dir."""
        os.chdir(FileHandler.get_project_path())
        os.chdir(FileHandler.DirectoryNames.initial_orbs)
        FileHandler.set_current_dir()

    @staticmethod
    def ch_to_initial_dmrg_dir():
        """Change directory to initial dmrg dir."""
        os.chdir(FileHandler.get_project_path())
        os.chdir(FileHandler.DirectoryNames.initial_dmrg)
        FileHandler.set_current_dir()

    @staticmethod
    def ch_to_final_calc_dir():
        """Change directory to final calc dir."""
        os.chdir(FileHandler.get_project_path())
        os.chdir(FileHandler.DirectoryNames.final_calc)
        FileHandler.set_current_dir()

    @staticmethod
    def ch_to_current_dir():
        """Change directory current_dir."""
        os.chdir(FileHandler.current_dir)

    @staticmethod
    def ch_to_project_dir():
        """Change directory current_dir."""
        os.chdir(FileHandler.current_dir)

    # -------------------------------------- Getters --------------------------------------
    # get path and dir information

    @staticmethod
    def get_scratch_path() -> str:
        """Return Project dir path.

        Returns
        -------
        path : str
            Path to the project directory
        """
        print(FileHandler.get_project_path() + "/" + FileHandler.DirectoryNames.scratch)
        return FileHandler.get_project_path() + "/" + FileHandler.DirectoryNames.scratch

    @staticmethod
    def get_project_path() -> str:
        """Return Project dir path.

        Returns
        -------
        path : str
            Path to the project directory
        """
        return FileHandler.project_dir + "/" + FileHandler.DirectoryNames.project_name

    @staticmethod
    def get_all_dmrg_dirs() -> List[str]:
        """Get all dmrg dirs in current project.

        Returns
        -------
        dmrg_dirs : List[str]
            all dir names starting with 'dmrg' in the autocas project
        """
        dmrg_dirs = []
        for dirpath in os.listdir(FileHandler.get_project_path()):
            if dirpath.startswith("dmrg"):
                dmrg_dirs.append(dirpath)
        return dmrg_dirs
