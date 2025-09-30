# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import re
import shutil
import subprocess
from typing import List, Optional, Tuple, Union

import numpy as np

from scine_autocas.interfaces.interface import Interface
from scine_autocas.interfaces.qcmaquis_utils import QcmaquisUtils
from scine_autocas.io import FileHandler
from scine_autocas.utils import Molecule

from .environment import Environment
from .input_handler import InputHandler
from .molcas_hdf5_utils import MolcasHdf5Utils

# from scine_autocas.utils.exceptions import InterfaceError


class Molcas(Interface):
    """Interface to the electronic structure program Molcas.

    Since there exists no Python interface right now, 'pymolcas' is
    called via subprocess, input files have to be written and hdf5 files
    have to be parsed and modified.

    Attributes
    ----------
    settings : Molcas.Settings
        provides all settings for molcas
    environment : Environment
        handles molcas environment variables
    qc_maquis : Qcmaquis
        handles qcmaquis related files and build entropies
    hdf5_utils : MolcasHdf5Utils
        handles molcas hdf5 files
    orbital_file : str
        the path of the orbital file
    project_name : str
        the molcas project name
    dump : bool
        flag to enable dumping
    input_handler : InputHandler
        handles molcas input files
    caspt2_energies : List[float]
        stores all caspt2 energies

    Notes
    -----
    The settings object is already in the __slots__ of the base interface.

    See Also
    --------
    settings : Interface.Settings
    """

    class Settings(Interface.Settings):
        """Controll the input parameters and writes input file for Molcas.

        Attributes
        ----------
        point_group : str, default = "C1"
            point group of the molecule in molcas notation
        cholesky : bool, default = True
            enable cholesky decomposition
        ci_root_string : str
            root string for excited states
        alpha_or_beta_string : int, default = 1
            UNUSED
        ipea : float, default = 0.00
            IPEA shift for CASPT2
        ci_size : List[int]
            size of the davidson
        weights : List[float]
            weight for each state
        roots : List[int]
            molcas nroots

        Notes
        -----
        For most attributes defined here, please see also the MOLCAS manual.
        """

        __slots__ = (
            "point_group",
            "cholesky",
            "ci_root_string",
            "alpha_or_beta_string",
            "ipea",
            "ci_size",
            "weights",
            "roots",
            "initial_orbitals",
            "initial_electrons",
            "active_orbitals",
            "active_electrons",
            "orbital_localisation",
            "localisation_space",
            "localisation_method",
            # "skip_scf",
        )

        def __init__(self, molecule: Molecule):  # , settings_dict: Optional[Dict] = None):
            """Construct class."""
            super().__init__(molecule=molecule)  # , settings_dict=settings_dict)
            self.point_group: str = "C1"
            """
            point group of the molecule in Molcas notation, e.g.
            "":          C1
            "XY":        C2
            "X":         Cs
            "XYZ":       Ci
            "XY Y":      C2v
            "XY XYZ":    C2h
            "XY XZ":     D2
            "XY XZ XYZ": C2h
            """
            self.cholesky: bool = True
            """toggle Cholesky decomposition"""
            self.ci_root_string: str = ""
            """string of the reference determinant"""
            self.alpha_or_beta_string: int = 1
            """string of determinants, NOT used"""
            self.ipea: float = 0.0
            """IPEA shift for CASPT2, only used if CASPT2 is in method"""
            self.ci_size: Optional[int] = None
            """size of the davidson"""
            self.roots: Optional[List[int]] = None
            """molcas nroots"""
            self.weights: Optional[List[int]] = None
            """weight for each state"""
            # TODO remove this
            self.active_orbitals: int = 0
            self.active_electrons: int = 0
            self.initial_orbitals: int = True
            self.initial_electrons: int = 0
            self.orbital_localisation = False
            self.localisation_space = "OCCUpied"
            self.localisation_method = "PIPEk-Mezey"
            # self.skip_scf = True

            # available method for this interface
            self._add_cas_method("casci", "casscf", "dmrgscf")
            self._add_post_cas_method("caspt2", "nevpt2")
            # if settings_dict:
            #     self.apply_settings(settings_dict)

    # Molcas
    __slots__ = (
        "environment",
        "qc_maquis",
        "hdf5_utils",
        "orbital_file",
        "project_name",
        "input_handler",
        "caspt2_energies",
    )

    def __init__(self, molecule: Molecule):
        """Construct a molcas interface.

        Parameters
        ----------
        molecule : Molecule
            contains molecular information

        See Also
        --------
        settings_dict : InputHandler
        """
        self.settings: Molcas.Settings
        # """provides all settings for molcas"""
        # initialize base interface
        super().__init__(molecule=molecule)
        self.settings = self.Settings(molecule)
        self.environment: Environment = Environment()
        """coltrols the enviroment for Molcas"""
        self.qc_maquis: QcmaquisUtils = QcmaquisUtils()
        """builds s1, s2 and Ipq from QCMaquis HDF5 files"""
        self.hdf5_utils: MolcasHdf5Utils = MolcasHdf5Utils()
        """handles Molcas HDF5 files"""
        self.orbital_file: str = ""
        """name of the orbital file, which has to be modified to enable Molcas "Typeindex" """
        self.project_name: str = "autocas_project"
        """name of the project"""
        self.input_handler: InputHandler = InputHandler()
        """handles molcas input files"""
        self.caspt2_energies: Optional[List[float]] = None
        """stores all caspt2 energies"""

        # always set up filehandler, to run molcas in specfific dir
        if not FileHandler.check_project_dir_exists():
            FileHandler.setup_project()
        else:
            FileHandler.ch_to_project_dir()

        # make ci string
        #

        # apply settings to subclasses
        # if settings_dict is not None:
        #     for key in settings_dict:
        #         if key == "environment":
        #             self.environment = Environment(settings_dict[key])
        #         elif key == "qc_maquis":
        #             pass
        #             # self.qc_maquis = Qcmaquis(settings_dict[key])
        #         elif key == "hdf5_utils":
        #             pass
        #             # self.hdf5_utils = MolcasHdf5Utils(settings_dict[key])
        #         elif key == "input_handler":
        #             pass
        #             # self.input_handler = InputHandler(settings_dict[key])
        #         elif key == "settings":
        #             pass
        #         elif hasattr(self, key):
        #             setattr(self, key, settings_dict[key])
    def make_ci_root(
        self,
        n_roots: Optional[int] = None,
        ci_size: Optional[int] = None,
        roots: Optional[List[int]] = None,
        weights: Optional[List[int]] = None,
    ):
        """Make the ci root string for excited state calculations.

        Parameters
        ----------
        n_roots : int
            number of roots
        ci_size : List[int]
            size of ci matrix in davidson
        roots : List[int]
            roots to use for calculation
        weights : List[float]
            weight of each root

        Returns
        -------
        root_string : str
            The full ci_root string, required for excited state calculations
        """
        if n_roots is None:
            n_roots = self.settings.n_excited_states
        root_string = str(n_roots) + " "
        if roots is not None and ci_size is not None:
            assert max(roots) >= ci_size
            root_string += str(ci_size) + "; "
            for root in roots:
                root_string += str(root) + " "
            root_string += "; "
            if weights is not None:
                assert len(weights) == len(roots)
                for weight in weights:
                    root_string += str(weight) + " "
            else:
                for _ in range(len(roots)):
                    root_string += str(1) + " "

        elif ci_size is not None:
            root_string += str(ci_size) + " " + str(1)

        else:
            root_string += str(n_roots) + " " + str(1)

        self.settings.ci_root_string = root_string
        return root_string

    def check_calculate_input(self, cas_occupation: Optional[List[int]] = None,
                              cas_indices: Optional[List[int]] = None) -> Tuple[List[int], List[int]]:
        """Check the input for the calculate function.

        Parameters
        ----------
        cas_occupation: : List[int], optional
            a list which contains the occupation of each orbital in the active space
        cas_indices: : List[int], optional
            a list which contains the indices of each orbital in the active space

        Returns
        -------
        cas_occupation: : List[int]
            a list which contains the occupation of each orbital in the active space
        cas_indices: : List[int]
            a list which contains the indices of each orbital in the active space

        Raises
        ------
        ValueError
            if different number of occupations and indices
        """
        if cas_occupation is None or cas_indices is None:
            return [], []
        if len(cas_occupation) is not len(cas_indices):
            raise ValueError(
                "Not the same number of orbital indices and orbital occupations"
            )
        return cas_occupation, cas_indices

    def run_molcas(self):
        """Run molcas, from current state.

        The function calls the molcas binary through a subprocess
        """
        # write input
        input_file = self.project_name + ".input"

        self.input_handler.write_input(self.settings, input_file, self.orbital_file)
        # Ensure formatted SCF orbitals exist in cwd so LUMORB/INPORB resolves
        try:
            scf_local = f"{self.project_name}.ScfOrb"
            if not os.path.exists(scf_local):
                cand = os.path.join(os.path.dirname(os.getcwd()), "initial", scf_local)
                if os.path.exists(cand):
                    shutil.copy(cand, scf_local)
        except Exception:
            pass
        # setup environment
        self.environment.project_name = self.project_name
        environment = self.environment.make_environment()
        # save current location
        pymolcas_string = self.environment.molcas_binary + " " + self.environment.molcas_flags
        print(self.environment.molcas_binary, self.environment.molcas_flags)
        print(self.environment)
        log_file = open("molcas.out", "w+")

        calculation_process = subprocess.Popen(
            [
                f"{pymolcas_string} {input_file} -nt {self.environment.get_nthreads()}"
            ],
            env=environment,
            shell=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
        calculation_process.communicate()

        molcas_internal_log_file = self.project_name + ".log"
        success_string = "Happy landing!"
        with open(molcas_internal_log_file, "r") as molcas_log:
            success = success_string in molcas_log.read()
        if not success:
            full_path = os.path.join(os.getcwd(), molcas_internal_log_file)
            err_file = molcas_internal_log_file.replace(".log", ".err")
            try:
                with open(molcas_internal_log_file, "r") as molcas_log:
                    print(f"Molcas log:\n{molcas_log.read()}")
                with open(err_file, "r") as err_log:
                    print(f"Molcas error log:\n{err_log.read()}")
            except FileNotFoundError:
                pass
            raise RuntimeError(f"Molcas calculation failed. Please check the log file {full_path} for details.")

    def reorder_indices_from_symmetry(self, cas_indices: List[int]) -> List[int]:
        """Reorder cas indices based on mo energies.

        Molcas changes the index order for calculations exploiting point group symmetry.
        Hence, these indices have to be reordered to fit entropies and other
        variables.

        Parameters
        ----------
        cas_indices : List[int], optional
            a list which contains the indice of each orbital in the active space

        Returns
        -------
        cas_indices_sym : List[int], optional
            a list which contains the indices of each orbital in the active space

        Raises
        ------
        ValueError
            if no MO energies can be read by the MolcasHdf5Utils instance
        """
        # check for mo energies
        if self.hdf5_utils.mo_energies:
            pass
        elif self.orbital_file != "":
            self.hdf5_utils.read_hdf5(self.orbital_file)
        else:
            raise ValueError("Please provide a scf.h5 file")

        # sort orbitals for correct indices
        tmp_mo_energies = np.array(self.hdf5_utils.mo_energies)
        sort_key = tmp_mo_energies.argsort()
        sort_key = sort_key.astype(int, copy=False)
        mo_energies_sorted = np.array(self.hdf5_utils.mo_energies)[sort_key].tolist()

        # get correct orbital indices
        cas_indices_sym = []
        for i in cas_indices:
            cas_indices_sym.append(
                self.hdf5_utils.mo_energies.index(mo_energies_sorted[i])
            )
        return cas_indices_sym

    def _initial_orbitals_impl(self):
        """Generate initial orbitals."""
        # check for dumping
        FileHandler.make_initial_orbital_dir()
        # start calculation
        print("Starting initial orbital calculation.")
        self.run_molcas()
        print("Initial orbitals finished.")

        # set class variables
        self.orbital_file = os.getcwd() + "/" + self.project_name + ".scf.h5"

        # get typeindice (and other variables)
        self.hdf5_utils.read_hdf5(self.orbital_file)

    def reorder_measurements_for_symmetry(self, cas_indices: List[int]):
        """Reorder dmrg measurements with respect to cas indices.

        If point group symmetry is used, the measurements usually are in an different order
        than the indices. Hence, this function is used to reorder them.

        Parameters
        ----------
        cas_indices : List[int], optional
            a list which contains the indices of each orbital in the active space
        """
        sort_key = (np.array(cas_indices).argsort()).argsort()
        self.qc_maquis.s1_entropy = self.qc_maquis.s1_entropy[sort_key]
        self.qc_maquis.s2_entropy = self.qc_maquis.s2_entropy[:, sort_key][sort_key, :]
        self.qc_maquis.mutual_information = self.qc_maquis.mutual_information[sort_key][:, sort_key]

    def _final_cas_impl(
        self, cas_occupation: List[int], cas_indices: List[int]
    ) -> List[float]:
        if self.settings.n_excited_states > 0:
            self.make_ci_root(
                self.settings.n_excited_states + 1,
                self.settings.ci_size,
                self.settings.roots,
                self.settings.weights,
            )
        self.settings.initial_orbitals = False
        FileHandler.make_final_calc_dir()
        orb_file_name = self.orbital_file.split('/')[-1]
        if not os.path.exists(orb_file_name):
            new_orb_file = os.getcwd() + f"/{orb_file_name}_sel"
            shutil.copy(self.orbital_file, new_orb_file)
            self.orbital_file = new_orb_file

        self.settings.active_electrons = sum(cas_occupation)
        self.settings.active_orbitals = len(cas_occupation)

        # make everything safe for symmetry
        if self.settings.point_group != "C1":
            cas_indices = self.reorder_indices_from_symmetry(cas_indices)
        self.hdf5_utils.modify_hdf5(self.orbital_file, cas_indices)

        # start calculation
        self.run_molcas()

        # make measurements
        try:
            self.make_dmrg_measurements()
            # sort entropies
            if self.settings.point_group != "C1":
                # self.reorder_indices_from_symmetry(cas_indices)
                self.reorder_measurements_for_symmetry(cas_indices)

            energy = [self.hdf5_utils.energy]
            if self.caspt2_energies is not None:
                energy = self.caspt2_energies
            if self.settings.n_excited_states > 0:
                return energy

            return energy

        except FileNotFoundError as exc:
            raise FileNotFoundError("Calculation failed") from exc

    def _initial_cas_impl(
        self, cas_occupation: List[int], cas_indices: List[int]
    ) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
        if self.settings.n_excited_states > 0:
            self.make_ci_root(
                # because molcas requires number of states
                self.settings.n_excited_states + 1,
                self.settings.ci_size,
                self.settings.roots,
                self.settings.weights,
            )

        self.settings.initial_orbitals = False
        if "dmrg" not in os.getcwd().split("/")[-1]:
            FileHandler.make_initial_dmrg_dir()
        orb_file_name = self.orbital_file.split('/')[-1]
        if not os.path.exists(orb_file_name):
            new_orb_file = os.getcwd() + f"/{orb_file_name}_sel"
            shutil.copy(self.orbital_file, new_orb_file)
            self.orbital_file = new_orb_file

        self.settings.active_electrons = sum(cas_occupation)
        self.settings.active_orbitals = len(cas_occupation)

        # make everything safe for symmetry
        if self.settings.point_group != "C1":
            cas_indices = self.reorder_indices_from_symmetry(cas_indices)
        self.hdf5_utils.modify_hdf5(self.orbital_file, cas_indices)

        # start calculation
        self.run_molcas()

        # make measurements
        try:
            s1_list, s2_list, mutual_information_list = self.make_dmrg_measurements()
            # sort entropies
            if self.settings.point_group != "C1":
                # self.reorder_indices_from_symmetry(cas_indices)
                self.reorder_measurements_for_symmetry(cas_indices)

            energy = [self.hdf5_utils.energy]
            if self.caspt2_energies is not None:
                energy = self.caspt2_energies
            if self.settings.n_excited_states > 0:
                return (energy, s1_list, s2_list, mutual_information_list)

            return energy, self.qc_maquis.s1_entropy, self.qc_maquis.s2_entropy, self.qc_maquis.mutual_information

        except FileNotFoundError as exc:
            raise FileNotFoundError("Calculation failed") from exc

    def _check_interface_exists(self):
        """Check that pymolcas is found."""
        self.environment.make_environment()
        if not self.environment.check_molcas_exists():
            raise EnvironmentError()

    def make_dmrg_measurements(self):
        """Build all entropies from QCMaquis HDF5 file."""
        try:
            s1_list = []
            s2_list = []
            mutual_information_list = []
            if self.settings.n_excited_states > 0:
                if os.path.isfile(os.getcwd() + "/" + self.project_name + ".results_state.0.h5"):
                    for i in range(self.settings.n_excited_states + 1):
                        self.qc_maquis.read_hdf5(
                            os.getcwd()
                            + "/"
                            + self.project_name
                            + ".results_state."
                            + str(i)
                            + ".h5"
                        )
                        self.qc_maquis.make_diagnostics()
                        s1_list.append(self.qc_maquis.s1_entropy)
                        s2_list.append(self.qc_maquis.s2_entropy)
                        mutual_information_list.append(self.qc_maquis.mutual_information)
                else:
                    # no dmrg
                    pass
            else:
                if os.path.isfile(os.getcwd() + "/" + self.project_name + ".results_state.0.h5"):
                    self.qc_maquis.read_hdf5(os.getcwd() + "/" + self.project_name + ".results_state.0.h5")
                    self.qc_maquis.make_diagnostics()
                else:
                    # no dmrg
                    pass

                s1_list = self.qc_maquis.s1_entropy
                s2_list = self.qc_maquis.s2_entropy
                mutual_information_list = self.qc_maquis.mutual_information
        # no DMRG
        except AttributeError:
            pass
        if os.path.isfile(os.getcwd() + "/" + self.project_name + ".dmrgscf.h5"):
            self.hdf5_utils.get_energy(os.getcwd() + "/" + self.project_name + ".dmrgscf.h5")
        elif os.path.isfile(os.getcwd() + "/" + self.project_name + ".rasscf.h5"):
            self.hdf5_utils.get_energy(os.getcwd() + "/" + self.project_name + ".rasscf.h5")
        else:
            pass

        self.get_caspt2_from_log(os.getcwd() + "/" + self.project_name + ".log")
        return s1_list, s2_list, mutual_information_list

    def _convert_type_indices(self, type_indices):
        """Convert type indices from hdf5 file to list"""
        indices = []
        for i, index in enumerate(type_indices):
            if (
                index.decode("UTF-8") == "2"
                or index.decode("UTF-8") == "1"
            ):
                indices.append(i)
        return indices

    def get_caspt2_from_log(self, file_name: str):
        """Get energy from molcas log.

        Parameters
        ----------
        file_name : str
            name of the file to read energies from
        """
        energies = []
        with open(file_name, "r", encoding="UTF-8") as molcas_log:
            caspt2_flag = False
            for line in molcas_log:
                if caspt2_flag is True:
                    if re.search("Total energy:", line):
                        energies.append(float(line.split()[-1]))
                if re.search("Total CASPT2 energies:", line):
                    caspt2_flag = True
        if caspt2_flag is True:
            self.caspt2_energies = energies

    def analyze(
        self, qcmaquis_result_file: Optional[str] = None, molcas_hdf5_file: Optional[str] = None
    ) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray, List[int], List[int]],
               Tuple[List[int], List[int]],
               Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Analyze output files."""
        orbital_indices: Optional[List[int]] = None
        if molcas_hdf5_file:
            self.hdf5_utils.read_hdf5(molcas_hdf5_file)
            orbital_indices = self._convert_type_indices(self.hdf5_utils.type_indices)
        if qcmaquis_result_file:
            self.qc_maquis.read_hdf5(qcmaquis_result_file)
            self.qc_maquis.make_diagnostics()

        if molcas_hdf5_file is not None and qcmaquis_result_file is not None and orbital_indices is not None:
            return (
                self.qc_maquis.s1_entropy,
                self.qc_maquis.s2_entropy,
                self.qc_maquis.mutual_information,
                orbital_indices,
                self.hdf5_utils.occupations,
            )
        if molcas_hdf5_file is not None and qcmaquis_result_file is None and orbital_indices is not None:
            return orbital_indices, self.hdf5_utils.occupations
        if molcas_hdf5_file is None and qcmaquis_result_file:
            return (
                self.qc_maquis.s1_entropy,
                self.qc_maquis.s2_entropy,
                self.qc_maquis.mutual_information,
            )
        raise ValueError("analyze function requires at least one argument!")

    def resource_estimate(self):
        """Estimate if DMRG is feasible."""
        raise NotImplementedError("Resource estimation is not implemented yet")
