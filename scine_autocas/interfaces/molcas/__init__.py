"""Provide interface to molcas.

Currently, no Python-Bindings exist to the Molcas electronic structure
programe, hence this whole folder exists. Only here implemented
functionalities are supported by this interface and "pymolcas" is called
via subprocess.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import re
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces import Interface
from scine_autocas.interfaces.molcas.dumper import Dumper
from scine_autocas.interfaces.molcas.environment import Environment
from scine_autocas.interfaces.molcas.input_handler import InputHandler
from scine_autocas.interfaces.molcas.molcas_hdf5_utils import MolcasHdf5Utils
from scine_autocas.interfaces.qcmaquis import Qcmaquis


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
    dumper : Dumper
        handles dumping of molcas calculations
    input_handler : InputHandler
        handles molcas input files
    s1_list : List[np.ndarray]
        stores all s1_entropies from large cas calculations
    s2_list : List[np.ndarra]
        stores all s2_entropies from large cas calculations
    mutual_information_list : List[np.ndarray]
        stores all mutual information from large cas calculations
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
        uhf : bool, default = False
            enable unrestricted hf
        active_electrons : int
            number of active electrons
        active_orbitals : int
            number of active orbitals
        ci_root_string : str
            root string for excited states
        alpha_or_beta_string : int, default = 1
            UNUSED
        fiedler : bool, default = False
            enable fiedler ordering for dmrg
        ipea : float, default = 0.00
            IPEA shift for CASPT2
        initial_orbitals : bool, default = True
            flag to turn off hf calculations for dmrg and final calculations
        only_hf : bool, default = False
            stop molcas from doing cas calculations
        orbital_localisation : bool, default = False
            flag for orbital localisation
        localisation_space : str, default = "OCCUpied"
            the orbital localisation space
        localisation_method : str, default = "PIPEk-Mezet"
            the localisation method
        n_excited_states : int, default = 0
            number of excited states
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
            "uhf",
            "active_electrons",
            "active_orbitals",
            "ci_root_string",
            "alpha_or_beta_string",
            "fiedler",
            "ipea",
            "initial_orbitals",
            "only_hf",
            "orbital_localisation",
            "localisation_space",
            "localisation_method",
            "ci_size",
            "weights",
            "roots",
            "skip_scf"
        )

        def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict] = None):
            """Construct class."""
            super().__init__(molecules=molecules, settings_dict=settings_dict)
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
            self.uhf: bool = False
            """enable unrestricted calculation, for an uneven number of electrons it is automatically enabled"""
            self.active_electrons: int = 0
            """number of active electrons"""
            self.active_orbitals: int = 0
            """number of active orbitals, only used if no orbital file with typeindices is provided"""
            self.ci_root_string: str = ""
            """string of the reference determinant"""
            self.alpha_or_beta_string: int = 1
            """string of determinants, NOT used"""
            self.fiedler: bool = False
            """enable fiedler ordering"""
            self.ipea: float = 0.0
            """IPEA shift for CASPT2, only used if CASPT2 is in method"""
            self.initial_orbitals: bool = True
            """flag to turn of subsequent HF calculations"""
            self.ci_size: Optional[int] = None
            """size of the davidson"""
            self.roots: Optional[List[int]] = None
            """molcas nroots"""
            self.weights: Optional[List[int]] = None
            """weight for each state"""
            self.only_hf: bool = False
            """stop molcas from doing cas calculations"""
            self.orbital_localisation: bool = False
            """enable orbital localisation"""
            self.localisation_space: str = "OCCUpied"
            """localisation space, e.g.
            OCCUpied, VIRTual
            """
            self.localisation_method: str = "PIPEk-Mezey"
            """localisation method, e.g.
            PIPEk-Mezey, BOYS, EDMIston-Ruedenberg, CHOLesky, PAO
            """
            self.skip_scf = False
            """If true, only one SCF iteration is performed."""
            if settings_dict:
                self.apply_settings(settings_dict)

    # Molcas
    __slots__ = (
        "environment",
        "qc_maquis",
        "hdf5_utils",
        "orbital_file",
        "project_name",
        "dump",
        "dumper",
        "input_handler",
        "s1_list",
        "s2_list",
        "mutual_information_list"
        "caspt2_energies",
    )

    def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict[str, Any]] = None):
        """Construct a molcas interface.

        Parameters
        ----------
        molecule : Molecule
            contains molecular information
        settings_dict : Dict[str, Any], optional
            holds all settings provided by a yaml input

        See Also
        --------
        settings_dict : InputHandler
        """
        self.settings: Molcas.Settings
        """provides all settings for molcas"""
        # initialize base interface
        super().__init__(molecules=molecules, settings_dict=settings_dict)
        self.environment: Environment = Environment()
        """coltrols the enviroment for Molcas"""
        self.qc_maquis: Qcmaquis = Qcmaquis()
        """builds s1, s2 and Ipq from QCMaquis HDF5 files"""
        self.hdf5_utils: MolcasHdf5Utils = MolcasHdf5Utils()
        """handles Molcas HDF5 files"""
        self.orbital_file: str = ""
        """name of the orbital file, which has to be modified to enable Molcas "Typeindex" """
        self.project_name: str = "autocas_project"
        """name of the project"""
        self.dump: bool = True
        """allows to dump all calculations in a separate directory"""
        self.dumper: Dumper = Dumper()
        """controls the dumping of calculations by creating a directory structure"""
        self.input_handler: InputHandler = InputHandler()
        """handles molcas input files"""
        self.s1_list: List[np.ndarray]
        """stores all s1_entropies from large cas calculations"""
        self.s2_list: List[np.ndarray]
        """stores all s2_entropies from large cas calculations"""
        self.mutual_information_list: List[np.ndarray]
        """stores all mutual information from large cas calculations"""
        self.caspt2_energies: Optional[List[float]] = None
        """stores all caspt2 energies"""

        # make ci string
        if self.settings.n_excited_states > 0:
            self.make_ci_root(
                self.settings.n_excited_states,
                self.settings.ci_size,
                self.settings.roots,
                self.settings.weights,
            )

        # apply settings to subclasses
        if settings_dict is not None:
            for key in settings_dict:
                if key == "environment":
                    self.environment = Environment(settings_dict[key])
                elif key == "qc_maquis":
                    pass
                    # self.qc_maquis = Qcmaquis(settings_dict[key])
                elif key == "hdf5_utils":
                    pass
                    # self.hdf5_utils = MolcasHdf5Utils(settings_dict[key])
                elif key == "dumper":
                    pass
                    # self.dumper = Dumper(settings_dict[key])
                elif key == "input_handler":
                    pass
                    # self.input_handler = InputHandler(settings_dict[key])
                elif key == "settings":
                    pass
                elif hasattr(self, key):
                    setattr(self, key, settings_dict[key])

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
            a list which contains the indice of each orbital in the active space

        Returns
        -------
        cas_occupation: : List[int]
            a list which contains the occupation of each orbital in the active space
        cas_indices: : List[int]
            a list which contains the indice of each orbital in the active space

        Raises
        ------
        ValueError
            if different number of occupations and indices
        """
        if cas_occupation is None or cas_indices is None:
            print("Calculating only mean field.")
            self.settings.only_hf = True
            return [], []
        if len(cas_occupation) is not len(cas_indices):
            raise ValueError(
                "Not the same number of orbital indices and orbital occupations"
            )
        return cas_occupation, cas_indices

    def set_cas_variables(self, cas_occupation: List[int]):
        """Set spin_multiplicity, active_electrons, active_orbitals.

        Parameters
        ----------
        cas_occupation : List[int]
            a list which contains the occupation of each orbital in the active space
        """
        # check multiplicity
        spin_mult = 1
        for i in cas_occupation:
            if i == 1:
                spin_mult += 1
        self.settings.active_orbitals = len(cas_occupation)
        self.settings.active_electrons = int(sum(cas_occupation))
        self.settings.spin_multiplicity = spin_mult

    def run_molcas(self):
        """Run molcas, from current state.

        The function calls the molcas binary through a subprocess
        """
        # write input
        input_file = self.settings.work_dir + "/" + self.project_name + ".input"
        self.input_handler.write_input(self.settings, input_file, self.orbital_file)
        # setup environment
        self.environment.project_name = self.project_name
        environment = self.environment.make_environment(self.settings.work_dir)
        # save current location
        current_dir = os.getcwd()
        os.chdir(self.settings.work_dir)
        pymolcas_string = self.environment.molcas_binary + " " + self.environment.molcas_flags
        print("Starting MOLCAS:")
        print(self.environment.molcas_binary, self.environment.molcas_flags)
        print("MOLCAS environment:")
        print(self.environment)
        calculation_process = subprocess.Popen(
            [
                f"{pymolcas_string} {input_file}"
            ],
            env=environment,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        calculation_process.communicate()
        # go back to current location
        os.chdir(current_dir)

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
            a list which contains the indice of each orbital in the active space

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

        # for i in cas_indices:
        #     cas_mos.append(mo_energies_sorted.index(self.hdf5_utils.mo_energies[i]))
        # cas_mos = np.array(cas_mos)
        return cas_indices_sym

    def run_initial_orbitals(self):
        """Generate initial orbitals."""
        # check for dumping
        if self.dump:
            self.dumper.setup_sub_dir()
            self.settings.work_dir = self.dumper.current_dir

        # start calculation
        print("Starting initial orbital calculation.")
        self.run_molcas()
        print("Initial orbitals finished.")

        # set class variables
        self.orbital_file = (self.settings.work_dir + "/" + self.project_name + ".scf.h5")
        self.dumper.orbital_file = self.orbital_file

        # get typeindice (and other variables)
        self.hdf5_utils.read_hdf5(self.orbital_file)
        self.settings.initial_orbitals = False

    def reorder_measurements_for_symmetry(self, cas_indices: List[int]):
        """Reorder dmrg measurements with respect to cas indices.

        If point group symmetry is used, the measurements usually are in an different order
        than the indices. Hence, this function is used to reorder them.

        Parameters
        ----------
        cas_indices : List[int], optional
            a list which contains the indice of each orbital in the active space
        """
        sort_key = (np.array(cas_indices).argsort()).argsort()
        self.qc_maquis.s1_entropy = self.qc_maquis.s1_entropy[sort_key]
        self.qc_maquis.s2_entropy = self.qc_maquis.s2_entropy[:, sort_key][sort_key, :]
        self.qc_maquis.mutual_information = self.qc_maquis.mutual_information[sort_key][:, sort_key]

    def calculate(self, cas_occupation: Optional[List[int]] = None, cas_indices: Optional[List[int]] = None
                  ) -> Tuple[Union[float, np.ndarray], Union[np.ndarray, List[np.ndarray]],
                             Union[np.ndarray, List[np.ndarray]], Union[np.ndarray, List[np.ndarray]]]:
        """Calculate a system with the correpsonding setting.

        If symmetry is enabled the cas_occupation and cas_indices are still expected to be ordered symmetry independent.
        Also the s1, s2 and mutual information matrices are reordered to be symmetry indipenedent.

        Parameters
        ----------
        cas_occupation : List[int], optional
            a list which contains the occupation of each orbital in the active space
        cas_indices : List[int], optional
            a list which contains the indice of each orbital in the active space

        Returns
        -------
        energy : Union[float, np.ndarray]
            energy, or enery for each state
        s1 : Union[np.ndarray, List[np.ndarray]]
            single orbital entropy, or s1 for each state
        s2 : Union[np.ndarray, List[np.ndarray]]
            two orbital entropy, or s2 for each state
        mutual_information : Union[np.ndarray, List[np.ndarray]]
            mutual information, or mutual information for each state
        """
        # check input
        cas_occupation, cas_indices = self.check_calculate_input(cas_occupation, cas_indices)

        # setup dumping
        if self.dump:
            self.dumper.project_dir = self.settings.work_dir
            self.dumper.project_name = self.project_name
            self.dumper.create_project_dir()

        # set internal variables
        if not self.settings.only_hf:
            self.set_cas_variables(cas_occupation)

        # check for initial orbitals
        if self.settings.initial_orbitals:
            self.run_initial_orbitals()

        # make everything safe for symmetry
        if self.settings.point_group != "C1":
            cas_indices = self.reorder_indices_from_symmetry(cas_indices)
        # stop after hartree fock
        if self.settings.only_hf:
            self.settings.work_dir = self.dumper.project_dir
            self.settings.only_hf = False
            return (0.0, np.array([]), np.array([]), np.array([]))

        # prepare dmrg calculations
        if self.dump:
            # ensure project dir exists and set corrects path
            self.dumper.create_project_dir()
            self.dumper.setup_sub_dir()
            self.settings.work_dir = self.dumper.current_dir

        # modify orbital file in dmrg or final
        self.orbital_file = (
            self.settings.work_dir + "/" + self.project_name + ".scf.h5_sel"
        )
        self.hdf5_utils.modify_hdf5(self.orbital_file, cas_indices)

        # start calculation
        print("Starting calculation.")
        self.run_molcas()
        print("Calculation finished.")

        # make measurements
        try:
            self.make_dmrg_measurements()
            # sort entropies
            if self.settings.point_group != "C1":
                # self.reorder_indices_from_symmetry(cas_indices)
                self.reorder_measurements_for_symmetry(cas_indices)

            self.settings.work_dir = self.dumper.project_dir

            energy = self.hdf5_utils.energy
            if self.caspt2_energies is not None:
                energy = np.array(self.caspt2_energies)
            if self.settings.n_excited_states > 0:
                return (energy, self.s1_list, self.s2_list, self.mutual_information_list)

            return (
                energy, self.qc_maquis.s1_entropy, self.qc_maquis.s2_entropy,
                self.qc_maquis.mutual_information
            )

        except FileNotFoundError as exc:
            raise FileNotFoundError("Calculation failed") from exc

    def make_dmrg_measurements(self):
        """Build all entropies from QCMaquis HDF5 file."""
        if self.settings.n_excited_states > 0:
            s1_list = []
            s2_list = []
            mutual_information_list = []
            if os.path.isfile(self.settings.work_dir + "/" + self.project_name + ".results_state.0.h5"):
                for i in range(self.settings.n_excited_states):
                    self.qc_maquis.read_hdf5(
                        self.settings.work_dir
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
                self.s1_list = s1_list
                self.s2_list = s2_list
                self.mutual_information_list = mutual_information_list
            else:
                # no dmrg
                pass
        else:
            if os.path.isfile(self.settings.work_dir + "/" + self.project_name + ".results_state.0.h5"):
                self.qc_maquis.read_hdf5(self.settings.work_dir + "/" + self.project_name + ".results_state.0.h5")
                self.qc_maquis.make_diagnostics()
            else:
                # no dmrg
                pass
        if os.path.isfile(self.settings.work_dir + "/" + self.project_name + ".dmrgscf.h5"):
            self.hdf5_utils.get_energy(self.settings.work_dir + "/" + self.project_name + ".dmrgscf.h5")
        elif os.path.isfile(self.settings.work_dir + "/" + self.project_name + ".rasscf.h5"):
            self.hdf5_utils.get_energy(self.settings.work_dir + "/" + self.project_name + ".rasscf.h5")
        else:
            pass
        self.get_caspt2_from_log(self.settings.work_dir + "/" + self.project_name + ".log")

    def _convert_type_indices(self, type_indices):
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

        Currently used for CASPT2.
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
