# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import math
import os
import shutil
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from scine_autocas.interfaces.interface import Interface
from scine_autocas.interfaces.qcmaquis_utils import QcmaquisUtils
from scine_autocas.utils import Molecule

from .cq_hdf5_reader import CQ_HDF5_Reader


class ChronusQ(Interface):

    class Settings(Interface.Settings):

        __slots__ = (
            "project_name",
            "scratch_dir",
            "input_file",
            "hdf5_file",
            "prev_hdf5_file",
            "out_file",
            "memory",
            "exe_file",
            "molecule",
            "max_iter",
            "broombridge_file",
            "broombridge_only",
            "frozen_scf",
            "memory_factor"
        )

        def __init__(self, molecule: Molecule, settings_dict: Optional[Dict] = None):
            """Construct Chronus Quantum Settings Class
            """
            super().__init__(molecule, settings_dict)

            # Set Default Settings
            self.project_name: str = "autocas"
            """name of the project, which is used for naming the generated files"""
            self.scratch_dir: str = ""
            """String that contains the scratch directory where the ChronusQ calculation is performed"""
            self.input_file: str = self.project_name + ".inp"
            """The name of the input file to be generated for ChronusQ"""
            self.hdf5_file: str = self.project_name + ".bin"
            """File name where ChronusQ will write its hdf5 file"""
            self.prev_hdf5_file: str = ""
            """String variable to store the path to the previous CQ HDF5 file."""
            self.out_file: str = self.project_name + ".out"
            """File name where ChronusQ will write its log file"""
            self.memory: str = "2GB"
            """String variable to store the amount of memory to be allocated during Chronus Quantum calculation."""
            self.exe_file: str = ""
            """String variable to store the executable path"""
            self.max_iter: int = 20
            """Variable that determines the number of GMRES iterations in CASPT2"""
            self.broombridge_file: str = ""
            """Variable to store the name of broombridge file to save."""
            self.broombridge_only: bool = False
            """Whether to exit the ChronusQ program after the broombridge file is saved."""
            self.frozen_scf: bool = False
            """Whether to skip the orbital optimization at the Hatree-Fock level (using a guess is greatly advised)"""
            self.memory_factor: int = 6
            """Memory factor for CASPT2 memory model"""

            if settings_dict:
                self.apply_settings(settings_dict)

    __slots__ = (
        "qc_maquis",
        "_job_iteration_cnt",
        "num_mo",
    )

    def __init__(self, molecule: Molecule, settings_dict: Optional[Dict[str, Any]] = None):
        """Construct a Chronus Quantum interface.

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

        self.settings = ChronusQ.Settings

        super().__init__(molecule=molecule, settings_dict=settings_dict)

        # Remove underscore from default method for use in ChronusQ
        self.settings.method = "dmrgci"

        self.qc_maquis = QcmaquisUtils()
        """builds s1, s2 and Ipq from QCMaquis HDF5 file"""
        self.hdf_reader = CQ_HDF5_Reader()
        """Object to read information from the ChronusQ HDF5 file"""
        self.num_mo: int = 0
        """Variable that stores the number of molecular orbitals upon reading the ChronusQ hdf5 file"""
        self._job_iteration_cnt: int = 0
        """Variable that stores the number of jobs that have been proformed so far in order to place
            the calculations in different directories"""
        self.molecule: Molecule = molecule
        """Object that stores the molecule information."""

    def update_file_names(self):
        self.settings.input_file = self.settings.project_name + ".inp"
        self.settings.hdf5_file = self.settings.project_name + ".bin"
        self.settings.out_file = self.settings.project_name + ".out"

    def create_post_hf_input_file(self, cas_occupation: List[int], cas_indices: List[int]):
        """Creates the CQ input file for the calculation"""
        self.update_file_names()
        with open(self.settings.scratch_dir+self.settings.input_file, 'w') as f:
            self.__write_geom(f)

            # Write Reference information
            f.write("\n[QM]\n")
            if self.settings.post_cas_method and self.settings.post_cas_method.upper() == "CASPT2":
                f.write("reference = GHF\n")
            elif self.settings.spin_multiplicity == 1:
                f.write("reference = RHF\n")
            else:
                f.write("reference = ROHF\n")

            # Write Job Type
            if self.settings.post_cas_method and self.settings.post_cas_method.upper() == "CASPT2":
                f.write("job = perturb\n")
            else:
                f.write("job = mcscf\n")

            # Write Basis set info
            f.write("\n[BASIS]\n")
            f.write("basis = " + self.settings.basis_set + "\n")

            # Add SCF terms to make it more likely to converge
            f.write("\n[SCF]\n")
            f.write("EneTol = 1E-8\n")
            f.write("DenTol = 1E-7\n")
            f.write("FDCTol = 1E-5\n")
            if len(self.settings.prev_hdf5_file) > 0:
                f.write("guess = readmo\n")
            if self.settings.post_cas_method and self.settings.post_cas_method.upper() == "CASPT2":
                f.write("alg = Skip\n")
            else:
                f.write("maxIter = 2000\n")
                f.write("nKeep = 20\n")
                f.write("diisAlg = CDIIS\n")
            if self.settings.frozen_scf:
                f.write("Alg = skip\n")

            # Write MCSCF Section
            f.write("\n[MCSCF]\n")
            if self.settings.method.lower() == "dmrg_ci":
                f.write("JobType = dmrgci\n")
            else:
                f.write("JobType = " + self.settings.method + "\n")
            f.write("NACTE = %i\n" % (sum(cas_occupation)))
            f.write("NRoots = 1\n")
            f.write("MaxSCFIter = 3000\n")
            if self.settings.method.upper() == "DMRG_CI" or self.settings.method == "DMRG_SCF":
                f.write("MaxM = " + str(self.settings.dmrg_bond_dimension) + "\n")
                f.write("MaxSweeps = " + str(self.settings.dmrg_sweeps) + "\n")
                f.write("ENTROPY = True\n")
            if (self.settings.post_cas_method and self.settings.post_cas_method.upper() == "CASPT2"):
                caspt2_indices = self.compute_ghf_indices(cas_indices)
                f.write("NACTO = " + str(len(caspt2_indices)) + "\n")
            else:
                f.write("NACTO = " + str(len(cas_indices)) + "\n")
                f.write("CASOrbital = " + self.compute_cas_orbital_string(cas_indices) + "\n")

            # Handle Broombridge options
            if self.settings.broombridge_only or len(self.settings.broombridge_file) > 0:
                if len(self.settings.broombridge_file) == 0:
                    self.settings.broombridge_file = self.project_name + ".brm"
                f.write("Broombridge = %s\n" % (self.settings.broombridge_file))
                if self.settings.broombridge_only:
                    f.write("saveBroombridgeOnly = True\n")

            # Handle Perturbation theory options
            if self.settings.post_cas_method and self.settings.post_cas_method.upper() == "CASPT2":
                f.write("\n[PERTURB]\n")
                f.write("MaxIter = %5i\n" % (self.settings.max_iter))
                frozen_core, frozen_virtual = self.caspt2_frozen_orbitals(len(caspt2_indices), int(sum(cas_occupation)))
                f.write("FrozenCore = " + str(frozen_core) + "\n")
                f.write("FrozenVirtual = " + str(frozen_virtual) + "\n")

            self.__write_misc(f)

    def create_hf_input_file(self, cas_occupation: List[int], cas_indices: List[int]):
        """Creates the CQ input file for the calculation"""
        self.update_file_names()
        with open(self.settings.scratch_dir+self.settings.input_file, 'w') as f:
            self.__write_geom(f)

            # Write Reference information
            f.write("\n[QM]\n")
            if self.settings.spin_multiplicity == 1:
                f.write("reference = RHF\n")
            else:
                f.write("reference = ROHF\n")

            # Write Job Type
            f.write("job = scf\n")

            # Write Basis set info
            f.write("\n[BASIS]\n")
            f.write("basis = " + self.settings.basis_set + "\n")

            # Add SCF terms to make it more likely to converge
            f.write("\n[SCF]\n")
            f.write("maxIter = 2000\n")
            f.write("nKeep = 20\n")
            f.write("diisAlg = CDIIS\n")

            self.__write_misc(f)

    def __write_geom(self, f):
        """Function to write the [MOLECULE] section of the ChronusQ input file"""
        f.write("[Molecule]\n")
        f.write("charge = " + str(self.settings.charge) + "\n")
        f.write("mult = " + str(self.settings.spin_multiplicity) + "\n")
        f.write("Geom:\n")
        lines = []
        with open(self.settings.xyz_file, 'r') as g:
            lines = g.readlines()
        lines = lines[2:]
        for line in lines:
            f.write("   " + line)

    def __write_misc(self, f):
        """Function to write the [MISC] section of the ChronusQ input file"""
        f.write("\n[MISC]\n")
        nThreads = "1"
        if "OMP_NUM_THREADS" in os.environ.keys():
            nThreads = os.environ["OMP_NUM_THREADS"]
        f.write("nSMP = " + nThreads + "\n")
        f.write("mem = " + self.settings.memory + "\n")

    def compute_cas_orbital_string(self, cas_indices: List[int]):
        """takes the cas_indices and creates a comma separated string to use as input to chronusq"""
        cas_string = ""
        for i in cas_indices:
            cas_string += str(i+1)+","
        cas_string = cas_string[:-1]
        return cas_string

    def compute_num_caspt2_dets(self, nV: int, nE: int, nAO: int, nAE: int):
        dets = 0
        if nAO-nAE > 1 and nAE > 1:
            for i in range(3):
                for j in range(3):
                    dets += math.comb(nE, i)*math.comb(nAO, nAE+i-j)*math.comb(nV, j)
        elif nAE < 1 or nAE == nAO:
            # No active space only singles and doubles
            dets += math.comb(nE, 1)*math.comb(nV, 1)
            dets += math.comb(nE, 2)*math.comb(nV, 2)
        else:
            raise RuntimeError("Invalid active space for caspt2")

        return dets

    def compute_caspt2_memory(self, nV: int, nE: int, nAO: int, nAE: int, maxIter: int):
        """Computes the amount of memory needed for caspt2"""

        # Compute the amount of memory from MO Integrals
        memory = 0
        num_orbitals = nV + nE + nAO + nAE
        if num_orbitals > 0:
            memory += 2*np.math.pow(num_orbitals, 4)

        # Compute total number of determinants
        memory += int(self.settings.memory_factor)*maxIter*self.compute_num_caspt2_dets(nV, nE, nAO, nAE)
        print("Memory Factor = " + str(self.settings.memory_factor))

        # Add buffer of 800MB
        # Note: right now memory is the number complex variables we are allocating
        memory += 50000000

        return int(16*memory)

    def caspt2_frozen_orbitals(self, nAO: int, nAE: int):
        """Computes the number of frozen core and frozen virtual orbitals to be used for
        CASPT2 based on the amount of memory requested
        """
        nCore = self.molecule.electrons - nAE
        raise RuntimeError("Fix above line")

        if self.num_mo == 0:
            raise RuntimeError("The number of MOs was never determined")
        nVirt = self.num_mo*2-nCore-nAO
        ratio = math.ceil(nVirt/nCore)
        memory_limit = self.memory_limit()
        for i in range(nCore):
            nV = nVirt - i*ratio
            nC = nCore - i
            if nV < 0:
                break
            mem = self.compute_caspt2_memory(nV, nC, nAO, nAE, self.max_iter)
            if mem < memory_limit:
                # return the number of frozen core and frozen virtual
                return i, i*ratio
        raise RuntimeError(
            "Could not find an optimal number of frozen orbitals for CASPT2 (perhaps you can request more memory)")

    def memory_limit(self):
        """Determines the integer number of bytes that have been requested by the memory member variable"""
        conv = 0
        substring = ""
        gb_pos = self.settings.memory.upper().find('GB')
        mb_pos = self.settings.memory.upper().find('MB')
        kb_pos = self.settings.memory.upper().find('KB')
        if gb_pos != -1:
            conv = 1E9
            substring = self.settings.memory[:gb_pos]
        elif mb_pos != -1:
            conv = 1E6
            substring = self.settings.memory[:mb_pos]
        elif kb_pos != -1:
            conv = 1E3
            substring = self.settings.memory[:kb_pos]
        else:
            raise RuntimeError("Could not determine the units used for the requested memory limit (e.g. GB,MB,KB)")

        return int(conv*float(substring))

    def compute_ghf_indices(self, cas_indices: List[int]):
        """Takes the orbital indices from a restricted calculation and converts them
        to the indices used when this is converted to GHF. THis is neccessary because
        ChronusQ's CASPT2 is only implemented for 2 component wave functions
        """
        caspt2_indices = []
        for i in cas_indices:
            caspt2_indices.append(2*i)
            caspt2_indices.append(2*i+1)
        return caspt2_indices

    def set_chronusq_env(self, cas_occupation: List[int], cas_indices: List[int]):
        """Creates the environment for ChronusQ calculation
        """
        # Check for CHRONUSQ environment variable
        self.settings.exe_file = os.environ["CHRONUSQ"]
        if not self.settings.exe_file:
            raise EnvironmentError(
                "Chronus Quantum is not installed or 'CHRONUSQ' system variable is not set"
            )

        # Create scratch directory
        if not os.path.isdir(self.settings.work_dir + '/scratch'):
            os.mkdir(self.settings.work_dir + '/scratch')

        # Create step folder in scratch directory
        step_dir = "/scratch/Step" + str(self._job_iteration_cnt) + "/"
        self.settings.scratch_dir = self.settings.work_dir + step_dir
        if not os.path.isdir(self.settings.scratch_dir):
            os.mkdir(self.settings.scratch_dir)

        # Create CQ Input file
        if cas_occupation is None and cas_indices is None:
            self.create_hf_input_file(cas_occupation, cas_indices)
        else:
            self.create_post_hf_input_file(cas_occupation, cas_indices)

    def run_chronusq(self):
        """Performs the ChronusQ calculation in the scratch directory
        """
        # Change working directory to scratch
        os.chdir(self.settings.scratch_dir)

        # Create ChronusQ command string
        exe_cmd = ""
        exe_cmd += self.settings.exe_file
        exe_cmd += " -i " + self.settings.input_file
        exe_cmd += " -o " + self.settings.out_file
        exe_cmd += " -b " + self.settings.hdf5_file

        # Use previous calculation as guess
        if len(self.settings.prev_hdf5_file) > 0:
            exe_cmd += " -s " + self.settings.prev_hdf5_file

        # Perform CHronusQ calculation
        if os.system(exe_cmd) != 0:
            raise RuntimeError("Chronus Quantum Job failed. See " + self.scratch_dir +
                               " for details.\n Command: %s" % (exe_cmd))

        # Save path to bin file
        self.settings.prev_hdf5_file = self.settings.scratch_dir + self.settings.hdf5_file

        # Move the broombridge file back to working directory.
        if len(self.settings.broombridge_file) > 0:
            print(self.settings.settings.work_dir + "/" + self.settings.broombridge_file)
            shutil.copyfile(self.settings.broombridge_file, self.settings.work_dir +
                            "/" + self.settings.broombridge_file)

        # Change back to previous working dir
        os.chdir(self.settings.work_dir)

        # Increment the job counter
        self._job_iteration_cnt += 1

    def calculate(self, cas_occupation: Optional[List[int]] = None, cas_indices: Optional[List[int]] = None
                  ) -> Tuple[Union[float, np.ndarray], Union[np.ndarray, List[np.ndarray]],
                             Union[np.ndarray, List[np.ndarray]], Union[np.ndarray, List[np.ndarray]]]:
        """Calculate a system with the correpsonding setting.

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

        is_hf = (cas_indices is None and cas_occupation is None)

        # Perform casscf calculation using RHF reference instead
        # GHF. Then perform CASPT2 using GHF reference.
        if (self.settings.method.lower() == "casscf" and self.settings.post_cas_method.lower() == "caspt2"):
            self.settings.post_cas_method = ""
            self.calculate(cas_occupation, cas_indices)
            self.settings.post_cas_method = "caspt2"
            self.settings.method = "casci"

        # Set ChronusQ environment
        self.set_chronusq_env(cas_occupation, cas_indices)

        # Generate input file
        self.run_chronusq()

        # Read ChronusQ HDF5 file
        energy = 0.
        if not self.settings.broombridge_only and not is_hf:
            self.hdf_reader.file_name = self.settings.scratch_dir + self.settings.hdf5_file
            self.hdf_reader.read_hdf5()
            energy = self.hdf_reader.energy
            self.num_mo = self.hdf_reader.num_mo

        # Read QCMaquis data
        if (
            (self.settings.method == "dmrg_ci" or self.settings.method ==
             "dmrg_scf") and not self.settings.broombridge_only and not is_hf
        ):
            self.qc_maquis.read_hdf5(self.settings.scratch_dir+"res.0.h5")
            self.qc_maquis.make_diagnostics()

        if cas_indices is not None and cas_occupation is not None:
            return energy, self.qc_maquis.s1_entropy, self.qc_maquis.s2_entropy, self.qc_maquis.mutual_information
        else:
            return
