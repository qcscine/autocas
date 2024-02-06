"""Provide interface to Serenity.

All calls are done through Serenity's Python interface.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import os

from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces import Interface


class Serenity(Interface):
    """Interface to the electronic structure program Serenity.

    Serenity is called through its Python interface.

    Notes
    -----
    The settings object is already in the __slots__ of the base interface.

    See Also
    --------
    settings : Interface.Settings
    """

    class Settings(Interface.Settings):
        """Control the input parameters and writes input file for Serenity."""

        __slots__ = (
            "uhf",
            "localisation_method",
            "alignment",
            "localize_virtuals",
            "read_and_write_to_molcas",
            "molcas_orbital_files",
            "partitioning_thresholds",
            "optimized_mapping",
            "xyz_files",
            "skip_localization",
            "score_start",
            "score_end",
            "system_names",
            "use_pi_bias",
            "write_orbital_map_file"
        )

        def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict] = None):
            """Construct class."""
            super().__init__(molecules=molecules, settings_dict=settings_dict)
            self.uhf: bool = False
            """ Run unrestricted """
            self.localisation_method: str = "IBO"
            """localisation method. Options are
            PIPEK_MEZEY, BOYS, EDMINSTON_RUEDENBERG, IBO.
            """
            self.alignment = True
            """Align orbitals for multiple structures"""
            self.localize_virtuals = True
            """Localize the virtual orbitals"""
            self.read_and_write_to_molcas = False
            """Replace the orbitals in a molcas orbital file"""
            self.molcas_orbital_files: List[str] = []
            """The path to the directory with the molcas orbital files"""
            self.partitioning_thresholds: List[float] = [3e-1, 5e-3, 1e-4]
            """Direct orbital selection mapping thresholds"""
            self.optimized_mapping = True
            """Optimize the first threshold of the direct orbital selection to give the best mapping"""
            self.xyz_files: List[str] = []
            """The list of xyz files"""
            self.skip_localization: bool = False
            """If true, the orbitals are not localized"""
            self.score_start = 5e-1
            """Start score for the threshold search"""
            self.score_end = 1e-2
            """Maximum score for the threshold search"""
            self.system_names: List[str] = []
            """The list of system names"""
            self.use_pi_bias: bool = False
            """If true the orbital comparison thresholds are scaled as presented in JCTC 16, 3607 (2020)."""
            self.write_orbital_map_file: bool = True
            """If true, Serenity will write a file containing the orbital mapping."""
            if settings_dict:
                self.apply_settings(settings_dict)

    # Serenity
    __slots__ = (
        "systems",
    )

    def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict[str, Any]] = None, load_paths=[]):
        """Construct a Serenity interface.

        All Serenity-System Controllers are constructed upon call.

        Parameters
        ----------
        molecule : Molecule
            contains molecular information
        settings_dict : Dict[str, Any], optional
            holds additional settings. Must contain the information about all xyz files.
        """
        import serenipy as spy  # type: ignore
        self.settings: Serenity.Settings
        # Build all systems
        self.systems: List[spy.System] = []
        if not load_paths:
            super().__init__(molecules=molecules, settings_dict=settings_dict)
            if not self.settings.xyz_files:
                raise ValueError("Serenity needs the XYZ files.")

            xyz_files = self.settings.xyz_files
            if len(xyz_files) != len(molecules) or len(self.settings.system_names) != len(molecules):
                raise ValueError("The number of XYZ files, molecules, and system names must the same.")
            for name, xyz_file in zip(self.settings.system_names, xyz_files):
                settings = spy.Settings()
                settings.basis.label = self.settings.basis_set
                settings.name = name
                settings.path = self.settings.work_dir
                assert self.settings.spin_multiplicity > 0
                settings.spin = self.settings.spin_multiplicity - 1
                settings.charge = self.settings.charge
                settings.geometry = xyz_file
                sys = spy.System(settings)
                self.systems.append(sys)
            self.molecules = molecules
        else:
            loaded_molecules = []
            xyz_files = []
            assert settings_dict
            for path, name in zip(load_paths, settings_dict["settings"]["system_names"]):
                settings = spy.Settings()
                settings.load = path
                settings.name = name
                settings.path = settings_dict["settings"]["work_dir"]
                sys = spy.System(settings)
                self.systems.append(sys)
                xyz_file_name = os.path.join(path, name, name + ".xyz")
                molecule = Molecule(xyz_file_name)
                loaded_molecules.append(molecule)
                xyz_files.append(xyz_file_name)
            settings_dict["settings"]["xyz_files"] = xyz_files
            self.molecules = loaded_molecules
            super().__init__(molecules=loaded_molecules, settings_dict=settings_dict)

    def calculate(
        self, cas_occupation: Optional[List[int]] = None, cas_indices: Optional[List[int]] = None
    ) -> Tuple[Union[float, np.ndarray], Union[np.ndarray, List[np.ndarray]],
               Union[np.ndarray, List[np.ndarray]], Union[np.ndarray, List[np.ndarray]]]:
        """DMRG calculations including orbitals corresponding to cas_indices
        and electrons corresponding to occupations.

        Parameters
        ----------
        cas_occupation: List[int], optional
            contains the occupation for each spatial orbital. 2: doubly occupied, 1: singly occupied, 0: virtual
        cas_indices: List[int], optional
            contains the indices of the orbitals for the CAS calculation

        Notes
        -----
        Either cas_occupations and cas_indices is provided and a CAS calculation is started,
        or if none are provided a plain HF calculation is started.
        This function must be implemented by the corresponding interface.
        """
        energies = []
        if cas_occupation is None and cas_indices is None:
            for sys in self.systems:
                energies.append(self.run_scf(sys, restricted=not self.settings.uhf))
        else:
            raise NotImplementedError
        return np.asarray(energies), [], [], []  # type: ignore

    def get_orbital_map(self) -> Tuple[List[List[List[int]]], List[List[List[int]]]]:
        """
        Getter for an orbital map in terms of orbital groups.
        See also:
            interfaces/__init__.py::get_orbital_map(self)
        """
        import serenipy as spy
        if not self.systems:
            return [], []

        if self.settings.read_and_write_to_molcas:
            self.load_or_write_molcas_orbitals()

        sys_zero = self.systems[0]
        loc_settings = spy.LocalizationTaskSettings()
        loc_settings.locType = self.get_localization_method(self.settings.localisation_method)
        loc_settings.splitValenceAndCore = True
        loc_settings.localizeVirtuals = self.settings.localize_virtuals
        loc_settings.replaceVirtuals = True
        if not self.settings.skip_localization:
            self.localize_orbitals(sys_zero, loc_settings)

        fragments = self.build_fragments(sys_zero)
        for i_sys in range(len(self.systems) - 1):
            sys = self.systems[i_sys + 1]
            if not self.settings.skip_localization:
                self.localize_orbitals(sys, loc_settings, sys_zero)
            fragments += self.build_fragments(sys)

        if self.settings.read_and_write_to_molcas:
            self.load_or_write_molcas_orbitals(True)

        return self.build_orbital_map(fragments)

    def write_molcas_orbitals(self):
        """
        Write the orbitals from Serenity to the Molcas files.
        """
        self.load_or_write_molcas_orbitals(True)

    def load_molcas_orbitals(self):
        """
        Load the orbitals from the Molcas files to Serenity.
        """
        self.load_or_write_molcas_orbitals(False)

    def load_or_write_molcas_orbitals(self, write=False):
        """
        Read or write orbitals from or to a molcas orbital file.
        Parameters:
        -----------
            write : bool
                If true, the orbitals in the molcas orbital file are replaced by Serenity's orbitals.
        """
        import serenipy as spy
        if len(self.settings.molcas_orbital_files) < len(self.systems):
            raise ValueError("Fewer loading paths for molcas orbitals then systems.")
        for i_sys, sys in enumerate(self.systems):
            loading_path = self.settings.molcas_orbital_files[i_sys]
            if self.settings.uhf:
                read = spy.ReadOrbitalsTask_U(sys)
            else:
                read = spy.ReadOrbitalsTask_R(sys)
            read.settings.fileFormat = spy.ORBITAL_FILE_TYPES.MOLCAS
            read.settings.resetCoreOrbitals = False
            read.settings.replaceInFile = write
            read.settings.path = loading_path
            read.run()

    @staticmethod
    def run_scf(sys, settings=None, restricted: bool = True) -> float:
        """
        Run a SCF calculation for the given system.#

        Parameters:
        -----------
            sys : spy.System
                The system to localize the orbitals for.
            settings : spy.ScfTaskSettings
                The Scf settings.
            restricted : bool
                Perform a restricted SCF calculation.
        """
        import serenipy as spy
        if restricted:
            scf = spy.ScfTask_R(sys)
            scf_mode = spy.SCF_MODES.RESTRICTED
        else:
            scf = spy.ScfTask_U(sys)
            scf_mode = spy.SCF_MODES.UNRESTRICTED

        if settings:
            scf.settings = settings
        scf.run()
        return sys.getEnergy(scf_mode, spy.ENERGY_CONTRIBUTIONS.HF_ENERGY)

    def localize_orbitals(self, sys, loc_settings, template=None) -> None:
        """
        Localize the orbitals of the given system. If a template system is given and orbital alignment is required,
        the orbitals are aligned first.

        Parameters:
        -----------
            sys : spy.System
                The system to localize the orbitals for.
            loc_settings : spy.LocalizationTaskSettings
                The settings for the orbital localization.
            template : spy.System
                The template system.
        """
        import serenipy as spy
        if self.settings.alignment and template:
            loc_settings.replaceVirtuals = True
            align = spy.LocalizationTask(sys, [template])
            align.settings = loc_settings
            align.settings.locType = spy.ORBITAL_LOCALIZATION_ALGORITHMS.ALIGN
            align.run()
            loc_settings.replaceVirtuals = False
        loc = spy.LocalizationTask(sys, [])
        loc.settings = loc_settings
        loc.settings.locType = self.get_localization_method(self.settings.localisation_method)
        loc.run()

    def build_fragments(self, sys):
        """
        Build fragments for the system. The number of fragments is determined by the number of partitioning thresholds
        plus one.

        Parameters:
            sys : spy.System
                The system.
        Returns:
            fragments : List[spy.System]
                The fragments.
        """
        import serenipy as spy
        name = sys.getSystemName()
        fragment_names = [name + "_" + str(i) for i in range(len(self.settings.partitioning_thresholds) + 1)]
        fragments = []
        for f_name in fragment_names:
            f_settings = sys.getSettings()
            f_settings.name = f_name
            f_settings.load = ""
            fragment = spy.System(f_settings)
            fragments.append(fragment)
        return fragments

    def build_orbital_map(self, fragments) -> Tuple[List[List[List[int]]], List[List[List[int]]]]:
        """
        Build the orbital group indices for all systems.
        Parameters:
            fragments : List[spy.System]
                The total list of fragments.
        Returns:
            Tuple[List[List[List[int]]], List[List[List[int]]]}
                The list of mappable orbital groups (vide supra), the list of unmappable orbital groups.
        """
        import serenipy as spy
        gdos = spy.GeneralizedDOSTask_R(self.systems, fragments)
        gdos.settings.similarityLocThreshold = self.settings.partitioning_thresholds
        gdos.settings.similarityKinEnergyThreshold = self.settings.partitioning_thresholds
        gdos.settings.mapVirtuals = self.settings.localize_virtuals
        gdos.settings.bestMatchMapping = self.settings.optimized_mapping
        gdos.settings.scoreEnd = self.settings.score_end
        gdos.settings.scoreStart = self.settings.score_start
        gdos.settings.usePiBias = self.settings.use_pi_bias
        gdos.settings.checkDegeneracies = True
        gdos.settings.writeGroupsToFile = self.settings.write_orbital_map_file
        gdos.run()
        return gdos.getOrbitalGroupIndices(), gdos.getUnmappableOrbitalGroupIndices()

    def get_localization_method(self, label: str):
        """
        Get the Serenity enum for the given keyword.
        """
        import serenipy as spy
        known_methods: dict = {
            "IBO": spy.ORBITAL_LOCALIZATION_ALGORITHMS.IBO,
            "PIPEK_MEZEY": spy.ORBITAL_LOCALIZATION_ALGORITHMS.PM,
            "EDMINSTON_RUEDENBERG": spy.ORBITAL_LOCALIZATION_ALGORITHMS.ER,
            "BOYS": spy.ORBITAL_LOCALIZATION_ALGORITHMS.FB
        }
        if label in known_methods:
            return known_methods[label]

        raise ValueError("Unknown orbital localization method. Known orbital localization methods are: "
                         + str(known_methods.keys()))
