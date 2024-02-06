"""Module to provide the Interface base class, required for each interface.

Each interface is required to inherite from this class, to be able to be
used for an autocas calculation.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from scine_autocas.autocas_utils.molecule import Molecule


class Interface:
    """The base class for all Interfaces.

    An Interface is an interface to a corresponding program, which
    creates inputs, runs the program, reads outputs and converts these
    into autoCAS compatible containers.

    Attributes
    ----------
    settings : Settings
        handles settings of an interface and defines required settings
    dumper : Union[Dumper, Any]
        provides the ability to dump calculation
    """

    class Settings:
        """The settings class for all Interfaces.

        An instance of Settings controls the input parameters for the
        underlying electronic structure program.

        Attributes
        ----------
        basis_set : str, default = "cc-pvdz"
            basis set for a calculation
        spin_multiplicity : int
            the spin multiplicity of the molecule
        charge : int
            the total charge of the molecule
        dmrg_sweeps : int
            the number of dmrg sweeps
        dmrg_bond_dimension : int
            the dmrg bond dimension
        xyz_file : str
            path to the xyz file
        method : str
            the active space method to use
        post_cas_method : str
            the post cas method to use
        work_dir : str
            the work dir of the electronic structure program

        Notes
        -----
        The string for method and post_cas_method are stripped from characters,
        to prevent misspellings. Removed characters are:
        " ", "-", "_", "/", "."
        Additionally the string are cast to uppercase, to further prevent misspellings.
        """

        class Methods(Enum):
            """Enum class for all methods.

            This class is a base class for CasMethods and
            PostCasMethods.
            """

            @classmethod
            def has(cls, key: str) -> str:
                """Check if method is part of class.

                A method is part of a class, if it is defined in the corresponding enum.
                This function also stripps the input key from some special characters
                and casts it to uppercase to prevent misspellings.

                Parameters
                ----------
                key : str
                    method string

                Returns
                -------
                str
                    if method is in enum class, return the key, else an empty string
                """
                chars_to_remove = set([" ", "-", "_", "/", "."])
                if key is not None:
                    key = key.upper()
                    key = "".join([c for c in key if c not in chars_to_remove])
                    if key not in cls.__members__:
                        key = ""
                return key

            @classmethod
            def key_value(cls, key: str) -> Union[int, bool]:
                """Return the value of the method if key exists in enum.

                Parameters
                ----------
                key : str
                    method string

                Returns
                -------
                Union[int, bool]
                    The value if method is in enum class, False otherwise
                """
                key = cls.has(key)
                for method in cls:
                    if str(method)[11:] == key:
                        return method.value
                    if str(method)[15:] == key:
                        return method.value
                return False

        class CasMethods(Methods):
            """Stores all possible CAS methods.

            Even numbers correspond to methods with orbital
            optimization, e.g. DMRGSCF and CASSCF and odd numbers to
            methods without orbital optimization like DMRGCI and CASCI.
            Numbers lower than 100 correspond to CI methods, e.g. CASCI
            and CASSCF and higher numbers to DMRG methods like DMRGCI
            and DMRGSCF.
            """

            CASCI = 1
            CASSCF = 2
            DMRGCI = 101
            DMRGSCF = 102
            SRDFTLRDMRGCI = 1001
            SRDFTLRDMRGSCF = 1002

        class PostCasMethods(Methods):
            """Store all possible post-CAS methods."""

            CASPT2 = 1
            NEVPT2 = 2

        # settings
        __slots__ = (
            "basis_set",
            "spin_multiplicity",
            "charge",
            "dmrg_sweeps",
            "dmrg_bond_dimension",
            "xyz_file",
            "method",
            "post_cas_method",
            "work_dir",
            "n_excited_states"
        )

        def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict[str, Any]] = None):
            """Construct a settings object.

            This is just the constructor of the base class. To use it, it needs to be called by
            super().__init__(molecules=molecules, settings_dict=settings_dict)

            Parameters
            ----------
            molecule : Molecule
                the molecule object provides all required information of the molecular system.
            settings_dict : Dict[str, Any], optional
                a dict, usually provided by the input_handler, which stores attributes and corresponding values

            See Also
            --------
            settings_dict : InputHandler
            """
            self.basis_set: str = "cc-pvdz"
            """basis set for the calculation"""
            self.spin_multiplicity: int = molecules[0].spin_multiplicity
            """spin multiplicity of the system"""
            self.charge: int = molecules[0].charge
            """total charge of the system"""
            self.dmrg_sweeps: int = 5
            """number of sweeps in a DMRG calculation"""
            self.dmrg_bond_dimension: int = 250
            """maximal bond dimension (m) in a DMRG calculation"""
            self.xyz_file: str = ""
            """path to the XYZ file"""
            self.method: str = "dmrg_ci"
            """defines the method used to evaluate the initial active space"""
            self.post_cas_method: str = ""
            """defines post-CAS methods like CASPT2 and NEVPT2"""
            self.work_dir: str = os.getcwd()
            """defines the working directory where produced data is stored"""
            self.n_excited_states: int = 0
            """the number of excited states to calculate"""
            if settings_dict:
                self.apply_settings(settings_dict)

        def apply_settings(self, settings: Dict[str, Any]):
            """Apply settings from a dict.

            The dict can come from an Inputhandler instance.

            Parameters
            ----------
            settings: Dict[str, Any]
                stores class attributes as string and the corresponding values
            """
            for key in settings:
                if hasattr(self, key):
                    setattr(self, key, settings[key])

    # interface
    def __init__(self, molecules: List[Molecule], settings_dict: Optional[Dict[str, Any]] = None):
        """Construct Interface from a molecule.

        This is just the constructor of the base class. To use it, it
        needs to be called by
        super().__init__(molecules=molecules, settings_dict=settings_dict)

        Parameters
        ----------
        dumper: Dumper, optional
            if a dumper is provided in the main class it can be used
        settings: Settings
            the settings of the interface.
        """
        self.dumper: Optional[Any] = None
        """if a dumper is provided in the main class it can be used"""
        self.settings: Interface.Settings
        """contains information on the input of the electronic structure program"""
        if settings_dict:
            self.settings = self.Settings(
                molecules=molecules, settings_dict=settings_dict["settings"]
            )
        else:
            self.settings = self.Settings(molecules=molecules)

    def calculate(
        self, cas_occupation: Optional[List[int]] = None, cas_indices: Optional[List[int]] = None
    ) -> Tuple[Union[float, np.ndarray], Union[np.ndarray, List[np.ndarray]],
               Union[np.ndarray, List[np.ndarray]], Union[np.ndarray, List[np.ndarray]]]:
        """DMRG calculations including orbitals corresponding to cas_indices
        and electrons corresponding to occupations.

        Parameters
        ----------
        cas_occupations: List[int], optional
            contains the occupation for each spatial orbital. 2: doubly occupied, 1: singly occupied, 0: virtual
        cas_indices: List[int], optional
            contains the indices of the orbitals for the CAS calculation

        Notes
        -----
        Either cas_occupations and cas_indices is provided and a CAS calculation is started,
        or if none are provided a plain HF calculation is started.
        This function must be implemented by the corresponding interface.
        """
        raise NotImplementedError

    def resource_estimate(self):
        """Estimate resources for the calculation."""
        return NotImplementedError

    def get_orbital_map(self):
        """
        Getter for an orbital map in terms of orbital groups, e.g,
        [
            [[3, 4, 5], [3, 4, 6]],
            [[6], [5]],
            [[7], [7]],
            ...
        ]
        This list means that the orbitals 3, 4, and 5 of the first system are mapped to the orbitals 3, 4, and 6 of
        the second system. The orbital 6 of system 1 is mapped to orbital 5 of system 2, and the orbital 7 of system
        1 is mapped to the orbital 7 of system 2.

        Returns:
        --------
            The list of orbital groups / the orbital mpa.
        """
        return NotImplementedError
