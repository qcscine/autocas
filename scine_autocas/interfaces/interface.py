# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from scine_autocas.interfaces.qcmaquis_utils import QcmaquisUtils
from scine_autocas.io import FileHandler
from scine_autocas.utils.defaults import CasMethods, Defaults, DmrgSolver, PostCasMethods
from scine_autocas.utils.exceptions import DumperError
from scine_autocas.utils.molecule import Molecule


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
    _state_handler : InterfaceState
        handle the state of the current calculation
    """

    class InterfaceDumper:
        """Handle dump files for interfaces and read qcmaquis data."""

        __slots__ = ('activated', 'qcmaquis_checkpoint_name', 'qcmaquis_result_name')

        def __init__(self) -> None:
            self.activated: bool = True
            """Enable/disable dumper"""
            self.qcmaquis_checkpoint_name = FileHandler.QcMaquisNames.qcmaquis_checkpoint_dir
            """Name of qcmaquis checkpoint dir"""
            self.qcmaquis_result_name = FileHandler.QcMaquisNames.qcmaquis_result_file
            """Name of qcmaquis result file"""

        def read_initial_orbitals(self, interface: Any):
            """Read initial orbital information from previous run"""
            if interface.settings.get_molecule().molden_file:
                print(f"Reading initial orbitals from {interface.settings.get_molecule().molden_file}")
                interface.init_from_molden()
            else:
                # TODO: this should be another error
                raise NotImplementedError

        def read_initial_dmrg(self):
            """Read initial dmrg information from previous run"""
            if not self.activated:
                raise DumperError("Dumper not enabled")

            if os.path.exists(self.qcmaquis_result_name):
                print(f"Found QCMaquis result file in {FileHandler.current_dir}")
                print(f"Reading results from file: {self.qcmaquis_result_name}")

                qcmaquis_autocas = QcmaquisUtils()
                qcmaquis_autocas.read_hdf5(self.qcmaquis_result_name)
                qcmaquis_autocas.make_diagnostics()
                energy = qcmaquis_autocas.energy
                s1 = qcmaquis_autocas.s1_entropy
                s2 = qcmaquis_autocas.s2_entropy
                mut_inf = qcmaquis_autocas.mutual_information
                return energy, s1, s2, mut_inf

            raise DumperError("No QCMaquis result file found")

        def read_final_calc(self):
            """Read initial dmrg information from previous run"""
            raise NotImplementedError("Implement this if you want to use it")

    class InterfaceState:
        """Handle the state of an interface

        Attributes
        ----------
        initial_orbitals_prepared : bool, default = False
            indicate if initial orbitals are prepared
        initial_cas_prepared : bool, default = False
            indicate if initial cas was evaluated
        """

        def __init__(self) -> None:
            self.initial_orbitals_prepared: bool = False
            """state of initial orbitals"""
            self.initial_cas_prepared: bool = False
            """state of initial cas calculation"""

    class Settings:
        """The settings class for all Interfaces.

        An instance of Settings controls the input parameters for the
        underlying electronic structure program.

        Attributes
        ----------
        _molecule : Molecule
            the molecule
        basis_set : str, default = "cc-pvdz"
            basis set for a calculation
        dmrg_sweeps : int
            the number of dmrg sweeps
        dmrg_bond_dimension : int
            the dmrg bond dimension
        uhf : bool
            if unrestricted orbitals are required
        n_excited_states : int
            number of excited states
        _cas_method : CasMethods
            the active space method to use
        _post_cas_method : PostCasMethods
            the post cas method to use
        _available_cas_methods : List[CasMethods]
            all available cas methods for an interface
        _available_post_cas_methods : List[PostCasMethods]
            all available post cas methods

        Notes
        -----
        The string for method and post_cas_method are stripped from characters,
        to prevent misspellings. Removed characters are:
        " ", "-", "_", "/", "."
        Additionally the string are cast to uppercase, to further prevent misspellings.
        """

        # settings
        __slots__ = (
            "_molecule",
            "basis_set",
            "dmrg_sweeps",
            "dmrg_bond_dimension",
            "init_dmrg_sweeps",
            "init_dmrg_bond_dimension",
            "uhf",
            "n_excited_states",
            "dmrg_solver",
            "init_cas_method",
            "cas_method",
            "post_cas_method",
            "_available_cas_methods",
            "_available_post_cas_methods",
            "_available_dmrg_solver",
            "init_fiedler",
            "fiedler",
            "init_orbital_order",
            "orbital_order",
            "large_cas",
        )

        def __init__(self, molecule: Molecule, settings_dict: Optional[Dict[str, Any]] = None):
            """Construct a settings object.

            This is just the constructor of the base class. To use it, it needs to be called by
            super().__init__(molecules=molecules)

            Parameters
            ----------
            molecule : Molecule
                the molecule object provides all required information of the molecular system.
            """

            self._molecule = molecule
            """the molecule with all molecular information"""
            self.basis_set: str = Defaults.Interface.basis_set
            """basis set for the calculation"""
            self.init_dmrg_sweeps: int = Defaults.Interface.init_dmrg_sweeps
            """number of sweeps in initial DMRG calculation"""
            self.init_dmrg_bond_dimension: int = Defaults.Interface.init_dmrg_bond_dimension
            """maximal bond dimension (m) in initial DMRG calculation"""
            self.dmrg_sweeps: int = Defaults.Interface.dmrg_sweeps
            """number of sweeps in a DMRG calculation"""
            self.dmrg_bond_dimension: int = Defaults.Interface.dmrg_bond_dimension
            """maximal bond dimension (m) in a DMRG calculation"""
            self.n_excited_states: int = Defaults.Interface.n_excited_states
            """the number of excited states to calculate"""
            self.uhf: bool = Defaults.Interface.uhf
            """the number of excited states to calculate"""
            self.init_cas_method: CasMethods = CasMethods.DMRGCI
            """defines the method used to evaluate the initial active space"""
            self.cas_method: CasMethods = CasMethods.get(Defaults.Interface.cas_method)
            """defines the method used to evaluate the active space"""
            self.post_cas_method: PostCasMethods
            """defines post-CAS methods like CASPT2 and NEVPT2"""
            self.dmrg_solver: DmrgSolver = DmrgSolver.get(Defaults.Interface.dmrg_solver)
            """defines the dmrg solver."""
            self.init_fiedler: bool = Defaults.Interface.init_fiedler
            """Enable fiedler ordering for DMRG in the active space selection."""
            self.fiedler: bool = Defaults.Interface.fiedler
            """Enable fiedler ordering for DMRG in the active space selection."""
            self._available_cas_methods: List[CasMethods] = [CasMethods.DMRGCI]
            """available CAS methods of an interface."""
            self._available_post_cas_methods: List[PostCasMethods] = []
            """available Post CAS methods of an interface."""
            self._available_dmrg_solver: List[DmrgSolver] = [DmrgSolver.QCMAQUIS]
            """available DMRG solver of an interface."""
            self.init_orbital_order: Optional[List[int]] = Defaults.Interface.init_orbital_order
            """orbital order for initial CAS calculation"""
            self.orbital_order: Optional[List[int]] = Defaults.Interface.orbital_order
            """orbital order for initial CAS calculation"""
            self.large_cas: bool = False
            self.apply_settings(settings_dict)

        def get_molecule(self) -> Molecule:
            """Getter for molecule"""
            return self._molecule

        def _add_dmrg_solver(self, *methods):
            """Add method to available dmrg solvers.

            Parameters
            ----------
            method : str
                The name of the method

            Notes
            -----
            see utils.defaults.DmrgSolver

            Raises
            ------
            NotImplementedError
                if method is not available
            """
            for method in methods:
                if DmrgSolver.has(method):
                    self._available_dmrg_solver.append(DmrgSolver.get(method))
                else:
                    raise NotImplementedError(f"The solver: {method} is not available")

        def _add_post_cas_method(self, *methods):
            """Add method to available methods.

            Parameters
            ----------
            method : str
                The name of the method

            Notes
            -----
            see utils.defaults.PostCasMethods

            Raises
            ------
            NotImplementedError
                if method is not available
            """
            for method in methods:
                if PostCasMethods.has(method):
                    self._available_post_cas_methods.append(PostCasMethods.get(method))
                else:
                    raise NotImplementedError(f"The method {method} is not available")

        def _add_cas_method(self, *methods):
            """Add method to available methods.

            Parameters
            ----------
            method: str
                The name of the method

            Notes
            -----
                see utils.defaults.CasMethods

            Raises
            ------
            NotImplementedError
                if method is not available
            """
            for method in methods:
                if CasMethods.has(method):
                    self._available_cas_methods.append(CasMethods.get(method))
                else:
                    raise NotImplementedError(f"The method {method} is not available")

        def get_cas_method(self) -> List[CasMethods]:
            """Get available methods."""
            return self._available_cas_methods

        def get_post_cas_method(self) -> List[PostCasMethods]:
            """Get available methods."""
            return self._available_post_cas_methods

        def apply_settings(self, settings: Optional[Dict[str, Any]] = None):
            """Apply settings from a dict.

            The dict can come from an Inputhandler instance.

            Parameters
            ----------
            settings: Dict[str, Any]
                stores class attributes as string and the corresponding values
            """
            if settings:
                for key in settings:
                    if hasattr(self, key):
                        setattr(self, key, settings[key])

    # interface
    __slots__ = (
        "settings",
        "_dumper",
        "_state_handler"
    )

    def __init__(self, molecule: Molecule, settings_dict: Optional[Dict[str, Any]] = None):
        """Construct Interface from a molecule.

        This is just the constructor of the base class. To use it, it
        needs to be called by
        super().__init__(molecules=molecules)

        Parameters
        ----------
        molecule : Molecule
            the molecule object
        """

        # super().__init__(molecule)
        self._dumper = self.InterfaceDumper()
        """if a dumper is provided in the main class it can be used"""
        self.settings = self.Settings(molecule=molecule)
        """contains information on the input of the electronic structure program"""
        self._state_handler = self.InterfaceState()
        """handle the state of the interface"""
        self.setup_from_dict(settings_dict)

    def set_orbital_state(self, state: bool):
        """Set state handler orbital state

        Parameters
        ----------
        state : bool
            True if orbitals are prepared, else False
        """
        self._state_handler.initial_orbitals_prepared = state

    def set_initial_cas_state(self, state: bool):
        """Set state handler cas state

        Parameters
        ----------
        state : bool
            True if initial cas is evaluated
        """
        self._state_handler.initial_cas_prepared = state

    # has to be implemented by new interface
    def _check_interface_exists(self):
        """Sanity check to validate backends exist."""
        raise NotImplementedError("This is a virtual class")

    # has to be implemented by new interface
    def _initial_orbitals_impl(self) -> List[float]:
        """Calculate initial orbitals (usually with HF)

        Returns
        -------
        energy : float
            mean field energy

        Raises
        ------
        NotImplementedError
            This method has to be overwritten by new interface
        """
        raise NotImplementedError("An interface has to implement this method")

    # has to be implemented by new interface
    def _initial_cas_impl(
        self, cas_occupation: List[int], cas_indices: List[int]
    ) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
        """Calculate initial active space with DMRG

        Parameters
        ----------
        cas_occupation : List[int]
            occupation for each orbital in cas
        cas_indices : List[int]
            orbital index (0-based) for each orbital in CAS

        Returns
        -------
        energy : float
            mean field energy
        s1 : np.array
            single orbital entropies
        s2 : np.array
            two orbital entropies
        mut_inf : np.array
            mutual information

        Raises
        ------
        NotImplementedError
            This method has to be overwritten by new interface
        """
        raise NotImplementedError("An interface has to implement this method")

    # has to be implemented by new interface
    def _final_cas_impl(
        self, cas_occupation: List[int], cas_indices: List[int]
    ) -> List[float]:
        """Calculate final active space

        Parameters
        ----------
        cas_occupation : List[int]
            occupation for each orbital in cas
        cas_indices : List[int]
            orbital index (0-based) for each orbital in CAS

        Returns
        -------
        energy : float
            mean field energy
        s1 : np.array
            single orbital entropies
        s2 : np.array
            two orbital entropies
        mut_inf : np.array
            mutual information

        Notes
        -----
        s1, s2 and mut_inf are optional and depend on the method used in the final calculation.

        Raises
        ------
        NotImplementedError
            This method has to be overwritten by new interface
        """
        raise NotImplementedError("An interface has to implement this method")

    # has to be implemented by new interface
    def _final_cc_impl(self) -> List[float]:
        """Calculate final energy for single reference systems (usually with CC)

        Returns
        -------
        energy : float
            mean field energy

        Raises
        ------
        NotImplementedError
            This method has to be overwritten by new interface
        """
        raise NotImplementedError("An interface has to implement this method")

    def init_from_molden(self, molden_file: Optional[str] = None):
        """Initialize molecule and mo coeffs from molden file.

        Parameters
        ----------
        molden_file: str
            molden file to read mo coeffs, etc.

        Raises
        ------
        InputError
            if no molden file is found
        """
        raise NotImplementedError("Interface has to implement this")

    def calculate(
        self, cas_occupation: Optional[List[int]] = None, cas_indices: Optional[List[int]] = None
    ) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
        """

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

        self._check_interface_exists()

        # check input
        cas_occupation, cas_indices = self._sanity_check(cas_occupation, cas_indices)

        # evaluate orbitals, or read orbitals from molden file
        energy: List[float] = [0.0]
        if not self._state_handler.initial_orbitals_prepared:
            print("Preparing initial orbitals")
            try:
                self._dumper.read_initial_orbitals(self)
            except DumperError:
                energy = self._initial_orbitals_impl()
            except NotImplementedError:
                energy = self._initial_orbitals_impl()
            self._state_handler.initial_orbitals_prepared = True

        if not cas_occupation and not cas_indices:
            return energy, np.array([]), np.array([]), np.array([])

        # initial CAS calculation
        if not self._state_handler.initial_cas_prepared and self._state_handler.initial_orbitals_prepared:
            print("Preparing single orbitals entropies from initial CAS")
            try:
                energy, s1, s2, mut_inf = self._dumper.read_initial_dmrg()
            except DumperError:
                energy, s1, s2, mut_inf = self._initial_cas_impl(cas_occupation, cas_indices)
            except NotImplementedError:
                energy, s1, s2, mut_inf = self._initial_cas_impl(cas_occupation, cas_indices)
            if not self.settings.large_cas:
                self._state_handler.initial_cas_prepared = True
            return energy, s1, s2, mut_inf

        print("Preparing final calculation")
        if len(cas_occupation) != 0:
            # Final CAS calculation
            try:
                energy = self._dumper.read_final_calc()
            except DumperError:
                energy = self._final_cas_impl(cas_occupation, cas_indices)
            except NotImplementedError:
                energy = self._final_cas_impl(cas_occupation, cas_indices)
            return energy, np.array([]), np.array([]), np.array([])

        # No CAS required, CC calculation
        try:
            energy = self._dumper.read_final_calc()
        except DumperError:
            energy = self._final_cc_impl()
        except NotImplementedError:
            energy = self._final_cc_impl()
        return energy, np.array([]), np.array([]), np.array([])

    def _sanity_check(
            self, cas_occupation: Optional[List[int]], cas_indices: Optional[List[int]]
    ) -> Tuple[List[int], List[int]]:
        """Checks length of occupation and indices list

        Parameters
        ----------
        cas_occupation : List[int]
            occupation of each orbital
        cas_indices : List[int]
            orbital index (0-based) for each orbital

        Raises
        ------
        ValueError
            if length of lists do not match
        """
        if cas_occupation and cas_indices:
            if len(cas_occupation) != len(cas_indices):
                print(f"""length of provided occupation vector {cas_occupation}
                             does not match length of provided indices {cas_indices}""")
                raise ValueError("input list have different length")
            return cas_occupation, cas_indices
        return [], []

    def initial_cas_prepared(self):
        """Change the state of the interface"""
        self._state_handler.initial_cas_prepared = True

    def set_post_cas_method(self, new_method: str):
        """Set a CAS method to evaluate the active space with.

        Parameters
        ----------
        new_method : str
            name of the method
        """
        has_method = False

        for method in self.settings.get_post_cas_method():
            if PostCasMethods.has(method.name) == new_method:
                has_method = True

        if has_method:
            self.settings.post_cas_method = PostCasMethods.get(new_method)
        else:
            print(f"This interface only supports the post cas methods: {self.settings.get_post_cas_method()}")
            raise NotImplementedError(f"This interface does not support {new_method}")

    def set_cas_method(self, new_method: str):
        """Set a CAS method to evaluate the active space with.

        Parameters
        ----------
        new_method : str
            name of the method
        """
        has_method = False
        for method in self.settings.get_cas_method():
            if CasMethods.has(method.name) == new_method:
                has_method = True

        if has_method:
            self.settings.cas_method = CasMethods.get(new_method)
        else:
            # logging.error(f"This interface only supports the cas methods: {self.settings.get_cas_method()}")
            print(f"This interface only supports the cas methods: {self.settings.get_cas_method()}")
            raise NotImplementedError(f"This interface does not support {new_method}")

    def setup_from_dict(self, settings: Optional[Dict[str, Any]] = None):
        """Apply settings from a dict to an interface object.

        Parameters
        ----------
        settings : Dict[str, Any]
            settings dict

        Notes
        -----
        The keys have to be the same name as the attributes.
        if a setting is an attribute from the interface.settings object,
        it is also set here.

        Examples
        --------
        >>>from scine_autocas import Molecule
        >>>from scine_autocas.interfaces import PyscfInterface
        >>>molecule = Molecule("mol.xyz")
        >>>settings = {"init_bond_dimension": 400, "basis" : "cc-pvdz"}
        >>>interface = PyscfInterface(molecule)
        >>>pyscf.setup_from_dict(settings)
        """
        if settings:
            interface_settings = settings["Interface"]

            for key in interface_settings:
                if hasattr(self, key):
                    setattr(self, key, interface_settings[key])
                elif hasattr(self.settings, key):
                    if key == "init_cas_method":
                        setattr(self.settings, key, CasMethods.get(interface_settings[key]))
                    elif key == "cas_method":
                        self.set_cas_method(interface_settings[key])
                    elif key == "post_cas_method":
                        self.set_post_cas_method(interface_settings[key])
                    elif key == "dmrg_solver":
                        setattr(self.settings, key, DmrgSolver.get(interface_settings[key]))
                    else:
                        setattr(self.settings, key, interface_settings[key])

    def get_orbital_map(self):
        """Getter for an orbital map in terms of orbital groups
        [[[3, 4, 5], [3, 4, 6]],
        [[6], [5]],
        [[7], [7]],
        ...]
        This list means that the orbitals 3, 4, and 5 of the first system are mapped to the orbitals 3, 4, and 6 of
        the second system. The orbital 6 of system 1 is mapped to orbital 5 of system 2, and the orbital 7 of system
        1 is mapped to the orbital 7 of system 2.

        Returns
        -------
        The list of orbital groups / the orbital mpa.
        """
        return NotImplementedError
