# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
from typing import Any, List, Optional, Tuple

import numpy as np

from scine_autocas.interfaces.interface import Interface
from scine_autocas.io import FileHandler
from scine_autocas.utils.defaults import CasMethods, DmrgSolver
from scine_autocas.utils.exceptions import InputError, InterfaceError
from scine_autocas.utils.molecule import Molecule

# PySCF is not always required
try:
    from pyscf import cc, gto, lib, mcscf, scf
    from pyscf.tools.molden import load
except ModuleNotFoundError:
    pass

try:
    from scine_qcmaquis.pyscf_interface import pyscf_interface
except ModuleNotFoundError:
    pass


class PyscfInterface(Interface):
    """Interface to the electronic structure program pyscf.

    Attributes
    ----------
    settings : Settings
        The interface settings
    pyscf_wave_function : Any
        The pyscf wave function
    pyscf_mol : pyscf.gto.Mole
        the pyscf molecule

    Notes
    -----
    Available CAS methods:
        CASCI, CASSCF, DMRGCI, DMRGSCF
    Available Post-CAS methods:
        NEVPT2

    See Also
    --------
    settings : Interface.Settings
    """

    class Settings(Interface.Settings):
        """Controll the input parameters for PySCF."""

        def __init__(self, molecule: Molecule):
            """Construct class."""
            super().__init__(molecule=molecule)
            # available method for this interface
            self._add_cas_method("casci", "casscf", "dmrgscf")
            self._add_post_cas_method("nevpt2")

    # Interface
    __slots__ = (
        "pyscf_hf",
        "pyscf_mol",
        "_pyscf_mo_coeffs",
        "_pyscf_cas_solver",
    )

    def __init__(self, molecule: Molecule):
        """Construct a Pyscf interface.

        Parameters
        ----------
        molecule : Molecule
            contains molecular information
        """

        self.settings: PyscfInterface.Settings
        """provides all settings for pyscf"""

        super().__init__(molecule=molecule)

        self.pyscf_hf: Any = None
        """The current pyscf object to pass around."""
        self.pyscf_mol: Any = None
        """The current pyscf molecule to pass around."""
        self._pyscf_mo_coeffs: Optional[Any] = None
        """Store pyscf mo_coeffs, especially important for molden input."""
        self._pyscf_cas_solver: Optional[Any] = None

    def _get_mo_coeffs(self) -> np.ndarray:
        """Get MO coeffs. First search for custom coeffs from moleden file, then
        use the orbtials from hf object.

        Returns
        -------
        mo_coeffs : np.ndarray
            the molecular orbital coefficient matrix.

        Raises
        ------
        ValueError
            if neither _pyscf_mo_coeffs nor pyscf_hf are set
        """
        if self._pyscf_mo_coeffs:
            return self._pyscf_mo_coeffs
        if self.pyscf_hf:
            return self.pyscf_hf.mo_coeff

        print(
            """No MO coefficients available. Either load them from a molden file
                <pyscf_interface.init_from_molden(path_to_moldenfile)>
                or evaluate the initial orbitals."""
        )
        raise ValueError("No Mo coeffs available")

    def _final_cc_impl(self) -> List[float]:
        """abs"""
        ccsd = cc.CCSD(self.pyscf_hf)
        ccsd.run()
        trip_e = ccsd.ccsd_t()
        return [ccsd.e_tot + trip_e]

    def _final_cas_impl(self, cas_occupation: List[int], cas_indices: List[int]) -> List[float]:
        """abs"""
        norbs = len(cas_occupation)
        nelec = int(sum(cas_occupation))

        # Start from coeffs (init from molden)
        initialize = None
        if self._pyscf_mo_coeffs is not None:
            print("Initialize from mo coefficients")
            initialize = self.pyscf_mol
        else:
            print("Initialize from pyscf HF wave function")
            initialize = self.pyscf_hf

        if self.settings.cas_method in (CasMethods.CASCI, CasMethods.DMRGCI):
            cas = mcscf.CASCI(initialize, norbs, nelec)
        elif self.settings.cas_method in (CasMethods.CASSCF, CasMethods.DMRGSCF):
            cas = mcscf.CASSCF(self.pyscf_hf, norbs, nelec)

        # sort_mo by default take the 1-based orbital indices.
        cas_list = [i+1 for i in cas_indices]
        print(f"CAS orbital indices (1-based): {cas_list}")
        if self._pyscf_mo_coeffs is not None:
            orbs = mcscf.sort_mo(cas, self._pyscf_mo_coeffs, cas_list)
        else:
            orbs = cas.sort_mo(cas_list)

        if self.settings.cas_method in (CasMethods.DMRGCI, CasMethods.DMRGSCF):
            # Set dmrg as fcisolver
            if self.settings.dmrg_solver == DmrgSolver.BLOCK2:
                print("Using Block2 as DMRG solver")
                cas.fcisolver = self._setup_block2(self.settings.init_orbital_order)
            else:
                print("Using QCMaquis as DMRG solver")
                cas.fcisolver = self._setup_qcmaquis(self.settings.init_orbital_order)

        print("Start final calculation")
        with lib.capture_stdout() as out:
            energy = cas.kernel(orbs)[0]
        print("Pyscf output:")
        print(out.read())

        self._pyscf_cas_solver = cas.fcisolver

        return [energy]

    def _sanity_check_pyscf_mol_input(self):
        if self.pyscf_mol.spin + 1 != self.settings.get_molecule().spin_multiplicity:
            print("Read molden file with different spin multiplicity than the specified one")
            print(f"Spin mutliplicity from molden file: {self.pyscf_mol.spin + 1}")
            print(f"Spin mutliplicity from input:       {self.settings.get_molecule().spin_multiplicity}")

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
        if molden_file is None:
            molden_file = self.settings.get_molecule().molden_file

        print(f"Reading initial orbitals from file: {molden_file}")

        # check if molden file exists
        if not os.path.exists(molden_file):
            print(f"Molden file: {molden_file} does not exist!")
            raise InputError("Molden file not found!")

        # alternative initialization of initial orbitals.
        self.set_orbital_state(True)

        # read molden file
        self.pyscf_mol, _, mo_coeffs, = load(molden_file)[:3]

        self._sanity_check_pyscf_mol_input()

        if self.pyscf_mol.spin > 0 or len(mo_coeffs) == 2:
            print("Using only Alpha orbitals")
            self._pyscf_mo_coeffs = mo_coeffs[0]
        else:
            self._pyscf_mo_coeffs = np.asarray(mo_coeffs)  # type: ignore

        print(f"MO coefficients from {molden_file}")
        print(f"{self._pyscf_mo_coeffs}")
        if self.settings.get_molecule().ecp_electrons != 0:
            # TODO make this easier to use
            print("Assuming ecp for Mo with def2-svp")
            print("Assuming basis def2-svp")
            self.pyscf_mol.build(
                spin=self.settings.get_molecule().spin_multiplicity - 1,
                charge=self.settings.get_molecule().charge,
                basis="def2-svp",
                ecp={"Mo": "def2-svp"},
            )
        else:
            print(
                f"""Build molecule:
                xyz-file:     {self.settings.get_molecule().xyz_file}
                basis:        {self.settings.basis_set}
                spin (2S+1):  {self.settings.get_molecule().spin_multiplicity}
                xyz unit:     {self.settings.get_molecule().unit}
                total charge: {self.settings.get_molecule().charge}""")
            self.pyscf_mol.basis = self.settings.basis_set
            self.pyscf_mol.spin = self.settings.get_molecule().spin_multiplicity - 1
            self.pyscf_mol.charge = self.settings.get_molecule().charge
            self.pyscf_mol.build()
        print(f"PySCF molecule: {self.pyscf_mol}")

    def _build_molecule(self):
        """Build pyscf molecule from settings."""
        print(
            f"""Build molecule:
            xyz-file:     {self.settings.get_molecule().xyz_file}
            basis:        {self.settings.basis_set}
            spin (2S+1):  {self.settings.get_molecule().spin_multiplicity}
            xyz unit:     {self.settings.get_molecule().unit}
            total charge: {self.settings.get_molecule().charge}""")
        self.pyscf_mol = gto.Mole()
        self.pyscf_mol.build(
            atom=self.settings.get_molecule().xyz_file,
            basis=self.settings.basis_set,
            symmetry=False,
            spin=self.settings.get_molecule().spin_multiplicity - 1,
            unit=self.settings.get_molecule().unit,
            charge=self.settings.get_molecule().charge
        )
        self.pyscf_mol.verbose = 3

    def _setup_block2(self, orbital_order: Optional[List[int]] = None):
        """Settings for initial DMRG calculations."""
        raise NotImplementedError("Block2 is not implemented yet")

    def _setup_qcmaquis(self, orbital_order: Optional[List[int]] = None):
        """Settings for initial DMRG calculations."""
        qcmaquis = pyscf_interface.QcMaquis(self.pyscf_mol)
        # has to be always on
        qcmaquis.parameters.set_entropies()
        qcmaquis.parameters.set("nsweeps", self.settings.init_dmrg_sweeps)
        qcmaquis.parameters.set("max_bond_dimension", self.settings.init_dmrg_bond_dimension)
        qcmaquis.fiedler = self.settings.fiedler

        if orbital_order:
            print(f"Using user defined orbital order: {orbital_order}")
            print("Disabling Fiedler order")
            qcmaquis.parameters.set("orbital_order", orbital_order)
            qcmaquis.fiedler = False

        if self._dumper.activated:
            qcmaquis.file_path = FileHandler.current_dir
            qcmaquis.parameters.set_result_path(FileHandler.QcMaquisNames.qcmaquis_result_file)
            qcmaquis.parameters.set_checkpoint_path(FileHandler.QcMaquisNames.qcmaquis_checkpoint_dir)
        else:
            qcmaquis.file_path = None
        return qcmaquis

    def _get_s1(self, cas_object: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get entropies from dmrg based ci solver.

        Parameters
        ----------
        cas_object : Any
            the pyscf cas object.

        Returns
        -------
        s1 : np.ndarray
            single orbital entropies
        s2 : np.ndarray
            two orbital entropies
        mut_inf : np.ndarray
            mutual information
        """
        if self.settings.dmrg_solver == DmrgSolver.QCMAQUIS:
            return cas_object.fcisolver.dmrg.get_entropies()
        print(f"Unable to get entropies from {self.settings.dmrg_solver}")
        raise NotImplementedError("DmrgSolver not implemented.""")
        # TODO: combine converge_dmrg and _run_initial_dmrg

    def _initial_orbitals_impl(self) -> List[float]:
        """Calculate initial orbitals.

        Returns
        -------
        energy : float
            mean field energy
        """
        if self.pyscf_mol is None:
            self._build_molecule()

        if self.pyscf_hf is not None:
            print("Using existing hartree fock run")
            return self.pyscf_hf.e_tot

        print("Starting Hartree Fock")
        with lib.capture_stdout() as out:
            self.pyscf_hf = scf.RHF(self.pyscf_mol)
            energy = self.pyscf_hf.scf()
        print("Pyscf output:")
        print(out.read())
        return [energy]

    def _initial_cas_impl(
        self, cas_occupation: List[int], cas_indices: List[int]
    ) -> Tuple[List[float], np.ndarray, np.ndarray, np.ndarray]:
        """Run initial DMRG calculation.

        Parameters
        ----------
        cas_occupation : List[int]
            list with occupations for each orbital, e.g. 2=doubly, 1=singly, 0=unoccupied
        cas_indices : List[int]
            list with orbital indices (0-based)
        """

        norbs = len(cas_occupation)
        nelec = int(sum(cas_occupation))

        # Start from coeffs (init from molden)
        if self._pyscf_mo_coeffs is not None:
            print("Initialize from mo coefficients")
            cas = mcscf.CASCI(self.pyscf_mol, norbs, nelec)
        else:
            print("Initialize from pyscf HF wave function")
            cas = mcscf.CASCI(self.pyscf_hf, norbs, nelec)

        # sort_mo by default take the 1-based orbital indices.
        cas_list = [i+1 for i in cas_indices]
        print(f"CAS orbital indices (1-based): {cas_list}")
        if self._pyscf_mo_coeffs is not None:
            orbs = mcscf.sort_mo(cas, self._pyscf_mo_coeffs, cas_list)
        else:
            orbs = cas.sort_mo(cas_list)

        # Set dmrg as fcisolver
        if self.settings.dmrg_solver == DmrgSolver.BLOCK2:
            print("Using Block2 as DMRG solver")
            cas.fcisolver = self._setup_block2(self.settings.init_orbital_order)
        else:
            print("Using QCMaquis as DMRG solver")
            cas.fcisolver = self._setup_qcmaquis(self.settings.init_orbital_order)

        print("Start initial DMRG calculation")

        # initial energy is useless
        energy: List[float] = [0.0]
        with lib.capture_stdout() as out:
            energy = [cas.kernel(orbs)[0]]
        print("Pyscf output:")
        print(out.read())
        self._print_orbital_order(cas.fcisolver)
        s1, s2, mut_inf = self._get_s1(cas)

        self._pyscf_cas_solver = cas.fcisolver

        return energy, s1, s2, mut_inf

    def _print_orbital_order(self, dmrg_object: Any):
        if isinstance(dmrg_object, pyscf_interface.QcMaquis):
            print(f"orbital order: {dmrg_object.parameters.get('orbital_order')}")

    def _check_interface_exists(self):
        """Check that pyscf is importable."""
        try:
            # pylint: disable=unused-import,import-outside-toplevel
            import pyscf  # noqa: F401

            # pylint: enable=unused-import,import-outside-toplevel
        except ModuleNotFoundError as exc:
            error_string = "PySCF not found. Consider installing pyscf via\n\n"
            error_string += "'pip install pyscf'\n\n"
            error_string += "in order to use the pyscf interface in autocas."
            raise InterfaceError(error_string) from exc
