# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
from dataclasses import dataclass
from enum import Enum
from typing import Any, List, Optional, Union


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
            key = key.lower()
            key = "".join([c for c in key if c not in chars_to_remove])
            if key.upper() not in cls.__members__:
                return ""

        return key

    @classmethod
    def get(cls, key: str) -> Any:
        """Get enum value based on key.

        Parameters
        ----------
        key : str
            method string

        Returns
        -------
        enum.name
            the enum name, which corresponds to the string

        Raises
        ------
        NotImplementedError
            if string is not an enum value
        """
        key = cls.has(key)
        if key:
            return cls(key)
        raise NotImplementedError(f"class {cls.__name__} has no value {key}")

    @classmethod
    def key_value(cls, key: str) -> Union[str, bool]:
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

    @classmethod
    def values(cls):
        """Get all available interfaces"""
        return [member.value for member in cls]


class AvailableActions(Methods):
    """All available autocas actions"""
    RUN = "run"
    INITCAS = "initcas"
    ANALYZE = "analyze"
    CREATE = "create"


class DmrgSolver(Methods):
    """All available dmrg solver"""
    QCMAQUIS = "qcmaquis"
    BLOCK2 = "block2"


class AvailableInteraces(Methods):
    """All available interfaces"""
    OPENMOLCAS = "openmolcas"
    PYSCF = "pyscf"
    MOLCAS = "molcas"


class CasMethods(Methods):
    """Stores all possible CAS methods.

    Even numbers correspond to methods with orbital
    optimization, e.g. DMRGSCF and CASSCF and odd numbers to
    methods without orbital optimization like DMRGCI and CASCI.
    Numbers lower than 100 correspond to CI methods, e.g. CASCI
    and CASSCF and higher numbers to DMRG methods like DMRGCI
    and DMRGSCF.
    """
    CASCI = "casci"
    CASSCF = "casscf"
    DMRGCI = "dmrgci"
    DMRGSCF = "dmrgscf"
    SRDFTLRDMRGCI = "srdftlrdmrgci"
    SRDFTLRDMRGSCF = "srdftlrdmrgscf"


class PostCasMethods(Methods):
    """Store all possible post-CAS methods."""
    CASPT2 = "caspt2"
    NEVPT2 = "nevpt2"


class NotInstantiatableClass:
    """Prevent class from instantiation"""
    def __new__(cls, *args, **kwargs):
        raise TypeError(f"{cls} is not instantiable")


@dataclass
class Defaults(NotInstantiatableClass):
    """Store all default values for AutoCAS."""
    @dataclass
    class DirName(NotInstantiatableClass):
        """Store default dir names"""
        initial_orbs = "initial"
        initial_dmrg = "dmrg"
        final_calc = "final"
        project_name = "autocas_project"
        dumper_name = "autocas_dump"
        scratch = "scratch"

    @dataclass
    class QcMaquisNames(NotInstantiatableClass):
        """Store default qcmaquis file names"""
        qcmaquis_result_file = "qcmaquis_result_file.h5"
        """Name of qcmaquis result file (should have h5 suffix)"""
        qcmaquis_checkpoint_dir = "qcmaquis_checkpoint.h5"
        """Name of qcmaquis checkpoint dir (should have h5 suffix for some reason)"""

    @dataclass
    class PlotNames:
        """Handle all plot names"""
        entanglement_file: str = "entanglement.pdf"
        """Name of the entanglement diagram"""
        threshold_file: str = "threshold.pdf"
        """Name of the threshold diagram"""

    @dataclass
    class Molecule(NotInstantiatableClass):
        """Store default values for Molecule objects."""
        charge: int = 0
        """Total charge of the molecule."""
        spin_multiplicity: int = 1
        """Spin multiplicity of the molecule, e.g. 2S+1."""
        double_d_shell: bool = True
        """For 3d transition metals, additionally include 4d orbitals in the initial active space"""
        unit: str = "ang"
        """Unit for the xyz file"""
        ecp_electrons: int = 0
        """Number of electrons on potential"""

    @dataclass
    class AutoCAS(NotInstantiatableClass):
        """Store default values for AutoCAS objects."""
        plateau_values: int = 10
        """Required number of consecutive threshold steps to form a plateau to determine the active space."""
        threshold_step: float = 0.01
        """One threshold step corresponds to a fraction of the maximal single orbital entropy. In %/100"""
        weak_correlation_threshold: float = 0.02
        """Any orbital with a s1 value below that threshold is directly excluded from cas."""
        single_reference_threshold: float = 0.14
        """If maximum of s1 is below that threshold a single reference method might be better."""
        large_cas_max_orbitals: int = 30
        """Settings for large active space protocol. Only used if large_cas is enabled. Maximum number of orbital
s in a sub-cas."""
        large_cas: bool = False
        """Flag to enable large cas calculations."""
        large_cas_seed: int = 42
        """Seed for random number generation in large cas."""
        large_cas_average_entanglement: bool = False
        """Flag to enable average of entanglement from sub cas calculations."""

    @dataclass
    class Interface(NotInstantiatableClass):
        """Store default values for Interface objects."""
        # interface: AvailableInteraces = AvailableInteraces.OPENMOLCAS
        interface: str = "pyscf"
        """Which interface to use."""
        dump: bool = False
        """If interface should write output and all related files."""
        init_dmrg_bond_dimension: int = 250
        """DMRG bond dimension for the active space selection."""
        init_dmrg_sweeps: int = 5
        """DMRG sweeps for the active space selection."""
        dmrg_bond_dimension: int = 3000
        """DMRG bond dimension"""
        dmrg_sweeps: int = 100
        """DMRG sweeps for the"""
        basis_set: str = "cc-pvdz"
        """The basis set."""
        dmrg_solver: str = "QCMaquis"
        """The DMRG solver to use"""
        init_cas_method: str = "dmrgci"
        """Method to evaluate the final active space."""
        cas_method: str = "dmrgci"
        """Method to evaluate the final active space."""
        post_cas_method: str = "caspt2"
        """Method to evaluate the final active space."""
        uhf: bool = False
        """Enable unrestriced HF. If your provided system is open-shell with the corresponding
        charge and/or spin multiplicity uhf is enabled automatically. """
        init_fiedler: bool = True
        """Enable fiedler ordering for DMRG in the active space selection."""
        fiedler: bool = True
        """Enable fiedler ordering for DMRG in the active space selection."""
        n_excited_states: int = 0
        """Number of states. 0 means that this option is disabled. Hence
        0 and 1 have the same meaning, that onlt the ground state is evaluated."""
        init_orbital_order: Optional[List[int]] = None
        """initial orbital order"""
        orbital_order: Optional[List[int]] = None
        """orbital order for final CAS"""


if __name__ == "__main__":
    VAL = Defaults.Molecule.unit
    Defaults.Molecule.unit = "bohr"
    X = Defaults.Molecule.unit
    Defaults.Molecule.unit = "ang"
    Z = Defaults.Molecule.unit
    print(VAL, X, Z)
    try:
        obj = Defaults.Molecule()
    except TypeError:
        print("This behavior is on purpose")
    for subclass in dir(Defaults):
        if not subclass.startswith('__'):
            print(subclass)
            for val in dir(Defaults.__dict__[subclass]):
                if not val.startswith('__'):
                    print(val)
            print()

    print("****************")
    print(CasMethods.get("DMRGCI"))
