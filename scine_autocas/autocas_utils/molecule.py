"""A molecular data container.

This module implementes the Molecule class, which stores all molecular
system related data, which is required to set up calculations with
autoCAS.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


from typing import Any, Dict, List, Optional

from scine_autocas.chemical_elements import Elements


class Molecule:
    """A class to provide useful information of a molecule.

    A molecule object stores all system dependend information to set up
    a Hartree-Fock and post-HF calculations.

    Attributes
    ----------
    atoms : List[str]
        contains the chemical element string for each atom in a molecule
    core_orbitals : int
        the number of core orbitals
    valence_orbitals : int
        the number of valence orbitals (occupied and virtual)
    electrons : int
        the number of electrons
    charge : int
        the total charge of the molecule
    spin_multiplicity : int
        the total spin multiplicity (2S+1) of the molecule
    double_d_shell : bool
        flag to add 4d orbitals to the number of valence orbitals for
        3d transition metals
    occupation : List[int]
        stores the occupation of all core and valence orbitals,
        e.g. 2 for doubly, 1 for singly and 0 for unoccupied.
    """

    __slots__ = (
        "atoms",
        "core_orbitals",
        "valence_orbitals",
        "electrons",
        "occupation",
        "charge",
        "spin_multiplicity",
        "double_d_shell",
        "n_basis_functions",
    )

    def __init__(self, xyz_file: Optional[str] = None, atoms: Optional[List[str]] = None, charge: int = 0,
                 spin_multiplicity: int = 1, n_basis_functions: int = 0, double_d_shell: bool = True,
                 settings_dict: Optional[Dict[str, Any]] = None):
        """Construct a molecule from any source of information.

        At creation a molecule is first assumed to be in the lowest possible spin state
        with no charge, and later updated. See update function.

        Parameters
        ----------
        xyz_file : str, optional
            the location of an xyz-file storing the molecular geometry
        atoms : List[str], optional
            a list containing element symbols for each atom in the molecule
        charge : int, optional
            the total charge of the molecule
        spin_multiplicity : int, optional
            the spin_multiplicity (2s+1) of the molecule
        n_basis_functions : int, optional
            the number of basis functions if more virtual orbitals should
            be part of the valence space
        double_d_shell : bool, optional
            flag to add 4d orbitals to total number of valence orbitals for
            3rd row transition metals
        settings_dict : Dict[str, Any], optional
            holds all settings provided by a yaml input

        Notes
        -----
        At least one of "xyz_file", "atoms" or "settings_dict" (with similar information) has to be provided.
        For unrestriced calculations double occupations are still assumed in the occupation.

        See Also
        --------
        settings_dict : InputHandler

        Raises
        ------
        ValueError
            If no atomic information is provided.
        """
        self.atoms: List[str] = []
        """contains the chemical element string for each atom in a molecule"""
        self.core_orbitals: int = -1
        """the number of core orbitals"""
        self.valence_orbitals: int = 0
        """the number of valence orbitals (occupied and virtual)"""
        self.electrons: int = 0
        """the number of electrons"""
        self.charge: int = 0
        """the total charge of the molecule"""
        self.spin_multiplicity: int = 1
        """the total spin multiplicity (2S+1) of the molecule"""
        self.n_basis_functions: int = 0
        """number of basis functions"""
        self.double_d_shell: bool = False
        """flag to add 4d orbitals to the number of valence orbitals for 3d transition metals"""
        self.occupation: List[int]
        """
        stores the occupation of all core and valence orbitals, e.g. 2 for doubly, 1 for singly and 0 for unoccupied.
        """

        if not settings_dict:
            try:
                if atoms:
                    self.atoms = atoms
                elif xyz_file:
                    self.atoms = self.get_atoms(xyz_file)
            except TypeError as exc:
                raise ValueError("No xyz-file or List of atom names provided") from exc
            # self.update()
            self.__setup_orbital_and_electrons(self.atoms, double_d_shell)
            self.__correct_charge(charge)
            # automatically set lowest spin multiplicity
            if spin_multiplicity == 1 and self.electrons % 2 != (spin_multiplicity - 1) % 2:
                spin_multiplicity = 2

            self.__correct_spin_multiplicity(spin_multiplicity)
            if n_basis_functions > 0:
                self.__setup_occupation(n_basis_functions)
            else:
                self.__setup_occupation()
        else:
            if "atoms" not in settings_dict:
                self.atoms = self.get_atoms(settings_dict["xyz_file"])
            for key in settings_dict:
                if hasattr(self, key):
                    setattr(self, key, settings_dict[key])
            self.update()

    def update(self):
        """Update the molecular information.

        This means the number of electrons is updated with respect to
        charge, 4d orbtials are included if necessary, the occupation is
        corrected with respect to spin and number of basis functions.
        """
        self.__setup_orbital_and_electrons(self.atoms, self.double_d_shell)
        self.__correct_charge(self.charge)
        self.__correct_spin_multiplicity(self.spin_multiplicity)
        if self.n_basis_functions > 0:
            self.__setup_occupation(self.n_basis_functions)
        else:
            self.__setup_occupation()

    def __setup_orbital_and_electrons(self, atoms: List[str], double_d_shell: bool):
        """Set up orbitals and electrons.

        This means the number of core and valence orbitals, as well as the number of electrons
        for the molecule is generated, based on all elements in the atoms list.

        Parameters
        ----------
        atoms : List[str]
            every string has to be a symbol for an atom
        double_d_shell : bool
            only important for 3rd row transition metals, if 4d orbitals are part of CAS.
        """
        self.double_d_shell = double_d_shell

        elements = Elements(self.double_d_shell)
        self.core_orbitals = 0
        self.valence_orbitals = 0
        self.electrons = 0

        for atom in atoms:
            self.core_orbitals += elements.get_core_orbitals(atom)
            self.valence_orbitals += elements.get_valence_orbitals(
                atom,
            )
            self.electrons += elements.get_electrons(atom)

    def __correct_charge(self, charge: int):
        """Update electrons with respect charge.

        The number of electrons is subtracted by the provided charge.

        Parameters
        ----------
        charge : int
            total charge of the system
        """
        self.charge = charge
        self.electrons -= self.charge

    def __correct_spin_multiplicity(self, spin_multiplicity: int):
        """Check spin multiplicity.

        Parameters
        ----------
        spin_multiplicity : int
            the spin multiplicity, e.g. (2s+1)

        Raises
        ------
        ValueError
            if the number of electrons prevents the spin multiplicity.
        """
        if self.electrons % 2 != (spin_multiplicity - 1) % 2:
            raise ValueError(
                "Spin multiplicity does not match with number of electrons. Did you set the charge correctly?"
            )
        self.spin_multiplicity = spin_multiplicity

    def __setup_occupation(self, n_basis_functions: int = -1):
        """Set up correct occupations.

        Fills the occupation list by determining if an orbital is doubly, singly or not
        occupied, based on provided information.

        Parameters
        ----------
        n_basis_functions : int, optional
            the number of basis functions

        Raises
        ------
        ValueError
            if number of basis functions is <= 0 or smaller than valence + core orbitals
        """
        if 0 < n_basis_functions < self.valence_orbitals + self.core_orbitals:
            raise ValueError(
                """Number of basis functions is smaller than minimum number of orbitals.
                Did you activate double d orbitals?"""
            )
        if n_basis_functions > 1:
            self.n_basis_functions = n_basis_functions
            n_orbitals = self.n_basis_functions
            print(f"more orbitals than minimal number of orbitals. Number of orbitals: {n_basis_functions}")
        else:
            n_orbitals = self.core_orbitals + self.valence_orbitals

        n_electrons = self.electrons - (self.spin_multiplicity - 1)
        singly_electrons = self.spin_multiplicity - 1
        self.occupation = []
        for _ in range(n_orbitals):
            if n_electrons > 1:
                self.occupation.append(2)
                n_electrons -= 2
            elif n_electrons == 1:
                self.occupation.append(1)
                n_electrons -= 1
            else:
                if singly_electrons > 0:
                    self.occupation.append(1)
                    singly_electrons -= 1
                else:
                    self.occupation.append(0)

    def get_atoms(self, xyz_file: str) -> List[str]:
        """Read atoms from XYZ file.

        A list is filled with all element string occuring in the xyz file.

        Parameters
        ----------
        xyz_file : str
            path to the XYZ file

        Returns
        -------
        atoms : List[str]
            List of all atoms in the XYZ file
        """
        atoms = []
        with open(xyz_file, "r", encoding="UTF-8") as xyz:
            n_atoms = int(xyz.readline().strip().split()[0])
            # skip comment line (2. line)
            xyz.readline()
            for _ in range(n_atoms):
                line = xyz.readline().strip().split()
                atoms.append(str(line[0]))
        return atoms


if __name__ == "__main__":
    molecule = Molecule()
    # molecule.get_atoms()
