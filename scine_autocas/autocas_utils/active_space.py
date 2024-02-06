"""Determine the active space.

This module handles the active space search for a give set of
diagnostics and system specific values.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import Any, Dict, List, Optional

import numpy as np

from scine_autocas.autocas_utils.diagnostics import Diagnostics
from scine_autocas.autocas_utils.molecule import Molecule


class ActiveSpace:
    """A Class to store active space information and search for one.

    A active space is give for a give occupation and orbital indices.

    Attributes
    ----------
    n_orbitals : int
        the number of orbitals in the active space search
    n_electrons : int
        the number of electrons in the active space search
    n_inactive_orbitals : int
        the number of inactive orbitals, to d subtract them from input occupations
    orbitals_indices : List[int]
        the index of each orbital in the active space
    occupation : List[int]
        the occupation of each orbital, e.g. 2, 1, 0
    excluded : bool
        internal flag if orbitals were excluded from the active space
    Notes
    -----
    The orbital_indices attribute indices start with 1.
    """

    __slots__ = (
        "n_orbitals",
        "n_electrons",
        "n_inactive_orbitals",
        "orbital_indices",
        "occupation",
        "excluded"
    )

    def __init__(self, settings_dict: Optional[Dict[str, Any]] = None):
        """Construct an ActiveSpace object.

        At creation every list is assumend empty and will be filled by the provided methods.
        If a settings dir is provided, with attributes in it, these attributes will be overwritten.

        Parameters
        ----------
        settings_dict : Dict[str, Any], optional
            a dict, usually provided by the input_handler, which stores attributes and corresponding values

        See Also
        --------
        settings_dict : InputHandler
        """
        self.n_orbitals: int = 0
        """number of orbitals in active space"""
        self.n_electrons: int = 0
        """number of electrons in active space"""
        self.n_inactive_orbitals: int = 0
        """number of inactive occupied orbitals"""
        self.orbital_indices: List[int]
        """List with orbital indices for active space"""
        self.occupation: List[int]
        """List with orbital occupations for CAS, e.g. "2" doubly occupied, "1" singly occupied, "0" virtual" """
        self.excluded: bool = False
        """internal flag if orbital were excluded"""
        if settings_dict is not None:
            for key in settings_dict:
                if hasattr(self, key):
                    setattr(self, key, settings_dict[key])

    def make_valence_cas(self, molecule: Molecule):
        """Make valence active space.

        Creates the valence active space for the provided molecule object.

        Parameters
        ----------
        molecule : Molecule
            the molecule object to create an active space for
        """
        self._make_valence_indices(molecule)
        self._make_valence_occupation(molecule)

    def _make_valence_occupation(self, molecule: Molecule):
        """Determine the number of electrons and orbitals, as well as the
        occupation for the provided molecule.

        The electrons, orbtials and occupation here corresponds to the valence active space.
        The occupation is generated with respect to the corresponding spin multiplicity.

        Parameters
        ----------
        molecule : Molecule
            the molecule to determine the occupation for
        """
        self.n_electrons = molecule.electrons - 2 * molecule.core_orbitals
        tmp_cas_electrons = self.n_electrons
        # -1 to get 0 for singlett
        tmp_spin_multiplicity = molecule.spin_multiplicity
        tmp_occupation = []

        n_orbitals = molecule.valence_orbitals
        if molecule.n_basis_functions > 1:
            n_orbitals = molecule.n_basis_functions

        for _ in range(n_orbitals):
            if tmp_cas_electrons > 1 and tmp_spin_multiplicity == 1:
                tmp_occupation.append(2)
                tmp_cas_electrons -= 2
            elif tmp_cas_electrons == 1:
                tmp_occupation.append(1)
                tmp_cas_electrons -= 1
            elif tmp_spin_multiplicity > 1 and tmp_cas_electrons > 1:
                tmp_occupation.append(1)
                tmp_cas_electrons -= 1
                tmp_spin_multiplicity -= 1
            else:
                tmp_occupation.append(0)

        tmp_occupation.sort(reverse=True)
        self.occupation = tmp_occupation
        self.n_orbitals = len(self.occupation)

    def _make_valence_indices(self, molecule: Molecule):
        """Generate a list of indices for the valence active space.

        The indices are either created from the number of basis functions of the number of valence orbitals.
        This depends on the provided molecule object.

        Parameters
        ----------
        molecule : Molecule
            the molecule to determine the cas indices for

        Notes
        ----_
        The orbital indices start with 1 for the first orbital.
        """
        self.orbital_indices = []
        if molecule.n_basis_functions < 1:
            for i in range(molecule.valence_orbitals):
                self.orbital_indices.append(molecule.core_orbitals + i)
        else:
            for i in range(molecule.n_basis_functions - molecule.core_orbitals):
                self.orbital_indices.append(molecule.core_orbitals + i)

    def update(self, molecule: Molecule, occupation: List[int], orbital_indices: Optional[List[int]] = None):
        """Get the orbital occupation from aufbau principle.

        If orbital_indices are provided, the aufbau principle here, means that the occupation
        will still be [2,2,0] even if the orbital_indices are [2,4,3]

        Parameters
        ----------
        molecule : Molecule
            the molecule object, which holds molecular system specific variables
        occupation : List[int]
            the orbital occupation, for each orbital in the active space
        orbital_indices : Optional[List[int]]
            the orbital indices, for each orbital in the active space
        """
        self.n_inactive_orbitals = int(molecule.electrons / 2 - sum(occupation) / 2)
        self.n_electrons = sum(occupation)
        if orbital_indices:
            self.n_orbitals = len(orbital_indices)
            self.orbital_indices = orbital_indices
        else:
            self.n_orbitals = len(occupation)
            self.occupation = occupation

    def exclude_orbitals(self, diagnostics: Diagnostics) -> List[float]:
        """Exclude weak correlated orbitals from active space.

        If not active space can be found, autocas tries to exclude orbitals with low single orbital entropy.
        The functions removes the corresponding orbital indices and occupations from the class attributes.

        Parameters
        ----------
        diagnostics : Diagnostics
            object, which holds diffferent threshold values and provides functions for diagnostics.

        Returns
        -------
        s1_new : List[float]
            single orbital entropies, where low s1 values are removed

        Notes
        -----
        This function set the "excluded" attribute to True
        """
        s1_new = []
        if not self.orbital_indices:

            self.orbital_indices = list(np.array(np.arange(1, self.n_orbitals + 1)))
        max_s1 = max(diagnostics.s1_entropy) * diagnostics.weak_correlation_threshold
        orbital_index = len(diagnostics.s1_entropy)
        excluded_orbitals = []
        for s1_entropy in reversed(diagnostics.s1_entropy):
            orbital_index -= 1
            if s1_entropy > max_s1:
                s1_new.append(s1_entropy)
            else:
                # self.occupation.pop(i)
                self.orbital_indices.pop(orbital_index)
                excluded_orbitals.append(orbital_index)
        print(f"Successfully excluded the orbitals: {excluded_orbitals}")
        self.n_orbitals = len(self.orbital_indices)
        self.excluded = True
        return s1_new[::-1]

    def get_from_plateau(self, plateau_vector: List[int], orbitals_index: List[int]):
        """Get number of cas electrons and cas indices from plateau vector and
        orbitals_index.

        The plateau vector contains orbital indices for the corresponding orbital_index list,
        which represents orbitals for an active space. This mapping is provided by the autoCAS class.

        Parameters
        ----------
        plateau_vector : List[int]
            contains indices from a found plateau
        orbitals_index : List[int]
            sorted list of orbital indices
        """
        self.orbital_indices = []
        self.n_electrons = 0
        occupation = []
        for i in orbitals_index[: int(plateau_vector[0])]:
            self.orbital_indices.append(self.n_inactive_orbitals + i - 1)
            occupation.append(self.occupation[i - 1])
        sort_key = np.array(self.orbital_indices).argsort()
        self.orbital_indices = list(np.array(self.orbital_indices)[sort_key].tolist())
        self.occupation = list(np.array(occupation)[sort_key].tolist())
        self.n_orbitals = len(self.orbital_indices)
        self.n_electrons = sum(self.occupation)
