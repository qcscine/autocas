# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Optional

import numpy as np

from scine_autocas.utils.molecule import Molecule

from .active_space import ActiveSpace
from .diagnostics import Diagnostics


class ActiveSpaceHandler:
    """A Class to store active space information.

    Attributes
    ----------
    _valence_cas : ActiveSpace
        the initial active space
    _current_cas : ActiveSpace
        the current active space
    _excluded : bool
        if orbitals were excluded from the final active space
    _molecule : Molecule
        the molecule to evaluate
    """

    __slots__ = (
        "_valence_cas",
        "_current_cas",
        "_excluded",
        "_molecule",
        "_successfully_excluded",
    )

    def __init__(self, molecule: Molecule):  # , settings_dict: Optional[Dict[str, Any]] = None):
        """Initialize the active space handler from a molecule.

        Parameters
        ----------
        molecule : Molecule
            The molecule to evaluate
        """

        self._valence_cas: ActiveSpace
        """The valence CAS"""
        self._current_cas: ActiveSpace
        """The current CAS"""
        self._excluded: bool = False
        """True if orbitals were excluded from the selected CAS."""
        self._molecule: Molecule = molecule
        """Molecule to evaluate"""
        self._successfully_excluded: bool = False
        """True if excluding was successfull"""

        self._make_valence_cas()

    def _set_valence_cas(self, cas: ActiveSpace):
        """Set the valence CAS.

        If current cas does not exists, set current cas to valence cas.

        Parameters
        ----------
        cas : ActiveSpace
            the active space to set as valence CAS.
        """
        self._valence_cas = cas
        try:
            if not self._current_cas:
                self._current_cas = cas
        except AttributeError:
            self._current_cas = cas

    def get_valence_cas(self) -> ActiveSpace:
        """Getter for the valence CAS.

        Returns
        -------
        self._valence_cas : ActiveSpace
            the valence active space

        Raises
        ------
        ValueError
            If no valence cas exists.
        """
        if self._valence_cas is None:
            print("Valence CAS not set.")
            raise ValueError("No valence cas exists.\nMake valence cas, before calling this method.")
        return self._valence_cas

    def custom_valence_cas(self, occupation: List[int], indices: List[int]):
        """Create custom valence CAS.

        Parameters
        ----------
        occupation : List[int]
            list with occupation
        indices : List[int]
            list with indices (0-based)

        Notes
        -----
        This function resets the current CAS.

        Raises
        ------
        ValueError
            If unpaired electrons from the occupation list does not match the spin of the molecule.
        ValueError
            If the occupation has more electrons than the molecule.
        """
        unpaired_elecs = self._molecule.spin_multiplicity - 1
        occupation_unpaired_elec = 0
        for i in occupation:
            if i == 1:
                occupation_unpaired_elec += 1
        if unpaired_elecs != occupation_unpaired_elec:
            print("Occupation in custom CAS does not match spin multiplicity.")
            raise ValueError("New occupation has wrong spin multiplicity")
        if sum(occupation) > self._molecule.electrons:
            print("The custom CAS has more electrons than the molecule.")
            raise ValueError("More electrons than possible in new occupation")
        print(f"New cas {occupation} {indices}")

        del self._current_cas
        self._set_valence_cas(ActiveSpace(occupation, indices))

    def _make_valence_cas(self):
        """Make valence active space."""
        orbital_indices = self._make_valence_indices()
        occupation = self._make_valence_occupation()
        self._set_valence_cas(ActiveSpace(occupation, orbital_indices))

    def _make_valence_occupation(self) -> List[int]:
        """Determine the number of electrons and orbitals, as well as the
        occupation for the provided molecule.

        The electrons, orbtials and occupation here corresponds to the valence active space.
        The occupation is generated with respect to the corresponding spin multiplicity.

        Returns
        -------
        occupation : List[int]
            List with orbital occupations
        """
        n_electrons = self._molecule.electrons - 2 * self._molecule.core_orbitals
        tmp_spin_multiplicity = self._molecule.spin_multiplicity
        occupation = []

        n_orbitals = self._molecule.valence_orbitals

        for _ in range(n_orbitals):
            if n_electrons > 1 and tmp_spin_multiplicity == 1:
                occupation.append(2)
                n_electrons -= 2
            elif n_electrons == 1:
                occupation.append(1)
                n_electrons -= 1
            elif tmp_spin_multiplicity > 1 and n_electrons > 1:
                occupation.append(1)
                n_electrons -= 1
                tmp_spin_multiplicity -= 1
            else:
                occupation.append(0)

        occupation.sort(reverse=True)
        return occupation

    def _make_valence_indices(self) -> List[int]:
        """Generate a list of indices for the valence active space.

        Returns
        -------
        orbital_indices : List[int]
            List with all orbital indices of valence CAS (0-based).

        Notes
        ----_
        The orbital indices start with 0 for the first orbital.
        """
        orbital_indices = []
        for i in range(self._molecule.valence_orbitals):
            orbital_indices.append(self._molecule.core_orbitals + i)
        return orbital_indices

    def store_valence_s1_entropies(self, s1_entropies: np.ndarray):
        """Store the s1 entropies corresponding to the valence CAS.

        Parameters
        ----------
        s1_entropies : np.ndarray
            Numpy array with the s1 entropies for each orbital in valence CAS.
        """
        self._valence_cas.store_s1(s1_entropies)
        self._set_valence_cas(self._valence_cas)

        # To update current CAS
        try:
            self._current_cas.get_s1_entropies()
        except LookupError:
            try:
                self._current_cas.store_s1(s1_entropies)
            except ValueError:
                pass

    def store_s1_entropy(self, s1_entropies: np.ndarray):
        """Store the s1 entropies for the current CAS.

        Parameters
        ----------
        s1_entropies : np.ndarray
            Numpy array with the s1 entropies for each orbital in current active space
        """
        if self._current_cas is None:
            self._current_cas = self._valence_cas
        self._current_cas.store_s1(s1_entropies)

    def reset(self) -> None:
        """Return if exclusion was successful"""
        if self._valence_cas:
            self._current_cas = self._valence_cas

    def update_cas(
            self,
            occupation: List[int],
            indices: List[int],
            s1_entropies: Optional[np.ndarray] = None
    ):
        """Update the current active space.

        Parameters
        ----------
        occupation : List[int]
            occupation list for the new CAS.
        indices : List[int]
            orbital indices (0-based) for the new CAS.
        s1_entropies : np.ndarray, optional
            Numpy array with the s1 entropies for each orbital in the new active space.
        """
        self._current_cas = ActiveSpace(occupation, indices)
        if s1_entropies is not None:
            self._current_cas.store_s1(s1_entropies)

    def get_cas(self) -> ActiveSpace:
        """Getter for the current active space."""
        return self._current_cas

    def get_n_electrons(self) -> int:
        """Getter for number of electrons of current CAS."""
        return self._current_cas.get_n_electrons()

    def get_n_orbitals(self) -> int:
        """Getter for number of orbitals in current CAS."""
        return self._current_cas.get_n_orbitals()

    def get_occupation(self) -> List[int]:
        """Getter for the occupation list of current CAS"""
        return self._current_cas.get_occupation()

    def get_indices(self) -> List[int]:
        """Getter for the index list of current CAS"""
        return self._current_cas.get_indices()

    def get_s1(self) -> np.ndarray:
        """Getter for the entropies of current CAS"""
        return self._current_cas.get_s1_entropies()

    def excluded(self) -> bool:
        """Get state of exclusion"""
        return self._excluded

    def successfully_excluded(self) -> bool:
        """Return if exclusion was successful"""
        return self._successfully_excluded

    def exclude_orbitals(self, diagnostics: Diagnostics) -> List[int]:
        """Exclude weak correlated orbitals from active space.

        If not active space can be found, autocas tries to exclude orbitals with low single orbital entropy.
        The functions removes the corresponding orbital indices and occupations from the class attributes.

        Parameters
        ----------
        diagnostics : Diagnostics
            object, which holds diffferent threshold values and provides functions for diagnostics.

        Notes
        -----
        This function set the "excluded" attribute to True
        """
        if self._excluded:
            return []

        print(
            f"Before excluding Orbials: {self._current_cas.get_occupation()}, {self._current_cas.get_indices()}")

        s1_new = []
        occ_new = []
        indices_new = []
        excluded_orb_indices = []
        excluded_orb_occ = []
        excluded_s1 = []
        max_s1 = max(self._current_cas.get_s1_entropies()) * diagnostics.weak_correlation_threshold
        for i, s1_entropy in enumerate(self._current_cas.get_s1_entropies()):
            if s1_entropy >= max_s1:
                s1_new.append(s1_entropy)
                occ_new.append(self._current_cas.get_occupation()[i])
                indices_new.append(self._current_cas.get_indices()[i])
            else:
                excluded_orb_indices.append(self._current_cas.get_indices()[i])
                excluded_orb_occ.append(self._current_cas.get_occupation()[i])
                excluded_s1.append(s1_entropy)

        if len(s1_new) == len(self._current_cas.get_s1_entropies()):
            print("No orbital could be excluded")
        else:
            print("Excluded orbitals with:")
            print(f"indices:    {excluded_orb_indices}")
            print(f"occupation: {excluded_orb_indices}")
            print(f"entropies:  {excluded_orb_indices}")
            self._successfully_excluded = True
            print("Successfully excluded orbitals with low single orbital entropy")

        self._excluded = True
        self._current_cas = ActiveSpace(occ_new, indices_new)
        self._current_cas.store_s1(np.array(s1_new))
        print(f"New CAS: {self._current_cas.get_occupation()}, {self._current_cas.get_indices()}")
        return excluded_orb_indices

    def set_from_plateau(self, plateau_vector: List[int], plateau_indeces: List[int]):
        """Get number of cas electrons and cas indices from plateau vector and
        orbitals_index.

        The plateau vector contains orbital indices for the corresponding orbital_index list,
        which represents orbitals for an active space. This mapping is provided by the autoCAS class.

        Parameters
        ----------
        plateau_vector : List[int]
            contains indices from a found plateau
        plateau_indices : List[int]
            list of plateau indices
        """
        new_orbital_indices = []
        new_occupation = []
        for i in plateau_indeces[: int(plateau_vector[0])]:
            new_orbital_indices.append(self._current_cas.get_indices()[i])
            new_occupation.append(self._current_cas.get_occupation()[i])
        sort_key = np.array(new_orbital_indices).argsort()
        new_orbital_indices = list(np.array(new_orbital_indices)[sort_key].tolist())
        new_occupation = list(np.array(new_occupation)[sort_key].tolist())
        self._current_cas = ActiveSpace(new_occupation, new_orbital_indices)
