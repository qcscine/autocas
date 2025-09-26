# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List

import numpy as np


class ActiveSpace:
    """Internal representation of an active space.

    Attributes
    ----------
    _occupation : List[int]
        occupation list of the active space
    _indices : List[int]
        index list of the active space
    _s1_entropies : np.ndarray
        s1 entropy of the active space
    """
    __slots__ = ("_occupation", "_indices", "_s1_entropies")

    def __init__(self, occupation: List[int], indices: List[int]):
        """Initialize active space.

        Parameters
        ----------
        occupation : List[int]
            occupation list of the active space
        indices : List[int]
            index list of the active space
        """
        self._sanity_check(occupation, indices)

        self._occupation: List[int] = occupation
        """Occupation (2, 1, 0) of every orbital in active space"""
        self._indices: List[int] = indices
        """Orbital index (0-based) of every orbital in active space"""
        self._s1_entropies: np.ndarray
        """S1 entropy for every orbital in active space"""
        print("New active space:")
        print(f"\torbital occuation:         {self._occupation}")
        print(f"\torbital indices (0-based): {self._indices}")

    def _sanity_check(self, occupation: List[int], indices: List[int]):
        """Check every orbital has an index and an occupation and check occupation range.

        Parameters
        ----------
        occupation : List[int]
            occupation of every orbital
        indices : List[int]
            index for each orbital

        Raises
        ------
        ValueError
            In case length of list is difference, or occupation range is wrong.
        """
        if len(occupation) != len(indices):
            print(f"number of orbitals does not match: {occupation}, {indices}")
            raise ValueError("occupation and indices do not match in length")
        for i in occupation:
            if i > 3 or i < 0:
                print(f"Found occupation of: {i} in occupation: {occupation}")
                raise ValueError("Occupation can only be 0, 1 or 2")

    def get_occupation(self) -> List[int]:
        """Getter for the occupation.

        Returns
        -------
        self._occupation : List[int]
            list with the occupation of each orbital in CAS.
        """
        return self._occupation

    def get_indices(self) -> List[int]:
        """Getter for the indices (0-based).

        Returns
        -------
        self._indices: List[int]
            list with the index of each orbital in CAS.
        """
        return self._indices

    def get_n_electrons(self) -> int:
        """Getter for the number of electrons in CAS.

        Returns
        -------
        nelec : int
            The number of electrons.
        """
        return sum(self._occupation)

    def get_n_orbitals(self) -> int:
        """Getter for the number of orbitals in CAS.

        Returns
        -------
        norb : int
            The number of orbitals.
        """
        return len(self._occupation)

    def get_s1_entropies(self) -> np.ndarray:
        """Getter for the s1 entropy of each orbital in CAS.

        Returns
        -------
        self._s1_entropies : np.ndarray
            Numpy array with the s1 entropies.

        Raises
        ------
        LookupError
            In case s1 entropies are not set yet.
        """
        if self._s1_entropies is None:
            print(f"No entropies found in CAS: {self._occupation}, {self._indices}")
            raise LookupError("No s1 correspond to this active space.")
        return self._s1_entropies

    def store_s1(self, s1_entropies: np.ndarray):
        """Store s1 entropies for the active space.

        Parameters
        ----------
        s1_entropies : np.ndarray
            Numpy array with the s1 entropies

        Raises
        ------
        ValueError
            In case s1 entropies has more or less entries than orbitals in CAS.
        """
        print(f"s1 in active space: {s1_entropies}")
        if len(s1_entropies) != len(self._occupation):
            raise ValueError("number of entropy values does not match nummber of orbitals")
        self._s1_entropies = s1_entropies
