"""A handler for large active space protocol.

This module implements the LargeSpaces class, which handles all
variables and functionalities required for the large active space
protocol.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from scine_autocas.autocas_utils.active_space import ActiveSpace
from scine_autocas.autocas_utils.molecule import Molecule


class LargeSpaces:
    """A Class to store information and handle functionalities for the large
    active space protocol.

    The large active space protocol enables active space searches in initial actives space with more
    than 200 orbitals. It divides the orbital space into an occupied and virtual sub space.
    These spaces are separated into many smaller subspaces. Afterwards all occupied subspaces are recombined
    with all virtual subspaces, hence creating more, but smaller, active space. All active spaces are
    evaluated by an "initial" DMRG calculation, to calculate the single orbital entropies. These entropies
    are then recombined, to approximate the single orbital entropies from the full initial cas.
    The "final" active space is then constructed from the approximated s1.

    """

    __slots__ = (
        "seed",
        "max_orbitals",
        "n_orbitals",
        "n_electrons",
        "orbital_indices",
        "occupation",
        "average_entanglement",
    )

    def __init__(self, settings_dict: Optional[Dict[str, Any]] = None):
        """Construct the LargeSpaces object.

        A LargeSpaces object stores all relevant data and provides routines to
        divide spaces into occupied and virtual as well as into sub spaces and
        the recombination of these subspaces.

        Parameters
        ----------
        settings_dict : Dict[str, Any], optional
            a dict, usually provided by the input_handler, which stores attributes and corresponding values

        See Also
        --------
        settings_dict : InputHandler
        """
        self.seed: Optional[int] = 42
        """sets the seed for np.random. Should not be modified."""
        self.max_orbitals: int = 30
        """maximum number of orbitals per active space. Active spaces contain the same number of orbitals."""
        self.n_orbitals: List[int]
        """contains number of orbitals per sub-CAS"""
        self.n_electrons: List[int]
        """contains number of electrons per sub-CAS"""
        self.orbital_indices: List[List[int]]
        """a List that contains a List of orbital indices"""
        self.occupation: List[List[int]]
        """a List that contains a List of orbital occupations"""
        self.average_entanglement: bool = False
        """Flag to average the entropies from all sub-CASs instead of taking the max value"""
        if settings_dict is not None:
            for key in settings_dict:
                if hasattr(self, key):
                    setattr(self, key, settings_dict[key])

    def _partition_space(self, orbital_indices: List[int]) -> List[List[int]]:
        """Create small active spaces within the valence space.

        New active space are created by recombining all subspaces from occupied and virtual space.

        Parameters
        ----------
        orbital_indices : List[int]
            contains all indices which correspond to an orbital space, e.g. occupied or virtual

        Returns
        -------
        partial_orbital_indices : List[List[int]]
            contains Lists which contains orbital indices from an orbital space, e.g. occupied or virtual
        """
        partial_orbital_indices = []
        i = 0
        while i < len(orbital_indices):
            subvector = []
            for _ in range(int(self.max_orbitals / 2)):
                if i < len(orbital_indices):
                    subvector.append(orbital_indices[i])
                    i += 1
                # fill last sub vector with random orbitals
                else:
                    if self.seed is not None:
                        np.random.seed(self.seed)
                    random_i = np.random.randint(len(orbital_indices))
                    tmp_var = 0
                    while tmp_var < 1:
                        if orbital_indices[random_i] not in subvector:
                            subvector.append(orbital_indices[random_i])
                            break
                        random_i = np.random.randint(len(orbital_indices))
            partial_orbital_indices.append(subvector)
        return partial_orbital_indices

    def separate_space(self, cas: ActiveSpace, molecule: Molecule) -> Tuple[List[int], List[int]]:
        """Separate the valence active space into occupied and virtual
        orbitals.

        The information for the separation is provided by an ActiveSpace and Molecule object.

        Parameters
        ----------
        cas : ActiveSpace
            an ActiveSpace object, which stores all active space related information
        molecule : Molecule
            a Molecule object, which stores all molecular system related information

        Returns
        -------
        occupied_orbitals : List[int]
            contains all orbital indices which correspond to occupied orbtials
        virtual_orbitals : List[int]
            contains all orbital indices which correspond to virtual orbtials
        """
        occupied_orbitals = []
        virtual_orbitals = []
        for i in cas.orbital_indices:
            if molecule.occupation[i] == 0:
                virtual_orbitals.append(i)
            else:
                occupied_orbitals.append(i)
        return occupied_orbitals, virtual_orbitals

    def generate_spaces(self, occupied_orbitals: List[int], virtual_orbitals: List[int], molecule: Molecule):
        """Generate sub active space from a set of occupied and virtual orbital
        indices.

        The provided indices are further devided into smaller lists, which are later recombined to
        generate many, small active spaces, from a set of occupied and virtual orbital indices.

        Parameters
        ----------
        occupied_orbitals : List[int]
            stores all indices, which represent occupied orbitals
        virtual_orbitals: List[int]
            stores all indices, which represent virtual orbitals
        molecule : Molecule
            handles all molecular information
        """

        partial_occupied_orbitals = self._partition_space(occupied_orbitals)
        partial_virtual_orbitals = self._partition_space(virtual_orbitals)
        self.orbital_indices = []
        self.occupation = []
        self.n_orbitals = []
        self.n_electrons = []
        for i in partial_occupied_orbitals:
            for j in partial_virtual_orbitals:
                sub_cas = i + j
                sub_space = []
                sub_occupation = []
                for orb in sub_cas:
                    sub_space.append(orb)
                    sub_occupation.append(molecule.occupation[orb])
                sort_key = np.array(sub_space).argsort()
                np_sub_space = np.array(sub_space)[sort_key]
                np_sub_occupation = np.array(sub_occupation)[sort_key]
                self.orbital_indices.append(list(np_sub_space.tolist()))
                self.n_electrons.append(sum(sub_occupation))
                self.n_orbitals.append(len(np_sub_space))
                self.occupation.append(list(np_sub_occupation.tolist()))
