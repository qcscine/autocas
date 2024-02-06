"""Combine active spaces and occupation lists.

These functions allow the combination of multiple active spaces along a reaction coordinate through orbital maps.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Tuple
import numpy as np


def transform_orbital_groups(orbital_groups: List[List[List[int]]]):
    """
    Transform the orbital mapping in terms of groups to a map of orbitals to groups.
    Parameters:
    -----------
        orbital_groups : List[List[List[int]]]
            The list of orbital groups.
    Returns:
    --------
        A matrix (orbitals x systems) containing the group indices for each orbital.
    """
    n_systems = len(orbital_groups[0])
    n_orbitals = 0
    for group in orbital_groups:
        n_orbitals += len(group[0])
    orbital_to_group = np.zeros((n_orbitals, n_systems), dtype=int)  # type: ignore
    for i_group, group in enumerate(orbital_groups):
        for i_sys, sys in enumerate(group):
            for i_orb in sys:
                orbital_to_group[i_orb, i_sys] = i_group
    return orbital_to_group


def combine_active_spaces(occupations: List[List[float]], active_spaces: List[List[int]],
                          orbital_groups: List[List[List[int]]]) -> Tuple[List[List[float]], List[List[int]]]:
    """
    Combine multiple active spaces (e.g., from points along a reaction coordinate) with an orbital map.

    Parameters:
    -----------
        occupations : List[List[float]]
            The occupations for each active space and each orbital.
            E.g., [[2, 2, 0, 0], [2, 2, 0, 0]] for two systems.
        active_spaces : List[List[int]]
            The indices of the orbitals considered as active (starting from 0) for each system.
        orbital_groups : List[List[List[int]]]
            Group-wise mapping of the orbitals. The orbitals are grouped into sets that are indistinguishable.
            For each group, and system the orbital indices are given, e.g.,
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
        The combined active spaces according to mapping.
    """
    # orbital_groups index 0: orbital group
    #                      1: system index
    #                      2: orbital index
    orbital_to_group = transform_orbital_groups(orbital_groups)
    # Collect the orbital groups that are active in some active space.
    active_orbital_groups: List[int] = []
    occupation_in_group: List[float] = []
    for i_sys, (occupation, active_space) in enumerate(zip(occupations, active_spaces)):
        for i_orb, occ in zip(active_space, occupation):
            i_group = orbital_to_group[i_orb, i_sys]
            if i_group not in active_orbital_groups:
                active_orbital_groups.append(i_group)
                occupation_in_group.append(occ)
            elif abs(occ - occupation_in_group[active_orbital_groups.index(i_group)]) > 1e-9:
                raise ValueError("Inconsistent occupation along the active space. The orbital mapping is incorrect.")
    # Build the combined active spaces as the union of all orbital groups considered active.
    n_systems = len(orbital_groups[0])
    new_active_spaces: List[List[int]] = [[] for _ in range(n_systems)]
    new_occupations: List[List[float]] = [[] for _ in range(n_systems)]
    for i_group, occ in zip(active_orbital_groups, occupation_in_group):
        group = orbital_groups[i_group]
        n_orbitals = len(group[0])
        for i_sys, orbitals in enumerate(group):
            new_active_spaces[i_sys] += orbitals
            new_occupations[i_sys] += [occ for _ in range(n_orbitals)]

    return new_occupations, new_active_spaces
