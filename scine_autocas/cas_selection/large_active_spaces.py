# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Optional, Tuple

import numpy as np

from scine_autocas.utils.defaults import Defaults
from scine_autocas.utils.molecule import Molecule

from .active_space import ActiveSpace


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
        "orbital_indices",
        "occupations",
        "average_entanglement",
    )

    def __init__(self) -> None:
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
        self.seed: Optional[int] = Defaults.AutoCAS.large_cas_seed
        """sets the seed for np.random. Should not be modified."""
        self.max_orbitals: int = Defaults.AutoCAS.large_cas_max_orbitals
        """maximum number of orbitals per active space. Active spaces contain the same number of orbitals."""
        self.orbital_indices: List[List[int]]
        """a List that contains a List of orbital indices"""
        self.occupations: List[List[int]]
        """a List that contains a List of orbital occupations"""
        self.average_entanglement = Defaults.AutoCAS.large_cas_average_entanglement
        """enable average entanglement for large cas instead of max entanglement"""

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

    def separate_space(self, cas: ActiveSpace, molecule: Molecule) -> Tuple[List[int], List[int], List[int]]:
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
        singly_occupied_spaces = []
        virtual_orbitals = []
        for i in cas.get_indices():
            # singly occupied orbitals count as occupied
            if molecule.occupation[i] == 0:
                virtual_orbitals.append(i)
            elif molecule.occupation[i] == 1:
                singly_occupied_spaces.append(i)
            else:
                occupied_orbitals.append(i)
        if len(singly_occupied_spaces) > 0:
            print("""Found singly occupied orbitals""")
        return occupied_orbitals, singly_occupied_spaces, virtual_orbitals

    def generate_spaces(
            self,
            occupied_orbitals: List[int],
            singly_occupied_spaces: List[int],
            virtual_orbitals: List[int],
            molecule: Molecule
    ):
        """Generate sub active space from a set of occupied and virtual orbital indices.

        The provided indices are further devided into smaller lists, which are later recombined to
        generate many, small active spaces, from a set of occupied and virtual orbital indices.

        Parameters
        ----------
        occupied_orbitals : List[int]
            stores all indices, which represent doubly occupied orbitals
        singly_orbitals : List[int]
            stores all indices, which represent singly occupied orbitals
        virtual_orbitals : List[int]
            stores all indices, which represent virtual orbitals
        molecule : Molecule
            handles all molecular information
        """
        if len(singly_occupied_spaces):
            self.max_orbitals -= len(singly_occupied_spaces)
            if self.max_orbitals <= 3:
                print("""
    Max number of orbitals for large cas protocol is smaller than minimum number of singly
    occupied orbitals. Please increase the size of the sub-CAS (large_cas_max_orbitals)
    to at least
        m singly occupied orbs + n occupied orbs + n virtual orbs
    where n should be > 0 and m the number of unpaired electrons / singly occupied orbitals
                                """)
                print(f"""
    Increasing max number of orbitals per sub-CAS to
        {len(singly_occupied_spaces)} singly occ. orbs + 2 occ. orbs + 2 virt. orbs
                                """)
                self.max_orbitals = 4

            print(f"""Found singly occupied orbitals.
                         Seeting max orbitals from {self.max_orbitals + len(singly_occupied_spaces)}
                         to {self.max_orbitals} and including singly occupied orbitals in every sub-CAS.
                         The final size of the sub-CASs will be {self.max_orbitals + len(singly_occupied_spaces)}.
                         """)
        # from here on max orbitals is the number of fully occupied and fully virtual orbtals
        # max orbs should be an even number
        # However, singly occupied orbitals account to
        modified_max_orbs = False
        if self.max_orbitals % 2 != 0:
            modified_max_orbs = True
            # self.max_orbitals -= 1

        partial_occupied_orbitals = self._partition_space(occupied_orbitals)

        # give orbital back to virtual
        if modified_max_orbs:
            self.max_orbitals += 1

        partial_virtual_orbitals = self._partition_space(virtual_orbitals)
        self.orbital_indices = []
        self.occupations = []
        for i in partial_occupied_orbitals:
            for j in partial_virtual_orbitals:
                sub_cas = i + j + singly_occupied_spaces
                sub_space = []
                sub_occupation = []
                for orb in sub_cas:
                    sub_space.append(orb)
                    sub_occupation.append(molecule.occupation[orb])
                sort_key = np.array(sub_space).argsort()
                np_sub_space = np.array(sub_space)[sort_key]
                np_sub_occupation = np.array(sub_occupation)[sort_key]
                self.orbital_indices.append(list(np_sub_space.tolist()))
                self.occupations.append(list(np_sub_occupation.tolist()))

    def get_spaces(self, cas: ActiveSpace, molecule: Molecule) -> Tuple[List[List[int]], List[List[int]]]:
        """Convenience function to generate occupations and indices for large CAS.

        Calls separate_space and generate_spaces

        Parameters
        ----------
        cas: ActiveSpace
            the active space to generate sub spaces from
        molecule: Molecule
            the molecule

        Returns
        -------
        self.occupations: List[List[int]]
            List which holds lists of occupations
        self.orbital_indices: List[List[int]]
            List which holds lists of indices
        """
        if self.max_orbitals > len(cas.get_indices()):
            warning_string = "Large CAS protocol is not required here, because:\n"
            warning_string += f"Number of orbitals in full cas: {len(cas.get_indices())}\n"
            warning_string += f"Number of orbitals per sub cas: {self.max_orbitals}\n"
            warning_string += "However, as it is requested will be done anyways with "
            warning_string += f"max orbitals = number of orbitals/2 = {len(cas.get_indices())/2}"""
            print(warning_string)
            self.max_orbitals = int(len(cas.get_indices()) / 2)

        occupied_orbitals, singly_occupied_orbitals, virtual_orbitals = self.separate_space(cas, molecule)
        self.generate_spaces(occupied_orbitals, singly_occupied_orbitals, virtual_orbitals, molecule)
        # if 1 in cas.get_occupation():
        #     self._replace_occupations_with_open_shell()
        return self.occupations, self.orbital_indices
