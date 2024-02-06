"""Module to provide the Autocas class, which controlls all active space
related methods."""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import Any, Dict, List, Optional, Tuple, Union, cast

import numpy as np

from scine_autocas.autocas_utils.active_space import ActiveSpace
from scine_autocas.autocas_utils.diagnostics import Diagnostics
from scine_autocas.autocas_utils.large_active_spaces import LargeSpaces
from scine_autocas.autocas_utils.molecule import Molecule

from ._version import __version__  # noqa: F401


class SingleReferenceException(Exception):
    """
    Raised when autoCAS determines that the system is single reference
    instead of multireference. This allows the client to catch the exception
    in case they want to perform a single reference calculation afterward.
    """
    pass


class Autocas:
    """
    Main class to handle all autoCAS functionality.

    Attributes
    ----------
    plateau_values : int, default = 10
        determines the minimum length of a plateau
    threshold_step : float, default 0.01
        determines the threshold steps in the plateau search and the number of threshold, respectively
    cas : ActiveSpace
        stores information on the active space
    molecule : Molecule
        stores information about the current molecule
    diagnostics : Diagnostics
        stores information from the DMRG calculation
    large_spaces : LargeSpaces
        controls settings for the large active space protocol
    _excited_states_orbital_indices : List[int]
        stores excited states CAS indices
    _excited_states_mod : int
        handles excited states in combination with large cas
    """

    __slots__ = (
        "plateau_values",
        "threshold_step",
        "cas",
        "molecule",
        "diagnostics",
        "large_spaces",
        "_excited_states_orbital_indices",
        "_excited_states_mod",
    )

    def __init__(self, molecule: Molecule, settings_dict: Optional[Dict[str, Any]] = None):
        """Initialize an autocas object.

        The autocas class provides all functionalities to set up active spaces from molecular information,
        or to apply an active space search based on entopy information from a dmrg calculation.

        Parameters
        ----------
        molecule : Molecule
            stores molecular information
        settings_dict : Dict[str, Any], optional
            holds settings from an autocas yaml input file
        """
        self.plateau_values: int = 10
        """determines the minimum length of a plateau"""
        self.threshold_step: float = 0.01
        """determines the threshold steps in the plateau search and the number of threshold, respectively"""
        self.cas: ActiveSpace = ActiveSpace()
        """stores information on the active space"""
        self.molecule: Molecule = molecule
        """stores information about the current molecule"""
        self.diagnostics: Diagnostics = Diagnostics()
        """stores information from the DMRG calculation"""
        self.large_spaces: LargeSpaces = LargeSpaces()
        """controls settings for the large active space protocol"""
        self._excited_states_orbital_indices: List[int] = []
        """stores excited states CAS indices"""
        self._excited_states_mod: int = -1
        """handles excited states in combination with large CAS"""

        if settings_dict is not None:
            for key in settings_dict:
                if key == "diagnostics":
                    self.diagnostics = Diagnostics(settings_dict[key])
                elif key == "large_spaces":
                    self.large_spaces = LargeSpaces(settings_dict[key])
                elif key == "cas":
                    self.cas = ActiveSpace(settings_dict[key])
                elif hasattr(self, key):
                    setattr(self, key, settings_dict[key])

    def make_initial_active_space(self) -> Tuple[List[int], List[int]]:
        """Build valence active space.

        Returns
        -------
        cas_occupations : List[int]
            a list with int occupations, e.g. 2: doubly occupied, 1: singly occupied, 0: virtual.
        cas_indices : List[int]
            contains the indices of the orbitals in the active space
        """
        self.cas.make_valence_cas(self.molecule)
        return self.cas.occupation, self.cas.orbital_indices

    def _sort_orbitals_by_s1(
        self, thresholds_list: np.ndarray, orbitals_index: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Sort orbital indices with respect to their threshold value.

        Parameters
        ----------
        thresholds_list : List[float]
            contains the thresholds of s1 values e.g. s1_i / max(s1)
        orbitals_index : List[int]
            a list with orbitals indices

        Returns
        -------
        thresholds_list : List[float]
            the ordered tresholds_list
        orbitals_index : List[int]
            the ordered orbitals_index with respect to the thresholds
        """
        # sort arrays decreasing
        sortkey = np.argsort(-thresholds_list)
        thresholds_list = thresholds_list[sortkey]
        orbitals_index = orbitals_index[sortkey]
        return thresholds_list, orbitals_index

    def _get_plateau(self, thresholds_list: np.ndarray) -> List[int]:
        """Create a plateau vector from s1 values.

        A plateau is complete after not changing the number of orbitals over 10%, e.g. self.plateau_values.

        Parameters
        ----------
        thresholds_list : np.ndarray
            contains fractions of s1 values, e.g. s1_i / max(s1)

        Returns
        -------
        plateau_vector : List[int]
            contains the indices from the plateau
        """
        percent = np.arange(0, 1.0, self.threshold_step)
        cas_orbs_over_threshold = np.zeros(int(1 / self.threshold_step + 1))
        for i, _ in enumerate(percent):
            cas_orbs_over_threshold[i] = sum(thresholds_list > percent[i])
        number_cas_orbs = cas_orbs_over_threshold[0]
        plateau_vector = []
        thresh_count = 1
        for i, threshold in enumerate(cas_orbs_over_threshold):
            if threshold == number_cas_orbs:
                thresh_count = thresh_count + 1
            else:
                thresh_count = 1
                number_cas_orbs = threshold
            if thresh_count >= self.plateau_values:
                plateau_vector.append(int(threshold))
        return plateau_vector

    def _find_plateaus(self) -> Tuple[List[int], np.ndarray]:
        """Search for a plateau of s1 values.

        Returns
        -------
        plateau_vector : List[int]
            the indices from the plateau
        orbitals_index : np.ndarray
                indices of the orbitals in the active space
        """
        max_s1 = max(self.diagnostics.s1_entropy)
        print(f"Maximal single orbital entropy: {max_s1}")
        if self._excited_states_orbital_indices:
            orbitals_index = np.array(self._excited_states_orbital_indices)
        else:
            orbitals_index = np.array(np.arange(1, self.cas.n_orbitals + 1))

        thresholds_list = np.zeros((self.cas.n_orbitals))
        # print(self.cas.n_orbitals)

        for i in range(self.cas.n_orbitals):
            # print(i)
            thresholds_list[i] = self.diagnostics.s1_entropy[i] / max_s1
        thresholds_list, orbitals_index = self._sort_orbitals_by_s1(
            thresholds_list, orbitals_index
        )
        plateau_vector = self._get_plateau(thresholds_list)
        return plateau_vector, orbitals_index

    def _make_active_space(self) -> bool:
        """Make active space by searching plteaus in s1 and reordering
        orbitals.

        Returns
        -------
        bool
            True if a plateau was found, e.g. a smaller active space; False if no plateau was found
        """
        plateau_vector, orbitals_index = self._find_plateaus()
        # found plateau

        if len(plateau_vector) != 0:
            self.cas.get_from_plateau(
                plateau_vector, list(orbitals_index.tolist())
            )
            print(f"Found plateau, including {self.cas.n_orbitals} orbitals.")
            return True
        return False

    def get_active_space(
        self, occupation: List[int], s1_entropy: np.ndarray, s2_entropy: Optional[np.ndarray] = None,
            mutual_information: Optional[np.ndarray] = None, force_cas: Optional[bool] = False,
            xyz_file: Optional[str] = None, molecule: Optional[Molecule] = None
    ) -> Tuple[List[int], List[int]]:
        """Get active space.

        autoCAS uses single orbital entropies to decrease the active space size by finding plateaus and
        including only orbitals with large single orbital entropies

        Parameters
        ----------
        occupation : List[int]
            orbital occupation for each orbital in active space
        s1_entropy : np.ndarray
            single orbital entropy
        s2_entropy : np.ndarray, optional
            two orbital entropy, not required right now
        mutual_information : np.ndarray, optional
            mutual information, not required right now
        force_cas :: bool
            forces to find an active space, even if no orbital has a high single orbital entropy

        Returns
        -------
        cas_occupations : List[int]
            the occupation of each orbital in the active space,
            e.g. 2: doubly occupied, 1: singly occupied, 0: virtual orbital
        cas_indices : List[int]
            indices of orbitals for the CAS calculation

        Raises
        ------
        AttributeError
            if no xyz file was provided
        FileNotFoundError
            if provided xyz file cannot be found
        Exception
            if the entropies indicate, that no cas is required, and not forced
        NotImplementedError
            if a larger cas than valence is theoretically required
        """

        print("")
        print("*******************************************************************************************")
        print("*                                                                                         *")
        print("*                                 Active space search                                     *")
        print("*                                                                                         *")
        print("*******************************************************************************************")

        # check if molecule is set
        if self.molecule.core_orbitals != -1 and molecule is None:
            pass
        else:
            try:
                if molecule is not None:
                    self.molecule = molecule
                elif xyz_file is not None:
                    self.molecule = Molecule(xyz_file)
                else:
                    raise AttributeError("No xyz-file provided.")
            except FileNotFoundError as exc:
                raise FileNotFoundError("No xyz-file provided.") from exc

        # check if excited states and large cas is used
        if self._excited_states_mod < -19:
            self._excited_states_mod = self.cas.n_orbitals

        # check if excited states are used
        if self._excited_states_orbital_indices is None:
            self.cas.update(self.molecule, occupation)
        else:
            self.cas.update(self.molecule, occupation, self._excited_states_orbital_indices)

        # for excited states and large cas
        if self._excited_states_mod > 1:
            self.cas.n_orbitals = self._excited_states_mod

        # set diagnostics
        self.diagnostics.s1_entropy = s1_entropy
        if s2_entropy is not None:
            self.diagnostics.s2_entropy = s2_entropy
        else:
            self.diagnostics.s2_entropy = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
        if mutual_information is not None:
            self.diagnostics.mutual_information = mutual_information
        else:
            self.diagnostics.mutual_information = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))

        self.cas.n_orbitals = len(s1_entropy)
        output_string = f"Start active space search for provided set of {self.cas.n_orbitals} orbitals"
        output_string += f" and {int(sum(occupation))} electrons."
        print(output_string)

        # check for single reference
        if self.diagnostics.is_single_reference() and not force_cas:
            print("Single orbital entropies (s1) are too low. AutoCAS suggests no active space.")
            print("Note: if you want an active space calculation anyway, you can enforce this")
            print("behavior, by adding 'force_cas=True' as agrument.")
            raise SingleReferenceException("No active space method required here")

        # make active space
        if self._make_active_space():
            print("Found smaller active space.")
            return self._pretty_return()

        print("No plateau found: autoCAS tries to exclude orbitals with low single orbital entropy.")

        self.cas.orbital_indices = []
        self.diagnostics.s1_entropy = np.array(self.cas.exclude_orbitals(self.diagnostics))

        if len(self.diagnostics.s1_entropy) < len(s1_entropy):
            occupation = self.cas.occupation
            self._excited_states_orbital_indices = self.cas.orbital_indices
            print("Successfully excluded orbitals with low single orbital entropy.")
            return self.get_active_space(
                occupation,
                self.diagnostics.s1_entropy,
                s2_entropy=self.diagnostics.s2_entropy,
                mutual_information=self.diagnostics.mutual_information,
                force_cas=force_cas,
                xyz_file=xyz_file,
            )
        self._excited_states_orbital_indices = []
        if self.cas.excluded and force_cas:
            print("Still no plateau found after excluding orbitals with low entropy.")
            return self._pretty_return()
        if force_cas:
            print("Enforcing CAS, even though there was no plateau.")
            return self._pretty_return()
        print("Not able to reduce active space size.")
        raise NotImplementedError("Repeat analysis with large active space.")
        # TODO: repeat analysis with larger cas

    def _pretty_return(self):
        """Ensure that cas occupation is always list of int."""
        self.cas.occupation = [int(i) for i in self.cas.occupation]
        print(f"Selected active space CAS(e, o): ({sum(self.cas.occupation)}, {len(self.cas.occupation)})")
        print(f"Orbital indices:     {self.cas.orbital_indices}")
        print(f"Orbital occupations: {self.cas.occupation}")
        return self.cas.occupation, self.cas.orbital_indices

    def get_large_active_spaces(self):
        """Create a list of active spaces for the large active space protocol
        in autoCAS.

        Returns
        -------
        occupation : List[List[int]]
            List of occupations for large CAS protocol
        orbital_indices : List[List[int]]
            List of orbitals indices for large CAS protocol
        """
        if self.large_spaces.max_orbitals > len(self.cas.orbital_indices):
            print(
                f"""Large CAS protocol is not required here, but it will be done anyways with
                \nmax orbitals = number of orbitals/2 = {len(self.cas.orbital_indices)/2}"""
            )
            self.large_spaces.max_orbitals = len(self.cas.orbital_indices) / 2
        occupied_orbitals, virtual_orbitals = self.large_spaces.separate_space(
            self.cas, self.molecule
        )
        self.large_spaces.generate_spaces(
            occupied_orbitals, virtual_orbitals, self.molecule
        )
        return self.large_spaces.occupation, self.large_spaces.orbital_indices

    def collect_entropies(
        self,
        indices_list: List[List[int]],
        occupations_list: List[List[int]],
        s1_list: List[np.ndarray],
        s2_list: Optional[List[np.ndarray]] = None,
        mut_inf_list: Optional[List[np.ndarray]] = None,
    ) -> Union[
        Tuple[np.ndarray, np.ndarray],
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ]:
        """Fill entropy matrices from lists of entropies produced by the large
        active space protocol.

        Parameters
        ----------
        indices_list : List[List[int]]
            contains several orbital indices from each active space calculation of the large CAS protocol
        occupations_list : List[List[int]]
            contains several orbital occupations for each active space calculation of the large CAS protocol
        s1_list : List[np.ndarray]
            contains several s1 entropies from the large CAS protocol
        s2_list : List[np.ndarray]
            contains several s2 entropies from the large CAS protocol
        mut_inf_list : List[np.ndarray]
            contains several mutual information from the large CAS protocol

        Returns
        -------
        occupation : List[int]
            the occupation for the large active space
        s1_entropy : np.ndarray
            s1 entropy build and scaled from several entropies
        s2_entropy : np.ndarray, optional
            s2 entropy build and scaled from several entropies
        mut_inf : np.ndarray, optional
            mutual information build and scaled from several mutual information
        """
        if self._excited_states_mod > 1:
            self.cas.n_orbitals = self._excited_states_mod
        occupation = np.zeros(self.cas.n_orbitals)
        s1_entropy = np.zeros(self.cas.n_orbitals)
        scale_1 = np.zeros(self.cas.n_orbitals)

        for i, _ in enumerate(s1_list):
            for count, index in enumerate(indices_list[i]):
                occupation[index - self.molecule.core_orbitals] = occupations_list[i][count]
                if self.large_spaces.average_entanglement:
                    s1_entropy[index - self.molecule.core_orbitals] += s1_list[i][count]
                    scale_1[index - self.molecule.core_orbitals] += 1
                else:
                    if s1_entropy[index - self.molecule.core_orbitals] < s1_list[i][count]:
                        s1_entropy[index - self.molecule.core_orbitals] = s1_list[i][count]
        if self.large_spaces.average_entanglement:
            # normalization
            for i, _ in enumerate(s1_entropy):
                s1_entropy[i] = s1_entropy[i] / scale_1[i]

        occupation = occupation.astype(int, copy=False)
        self.diagnostics.s1_entropy = s1_entropy
        self.diagnostics.s2_entropy = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
        self.diagnostics.mutual_information = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
        self.cas.occupation = list(occupation.tolist())

        if s2_list and mut_inf_list:
            s2_entropy = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
            mut_inf = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
            scale_2 = np.zeros((self.cas.n_orbitals, self.cas.n_orbitals))
            for i in range(len(s1_list)):
                for count1, index_1 in enumerate(indices_list[i]):
                    for count2, index_2 in enumerate(indices_list[i]):
                        if self.large_spaces.average_entanglement:
                            s2_entropy[
                                index_1 - self.molecule.core_orbitals,
                                index_2 - self.molecule.core_orbitals,
                            ] += s2_list[i][count1, count2]
                            mut_inf[
                                index_1 - self.molecule.core_orbitals,
                                index_2 - self.molecule.core_orbitals,
                            ] += mut_inf_list[i][count1, count2]
                            scale_2[
                                index_1 - self.molecule.core_orbitals,
                                index_2 - self.molecule.core_orbitals,
                            ] += 1
                        else:
                            if s2_entropy[
                                index_1 - self.molecule.core_orbitals,
                                index_2 - self.molecule.core_orbitals,
                            ] < s2_list[i][count1, count2]:
                                s2_entropy[
                                    index_1 - self.molecule.core_orbitals,
                                    index_2 - self.molecule.core_orbitals,
                                ] = s2_list[i][count1, count2]
                            if mut_inf[
                                index_1 - self.molecule.core_orbitals,
                                index_2 - self.molecule.core_orbitals,
                            ] < mut_inf_list[i][count1, count2]:
                                mut_inf[
                                    index_1 - self.molecule.core_orbitals,
                                    index_2 - self.molecule.core_orbitals,
                                ] = mut_inf_list[i][count1, count2]
            if self.large_spaces.average_entanglement:
                # normalization
                for i in range(len(s1_entropy)):
                    for j in range(len(s1_entropy)):
                        if scale_2[i, j] != 0:
                            s2_entropy[i, j] = s2_entropy[i, j] / scale_2[i, j]
                            mut_inf[i, j] = mut_inf[i, j] / scale_2[i, j]
            self.diagnostics.s2_entropy = s2_entropy
            self.diagnostics.mutual_information = mut_inf
            return occupation, s1_entropy, s2_entropy, mut_inf

        return occupation, s1_entropy

    def get_cas_from_large_cas(
        self,
        indices_list: List[List[int]],
        occupations_list: List[List[int]],
        s1_list: List[np.ndarray],
        s2_list: Optional[List[np.ndarray]] = None,
        mut_inf_list: Optional[List[np.ndarray]] = None,
        force_cas: Optional[bool] = False,
    ) -> Tuple[List[int], List[int]]:
        """Shortcut function to directly get an active space from the large CAS
        protocol.

        Parameters
        ----------
        indices_list : List[List[int]]
            contains several orbital indices from each active space calculation of the large CAS protocol
        occupations_list : List[List[int]]
            contains several orbital occupations for each active space calculation of the large CAS protocol
        s1_list : List[np.ndarray]
            contains several s1 entropies from the large CAS protocol
        s2_list : List[np.ndarray], optional
            contains several s2 entropies from the large CAS protocol
        mut_inf_list : List[np.ndarray], optional
            contains several mutual information from the large CAS protocol
        force_cas : bool, optional, default = False
            forces to find an active space, even if no orbital has a high single orbital entropy

        Returns
        -------
        cas_occupations : List[int]
            the occupation of each orbital in the active space,
            e.g. 2: doubly occupied, 1: singly occupied, 0: virtual orbital
        cas_indices : List[int]
            indices of orbitals for the CAS calculation

        Raises
        ------
        ValueError
            if something strange happened
        """
        if s2_list and mut_inf_list:
            entropies = self.collect_entropies(
                indices_list, occupations_list, s1_list, s2_list, mut_inf_list
            )
            if len(entropies) == 4:
                entropies = cast(Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], entropies)
                occupation = list(entropies[0].tolist())
                s1_entropy = entropies[1]
                s2_entropy = entropies[2]
                mut_inf = entropies[3]
                return self.get_active_space(occupation, s1_entropy, s2_entropy, mut_inf, force_cas=force_cas)

        entropies = self.collect_entropies(
            indices_list, occupations_list, s1_list, s2_list, mut_inf_list
        )
        if len(entropies) == 2:
            occupation = list(entropies[0].tolist())
            s1_entropy = entropies[1]
            return self.get_active_space(occupation, s1_entropy, force_cas=force_cas)

        raise ValueError("Something went wrong in get_cas_from_large_cas.")

    def get_cas_from_excited_states(
        self, occupation: List[int], s1_list: List[np.ndarray], s2_list: Optional[List[np.ndarray]] = None,
        mutual_information_list: Optional[List[np.ndarray]] = None, force_cas: Optional[bool] = False,
    ) -> Tuple[List[int], List[int]]:
        """Shortcut function to directly get an active space from excited
        states.

        Parameters
        ----------
        occupation : List[int]
            contains the orbital occupation for the active space.
        s1_list : List[np.ndarray]
            contains the one orbital entropies for each provided occupation
        s2_list : List[np.ndarray], optional
            contains the two orbital entropies for each provided occupation
        mutual_information_list : List[np.ndarray], optional
            contains the mutual information for each provided occupation

        Returns
        -------
        final_occupation : List[int]
            contains the occupations of the found active space
        final_orbital_indices : List[int]
            contains the orbital indices of the found active space
        """
        final_occupation_list = []
        final_orbital_list = []
        for i, s1_entropy in enumerate(s1_list):
            if mutual_information_list:
                if not s2_list:
                    s2_list = mutual_information_list
                final_occupation, final_orbital_indices = self.get_active_space(
                    occupation,
                    s1_entropy,
                    s2_entropy=s2_list[i],
                    mutual_information=mutual_information_list[i],
                    force_cas=force_cas,
                )
            else:
                final_occupation, final_orbital_indices = self.get_active_space(
                    occupation,
                    s1_entropy,
                    force_cas=force_cas,
                )

            final_occupation_list.append(final_occupation)
            final_orbital_list.append(final_orbital_indices)
        final_occupation = []
        final_orbital_indices = []
        for i, final_orbitals in enumerate(final_orbital_list):
            for j, final_orbital in enumerate(final_orbitals):
                if final_orbital not in final_orbital_indices:
                    final_orbital_indices.append(final_orbital)
                    final_occupation.append(final_occupation_list[i][j])
        return final_occupation, final_orbital_indices

    def get_cas_from_large_cas_excited_states(
        self, indices: List[List[int]], occupations: List[List[int]], s1_entropy: List[List[np.ndarray]],
        s2_entropy: Optional[List[List[np.ndarray]]] = None, mut_inf: Optional[List[List[np.ndarray]]] = None,
        force_cas: bool = False
    ):
        """Get active space from a large CAS excited state calculation.

        Parameters
        ----------
        indices :
            contains a list of indices for each state
        """

        self._excited_states_mod = -20

        # re-sorting entropies by states
        partial_s1_per_state: List[List[np.ndarray]] = []
        for i in range(len(s1_entropy[0])):
            partial_s1_per_state.append([])

        for i, part_s1 in enumerate(s1_entropy):
            for j, part_s1_1 in enumerate(part_s1):
                partial_s1_per_state[j].append(part_s1_1)

        if s2_entropy is not None and mut_inf is not None:

            partial_s2_per_state: List[List[np.ndarray]] = []
            partial_mut_inf_per_state: List[List[np.ndarray]] = []

            for i in range(len(s1_entropy[0])):
                partial_s2_per_state.append([])
                partial_mut_inf_per_state.append([])

            for i, part_s1 in enumerate(s1_entropy):
                for j, part_s1_1 in enumerate(part_s1):
                    partial_s2_per_state[j].append(s2_entropy[i][j])
                    partial_mut_inf_per_state[j].append(mut_inf[i][j])

        # get best active space per state
        final_occupation_list = []
        final_orbital_list = []

        final_orbital_indices = []
        final_occupation = []
        for i, partial_s1_state in enumerate(partial_s1_per_state):
            print(f"Evaluating state: {i}")
            if mut_inf is not None and s2_entropy is not None:
                occupation_per_state, orbital_indices_per_state = self.get_cas_from_large_cas(
                    indices, occupations, partial_s1_state, partial_s2_per_state[i],
                    partial_mut_inf_per_state[i], force_cas=force_cas,
                )
            else:
                occupation_per_state, orbital_indices_per_state = self.get_cas_from_large_cas(
                    indices, occupations, partial_s1_state, force_cas=force_cas,
                )

            print_n_electrons = int(sum(occupation_per_state))
            print_n_orbitals = len(occupation_per_state)
            print(
                f"Optimal CAS for current state CAS(e, o): CAS({print_n_electrons}, {print_n_orbitals})"
            )
            print(f"Optimal Indices:     {orbital_indices_per_state}")
            print(f"Optimal Occupations: {occupation_per_state}")

            final_occupation_list.append(occupation_per_state)
            final_orbital_list.append(orbital_indices_per_state)

            # merge active spaces and occupations
            for j, orbital_index in enumerate(orbital_indices_per_state):
                if orbital_index not in final_orbital_indices:
                    final_orbital_indices.append(orbital_index)
                    final_occupation.append(int(occupation_per_state[j]))

            # for i, final_orbitals in enumerate(final_orbital_list):
            #     for j, final_orbital in enumerate(final_orbitals):
            #         if final_orbital not in final_orbital_indices:
            #             final_orbital_indices.append(final_orbital)
            #             final_occupation.append(final_occupation_list[i][j])

        # sort indices and occupations
        final_orbital_indices_tmp = np.array(final_orbital_indices)
        final_occupation_tmp = np.array(final_occupation)
        sort_key = final_orbital_indices_tmp.argsort()
        final_orbital_indices_tmp = final_orbital_indices_tmp[sort_key]
        final_occupation_tmp = final_occupation_tmp[sort_key]
        final_occupation = list(final_occupation_tmp.tolist())
        final_orbital_indices = list(final_orbital_indices_tmp.tolist())

        # output
        n_states = len(final_occupation_list)
        print_n_electrons = int(sum(occupation_per_state))
        print_n_orbitals = len(occupation_per_state)
        print(f"Used entropies from {n_states} for evaluation of the active space.")
        print(f"Best active space involving all states: CAS({print_n_electrons}, {print_n_orbitals})")
        print(f"Optimal Indices:     {final_orbital_indices}")
        print(f"Optimal Occupations: {final_occupation}")

        return final_occupation, final_orbital_indices
