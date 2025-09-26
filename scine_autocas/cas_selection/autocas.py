# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from scine_autocas.utils.defaults import Defaults
from scine_autocas.utils.exceptions import SingleReferenceException
from scine_autocas.utils.molecule import Molecule

from .active_space_handler import ActiveSpaceHandler
from .diagnostics import Diagnostics
from .large_active_spaces import LargeSpaces


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

    def __init__(self, molecule: Molecule):
        """Initialize an autocas object.

        The autocas class provides all functionalities to set up active spaces from molecular information,
        or to apply an active space search based on entopy information from a dmrg calculation.

        Parameters
        ----------
        molecule : Molecule
            stores molecular information
        """
        self.plateau_values: int = Defaults.AutoCAS.plateau_values
        """determines the minimum length of a plateau"""
        self.threshold_step: float = Defaults.AutoCAS.threshold_step
        """determines the threshold steps in the plateau search and the number of threshold, respectively"""
        self.molecule: Molecule = molecule
        """stores information about the current molecule"""
        self.cas: ActiveSpaceHandler = ActiveSpaceHandler(self.molecule)
        """stores information on the active space"""
        self.diagnostics: Diagnostics = Diagnostics()
        """stores information from the DMRG calculation"""
        self.large_spaces: LargeSpaces = LargeSpaces()
        """controls settings for the large active space protocol"""
        self._excited_states_orbital_indices: List[int]
        """stores excited states CAS indices"""
        self._excited_states_mod: int = -1
        """handles excited states in combination with large CAS"""

    def setup_from_dict(self, settings: Dict[str, Any]):
        """Apply settings from a dict to the autocas object.

        Parameters
        ----------
        settings : Dict[str, Any]
            settings dict

        Notes
        -----
        The keys have to be the same name as the attributes.
        if a setting is an attribute from the interface.settings object,
        it is also set here.
        """
        autocas_settings = settings["AutoCAS"]
        for key in autocas_settings:
            large_cas_key = key
            if key.startswith("large_cas_"):
                large_cas_key = key[10:]
            if hasattr(self, key):
                setattr(self, key, autocas_settings[key])
            elif hasattr(self.diagnostics, key):
                setattr(self.diagnostics, key, autocas_settings[key])
            elif hasattr(self.large_spaces, large_cas_key):
                setattr(self.large_spaces, large_cas_key, autocas_settings[key])

    def make_initial_active_space(self) -> Tuple[List[int], List[int]]:
        """Build valence active space.

        Returns
        -------
        cas_occupations : List[int]
            a list with int occupations, e.g. 2: doubly occupied, 1: singly occupied, 0: virtual.
        cas_indices : List[int]
            contains the indices of the orbitals in the active space
        """
        print("Make valence CAS")
        self.cas.get_valence_cas()
        init_cas = self.cas.get_cas()
        return init_cas.get_occupation(), init_cas.get_indices()

    def _sort_orbitals_by_s1(
        self,
        thresholds_list: np.ndarray,
        plateau_indices: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Sort orbital indices with respect to their threshold value.

        Parameters
        ----------
        thresholds_list : List[float]
            contains the thresholds of s1 values e.g. s1_i / max(s1)
        plateau_indices : List[int]
            a list with plateau indices

        Returns
        -------
        thresholds_list : List[float]
            the ordered tresholds_list
        plateau_indices : List[int]
            the ordered orbitals_index with respect to the thresholds
        """
        # sort arrays decreasing
        sortkey = np.argsort(-thresholds_list)
        thresholds_list = thresholds_list[sortkey]
        plateau_indices = plateau_indices[sortkey]

        return thresholds_list, plateau_indices

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

        # print(f"plateau_vector: {set(plateau_vector)}")
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
        cas = self.cas.get_cas()
        max_s1 = max(cas.get_s1_entropies())
        print(f"Maximal single orbital entropy: {max_s1}")
        plateau_indices = np.array(np.arange(0, cas.get_n_orbitals()))
        thresholds_list = np.zeros((cas.get_n_orbitals()))
        for i in range(cas.get_n_orbitals()):
            thresholds_list[i] = cas.get_s1_entropies()[i] / max_s1
        thresholds_list, plateau_indices = self._sort_orbitals_by_s1(
            thresholds_list, plateau_indices
        )
        plateau_vector = self._get_plateau(thresholds_list)
        return plateau_vector, plateau_indices

    def _make_active_space(self) -> bool:
        """Make active space by searching plteaus in s1 and reordering
        orbitals.

        Returns
        -------
        bool
            True if a plateau was found, e.g. a smaller active space; False if no plateau was found
        """
        plateau_vector, plateau_indices = self._find_plateaus()

        # found plateau
        if len(plateau_vector) != 0:
            self.cas.set_from_plateau(plateau_vector, list(plateau_indices.tolist()))
            print(f"Found plateau, including {self.cas.get_n_orbitals()} orbitals.")
            return True

        print("No plateau found")
        return False

    def get_cas_suggestions(
            self,
            occupation: List[int],
            indices: List[int],
            s1_entropy: np.ndarray,
            s2_entropy: Optional[np.ndarray] = None,
            mutual_information: Optional[np.ndarray] = None
    ):
        """Return all plateaus found"""
        _ = s2_entropy
        _ = mutual_information
        self.cas.update_cas(occupation, indices, s1_entropy)
        output_string = f"Start active space search for provided set of {self.cas.get_n_orbitals()} orbitals"
        output_string += f" and {int(sum(occupation))} electrons."
        print(output_string)
        plateau_vector, plateau_index = self._find_plateaus()
        print(f"Plateau vector: {plateau_vector}, plateau: {plateau_index}")

        suggestions = []
        for plateau in list(set(plateau_vector))[::-1]:
            print(f"plateau: {plateau} orbitals: {plateau_index[:plateau]}")
            tmp_orbitals = list(plateau_index[:plateau])
            new_occupations = []
            new_orbitals = []
            for orb_index in tmp_orbitals:
                new_occupations.append(self.cas.get_occupation()[orb_index])
                new_orbitals.append(self.cas.get_indices()[orb_index])
            new_orbitals.sort()
            new_occupations.sort()
            suggestions.append((new_occupations[::-1], new_orbitals))

        return suggestions

    def get_active_space(
            self,
            s1_entropy: np.ndarray,
            s2_entropy: Optional[np.ndarray] = None,
            mutual_information: Optional[np.ndarray] = None,
            occupation: Optional[List[int]] = None,
            indices: Optional[List[int]] = None,
            force_cas: bool = False
    ) -> Tuple[List[int], List[int]]:
        """Get active space.

        autoCAS uses single orbital entropies to decrease the active space size by finding plateaus and
        including only orbitals with large single orbital entropies.
        To find a CAS for a custom initial active space, use the occupation and indices lists, to provide
        autocas with the required information.

        Parameters
        ----------
        s1_entropy : np.ndarray
            single orbital entropy
        s2_entropy : np.ndarray, optional
            two orbital entropy, not required right now
        mutual_information : np.ndarray, optional
            mutual information, not required right now
        occupation : List[int], optional
            list with occupation of each orbital
        indices : List[int], optional
            list with orbital indices
        force_cas : bool
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
        SingleReferenceException
            if the entropies indicate, that no cas is required, and not forced
        TODO: implement this
        NotImplementedError
            if a larger cas than valence is theoretically required
        """

        if occupation and indices:
            self.cas.update_cas(occupation, indices, s1_entropy)

        self._check_single_reference(s1_entropy, force_cas)
        # check if excited states are used
        self.cas.store_s1_entropy(s1_entropy)
        output_string = f"Start active space search for provided set of {self.cas.get_n_orbitals()} orbitals"
        output_string += f" and {self.cas.get_n_electrons()} electrons."
        print(output_string)

        # make active space
        if self._make_active_space():
            print("Found smaller active space.")
            return self._pretty_return()

        print("No plateau found: autoCAS tries to exclude orbitals with low single orbital entropy.")
        self.cas.exclude_orbitals(self.diagnostics)
        if self.cas.successfully_excluded():
            return self.get_active_space(
                self.cas.get_s1(),
                s2_entropy=s2_entropy,
                mutual_information=mutual_information,
                occupation=self.cas.get_occupation(),
                indices=self.cas.get_indices(),
                force_cas=force_cas,
            )
        if self.cas.excluded() and force_cas:
            print("No plateau found after excluding orbitals with low entropy.")
            return self._pretty_return()

        print("Not able to reduce active space size.")
        raise NotImplementedError("Repeat analysis with large active space or use force_cas=True.")
        # TODO: repeat analysis with larger cas

    def _check_single_reference(self, s1: np.ndarray, force_cas: bool):
        """Check if system is single reference.

        Parameters
        ----------
        s1 : np.ndarray
            single orbital entropies
        force_cas : bool
            if force_cas is set to true, skip single reference check

        Raises
        ------
        SingleReferenceException
            if diagnostics result in single reference and no cas is forced.
        """
        # check for single reference
        if self.diagnostics.is_single_reference(s1) and not force_cas:
            output_string = "Single orbital entropies (s1) are too low. AutoCAS suggests no active space.\n"
            output_string += "Note: if you want an active space calculation anyway, you can enforce this\n"
            output_string += "behavior, by adding 'force_cas=True' as agrument."
            print(output_string)
            raise SingleReferenceException("No active space method required here")

    def _pretty_return(self):
        """Ensure that cas occupation is always list of int."""
        # print(f"Selected active space CAS(e, o): ({self.cas.get_n_electrons()}, {self.cas.get_n_orbitals()})")
        # print(f"Orbital indices:     {self.cas.get_indices()}")
        # print(f"Orbital occupations: {self.cas.get_occupation()}")
        return self.cas.get_occupation(), self.cas.get_indices()

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
        return self.large_spaces.get_spaces(self.cas, self.molecule)

    def collect_entropies(
        self,
        indices_list: List[List[int]],
        occupations_list: List[List[int]],
        s1_list: List[np.ndarray],
        s2_list: Optional[List[np.ndarray]] = None,
        mut_inf_list: Optional[List[np.ndarray]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        occupation = np.zeros(self.cas.get_n_orbitals())
        s1_entropy = np.zeros(self.cas.get_n_orbitals())
        scale_1 = np.zeros(self.cas.get_n_orbitals())

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

        # pylint: disable=R1702
        if s2_list and mut_inf_list:
            s2_entropy = np.zeros((self.cas.get_n_orbitals(), self.cas.get_n_orbitals()))
            mut_inf = np.zeros((self.cas.get_n_orbitals(), self.cas.get_n_orbitals()))
            scale_2 = np.zeros((self.cas.get_n_orbitals(), self.cas.get_n_orbitals()))
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
            return occupation, s1_entropy, s2_entropy, mut_inf

        return occupation, s1_entropy, np.array([]), np.array([])

    def get_cas_from_large_cas(
        self,
        indices_list: List[List[int]],
        occupations_list: List[List[int]],
        s1_list: List[np.ndarray],
        s2_list: Optional[List[np.ndarray]] = None,
        mut_inf_list: Optional[List[np.ndarray]] = None,
        force_cas: bool = False,
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
                occupation = list(entropies[0].tolist())
                s1_entropy = entropies[1]
                s2_entropy = entropies[2]
                mut_inf = entropies[3]
                return self.get_active_space(s1_entropy, s2_entropy, mut_inf, force_cas=force_cas)

        entropies = self.collect_entropies(
            indices_list, occupations_list, s1_list, s2_list, mut_inf_list
        )
        if len(entropies) == 2:
            occupation = list(entropies[0].tolist())
            s1_entropy = entropies[1]
            return self.get_active_space(occupation, s1_entropy, force_cas=force_cas)

        raise ValueError("Something went wrong in get_cas_from_large_cas.")

    def get_cas_from_excited_states(
        self, s1_list: List[np.ndarray], s2_list: Optional[List[np.ndarray]] = None,
        mutual_information_list: Optional[List[np.ndarray]] = None, force_cas: bool = False,
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
                    s1_entropy,
                    s2_entropy=s2_list[i],
                    mutual_information=mutual_information_list[i],
                    force_cas=force_cas,
                )
                self.cas.reset()
            else:
                final_occupation, final_orbital_indices = self.get_active_space(
                    s1_entropy,
                    force_cas=force_cas,
                )
                self.cas.reset()

            final_occupation_list.append(final_occupation)
            final_orbital_list.append(final_orbital_indices)
            print("")
        final_occupation = []
        final_orbital_indices = []
        for i, final_orbitals in enumerate(final_orbital_list):
            for j, final_orbital in enumerate(final_orbitals):
                if final_orbital not in final_orbital_indices:
                    final_orbital_indices.append(final_orbital)
                    final_occupation.append(final_occupation_list[i][j])
        # sort indices and corresponding occupations
        final_occupation = [x for _, x in sorted(zip(final_orbital_indices, final_occupation))]
        final_orbital_indices.sort()

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
            print("")
            print(f"Evaluating state: {i}")
            if mut_inf is not None and s2_entropy is not None:
                occupation_per_state, orbital_indices_per_state = self.get_cas_from_large_cas(
                    indices, occupations, partial_s1_state, partial_s2_per_state[i],
                    partial_mut_inf_per_state[i], force_cas=force_cas,
                )
                self.cas.reset()
            else:
                occupation_per_state, orbital_indices_per_state = self.get_cas_from_large_cas(
                    indices, occupations, partial_s1_state, force_cas=force_cas,
                )
                self.cas.reset()

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
        print("")
        print(f"Used entropies from {n_states} for evaluation of the active space.")
        print(f"Best active space involving all states: CAS({print_n_electrons}, {print_n_orbitals})")
        print(f"Optimal Indices:     {final_orbital_indices}")
        print(f"Optimal Occupations: {final_occupation}")
        print("")
        print("")

        return final_occupation, final_orbital_indices
