# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
from scine_autocas.io import FileHandler, logger

from .workflow import Workflow


class LargeCasExcitedStateWorkflow(Workflow):
    """The workflow"""

    def _initial_dmrg_impl(self):
        """Large CAS initial DMRG calculations."""
        self.interface.settings.large_cas = True
        print("Check large cas workflow")
        print("Sub-CAS in large cas protocol")
        print(f"Averaging of entropies is {self.autocas.large_spaces.average_entanglement}")

        large_cas_spaces = self.autocas.get_large_active_spaces()
        large_cas_occupations = large_cas_spaces[0]
        large_cas_indices = large_cas_spaces[1]

        orig_entang_name = FileHandler.PlotNames.entanglement_file
        orig_thresh_name = FileHandler.PlotNames.threshold_file
        print(f"Number of sub-CASs: {len(large_cas_occupations)}")
        for i, partial_occupation in enumerate(large_cas_occupations):
            print(f" Cas space: {i+1}/{len(large_cas_occupations)}")
            print(f"  occupation {partial_occupation}")
            print(f"  indices    {large_cas_indices[i]}")
            print("")

        # store partial entropies for each state
        partial_s1_list = []
        partial_s2_list = []
        partial_mut_inf_list = []
        logger.frame("Initial DMRG")
        for i, partial_occupation in enumerate(large_cas_occupations):
            dir_name = f"dmrg_{i}"
            FileHandler.make_directory_from_project_root(dir_name)
            print(f" Calculating Cas space: {i+1}/{len(large_cas_occupations)}")
            print(f"  occupation {partial_occupation}")
            print(f"  indices    {large_cas_indices[i]}")
            _, s1_partial, s2_partial, mut_inf_partial = self.interface.calculate(
                partial_occupation, large_cas_indices[i]
            )
            print("")
            partial_s1_list.append(s1_partial)
            partial_s2_list.append(s2_partial)
            partial_mut_inf_list.append(mut_inf_partial)
        self.interface.initial_cas_prepared()

        # shuffle lists into list of states of subcas
        partial_s1_per_state = []
        partial_s2_per_state = []
        partial_mut_inf_per_state = []

        for i in range(len(partial_s1_list[0])):
            partial_s1_per_state.append([])
            partial_s2_per_state.append([])
            partial_mut_inf_per_state.append([])

        for i, part_s1 in enumerate(partial_s1_list):
            for j, part_s1_j in enumerate(part_s1):
                partial_s1_per_state[j].append(part_s1_j)
                partial_s2_per_state[j].append(partial_s2_list[i][j])
                partial_mut_inf_per_state[j].append(partial_mut_inf_list[i][j])

        s1_per_state = []
        s2_per_state = []
        mut_inf_per_state = []

        # occ, s1, s2, mut_inf
        for i, tmp_s1 in enumerate(partial_s1_per_state):
            (
                _, initial_s1, initial_s2, initial_mut_inf
            ) = self.autocas.collect_entropies(
                large_cas_indices,
                large_cas_occupations,
                tmp_s1,
                partial_s2_per_state[i],
                partial_mut_inf_per_state[i],
            )
            s1_per_state.append(initial_s1)
            s2_per_state.append(initial_s2)
            mut_inf_per_state.append(initial_mut_inf)

            print(f"initial s1: {initial_s1}")
            print(f"max s1: {max(initial_s1)}")
            self.results["initial_s1"] = initial_s1
            self.results["initial_s2"] = initial_s2
            self.results["initial_mutual_information"] = initial_mut_inf
            logger.empty_line()

            # create entanglement and threshold diagrams
            FileHandler.PlotNames.entanglement_file = orig_entang_name[:-4] + f"_{i}.pdf"
            FileHandler.PlotNames.threshold_file = orig_thresh_name[:-4] + f"_{i}.pdf"
            self.plotting()
            logger.empty_line()
        self.results["initial_s1"] = s1_per_state
        self.results["initial_s2"] = s2_per_state
        self.results["initial_mutual_information"] = mut_inf_per_state
        FileHandler.PlotNames.entanglement_file = orig_entang_name
        FileHandler.PlotNames.threshold_file = orig_thresh_name

        print("Final CAS")
        final_occupation, final_orbital_indices = self.autocas.get_cas_from_large_cas_excited_states(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_list,  # type: ignore
            partial_s2_list,  # type: ignore
            partial_mut_inf_list,  # type: ignore
        )
        logger.final_cas(final_occupation, final_orbital_indices)
        logger.empty_line()
        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices
