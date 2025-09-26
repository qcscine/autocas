# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
from scine_autocas.io import FileHandler, logger

from .workflow import Workflow


class ClassicExcitedStatesWorkflow(Workflow):
    """The workflow"""

    def _initial_dmrg_impl(self):
        """Initial DMRG calculation"""
        logger.frame("Initial DMRG calculation")
        FileHandler.make_initial_dmrg_dir()
        _, initial_s1, initial_s2, init_mut_inf = self.interface.calculate(
            self.results["initial_occupation"], self.results["initial_orbital_indices"]
        )
        orig_entang_name = FileHandler.PlotNames.entanglement_file
        orig_thresh_name = FileHandler.PlotNames.threshold_file
        for i, s1 in enumerate(initial_s1):
            print(f"State: {i}")
            print(f"initial s1: {s1}")
            print(f"max s1: {max(s1)}")
            self.results["initial_s1"] = initial_s1[i]
            self.results["initial_s2"] = initial_s2[i]
            self.results["initial_mutual_information"] = init_mut_inf[i]
            # FileHandler.
            logger.empty_line()

            # create entanglement and threshold diagrams
            FileHandler.PlotNames.entanglement_file = orig_entang_name[:-4] + f"_{i}.pdf"
            FileHandler.PlotNames.threshold_file = orig_thresh_name[:-4] + f"_{i}.pdf"
            self.plotting()

        self.results["initial_s1"] = initial_s1
        self.results["initial_s2"] = initial_s2
        self.results["initial_mutual_information"] = init_mut_inf

        FileHandler.PlotNames.entanglement_file = orig_entang_name
        FileHandler.PlotNames.threshold_file = orig_thresh_name
        logger.empty_line()

        logger.frame("Active space search")
        final_occupation, final_orbital_indices = self.autocas.get_cas_from_excited_states(
            initial_s1,
            mutual_information_list=init_mut_inf,
            force_cas=False,
        )
        logger.final_cas(final_occupation, final_orbital_indices)
        logger.empty_line()
        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices
