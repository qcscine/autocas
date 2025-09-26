# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
# from scine_autocas.cas_selection import Autocas
# from scine_autocas.interfaces.interface import Interface
from scine_autocas.io import FileHandler

from .analyze import Analyze


class AnalyzeConventional(Analyze):
    """Anlyze"""
    # def __init__(self, autocas: Autocas = None, interface: Interface = None):
    #     super().__init__(autocas, interface)

    def run(self):
        """run"""
        FileHandler.ch_to_initial_dmrg_dir()
        # occ, index
        init_occ, indices = self.autocas.make_initial_active_space()
        s1, s2, mut_inf = self.analyze_dmrg_dir()

        self.results["initial_s1"] = s1
        self.results["initial_mutual_information"] = mut_inf
        self.results["initial_orbital_indices"] = indices
        self.plotting()

        print(self.autocas.get_active_space(s1, s2, mut_inf, init_occ, indices))

        print(self.autocas.get_cas_suggestions(init_occ, indices, s1, s2, mut_inf))
