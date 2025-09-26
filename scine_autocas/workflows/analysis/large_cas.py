# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
import re

from scine_autocas.io import FileHandler

from .analyze import Analyze


class AnalyzeLargeCAS(Analyze):
    """Analyze large cas dir"""

    def numerical_sort(self, string: str):
        """Sort list of strings with numerical values acoording to their value"""
        return list(map(int, re.findall(r'\d+', string)))[0]

    def run(self):
        """run"""
        dmrg_dirs = FileHandler.get_all_dmrg_dirs()
        s1_entropies = []
        s2_entropies = []
        mut_infs = []
        # occ, index
        occupation, indices = self.autocas.make_initial_active_space()

        dmrg_dirs.sort(key=self.numerical_sort)

        for dmrg_dir in dmrg_dirs:
            os.chdir(dmrg_dir)
            print(f"Analyze: {os.getcwd()}")
            s1, s2, mut_inf = self.analyze_dmrg_dir()
            s1_entropies.append(s1)
            n_orbitals = len(s1)
            s2_entropies.append(s2)
            mut_infs.append(mut_inf)
            os.chdir(FileHandler.get_project_path())

        self.autocas.large_spaces.max_orbitals = n_orbitals
        tmp_occupations, tmp_indices = self.autocas.get_large_active_spaces()

        # occ, s1, s2, mut_inf
        _, s1, s2, mut_inf = self.autocas.collect_entropies(
            tmp_indices, tmp_occupations, s1_entropies, s2_entropies, mut_infs)
        self.results["initial_s1"] = s1
        self.results["initial_mutual_information"] = mut_inf
        self.results["initial_orbital_indices"] = indices

        self.plotting()

        print(self.autocas.get_active_space(s1, s2, mut_inf, occupation, indices))
