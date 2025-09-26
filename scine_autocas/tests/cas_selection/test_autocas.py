# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np

from scine_autocas.cas_selection.autocas import Autocas
from scine_autocas.utils.molecule import Molecule


class TestAutocas(unittest.TestCase):
    def setUp(self):
        path = os.path.dirname(os.path.abspath(__file__))
        self._xyz_name = f"{path}/n2_molecule_test.xyz"
        n2_xyz = open(self._xyz_name, "w")
        n2_xyz.write("2\n\n")
        n2_xyz.write("N 0 0 0\n")
        n2_xyz.write("N 0 0 2\n")
        n2_xyz.close()
        self.molecule = Molecule(self._xyz_name)

    def tearDown(self):
        try:
            os.remove(self._xyz_name)
        except FileNotFoundError:
            pass

    def test_init(self):
        autocas = Autocas(self.molecule)
        init_occ, init_ind = autocas.make_initial_active_space()
        self.assertEqual(init_occ, [2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(init_ind, [2, 3, 4, 5, 6, 7, 8, 9])

    def test_get_cas_1(self):
        autocas = Autocas(self.molecule)
        # works because autocas assumes valence CAS occupation and indices
        s1 = np.array([0.01, 0.02, 0.65, 0.7, 0.6, 0.61, 0.62, 0.69])
        occupation, indices = autocas.get_active_space(s1)
        self.assertEqual(occupation, [2, 2, 2, 0, 0, 0])
        self.assertEqual(indices, [4, 5, 6, 7, 8, 9])

    def test_get_cas_2(self):
        autocas = Autocas(self.molecule)
        # use one more orbital than valence
        init_occ = [2, 2, 2, 2, 2, 2, 0, 0, 0]
        init_ind = [1, 2, 3, 8, 15, 16, 27, 38, 39]
        s1 = np.array([0.03, 0.01, 0.02, 0.65, 0.7, 0.6, 0.61, 0.62, 0.69])
        occupation, indices = autocas.get_active_space(s1, occupation=init_occ, indices=init_ind)
        self.assertEqual(occupation, [2, 2, 2, 0, 0, 0])
        self.assertEqual(indices, [8, 15, 16, 27, 38, 39])

    def test_get_cas_suggestion(self):
        autocas = Autocas(self.molecule)
        init_occ = [2, 2, 2, 2, 1, 1, 0, 0, 0]
        init_ind = [1, 2, 3, 8, 15, 16, 27, 38, 39]
        s1 = np.array([0.03, 0.01, 0.02, 0.65, 0.9, 0.94, 0.61, 0.62, 0.69])
        suggestions = autocas.get_cas_suggestions(occupation=init_occ, indices=init_ind, s1_entropy=s1)
        self.assertEqual([([2, 1, 1, 0, 0, 0], [8, 15, 16, 27, 38, 39]), ([1, 1], [15, 16])], suggestions)


if __name__ == "__main__":
    unittest.main()
