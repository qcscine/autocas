# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import unittest

import numpy as np

from scine_autocas.cas_selection.active_space import ActiveSpace


class TestActiveSpace(unittest.TestCase):

    def test_init(self):
        occupation = [2, 2, 2, 2, 1, 1, 0, 0]
        indices = [2, 4, 9, 10, 11, 14, 16, 18]
        cas = ActiveSpace(occupation, indices)
        self.assertEqual(occupation, cas.get_occupation())
        self.assertEqual(indices, cas.get_indices())
        self.assertEqual(8, cas.get_n_orbitals())
        self.assertEqual(10, cas.get_n_electrons())

    def test_s1(self):
        occupation = [2, 2, 2, 2, 1, 1, 0, 0]
        indices = [2, 4, 9, 10, 11, 14, 16, 18]
        s1 = np.array([0.1, 0.2, 1.1, 0.9, 0.8, 0.6, 0.2, 0.3])
        cas = ActiveSpace(occupation, indices)
        cas.store_s1(s1)
        self.assertEqual(occupation, cas.get_occupation())
        self.assertEqual(indices, cas.get_indices())
        self.assertEqual(8, cas.get_n_orbitals())
        self.assertEqual(10, cas.get_n_electrons())
        self.assertTrue(np.allclose(s1, cas.get_s1_entropies()))


if __name__ == "__main__":
    unittest.main()
