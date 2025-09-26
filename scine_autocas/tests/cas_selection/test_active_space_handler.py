# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np

from scine_autocas.cas_selection.active_space_handler import ActiveSpaceHandler
from scine_autocas.cas_selection.diagnostics import Diagnostics
from scine_autocas.utils.molecule import Molecule


class TestActiveSpaceHandler(unittest.TestCase):
    def setUp(self):
        path = os.path.dirname(os.path.abspath(__file__))
        xyz_file = f"{path}/files/n2.xyz"
        self.molecule_1 = Molecule(xyz_file)
        self.molecule_2 = Molecule(xyz_file)
        self.molecule_2.charge = -1
        self.molecule_2.spin_multiplicity = 4
        self.molecule_2.ecp_electrons = 2
        self.molecule_2.update()
        self.diagnostics = Diagnostics()
        self.diagnostics.weak_correlation_threshold = 0.1
        self.custom_s1 = np.array([0.0, 0.1, 0.1, 0.6, 0.7, 0.9, 0.8, 0.7, 0.4, 1.1, 0.0, 0.9, 0.8, 0.7])
        self.custom_occ = [2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        self.custom_inidices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17]
        self.expected_occ = [2, 2, 1, 1, 1, 0, 0, 0, 0, 0]
        self.expected_inidices = [3, 4, 5, 6, 7, 8, 9, 12, 13, 17]
        self.expected_s1 = np.array([0.6, 0.7, 0.9, 0.8, 0.7, 0.4, 1.1, 0.9, 0.8, 0.7])

    def test_valence_cas_closed_shell(self):
        cas_handler = ActiveSpaceHandler(self.molecule_1)
        valence_cas = cas_handler.get_valence_cas()
        occupation = valence_cas.get_occupation()
        indices = valence_cas.get_indices()
        self.assertEqual(occupation, [2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(indices, [2, 3, 4, 5, 6, 7, 8, 9])

    def test_valence_cas_open_shell(self):
        cas_handler = ActiveSpaceHandler(self.molecule_2)
        valence_cas = cas_handler.get_valence_cas()
        occupation = valence_cas.get_occupation()
        indices = valence_cas.get_indices()
        self.assertEqual(occupation, [2, 2, 2, 2, 1, 1, 1, 0])
        self.assertEqual(indices, [1, 2, 3, 4, 5, 6, 7, 8])

    def test_custom_valence_cas(self):
        custom_occ = [2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        custom_inidices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17]
        cas_handler = ActiveSpaceHandler(self.molecule_2)
        cas_handler.custom_valence_cas(custom_occ, custom_inidices)
        valence_cas = cas_handler.get_valence_cas()
        occupation = valence_cas.get_occupation()
        indices = valence_cas.get_indices()
        self.assertEqual(occupation, custom_occ)
        self.assertEqual(indices, custom_inidices)

    def test_exclude_orbitals(self):
        cas_handler = ActiveSpaceHandler(self.molecule_2)
        cas_handler.custom_valence_cas(self.custom_occ, self.custom_inidices)
        cas_handler.store_valence_s1_entropies(self.custom_s1)

        cas_handler.exclude_orbitals(self.diagnostics)

        valence_cas = cas_handler.get_valence_cas()
        occupation = valence_cas.get_occupation()
        indices = valence_cas.get_indices()
        s1 = valence_cas.get_s1_entropies()
        self.assertEqual(occupation, self.custom_occ)
        self.assertEqual(indices, self.custom_inidices)
        self.assertTrue(np.allclose(s1, self.custom_s1))

        cas = cas_handler.get_cas()
        occupation = cas.get_occupation()
        indices = cas.get_indices()
        s1 = cas.get_s1_entropies()
        self.assertEqual(occupation, self.expected_occ)
        self.assertEqual(indices, self.expected_inidices)
        self.assertTrue(np.allclose(s1, self.expected_s1))

    def test_getter(self):
        cas_handler = ActiveSpaceHandler(self.molecule_2)
        cas_handler.custom_valence_cas(self.custom_occ, self.custom_inidices)
        cas_handler.store_valence_s1_entropies(self.custom_s1)

        cas = cas_handler.get_cas()
        nelec = cas_handler.get_n_electrons()
        indices = cas_handler.get_indices()
        norbs = cas_handler.get_n_orbitals()
        occ = cas_handler.get_occupation()
        valence_cas = cas_handler.get_valence_cas()
        self.assertEqual(cas.get_occupation(), self.custom_occ)
        self.assertEqual(cas.get_indices(), self.custom_inidices)
        self.assertEqual(cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(valence_cas.get_occupation(), self.custom_occ)
        self.assertEqual(valence_cas.get_indices(), self.custom_inidices)
        self.assertEqual(valence_cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(valence_cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(valence_cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(norbs, len(self.custom_inidices))
        self.assertEqual(nelec, sum(self.custom_occ))
        self.assertEqual(indices, self.custom_inidices)
        self.assertEqual(occ, self.custom_occ)

        cas_handler.exclude_orbitals(self.diagnostics)
        cas = cas_handler.get_cas()
        nelec = cas_handler.get_n_electrons()
        indices = cas_handler.get_indices()
        norbs = cas_handler.get_n_orbitals()
        occ = cas_handler.get_occupation()
        valence_cas = cas_handler.get_valence_cas()
        excluded = cas_handler.excluded()
        successfully_excluded = cas_handler.successfully_excluded()

        self.assertEqual(cas.get_occupation(), self.expected_occ)
        self.assertEqual(cas.get_indices(), self.expected_inidices)
        self.assertEqual(cas.get_n_orbitals(), len(self.expected_inidices))
        self.assertEqual(cas.get_n_electrons(), sum(self.expected_occ))
        self.assertTrue(np.allclose(cas.get_s1_entropies(), self.expected_s1))
        self.assertEqual(valence_cas.get_occupation(), self.custom_occ)
        self.assertEqual(valence_cas.get_indices(), self.custom_inidices)
        self.assertEqual(valence_cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(valence_cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(valence_cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(norbs, len(self.expected_occ))
        self.assertEqual(nelec, sum(self.expected_occ))
        self.assertEqual(indices, self.expected_inidices)
        self.assertEqual(occ, self.expected_occ)
        self.assertEqual(excluded, True)
        self.assertEqual(successfully_excluded, True)

    def test_update_cas(self):
        cas_handler = ActiveSpaceHandler(self.molecule_2)
        cas_handler.custom_valence_cas(self.custom_occ, self.custom_inidices)
        cas_handler.store_valence_s1_entropies(self.custom_s1)

        cas = cas_handler.get_cas()
        nelec = cas_handler.get_n_electrons()
        indices = cas_handler.get_indices()
        norbs = cas_handler.get_n_orbitals()
        occ = cas_handler.get_occupation()
        valence_cas = cas_handler.get_valence_cas()
        self.assertEqual(cas.get_occupation(), self.custom_occ)
        self.assertEqual(cas.get_indices(), self.custom_inidices)
        self.assertEqual(cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(valence_cas.get_occupation(), self.custom_occ)
        self.assertEqual(valence_cas.get_indices(), self.custom_inidices)
        self.assertEqual(valence_cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(valence_cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(valence_cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(norbs, len(self.custom_inidices))
        self.assertEqual(nelec, sum(self.custom_occ))
        self.assertEqual(indices, self.custom_inidices)
        self.assertEqual(occ, self.custom_occ)

        cas_handler.update_cas(self.expected_occ, self.expected_inidices, self.expected_s1)
        cas = cas_handler.get_cas()
        nelec = cas_handler.get_n_electrons()
        indices = cas_handler.get_indices()
        norbs = cas_handler.get_n_orbitals()
        occ = cas_handler.get_occupation()
        valence_cas = cas_handler.get_valence_cas()
        self.assertEqual(cas.get_occupation(), self.expected_occ)
        self.assertEqual(cas.get_indices(), self.expected_inidices)
        self.assertEqual(cas.get_n_orbitals(), len(self.expected_inidices))
        self.assertEqual(cas.get_n_electrons(), sum(self.expected_occ))
        self.assertTrue(np.allclose(cas.get_s1_entropies(), self.expected_s1))
        self.assertEqual(valence_cas.get_occupation(), self.custom_occ)
        self.assertEqual(valence_cas.get_indices(), self.custom_inidices)
        self.assertEqual(valence_cas.get_n_orbitals(), len(self.custom_inidices))
        self.assertEqual(valence_cas.get_n_electrons(), sum(self.custom_occ))
        self.assertTrue(np.allclose(valence_cas.get_s1_entropies(), self.custom_s1))
        self.assertEqual(norbs, len(self.expected_occ))
        self.assertEqual(nelec, sum(self.expected_occ))
        self.assertEqual(indices, self.expected_inidices)
        self.assertEqual(occ, self.expected_occ)

    def test_get_from_plateau(self):
        plateau_vector = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        orbital_indices = [2, 6, 7, 5, 4, 3, 0, 1]
        cas_handler = ActiveSpaceHandler(self.molecule_1)
        cas_handler.set_from_plateau(plateau_vector, orbital_indices)
        self.assertEqual(cas_handler._current_cas.get_occupation(), [2, 2, 2, 0, 0, 0])
        self.assertEqual(cas_handler._current_cas.get_indices(), [4, 5, 6, 7, 8, 9])


if __name__ == "__main__":
    unittest.main()
