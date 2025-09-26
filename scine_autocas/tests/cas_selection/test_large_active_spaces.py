# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import unittest

from scine_autocas.cas_selection.active_space import ActiveSpace
from scine_autocas.cas_selection.large_active_spaces import LargeSpaces
from scine_autocas.utils.molecule import Molecule


class TestLargeSpaces(unittest.TestCase):

    def test_partition_space(self):
        large_cas = LargeSpaces()
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        orb_indices = [3, 4, 5, 6, 7, 8, 9, 10]
        partitioned_indices = large_cas._partition_space(orb_indices)
        self.assertEqual(partitioned_indices, [[3, 4], [5, 6], [7, 8], [9, 10]])

    def test_generate_spaces(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 14
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 2, 0, 0, 0]

        large_cas = LargeSpaces()
        occupied_orbital_indices = [0, 1, 2, 3, 4]
        virtual_orbital_indices = [5, 6, 7]
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        large_cas.generate_spaces(occupied_orbital_indices, [], virtual_orbital_indices, molecule)
        self.assertEqual(large_cas.occupations,
                         [[2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0]])
        self.assertEqual(large_cas.orbital_indices,
                         [[0, 1, 5, 6], [0, 1, 5, 7], [2, 3, 5, 6], [2, 3, 5, 7], [3, 4, 5, 6], [3, 4, 5, 7]])

    def test_generate_spaces_hard(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 15
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 1, 0, 0, 0]
        molecule.spin_multiplicity = 2
        molecule.charge = 1

        large_cas = LargeSpaces()
        occupied_orbital_indices = [0, 1, 2, 3]
        singly_occupied_orbital_indices = [4]
        virtual_orbital_indices = [5, 6, 7]
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        large_cas.generate_spaces(occupied_orbital_indices, singly_occupied_orbital_indices,
                                  virtual_orbital_indices, molecule)
        self.assertEqual(
            large_cas.occupations,
            [[2, 2, 1, 0, 0], [2, 2, 1, 0, 0], [2, 2, 1, 0, 0], [2, 2, 1, 0, 0]]
        )
        self.assertEqual(
            large_cas.orbital_indices,
            [[0, 1, 4, 5, 6], [0, 1, 4, 5, 7], [2, 3, 4, 5, 6], [2, 3, 4, 5, 7]]
        )

    def test_generate_spaces_hard_2(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 15
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 1, 0, 0, 0]
        molecule.spin_multiplicity = 2
        molecule.charge = 1

        large_cas = LargeSpaces()
        occupied_orbital_indices = [0, 1, 2, 3]
        singly_occupied_orbital_indices = [4]
        virtual_orbital_indices = [5, 6, 7]
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 6
        large_cas.generate_spaces(occupied_orbital_indices, singly_occupied_orbital_indices,
                                  virtual_orbital_indices, molecule)
        self.assertEqual(
            large_cas.occupations,
            [[2, 2, 1, 0, 0, 0], [2, 2, 1, 0, 0, 0]]
        )
        self.assertEqual(
            large_cas.orbital_indices,
            [[0, 1, 4, 5, 6, 7], [2, 3, 4, 5, 6, 7]]
        )

    def test_separate_spaces(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 14
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 2, 2, 2, 0, 0, 0]

        cas_indices = [2, 3, 4, 5, 6, 7, 8, 9]
        cas_occupation = [2, 2, 2, 2, 2, 0, 0, 0]
        # indices = [0, 1, 2, 3, 4, 5, 6, 7]
        cas = ActiveSpace(cas_occupation, cas_indices)

        large_cas = LargeSpaces()
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        occupied_orbs, singly_occupied_orbs, virtual_orbs = large_cas.separate_space(cas, molecule)
        self.assertEqual(occupied_orbs, [2, 3, 4, 5, 6])
        self.assertEqual(singly_occupied_orbs, [])
        self.assertEqual(virtual_orbs, [7, 8, 9])

    def test_separate_spaces_hard(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 14
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 2, 2, 1, 0, 0, 0]
        molecule.spin_multiplicity = 2
        molecule.charge = 1

        cas_indices = [2, 3, 4, 5, 6, 7, 8, 9]
        cas_occupation = [2, 2, 2, 2, 1, 0, 0, 0]
        # indices = [0, 1, 2, 3, 4, 5, 6, 7]
        cas = ActiveSpace(cas_occupation, cas_indices)

        large_cas = LargeSpaces()
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        occupied_orbs, singly_occupied_orbs, virtual_orbs = large_cas.separate_space(cas, molecule)
        self.assertEqual(occupied_orbs, [2, 3, 4, 5])
        self.assertEqual(singly_occupied_orbs, [6])
        self.assertEqual(virtual_orbs, [7, 8, 9])

    def test_get_spaces(self):
        molecule = Molecule(atom_list=["N", "N"])
        molecule.electrons = 14
        molecule.valence_orbitals = 8
        molecule.core_orbitals = 2
        molecule.occupation = [2, 2, 2, 2, 2, 2, 2, 0, 0, 0]

        cas_indices = [2, 3, 4, 5, 6, 7, 8, 9]
        cas_occupation = [2, 2, 2, 2, 2, 0, 0, 0]
        # indices = [0, 1, 2, 3, 4, 5, 6, 7]
        cas = ActiveSpace(cas_occupation, cas_indices)

        large_cas = LargeSpaces()
        # max_orb = 4 only for testing
        large_cas.max_orbitals = 4
        occupations, indices = large_cas.get_spaces(cas, molecule)
        self.assertEqual(occupations,
                         [[2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0], [2, 2, 0, 0]])
        self.assertEqual(indices,
                         [[2, 3, 7, 8], [2, 3, 7, 9], [4, 5, 7, 8], [4, 5, 7, 9], [5, 6, 7, 8], [5, 6, 7, 9]])


if __name__ == "__main__":
    unittest.main()
