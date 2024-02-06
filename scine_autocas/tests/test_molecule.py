# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

from scine_autocas.autocas_utils.molecule import Molecule


class TestAutocasClasses(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))

    def test_xyz_file(self):
        molecule = Molecule(self.path + "/files/n2.xyz")
        self.assertEqual(molecule.core_orbitals, 2)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 14)
        self.assertEqual(molecule.spin_multiplicity, 1)
        # fmt: off
        self.assertEqual(molecule.occupation, [2, 2, 2, 2, 2, 2, 2, 0, 0, 0])
        # fmt: on

    def test_charge(self):
        molecule = Molecule(self.path + "/files/n2.xyz", charge=2)
        self.assertEqual(molecule.core_orbitals, 2)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 12)
        self.assertEqual(molecule.spin_multiplicity, 1)
        # fmt: off
        self.assertEqual(molecule.occupation, [2, 2, 2, 2, 2, 2, 0, 0, 0, 0])
        # fmt: on

    def test_multiplicity(self):
        molecule = Molecule(self.path + "/files/n2.xyz", charge=-1, spin_multiplicity=4)
        self.assertEqual(molecule.core_orbitals, 2)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.spin_multiplicity, 4)
        # fmt: off
        self.assertEqual(molecule.occupation, [2, 2, 2, 2, 2, 2, 1, 1, 1, 0])
        # fmt: on

    def test_additional_basis_functions(self):
        molecule = Molecule(self.path + "/files/n2.xyz", charge=-1, spin_multiplicity=4, n_basis_functions=15)
        self.assertEqual(molecule.core_orbitals, 2)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.spin_multiplicity, 4)
        # fmt: off
        self.assertEqual(
            molecule.occupation, [2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        )
        # fmt: on

    def test_update(self):
        molecule = Molecule(self.path + "/files/n2.xyz", charge=-1, spin_multiplicity=4, n_basis_functions=15)
        molecule.n_basis_functions = 16
        molecule.update()
        self.assertEqual(molecule.core_orbitals, 2)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.spin_multiplicity, 4)
        # fmt: off
        self.assertEqual(molecule.occupation, [2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
        # fmt: on


if __name__ == "__main__":
    unittest.main()
