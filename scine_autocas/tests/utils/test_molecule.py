# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

from scine_autocas.utils.molecule import Molecule


class TestMolecule(unittest.TestCase):
    def setUp(self):
        path = os.path.dirname(os.path.abspath(__file__))
        self.xyz_name = f"{path}/n2_molecule_test.xyz"
        n2_xyz = open(self.xyz_name, "w")
        n2_xyz.write("2\n\n")
        n2_xyz.write("N 0 0 0\n")
        n2_xyz.write("N 0 0 2\n")
        n2_xyz.close()

    def tearDown(self):
        try:
            os.remove(self.xyz_name)
        except FileNotFoundError:
            pass

    def test_update(self):
        molecule = Molecule(self.xyz_name)
        self.assertEqual(molecule.electrons, 14)
        self.assertEqual(molecule.charge, 0)
        self.assertEqual(molecule.spin_multiplicity, 1)

        molecule.charge = -1
        molecule.spin_multiplicity = 4
        molecule.update()
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.core_orbitals, 2)

        molecule.ecp_electrons = 2
        molecule.update()
        self.assertEqual(molecule.atoms, ["N", "N"])
        self.assertEqual(molecule.core_orbitals, 1)
        self.assertEqual(molecule.valence_orbitals, 8)
        self.assertEqual(molecule.electrons, 13)
        self.assertEqual(molecule.spin_multiplicity, 4)
        self.assertEqual(molecule.double_d_shell, True)
        self.assertEqual(molecule.unit, "ang")
        self.assertEqual(molecule.charge, -1)
        self.assertEqual(molecule.occupation, [2, 2, 2, 2, 2, 1, 1, 1, 0])

    def test_initialize(self):
        molecule = Molecule(self.xyz_name)
        self.assertEqual(molecule.charge, 0)
        self.assertEqual(molecule.spin_multiplicity, 1)
        molecule.charge = -1
        molecule.spin_multiplicity = 4
        molecule.update()
        self.assertEqual(molecule.charge, -1)
        self.assertEqual(molecule.spin_multiplicity, 4)
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.core_orbitals, 2)
        molecule.ecp_electrons = 2
        molecule.update()
        self.assertEqual(molecule.charge, -1)
        self.assertEqual(molecule.spin_multiplicity, 4)
        self.assertEqual(molecule.electrons, 13)
        self.assertEqual(molecule.core_orbitals, 1)
        molecule.ecp_electrons = 0
        molecule.update()
        self.assertEqual(molecule.charge, -1)
        self.assertEqual(molecule.spin_multiplicity, 4)
        self.assertEqual(molecule.electrons, 15)
        self.assertEqual(molecule.core_orbitals, 2)

        molecule2 = Molecule(atom_list=["N", "N"], charge=-1, spin_multiplicity=4)
        self.assertEqual(molecule.charge, molecule2.charge)
        self.assertEqual(molecule.spin_multiplicity, molecule2.spin_multiplicity)
        self.assertEqual(molecule.electrons, molecule2.electrons)
        self.assertEqual(molecule.core_orbitals, molecule2.core_orbitals)

        molecule3 = Molecule(atom_list=["N", "N"], charge=-1, spin_multiplicity=4, ecp_electrons=2)
        self.assertEqual(molecule3.charge, molecule2.charge)
        self.assertEqual(molecule3.spin_multiplicity, molecule2.spin_multiplicity)
        self.assertEqual(molecule3.electrons, molecule2.electrons - 2)
        self.assertEqual(molecule3.core_orbitals, molecule2.core_orbitals - 1)

        molecule4 = Molecule(atom_list=["N", "N"], charge=-1, spin_multiplicity=4, ecp_electrons=2)
        molecule4.update()
        self.assertEqual(molecule4.charge, molecule2.charge)
        self.assertEqual(molecule4.spin_multiplicity, molecule2.spin_multiplicity)
        self.assertEqual(molecule4.electrons, molecule2.electrons - 2)
        self.assertEqual(molecule4.core_orbitals, molecule2.core_orbitals - 1)

        molecule5 = Molecule(atom_list=["N", "N"], charge=-1, spin_multiplicity=4)
        molecule5.ecp_electrons = 2
        molecule5.update()
        self.assertEqual(molecule5.charge, molecule2.charge)
        self.assertEqual(molecule5.spin_multiplicity, molecule2.spin_multiplicity)
        self.assertEqual(molecule5.electrons, molecule2.electrons - 2)
        self.assertEqual(molecule5.core_orbitals, molecule2.core_orbitals - 1)


if __name__ == "__main__":
    unittest.main()
