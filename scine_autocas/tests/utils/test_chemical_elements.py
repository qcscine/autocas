# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import unittest

from scine_autocas.utils.chemical_elements import Elements


class TestMolecule(unittest.TestCase):
    def test_get_core_orb_labels(self):
        chem_ele = Elements()
        self.assertEqual(chem_ele.get_core_orb_labels("Fe"), ["1s", "2s", "2p", "3s", "3p"])
        self.assertEqual(chem_ele.get_core_orb_labels("H"), [])

    def test_get_core_orbitals(self):
        chem_ele = Elements()
        self.assertEqual(chem_ele.get_core_orbitals("Fe"), 9)
        self.assertEqual(chem_ele.get_core_orbitals("H"), 0)

    def test_valence_orb_labels(self):
        chem_ele = Elements(double_d_shell=False)
        self.assertEqual(chem_ele.get_valence_orb_labels("Fe"), ["4s", "3d", "4p"])
        self.assertEqual(chem_ele.get_valence_orb_labels("H"), ["1s"])
        chem_ele = Elements(double_d_shell=True)
        self.assertEqual(chem_ele.get_valence_orb_labels("Fe"), ["4s", "3d", "4p", "4d"])
        self.assertEqual(chem_ele.get_valence_orb_labels("H"), ["1s"])

    def test_get_valence_orbitals(self):
        chem_ele = Elements(double_d_shell=False)
        self.assertEqual(chem_ele.get_valence_orbitals("Fe"), 9)
        self.assertEqual(chem_ele.get_valence_orbitals("H"), 1)
        chem_ele = Elements(double_d_shell=True)
        self.assertEqual(chem_ele.get_valence_orbitals("Fe"), 14)
        self.assertEqual(chem_ele.get_valence_orbitals("H"), 1)

    def test_get_electrons(self):
        chem_ele = Elements()
        self.assertEqual(chem_ele.get_electrons("Fe"), 26)
        self.assertEqual(chem_ele.get_electrons("H"), 1)

    def test_get(self):
        chem_ele = Elements(double_d_shell=False)
        self.assertEqual(chem_ele.get("Fe")["name"], "Fe")
        self.assertEqual(chem_ele.get("Fe")["atomic number"], 26)
        self.assertEqual(chem_ele.get("Fe")["number of core orbitals"], 9)
        self.assertEqual(chem_ele.get("Fe")["number of valence orbitals"], 9)
        self.assertEqual(chem_ele.get("Fe")["core orbitals"], ["1s", "2s", "2p", "3s", "3p"])
        self.assertEqual(chem_ele.get("Fe")["valence orbitals"], ["4s", "3d", "4p"])
        self.assertEqual(chem_ele.get("H")["name"], "H")
        self.assertEqual(chem_ele.get("H")["atomic number"], 1)
        self.assertEqual(chem_ele.get("H")["number of core orbitals"], 0)
        self.assertEqual(chem_ele.get("H")["number of valence orbitals"], 1)
        self.assertEqual(chem_ele.get("H")["core orbitals"], [])
        self.assertEqual(chem_ele.get("H")["valence orbitals"], ["1s"])
        chem_ele = Elements(double_d_shell=True)
        self.assertEqual(chem_ele.get("Fe")["name"], "Fe")
        self.assertEqual(chem_ele.get("Fe")["atomic number"], 26)
        self.assertEqual(chem_ele.get("Fe")["number of core orbitals"], 9)
        self.assertEqual(chem_ele.get("Fe")["number of valence orbitals"], 14)
        self.assertEqual(chem_ele.get("Fe")["core orbitals"], ["1s", "2s", "2p", "3s", "3p"])
        self.assertEqual(chem_ele.get("Fe")["valence orbitals"], ["4s", "3d", "4p", "4d"])
        self.assertEqual(chem_ele.get("H")["name"], "H")
        self.assertEqual(chem_ele.get("H")["atomic number"], 1)
        self.assertEqual(chem_ele.get("H")["number of core orbitals"], 0)
        self.assertEqual(chem_ele.get("H")["number of valence orbitals"], 1)
        self.assertEqual(chem_ele.get("H")["core orbitals"], [])
        self.assertEqual(chem_ele.get("H")["valence orbitals"], ["1s"])


if __name__ == "__main__":
    unittest.main()
