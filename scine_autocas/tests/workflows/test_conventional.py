# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np
import pytest

from scine_autocas import Autocas, Molecule
from scine_autocas.interfaces import PyscfInterface
from scine_autocas.workflows.conventional import ClassicWorkflow


class TestConventionalWorkflow(unittest.TestCase):
    def setUp(self):

        path = os.path.dirname(os.path.abspath(__file__))
        self._xyz_name = f"{path}/n2_molecule_test.xyz"
        n2_xyz = open(self._xyz_name, "w")
        n2_xyz.write("2\n\n")
        n2_xyz.write("N 0 0 0\n")
        n2_xyz.write("N 0 0 2\n")
        n2_xyz.close()
        self.n2 = Molecule(self._xyz_name)

    def tearDown(self):
        try:
            os.remove(self._xyz_name)
        except FileNotFoundError:
            pass

    def test_n2(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")
        autocas = Autocas(self.n2)
        interface = PyscfInterface(self.n2)
        interface._dumper.activated = False
        workflow = ClassicWorkflow(autocas, interface)
        workflow.run()
        # for i in workflow.results:
        #     print(i, workflow.results[i])

        self.assertEqual(workflow.results["initial_occupation"], [2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(workflow.results["initial_orbital_indices"], [2, 3, 4, 5, 6, 7, 8, 9])
        self.assertTrue(np.allclose(workflow.results["initial_energy"], -108.75697657579984))
        print(workflow.results["initial_s1"])
        self.assertTrue(np.allclose(
            workflow.results["initial_s1"],
            [0.04850772, 0.04476995, 0.76259808, 1.1422812, 1.14228117, 1.14299842, 1.1429984, 0.78126309],
            atol=1e-5, rtol=1e-1))
        self.assertEqual(workflow.results["final_occupation"],  [2, 2, 2, 0, 0, 0])
        self.assertEqual(workflow.results["final_orbital_indices"],  [4, 5, 6, 7, 8, 9])
        self.assertTrue(np.allclose(workflow.results["final_energy"],  [-108.74113359189914]))
        # print(workflow.results["final_s1"])
        # self.assertTrue(np.allclose(
        #     workflow.results["final_s1"],
        #     [0.72811438, 1.13636396, 1.13636398, 1.13650361, 1.13650359, 0.72746021]))


if __name__ == "__main__":
    unittest.main()
