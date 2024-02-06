# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import shutil
import unittest

from scine_autocas import Autocas
from scine_autocas.__main__ import MainFunctions
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas


class TestMainClass(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.dummy_project_dir = self.path + "/dummy_dir"
        os.makedirs(self.dummy_project_dir, exist_ok=True)

    def tearDown(self):
        try:
            shutil.rmtree(self.dummy_project_dir)
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(self.path + "/autocas_project")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(self.path + "/../../autocas_project")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(self.path + "/molcas_scratch")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(self.path + "/../../molcas_scratch")
        except FileNotFoundError:
            pass

    def test_conventional(self):
        mainfunctions = MainFunctions()
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule=molecule)
        molcas = Molcas(molecules=[molecule])
        molcas.project_name = "test_conventional"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        # check if molcas is available
        try:
            occupation, indices = mainfunctions.conventional(autocas, molcas)
        except OSError:
            print("No MOLCAS binary found.")
        else:
            # fmt: off
            self.assertEqual(indices, [4, 5, 6, 7, 8, 9])
            self.assertEqual(occupation, [2, 2, 2, 0, 0, 0])
            # fmt: on

    def test_large_cas(self):
        mainfunctions = MainFunctions()
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule=molecule)
        autocas.large_spaces.seed = 1
        autocas.large_spaces.max_orbitals = 4
        molcas = Molcas(molecules=[molecule])
        molcas.project_name = "test_large_cas"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        # check if molcas is available
        try:
            occupation, indices = mainfunctions.large_cas(autocas, molcas)
        except OSError:
            print("No MOLCAS binary found.")
        else:
            # fmt: off
            self.assertEqual(indices, [4, 5, 6, 7, 8, 9])
            self.assertEqual(occupation, [2, 2, 2, 0, 0, 0])
            # fmt: on

    def test_excited_states(self):
        mainfunctions = MainFunctions()
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule=molecule)
        molcas = Molcas(molecules=[molecule])
        molcas.project_name = "test_excited_states"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.settings.n_excited_states = 4
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.make_ci_root()
        # check if molcas is available
        try:
            occupation, indices = mainfunctions.excited_states(autocas, molcas)
        except OSError:
            print("No MOLCAS binary found.")
        else:
            # fmt: off
            self.assertEqual(indices, [4, 5, 6, 7, 8, 9, 2, 3])
            self.assertEqual(occupation, [2, 2, 2, 0, 0, 0, 2, 2])
            # fmt: on

    def test_large_cas_excited_states(self):
        mainfunctions = MainFunctions()
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule=molecule)
        autocas.large_spaces.max_orbitals = 4
        autocas.large_spaces.seed = 1
        molcas = Molcas(molecules=[molecule])
        molcas.settings.n_excited_states = 2
        molcas.project_name = "test_large_cas_excited_states"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.make_ci_root()
        # check if molcas is available
        try:
            occupation, indices = mainfunctions.large_cas_excited_states(autocas, molcas)
        except OSError:
            print("No MOLCAS binary found.")
        else:
            # fmt: off
            self.assertEqual(sorted(indices), [3, 4, 5, 6, 7, 8, 9])
            self.assertEqual(occupation, [2, 2, 2, 2, 0, 0, 0])
            # fmt: on

    def test_main_function_with_yaml(self):
        path = os.path.dirname(os.path.abspath(__file__))
        test_yaml = path + "/files/working.yml"
        settings = {"yaml_input": test_yaml}
        print(settings)
        main_functions = MainFunctions()
        try:
            main_functions.main(settings)
        except FileNotFoundError:
            print("YAML file includes absolute path, hence this test is")
            print("not possible to work.")
        except OSError:
            print("No MOLCAS binary found.")
        else:
            self.assertTrue(abs(main_functions.results["energy"][0] - -108.96105402) < 1e-6)

    def test_main_function_without_yaml(self):
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = path + "/files/n2.xyz"
        # to make test possible

        settings = {"method": "dmrgci", "xyz_file": xyz, "yaml_input": None, "basis_set": "cc-pvdz"}
        main_functions = MainFunctions()
        try:
            main_functions.main(settings)
        except OSError:
            print("No MOLCAS binary found.")
        else:
            self.assertTrue(abs(main_functions.results["energy"] - -108.73628442571527) < 1e-6)


if __name__ == "__main__":
    unittest.main()
