# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

from scine_autocas.autocas_utils.input_handler import InputHandler


class TestInputHandler(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.test_yaml = self.path + "/files/test.yaml"

    def test_get_objects_from_input(self):
        input_handler = InputHandler(self.test_yaml)
        input_handler.settings_dir["molecule"]["xyz_file"] = self.path + "/" + \
            input_handler.settings_dir["molecule"]["xyz_file"]
        input_handler.settings_dir["interface"]["settings"]["xyz_file"] = self.path + "/" + \
            input_handler.settings_dir["interface"]["settings"]["xyz_file"]

        test_molecule = input_handler.get_molecule()
        self.assertEqual(test_molecule.charge, 2)
        self.assertEqual(test_molecule.spin_multiplicity, 3)
        self.assertEqual(test_molecule.core_orbitals, 2)
        self.assertEqual(test_molecule.electrons, 12)
        self.assertEqual(test_molecule.atoms, ["N", "N"])
        # fmt: off
        self.assertEqual(
            test_molecule.occupation,
            [2, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        )
        # fmt: on

        test_autocas = input_handler.get_autocas()
        self.assertEqual(test_autocas.plateau_values, 20)
        self.assertEqual(test_autocas.threshold_step, 0.005)
        self.assertEqual(test_autocas.diagnostics.weak_correlation_threshold, 0.03)
        self.assertEqual(test_autocas.diagnostics.single_reference_threshold, 0.15)
        self.assertEqual(test_autocas.large_spaces.max_orbitals, 25)

        test_molcas = input_handler.get_interface()
        self.assertEqual(test_molcas.dump, False)
        self.assertEqual(test_molcas.project_name, "blub")
        self.assertEqual(test_molcas.settings.dmrg_bond_dimension, 350)
        self.assertEqual(test_molcas.settings.dmrg_sweeps, 4)
        self.assertEqual(test_molcas.settings.basis_set, "6-31g")
        self.assertEqual(test_molcas.settings.method, "dmrg-scf")
        self.assertEqual(test_molcas.settings.post_cas_method, "caspt2")
        self.assertEqual(
            test_molcas.settings.work_dir, "/random/path/for/test"
        )
        self.assertEqual(
            test_molcas.settings.xyz_file,
            self.path + "/files/n2.xyz",
        )
        self.assertEqual(test_molcas.settings.point_group, "C2h")
        self.assertEqual(test_molcas.settings.ipea, 1.1)
        self.assertEqual(test_molcas.settings.orbital_localisation, True)
        self.assertEqual(test_molcas.settings.localisation_space, "Virtual")
        self.assertEqual(test_molcas.settings.localisation_method, "BOYS")
        self.assertEqual(test_molcas.settings.cholesky, True)
        self.assertEqual(test_molcas.settings.uhf, True)
        self.assertEqual(test_molcas.settings.fiedler, True)
        self.assertEqual(test_molcas.settings.n_excited_states, 4)
        self.assertEqual(test_molcas.settings.ci_root_string, "4 4 1")


if __name__ == "__main__":
    unittest.main()
