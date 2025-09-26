# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import yaml

from scine_autocas.io.input_handler import InputHandler


class TestInputHandler(unittest.TestCase):
    def setUp(self):
        self.yaml_name = "test.yaml"

        self.molecule_dict = {
            "Molecule": {
                "xyz_file": "/abc.xyz",
                "unit": "bohr",
                "ecp_electrons": 0,
                "charge": -1,
                "spin_multiplicity": 2,
                "double_d_shell": True
            }
        }

        self.autocas_dict = {
            "AutoCAS": {
                "large_cas": False,
                "large_cas_average_entanglement": False,
                "large_cas_seed": 42,
                "large_cas_max_orbitals": 4,
                "plateau_values": 10,
                "single_reference_threshold": 0.14,
                "threshold_step": 0.01,
                "weak_correlation_threshold": 0.02,
            }
        }

        self.interface_dict = {
            "Interface": {
                "basis_set": "cc-pvtz",
                "interface": "molcas",
                "dmrg_sweeps": 100,
                "dmrg_bond_dimension": 400,
                "init_dmrg_sweeps": 5,
                "init_dmrg_bond_dimension": 250,
                "uhf": False,
                "n_excited_states": 0,
                "dmrg_solver": "QCMaquis",
                "init_cas_method": "dmrgscf",
                "cas_method": "dmrgscf",
                "post_cas_method": "nevpt2",
                "dump": False,
                "init_fiedler": False,
                "fiedler": True,
                "init_orbital_order": None,
                "orbital_order": "3,2,5,4,1",
            }

        }

        self.cli_dict = {
            "large_cas": False,
            "large_cas_average_entanglement": False,
            "large_cas_seed": 42,
            "large_cas_max_orbitals": 4,
            "cas_method": "dmrgscf",
            'yaml': None,
            'action': 'run',
            'interface': 'molcas',
            'dmrg_bond_dimension': 400,
            'dmrg_sweeps': 100,
            "init_fiedler": False,
            'dump': False,
            "post_cas_method": "nevpt2",
            "init_cas_method": "dmrgscf",
            # '/' prevents xyz file to be correct to absolute path
            'xyz_file': "/abc.xyz",
            'spin_multiplicity': 2,
            'charge': -1,
            'basis_set': 'cc-pvtz',
            'unit': 'bohr',
            "orbital_order": "3,2,5,4,1",
        }

        self.yaml_dict = {
            "AutoCAS": {
                "large_cas": False,
                "large_cas_average_entanglement": False,
                "large_cas_seed": 42,
                "plateau_values": 10,
                "large_cas_max_orbitals": 4,
            },
            "Defaults": {
                "Molecule": {
                    "unit": "ang",
                    "charge": -1,
                    "spin_multiplicity": 2,
                    "double_d_shell": True
                },
                "AutoCAS": {
                    "plateau_values": 9,
                    "single_reference_threshold": 0.14,
                    "threshold_step": 0.01,
                    "weak_correlation_threshold": 0.02,
                },
                "Interface": {
                    "basis_set": "cc-pvtz",
                    "interface": "molcas",
                    "dmrg_sweeps": 20,
                    "dmrg_bond_dimension": 400,
                    "init_dmrg_sweeps": 5,
                    "init_dmrg_bond_dimension": 250,
                    "cas_method": "dmrgci",
                    "init_cas_method": "dmrgscf",
                    "orbital_order": "3,2,5,4,1"
                }
            },
            "Interface": {
                "basis_set": "cc-pvtz",
                "interface": "molcas",
                "dmrg_sweeps": 100,
                "init_fiedler": False,
                "cas_method": "dmrgscf",
                "post_cas_method": "nevpt2",
            },
            "Molecule": {
                "xyz_file": "/abc.xyz",
                "unit": "bohr",
            },
        }

    def tearDown(self):
        try:
            os.remove(self.yaml_name)
        except FileNotFoundError:
            pass

    def test_molecule_from_cli(self):
        inp = InputHandler(self.cli_dict)
        test_mol_dict = inp.get_molecule_options()
        self.assertEqual(test_mol_dict, self.molecule_dict)

    def test_autocas_from_cli(self):
        inp = InputHandler(self.cli_dict)
        test_autocas_dict = inp.get_autocas_options()
        self.assertEqual(test_autocas_dict, self.autocas_dict)

    def test_interface_from_cli(self):
        inp = InputHandler(self.cli_dict)
        test_interface_dict = inp.get_interface_options()
        self.assertEqual(test_interface_dict, self.interface_dict)

    def test_yaml_only(self):
        with open(self.yaml_name, 'w') as outfile:
            yaml.dump(self.yaml_dict, outfile)
        inp = InputHandler({"yaml": self.yaml_name, "action": "run", "xyz_file": "/abc.xyz"})

        test_interface_dict = inp.get_interface_options()
        self.assertEqual(test_interface_dict, self.interface_dict)
        test_autocas_dict = inp.get_autocas_options()
        self.assertEqual(test_autocas_dict, self.autocas_dict)
        test_mol_dict = inp.get_molecule_options()
        self.assertEqual(test_mol_dict, self.molecule_dict)

    def test_cli_plus_yaml(self):
        self.yaml_dict["Molecule"]["xyz_file"] = "haha.xyz"
        self.yaml_dict["Interface"]["basis_set"] = "6-31g"
        self.yaml_dict["AutoCAS"]["large_cas_average_entanglement"] = True
        self.yaml_dict["Defaults"]["Molecule"]["charge"] = 9
        self.yaml_dict["Defaults"]["Interface"]["cas_method"] = "casscf"
        self.yaml_dict["Defaults"]["AutoCAS"]["large_cas_max_orbitals"] = 100
        self.yaml_dict["Defaults"]["Interface"]["orbital_order"] = "3,2,5,4,1"
        with open(self.yaml_name, 'w') as outfile:
            yaml.dump(self.yaml_dict, outfile)
        self.cli_dict["yaml"] = self.yaml_name
        self.cli_dict["action"] = "run"
        inp = InputHandler(self.cli_dict)
        test_interface_dict = inp.get_interface_options()
        self.assertEqual(test_interface_dict, self.interface_dict)
        test_autocas_dict = inp.get_autocas_options()
        self.assertEqual(test_autocas_dict, self.autocas_dict)
        test_mol_dict = inp.get_molecule_options()
        self.assertEqual(test_mol_dict, self.molecule_dict)


if __name__ == "__main__":
    unittest.main()
