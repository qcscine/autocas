# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
import unittest

import pytest

from scine_autocas.io.actions import run_action
from scine_autocas.io.input_handler import InputHandler


class TestRunActions(unittest.TestCase):
    def setUp(self):
        path = os.path.dirname(os.path.abspath(__file__))
        os.chdir(path)
        self.cli_dict = {
            'action': 'run',
            "large_cas": False,
            "large_cas_seed": 42,
            "large_cas_max_orbitals": 5,
            "cas_method": "dmrgscf",
            'yaml': None,
            'interface': 'pyscf',
            'dmrg_bond_dimension': 100,
            'dmrg_sweeps': 100,
            "init_fiedler": False,
            'dump': False,
            "init_cas_method": "dmrgci",
            # '/' prevents xyz file to be correct to absolute path
            'xyz_file': "files/n2.xyz",
            'spin_multiplicity': 2,
            'charge': -1,
            'basis_set': 'cc-pvtz',
            'unit': 'bohr'
        }

    # just test if it runs through
    def test_conventional(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")
        inp = InputHandler(self.cli_dict)
        run_action(inp)

    # just test if it runs through
    def test_large_cas(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")
        self.cli_dict["large_cas"] = True
        inp = InputHandler(self.cli_dict)
        run_action(inp)


if __name__ == "__main__":
    unittest.main()
