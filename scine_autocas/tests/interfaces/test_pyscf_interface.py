# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import os
import shutil
import unittest

import pytest

from scine_autocas.interfaces.pyscf import PyscfInterface
from scine_autocas.utils.molecule import Molecule

try:
    from pyscf import gto, mcscf, scf
except ModuleNotFoundError:
    pass


class TestPyscfInterface(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.dummy_dir = self.path + "/dummy_dir"
        os.makedirs(self.dummy_dir)
        os.chdir(self.dummy_dir)
        self.xyz_name = f"{self.dummy_dir}/n2_molecule_test.xyz"
        n2_xyz = open(self.xyz_name, "w")
        n2_xyz.write("2\n\n")
        n2_xyz.write("N 0 0 0\n")
        n2_xyz.write("N 0 0 2\n")
        n2_xyz.close()
        self.charge = 1
        self.spin_mult = 4

        self.autocas_mol = Molecule(self.xyz_name)
        self.autocas_mol.charge = self.charge
        self.autocas_mol.spin_multiplicity = self.spin_mult
        self.autocas_mol.update()

    def tearDown(self):
        os.chdir(self.path)
        try:
            shutil.rmtree(self.dummy_dir)
        except FileNotFoundError:
            pass

    def test_molecule(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")
        mol = gto.Mole()
        mol.atom = self.xyz_name
        mol.build(charge=self.charge, spin=self.spin_mult-1, basis="cc-pvdz")

        pyscf_interface = PyscfInterface(self.autocas_mol)
        pyscf_interface._build_molecule()
        self.assertEqual(self.autocas_mol.charge, mol.charge)
        self.assertEqual(self.autocas_mol.spin_multiplicity-1, mol.spin)

    def test_hartree_fock(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")
        mol = gto.Mole()
        mol.atom = self.xyz_name
        mol.build(charge=self.charge, spin=self.spin_mult-1, basis="cc-pvdz")
        hf = scf.RHF(mol)
        e = hf.scf()

        pyscf_interface = PyscfInterface(self.autocas_mol)
        pyscf_interface._build_molecule()
        pyscf_interface._initial_orbitals_impl()
        energy = pyscf_interface.pyscf_hf.e_tot
        self.assertAlmostEqual(e, energy)

    def test_cas(self):
        try:
            import pyscf  # noqa: F401
        except ModuleNotFoundError:
            pytest.skip("PySCF not found")

        mol = gto.Mole()
        mol.atom = self.xyz_name
        mol.build(charge=0, spin=0, basis="cc-pvdz")
        hf = scf.RHF(mol)
        _ = hf.scf()
        mc = mcscf.CASCI(hf, 6, 6)
        e = mc.kernel()[0]

        pyscf_interface = PyscfInterface(self.autocas_mol)
        pyscf_interface._dumper.activated = False
        pyscf_interface.pyscf_mol = mol
        pyscf_interface.pyscf_hf = hf
        energy = pyscf_interface.calculate()[0]
        energy = pyscf_interface.calculate([2, 2, 2, 0, 0, 0], [4, 5, 6, 7, 8, 9])[0][0]
        self.assertAlmostEqual(e, energy)


if __name__ == "__main__":
    unittest.main()
