# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import os
import shutil
import unittest

try:
    from pyscf import gto, mcscf, scf
    pyscf_avail = True
except ModuleNotFoundError:
    pyscf_avail = False

try:
    from scine_qcmaquis.pyscf_interface import QcMaquis
    qcmaquis_avail = True
except ModuleNotFoundError:
    qcmaquis_avail = False


class TestQCMaquisBasics(unittest.TestCase):
    def setUp(self):
        self.path = os.path.join(os.path.dirname(__file__), "qcmaquis_test")
        os.makedirs(self.path, exist_ok=True)
        os.chdir(self.path)

    def tearDown(self):
        os.chdir(self.path)
        try:
            shutil.rmtree(self.path)
        except FileNotFoundError:
            pass

    @unittest.skipIf(not pyscf_avail, "pyscf not found")
    @unittest.skipIf(not qcmaquis_avail, "qcmaquis not found")
    def test_single_sweep(self):
        mol = gto.M(atom="""N 0 0 0; N 0 0 2""", basis="cc-pvdz", charge=0, spin=0)
        mol.build()
        rhf = scf.RHF(mol)
        rhf.kernel()

        # only one sweep DMRG
        maquis = QcMaquis(mol)
        maquis.parameters.set("nsweeps", 1)
        maquis.parameters.set("max_bond_dimension", 2000)
        maquis.file_path = None
        cas_unco = mcscf.CASCI(rhf, 10, 8)
        cas_unco.fcisolver = maquis
        e_unco = cas_unco.kernel()[0]

        # converge DMRG
        maquis2 = QcMaquis(mol)
        maquis2.parameters.set("nsweeps", 100)
        maquis2.parameters.set("max_bond_dimension", 2000)
        maquis2.file_path = None
        cas_conv = mcscf.CASCI(rhf, 10, 8)
        cas_conv.fcisolver = maquis2
        e_conv = cas_conv.kernel()[0]

        cas = mcscf.CASCI(rhf, 10, 8)
        e_ref = cas.kernel()[0]

        # assert difference between converged and unconverged DMRG
        assert e_conv < e_unco
        # assert that DMRG converges to CASCI energy
        self.assertAlmostEqual(e_ref, e_conv)
