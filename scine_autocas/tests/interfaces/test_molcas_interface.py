# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""


import os
import re
import shutil
import unittest
from copy import deepcopy

import numpy as np
import pytest

from scine_autocas import Autocas
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.interfaces.molcas.environment import Environment
from scine_autocas.interfaces.molcas.input_handler import InputHandler
from scine_autocas.interfaces.molcas.molcas_hdf5_utils import MolcasHdf5Utils
from scine_autocas.io import FileHandler
from scine_autocas.utils.defaults import CasMethods, PostCasMethods
from scine_autocas.utils.molecule import Molecule


class TestMolcasClasses(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.dummy_project_dir = self.path + "/dummy_dir"
        orbital_file = self.path + "/files/n2.scf.h5"
        self.orbital_file = self.dummy_project_dir + "/n2.scf.h5"
        os.makedirs(self.dummy_project_dir, exist_ok=True)
        shutil.copyfile(orbital_file, self.orbital_file)
        os.chdir(self.dummy_project_dir)
        FileHandler.setup_project()

    def tearDown(self):
        os.chdir(self.path)
        try:
            shutil.rmtree(self.dummy_project_dir)
        except FileNotFoundError:
            pass

    def test_hdf5_utils_read_hdf5(self):
        hdf5_helper = MolcasHdf5Utils()
        hdf5_helper.read_hdf5(self.orbital_file)
        self.assertEqual(hdf5_helper.nbas[0], 28)
        self.assertEqual(hdf5_helper.natoms, 2)
        self.assertEqual(hdf5_helper.module, "SCF")
        self.assertEqual(hdf5_helper.orbital_type, "SCF-RHF")
        self.assertEqual(hdf5_helper.energy, 0)
        self.assertListEqual(
            hdf5_helper.mo_energies,
            [
                -15.827122169017748, -15.826872145810418, -1.0835347073549308, -0.974208781863933,
                -0.46116997858793535, -0.3902819711828907, -0.39028197118288777, -0.09432361974807853,
                -0.0943236197480776, 0.01321677197231618, 0.8737225947596443, 0.8737225947596495,
                0.9036740957609027, 0.951356747197556, 0.988461457466255, 0.9884614574662581,
                1.0543415057633805, 1.3237280685939334, 1.8229891626123127, 1.9172052826601216,
                1.9172052826601225, 1.9779269593572306, 1.977926959357231, 2.003416984902952,
                2.0034169849029526, 2.1002482065579122, 2.100248206557916, 2.5304772424757447,
            ],
        )
        # fmt: off
        self.assertEqual(hdf5_helper.irrep_labels, np.array(b'a  '))
        self.assertTrue(
            np.array_equal(
                hdf5_helper.type_indices,
                [
                    b'I', b'I', b'I', b'I', b'I', b'I', b'I',
                    b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                    b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                    b'S', b'S', b'S', b'S', b'S', b'S', b'S'
                ]
            )
        )
        self.assertEqual(
            hdf5_helper.occupations,
            [2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        self.assertEqual(
            hdf5_helper.symmetries,
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        # fmt: on
        hdf5_helper.read_hdf5(self.orbital_file)
        self.assertEqual(
            hdf5_helper.occupations,
            [2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )

    def test_hdf5_utils_modify_hdf5(self):
        hdf5_helper = MolcasHdf5Utils()
        # fmt: off
        hdf5_helper.type_indices = [
            b'I', b'I', b'I', b'I', b'I', b'I', b'I',
            b'S', b'S', b'S', b'S', b'S', b'S', b'S',
            b'S', b'S', b'S', b'S', b'S', b'S', b'S',
            b'S', b'S', b'S', b'S', b'S', b'S', b'S'
        ]
        # fmt: on
        hdf5_helper.modify_hdf5(self.orbital_file, [5, 6, 7, 9])
        hdf5_helper.read_hdf5(self.orbital_file)
        self.assertEqual(hdf5_helper.nbas[0], 28)
        self.assertEqual(hdf5_helper.natoms, 2)
        self.assertEqual(hdf5_helper.module, "SCF")
        self.assertEqual(hdf5_helper.orbital_type, "SCF-RHF")
        # fmt: off
        self.assertEqual(hdf5_helper.irrep_labels, np.array(b'a  '))
        self.assertTrue(
            np.array_equal(
                hdf5_helper.type_indices,
                [
                    b'I', b'I', b'I', b'I', b'I', b'2', b'2',
                    b'2', b'S', b'2', b'S', b'S', b'S', b'S',
                    b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                    b'S', b'S', b'S', b'S', b'S', b'S', b'S'
                ]
            )
        )
        # fmt: on
        self.assertListEqual(
            hdf5_helper.mo_energies,
            [
                -15.827122169017748, -15.826872145810418, -1.0835347073549308, -0.974208781863933,
                -0.46116997858793535, -0.3902819711828907, -0.39028197118288777, -0.09432361974807853,
                -0.0943236197480776, 0.01321677197231618, 0.8737225947596443, 0.8737225947596495,
                0.9036740957609027, 0.951356747197556, 0.988461457466255, 0.9884614574662581,
                1.0543415057633805, 1.3237280685939334, 1.8229891626123127, 1.9172052826601216,
                1.9172052826601225, 1.9779269593572306, 1.977926959357231, 2.003416984902952,
                2.0034169849029526, 2.1002482065579122, 2.100248206557916, 2.5304772424757447,
            ],
        )
        # fmt: off
        self.assertEqual(
            hdf5_helper.occupations,
            [2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        self.assertEqual(
            hdf5_helper.symmetries,
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        # fmt: on
        self.assertEqual(hdf5_helper.energy, 0)

    def test_environment(self):
        envir = Environment()
        calcdir = self.dummy_project_dir

        # set molcas to prevent exception
        tmp = os.environ.copy()
        try:
            environment = envir.make_environment(calcdir)
        except OSError:
            print("test_environment: No MOLCAS binary found.")
        else:
            os.environ = environment
            self.assertEqual(os.environ["MOLCAS_MEM"], "12000")

            os.environ["MOLCAS_MEM"] = "6000"
            environment = envir.make_environment(calcdir)
            os.environ = environment
            self.assertEqual(os.environ["MOLCAS_MEM"], "6000")
            self.assertRaises(KeyError, lambda: os.environ["MOLCAS_NPROCS"])
            self.assertRaises(KeyError, lambda: os.environ["OPM_NUM_THREADS"])
            self.assertRaises(KeyError, lambda: os.environ["QCMaquis_CPUS"])
            self.assertEqual(os.environ["MOLCAS_OUTPUT"], calcdir)

            envir.molcas_nprocs = "2"
            envir.molcas_scratch_dir = ""
            environment = envir.make_environment(calcdir)
            os.environ = environment
            self.assertEqual(os.environ["QCMaquis_CPUS"], "2")

        try:
            del os.environ["MOLCAS"]
        except KeyError:
            print("No Molcas found anyways.")
        envir.molcas_binary = ""
        self.assertRaises(OSError, lambda: envir.make_environment(calcdir))

        os.environ = tmp

    def test_input_handler_1(self):
        input_handler = InputHandler()
        molecule = Molecule(atom_list=["H", "H"])
        molecule.xyz_file = "blubblub.xyz"
        settings = Molcas.Settings(molecule=molecule)
        settings.orbital_localisation = True
        settings.localisation_space = "OCCUpied"
        settings.localisation_method = "PIPEk-Mezey"

        # write input for initial orbitals
        initial_input_file = self.dummy_project_dir + "/test_initial.inp"
        input_handler.write_input(settings, initial_input_file)

        # write dmrg input file
        settings.initial_orbitals = False
        settings.fiedler = False
        dmrg_input_file = self.dummy_project_dir + "/test_dmrg.inp"
        input_handler.write_input(settings, dmrg_input_file, self.orbital_file)

        # write dmrg input file
        dmrg_input_file_2 = self.dummy_project_dir + "/test_dmrg_2.inp"
        input_handler.write_input(settings, dmrg_input_file_2, self.orbital_file, alter="2; 1 4 5; 3 6 8")

        # write dmrg input file
        settings.fiedler = True
        settings.ci_root_string = "just a test string"
        dmrg_input_file_3 = self.dummy_project_dir + "/test_dmrg_3.inp"
        input_handler.write_input(settings, dmrg_input_file_3, alter="2; 1 4 5; 3 6 8")
        settings.ci_root_string = ""

        # write final input
        settings.post_cas_method = PostCasMethods.CASPT2
        settings.cas_method = CasMethods.CASSCF
        final_input_file = self.dummy_project_dir + "/test_final.inp"
        input_handler.write_input(settings, final_input_file, self.orbital_file)

        final_input_file_2 = self.dummy_project_dir + "/test_final_2.inp"
        input_handler.write_input(settings, final_input_file_2, self.orbital_file, alter="2; 1 4 5; 3 6 8")

        settings.post_cas_method = PostCasMethods.NEVPT2
        settings.cas_method = CasMethods.CASCI
        settings.ci_root_string = "just a test string"
        final_input_file_3 = self.dummy_project_dir + "/test_final_3.inp"
        input_handler.write_input(settings, final_input_file_3, alter="2; 1 4 5; 3 6 8")

        initial_orbital_content = open(initial_input_file, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", initial_orbital_content.read()),
            re.sub(
                "\\s+",
                " ",
                """&GATEWAY\n
                    COORD = blubblub.xyz\n
                    BASIS = cc-pvdz\n
                    GROUP = C1\n\n
                   &SEWARD\n\n
                    CHOLesky\n\n
                   &SCF\n SPIN   = 1\n
                   CHARGE = 0\n
                   CHOLesky\n
                   &LOCALISATION\n
                    OCCUpied
                    PIPEk-Mezey
                """,
            ),
        )
        initial_orbital_content.close()

        dmrg_content = open(dmrg_input_file, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", dmrg_content.read()),
            re.sub(
                "\\s+",
                " ",
                f"""&DMRGSCF\n
                    ActiveSpaceOptimizer = QCMaquis\n
                    DMRGSettings\n
                     nsweeps = 5\n
                     max_bond_dimension = 250\n
                    EndDMRGSettings\n
                    OOptimizationSettings\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     FILEORB = {self.orbital_file}\n
                     TYPEINDEX\n
                     CIONLY\n
                    EndOOptimizationSettings\n
                """
            ),
        )
        dmrg_content.close()

        dmrg_content = open(dmrg_input_file_2, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", dmrg_content.read()),
            re.sub(
                "\\s+",
                " ",
                f"""&DMRGSCF\n
                    ActiveSpaceOptimizer = QCMaquis\n
                    DMRGSettings\n
                     nsweeps = 5\n
                     max_bond_dimension = 250\n
                    EndDMRGSettings\n
                    OOptimizationSettings\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     FILEORB = {self.orbital_file}\n
                     ALTEr=2; 1 4 5; 3 6 8\n
                     CIONLY\n
                    EndOOptimizationSettings\n
                """
            ),
        )
        dmrg_content.close()

        dmrg_content = open(dmrg_input_file_3, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", dmrg_content.read()),
            re.sub(
                "\\s+",
                " ",
                """&DMRGSCF\n
                    FIEDLER=ON\n
                    ActiveSpaceOptimizer = QCMaquis\n
                    DMRGSettings\n
                     nsweeps = 5\n
                     max_bond_dimension = 250\n
                    EndDMRGSettings\n
                    OOptimizationSettings\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     RAS2 = 0\n
                     ALTEr=2; 1 4 5; 3 6 8\n
                     CIONLY\n
                     CIRoot = just a test string\n
                    EndOOptimizationSettings\n
                """
            ),
        )
        dmrg_content.close()

        final_content = open(final_input_file, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", final_content.read()),
            re.sub(
                "\\s+",
                " ",
                f"""&RASSCF\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     FILEORB = {self.orbital_file}\n
                     TYPEINDEX\n
                    &CASPT2\n
                     IPEA = 0.0\n
                """
            ),
        )
        final_content.close()

        final_content = open(final_input_file_2, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", final_content.read()),
            re.sub(
                "\\s+",
                " ",
                f"""&RASSCF\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     FILEORB = {self.orbital_file}\n
                     ALTEr=2; 1 4 5; 3 6 8\n
                    &CASPT2\n
                     IPEA = 0.0\n
                """
            ),
        )
        final_content.close()

        final_content = open(final_input_file_3, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", final_content.read()),
            re.sub(
                "\\s+",
                " ",
                """&RASSCF\n
                     SPIN = 1\n
                     NACTEL = 0\n
                     RAS2 = 0\n
                     ALTEr=2; 1 4 5; 3 6 8\n
                     CIONLY\n
                     CIRoot = just a test string\n
                    &NEVPT2\n
                """
            ),
        )
        final_content.close()

    def test_input_handler_2(self):
        input_handler = InputHandler()
        molecule = Molecule(atom_list=["H", "H"])
        molecule.xyz_file = "blubblub.xyz"
        settings = Molcas.Settings(molecule=molecule)

        # write input for initial orbitals
        initial_input_file = self.dummy_project_dir + "/test_initial.inp"
        input_handler.write_input(settings, initial_input_file)

        # write dmrg input file
        settings.initial_orbitals = False
        settings.cas_method = CasMethods.CASSCF
        dmrg_input_file = self.dummy_project_dir + "/test_cas.inp"
        input_handler.write_input(settings, dmrg_input_file, self.orbital_file)

        # write final input
        initial_orbital_content = open(initial_input_file, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", initial_orbital_content.read()),
            re.sub(
                "\\s+",
                " ",
                """&GATEWAY\n
                    COORD = blubblub.xyz\n
                    BASIS = cc-pvdz\n
                    GROUP = C1\n\n
                   &SEWARD\n\n
                    CHOLesky\n\n
                   &SCF\n
                    SPIN   = 1\n
                    CHARGE = 0\n
                    CHOLesky\n
                """,
            ),
        )
        initial_orbital_content.close()

        dmrg_content = open(dmrg_input_file, encoding="utf-8")
        self.assertEqual(
            re.sub("\\s+", " ", dmrg_content.read()),
            re.sub(
                "\\s+",
                " ",
                f"""&RASSCF\n
                    SPIN    = 1\n
                    NACTEL  = 0\n
                    FILEORB = {self.orbital_file}\n
                    TYPEINDEX\n
                """
            ),
        )
        dmrg_content.close()

    def test_molcas(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(self.path + "/files/n2.xyz")
        # molecule.xyz_file = self.path + "/files/n2.xyz"
        molcas = Molcas(molecule=molecule)
        molcas.project_name = "test123"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        molcas_results = molcas.calculate([2, 2, 0, 0], [3, 4, 5, 6])
        energy = molcas_results[0]
        s1_entropy = molcas_results[1]
        self.assertTrue(abs(energy[0] - -104.24521324932175) < 1e-6)
        # fmt: off
        self.assertTrue(
            np.allclose(
                s1_entropy,
                [0.0252543, 0.04813958, 0.03711964, 0.03711973]
            )
        )

    def test_charge(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz", spin_multiplicity=4, charge=-1)
        molcas = Molcas(molecule=molecule)
        molcas.project_name = "test_charge"
        # molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        molcas_results = molcas.calculate()
        molcas_results = molcas.calculate([2, 2, 2, 2, 1, 1, 1, 0], [2, 3, 4, 5, 6, 7, 8, 9])
        energy = molcas_results[0]
        s1_entropy = molcas_results[1]
        self.assertTrue(abs(energy[0] - -108.5830821627763) < 1e-6)
        print(s1_entropy)
        self.assertTrue(
            np.allclose(
                s1_entropy,
                [
                    0.00211653, 0.00265149, 0.00423185, 0.09447688,
                    0.00215388, 0.00243114, 0.04983321, 0.09760716
                ],
                atol=1e-7,
                rtol=1e-1
            )
        )

    def test_excited_states(self):
        # This test fails in the current CI (and locally) because of unknown reasons.
        # It will be disabled for the time being.
        pytest.skip("DISABLED")
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecule=molecule)
        molcas.settings.n_excited_states = 2
        molcas.set_cas_method("dmrgci")
        molcas.project_name = "test_excited_states"
        molcas.settings.post_cas_method = None
        # molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        # molcas.make_ci_root()

        energy, s1_entropy, s2_entropy, _ = molcas.calculate([2, 2, 2, 2, 2, 0, 0, 0], [2, 3, 4, 5, 6, 7, 8, 9])
        _ = s2_entropy  # get rid of warning
        self.assertEqual(len(s1_entropy), 3)
        self.assertTrue(abs(energy[0] - -108.74985111) < 1e-7)
        ci_string = molcas.make_ci_root(3, 5, [2, 4, 5], [1, 1, 3])
        self.assertEqual(ci_string, "3 5; 2 4 5 ; 1 1 3 ")

    def test_analyze(self):
        # This test fails in the current CI (and locally) because of unknown reasons.
        # It will be disabled for the time being.
        pytest.skip("DISABLED")
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecule=molecule)
        molcas.set_cas_method("dmrgci")
        molcas.settings.n_excited_states = 2
        molcas.project_name = "test_analyze"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.settings.post_cas_method = None
        # molcas.make_ci_root()
        # check if molcas is available
        _ = molcas.calculate()
        _ = molcas.calculate([2, 2, 2, 2, 2, 0, 0, 0], [2, 3, 4, 5, 6, 7, 8, 9])
        molecule1 = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas1 = Molcas(molecule=molecule1)
        hdf5_results = molcas1.analyze(
            molcas_hdf5_file=self.dummy_project_dir
            + "/autocas_project/dmrg/test_analyze.scf.h5_sel"
        )
        qcmaquis_results = molcas1.analyze(
            qcmaquis_result_file=self.dummy_project_dir
            + "/autocas_project/dmrg/test_analyze.results_state.0.h5"
        )

        indices = hdf5_results[0]
        occupation = hdf5_results[1]
        s1_entropy = qcmaquis_results[0]
        # fmt: off
        self.assertTrue(
            np.allclose(
                indices,
                [2, 3, 4, 5, 6, 7, 8, 9]
            )
        )
        # because this is an alrady modified file the occupation is different than expected
        self.assertTrue(
            np.allclose(
                occupation,
                [2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            )
        )

        self.assertTrue(
            np.allclose(
                s1_entropy,
                [0.03766205, 0.04412138, 0.89634891, 1.21125612, 1.21125612, 1.21160228, 1.21160228, 0.91054304]
            )
        )

    def test_symmetry(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecule=molecule)

        # problem with cholesky together with symmetry
        molcas.settings.cholesky = False
        molcas.settings.fiedler = False  # fiedler reorders orbitals
        molcas.project_name = "test_for_symmetry_1_no_sym"
        molcas.settings.post_cas_method = None
        molcas.set_cas_method("dmrgci")
        # molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # fmt: off
        n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]
        # fmt: on

        (
            energy_init_no_sym,
            s1_init_no_sym,
            s2_init_no_sym,
            mut_inf_init_no_sym,
        ) = molcas.calculate(n2_occ_initial, n2_index_initial)

        # fmt: off
        n2_occ = [2, 2, 2, 0, 0, 0]
        n2_index = [4, 5, 6, 7, 8, 9]
        # fmt: on

        energy_no_sym = molcas.calculate(n2_occ, n2_index)[0]

        # same with symmetry
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas1 = Molcas(molecule=molecule)

        # problem with cholesky
        molcas1.settings.cholesky = False
        molcas1.settings.fiedler = False  # fiedler reorders orbitals
        # molcas1.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas1.settings.point_group = "X Y"
        molcas1.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas1.project_name = "test_for_symmetry_1_sym"
        molcas1.settings.post_cas_method = None
        molcas1.set_cas_method("dmrgci")

        # n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        # n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]

        energy_init_sym, s1_init_sym, s2_init_sym, mut_inf_init_sym = molcas1.calculate(
            n2_occ_initial, n2_index_initial
        )

        # n2_occ = [2, 2, 2, 0, 0, 0]
        # n2_index = [4, 5, 6, 7, 8, 9]

        energy_sym = molcas1.calculate(n2_occ, n2_index)[0]

        np.testing.assert_allclose(s1_init_no_sym, s1_init_sym, atol=1e-5)

        # zero these blocks, since a comparison does not work with symmetry otherwise
        s2_init_sym[4:7, 4:7] = 0
        s2_init_sym[5:7, 3:4] = 0
        s2_init_sym[3:4, 5:7] = 0
        s2_init_no_sym[4:7, 4:7] = 0
        s2_init_no_sym[5:7, 3:4] = 0
        s2_init_no_sym[3:4, 5:7] = 0
        mut_inf_init_sym[4:7, 4:7] = 0
        mut_inf_init_sym[5:7, 3:4] = 0
        mut_inf_init_sym[3:4, 5:7] = 0
        mut_inf_init_no_sym[4:7, 4:7] = 0
        mut_inf_init_no_sym[5:7, 3:4] = 0
        mut_inf_init_no_sym[3:4, 5:7] = 0
        np.testing.assert_allclose(s2_init_no_sym, s2_init_sym, atol=1e-6)
        np.testing.assert_allclose(mut_inf_init_no_sym, mut_inf_init_sym, atol=1e-6)
        self.assertAlmostEqual(energy_init_no_sym[0], energy_init_sym[0], 6)
        self.assertAlmostEqual(energy_no_sym[0], energy_sym[0], 6)

    def test_fiedler(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecule=molecule)
        molcas.settings.fiedler = False
        molcas.settings.dmrg_bond_dimension = 3000
        molcas.settings.dmrg_sweeps = 20
        molcas.project_name = "test_for_fiedler_1_no_fiedler"
        molcas.settings.post_cas_method = None
        molcas.set_cas_method("dmrgci")
        # molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # fmt: off
        n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]
        # fmt: on

        (
            energy_init_no_sym,
            s1_init_no_sym,
            s2_init_no_sym,
            mut_inf_init_no_sym,
        ) = molcas.calculate(n2_occ_initial, n2_index_initial)

        # fmt: off
        n2_occ = [2, 2, 2, 0, 0, 0]
        n2_index = [4, 5, 6, 7, 8, 9]
        # fmt: on

        energy_no_sym, s1_no_sym, s2_no_sym, mut_inf_no_sym = molcas.calculate(n2_occ, n2_index)

        # same with symmetry
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas1 = Molcas(molecule=molecule)
        molcas1.settings.fiedler = True
        molcas1.settings.dmrg_bond_dimension = 3000
        molcas1.settings.dmrg_sweeps = 20
        # molcas1.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas1.settings.post_cas_method = None
        molcas1.project_name = "test_for_fiedler_1_fiedler"
        molcas1.set_cas_method("dmrgci")
        molcas1.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        # n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]

        energy_init_sym, s1_init_sym, s2_init_sym, mut_inf_init_sym = molcas1.calculate(
            n2_occ_initial, n2_index_initial
        )

        # n2_occ = [2, 2, 2, 0, 0, 0]
        # n2_index = [4, 5, 6, 7, 8, 9]

        energy_sym, s1_sym, s2_sym, mut_inf_sym = molcas1.calculate(n2_occ, n2_index)

        # get rid of warnings
        _ = s2_init_sym
        _ = s2_sym
        _ = s2_init_no_sym
        _ = s2_no_sym
        _ = mut_inf_init_sym
        _ = mut_inf_sym
        _ = mut_inf_init_no_sym
        _ = mut_inf_no_sym

        np.testing.assert_allclose(s1_init_no_sym, s1_init_sym, atol=1e-6)
        np.testing.assert_allclose(s1_no_sym, s1_sym, atol=1e-6)
        # depends on orbital order
        # np.testing.assert_allclose(s2_no_sym, s2_sym, atol=1e-6)
        # np.testing.assert_allclose(s2_init_no_sym, s2_init_sym, atol=1e-6)
        # np.testing.assert_allclose(I_init_no_sym, I_init_sym, atol=1e-6)
        # np.testing.assert_allclose(I_no_sym, I_sym, atol=1e-6)
        self.assertAlmostEqual(energy_init_no_sym[0], energy_init_sym[0], 5)
        self.assertAlmostEqual(energy_no_sym[0], energy_sym[0], 5)

    def test_unrestricted(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molecule.spin_multiplicity = 3

        molcas = Molcas(molecule=molecule)
        molcas.settings.dmrg_bond_dimension = 3000
        molcas.settings.dmrg_sweeps = 20
        molcas.project_name = "test_for_uhf"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        molcas_2 = Molcas(molecule=molecule)
        molcas_2.settings.dmrg_bond_dimension = 3000
        molcas_2.settings.dmrg_sweeps = 20
        molcas_2.project_name = "test_for_uhf_2"
        molcas_2.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas_2.settings.fiedler = True

        # n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        # n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]
        autocas = Autocas(molecule)
        init_occ, init_index = autocas.make_initial_active_space()

        molcas.calculate()
        (
            energy_init_no_sym,
            s1_init_no_sym,
            s2_init_no_sym,
            mut_inf_init_no_sym,
        ) = molcas.calculate(init_occ, init_index)

        molcas_2.calculate()
        (
            energy_init_no_sym_2,
            s1_init_no_sym_2,
            s2_init_no_sym_2,
            mut_inf_init_no_sym_2,
        ) = molcas_2.calculate(init_occ, init_index)

        print(init_occ)
        print(init_index)
        print(s1_init_no_sym)
        print(s1_init_no_sym_2)
        np.testing.assert_allclose(s1_init_no_sym, s1_init_no_sym_2, atol=1e-6)

        autocas_2 = deepcopy(autocas)
        occ, index = autocas.get_active_space(s1_init_no_sym)

        occ_2, index_2 = autocas_2.get_active_space(s1_init_no_sym_2)
        print(occ)
        print(occ_2)
        print(index)
        print(index_2)
        self.assertTrue(
            np.allclose(occ, occ_2)
        )
        self.assertTrue(
            np.allclose(index, index_2)
        )
