# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import pathlib
import re
import shutil
import unittest

import numpy as np

from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.interfaces.molcas.dumper import Dumper
from scine_autocas.interfaces.molcas.environment import Environment
from scine_autocas.interfaces.molcas.input_handler import InputHandler
from scine_autocas.interfaces.molcas.molcas_hdf5_utils import MolcasHdf5Utils


class TestMolcasClasses(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.dummy_project_dir = self.path + "/dummy_dir"
        orbital_file = self.path + "/files/n2.scf.h5"
        self.orbital_file = self.dummy_project_dir + "/n2.scf.h5"
        os.makedirs(self.dummy_project_dir, exist_ok=True)
        shutil.copyfile(orbital_file, self.orbital_file)

    def tearDown(self):
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
                -15.827122169017748,
                -15.826872145810418,
                -1.0835347073549308,
                -0.974208781863933,
                -0.46116997858793535,
                -0.3902819711828907,
                -0.39028197118288777,
                -0.09432361974807853,
                -0.0943236197480776,
                0.01321677197231618,
                0.8737225947596443,
                0.8737225947596495,
                0.9036740957609027,
                0.951356747197556,
                0.988461457466255,
                0.9884614574662581,
                1.0543415057633805,
                1.3237280685939334,
                1.8229891626123127,
                1.9172052826601216,
                1.9172052826601225,
                1.9779269593572306,
                1.977926959357231,
                2.003416984902952,
                2.0034169849029526,
                2.1002482065579122,
                2.100248206557916,
                2.5304772424757447,
            ],
        )
        # fmt: off
        self.assertEqual(hdf5_helper.irrep_labels, np.array(b'a  '))
        self.assertTrue(
            np.array_equal(
                hdf5_helper.type_indices,
                [b'I', b'I', b'I', b'I', b'I', b'I', b'I',
                 b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                 b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                 b'S', b'S', b'S', b'S', b'S', b'S', b'S']
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
         b'S', b'S', b'S', b'S', b'S', b'S', b'S']
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
                [b'I', b'I', b'I', b'I', b'I', b'2', b'2',
                 b'2', b'S', b'2', b'S', b'S', b'S', b'S',
                 b'S', b'S', b'S', b'S', b'S', b'S', b'S',
                 b'S', b'S', b'S', b'S', b'S', b'S', b'S']
            )
        )
        # fmt: on
        self.assertListEqual(
            hdf5_helper.mo_energies,
            [
                -15.827122169017748,
                -15.826872145810418,
                -1.0835347073549308,
                -0.974208781863933,
                -0.46116997858793535,
                -0.3902819711828907,
                -0.39028197118288777,
                -0.09432361974807853,
                -0.0943236197480776,
                0.01321677197231618,
                0.8737225947596443,
                0.8737225947596495,
                0.9036740957609027,
                0.951356747197556,
                0.988461457466255,
                0.9884614574662581,
                1.0543415057633805,
                1.3237280685939334,
                1.8229891626123127,
                1.9172052826601216,
                1.9172052826601225,
                1.9779269593572306,
                1.977926959357231,
                2.003416984902952,
                2.0034169849029526,
                2.1002482065579122,
                2.100248206557916,
                2.5304772424757447,
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
            del os.environ["MOLCAS"]
        except KeyError:
            print("No Molcas found anyways.")
        self.assertRaises(EnvironmentError, lambda: envir.make_environment(calcdir))

        os.environ["MOLCAS"] = "dummy"
        environment = envir.make_environment(calcdir)
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
        self.assertEqual(os.environ["MOLCAS_NPROCS"], "2")
        self.assertEqual(os.environ["OMP_NUM_THREADS"], "2")
        self.assertEqual(os.environ["QCMaquis_CPUS"], "2")
        self.assertEqual(os.environ["WorkDir"], "/tmp")

        os.environ = tmp

    def test_input_handler_1(self):
        input_handler = InputHandler()
        molecule = Molecule(atoms=["H", "H"])
        settings = Molcas.Settings(molecules=[molecule])
        settings.xyz_file = "blubblub.xyz"
        settings.orbital_localisation = True
        settings.localisation_space = "OCCUpied"
        settings.localisation_method = "PIPEk-Mezey"

        # write input for initial orbitals
        initial_input_file = self.dummy_project_dir + "/test_initial.inp"
        input_handler.write_input(settings, initial_input_file)

        # write dmrg input file
        settings.initial_orbitals = False
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
        settings.post_cas_method = "caspt2"
        settings.method = "casscf"
        final_input_file = self.dummy_project_dir + "/test_final.inp"
        input_handler.write_input(settings, final_input_file, self.orbital_file)

        final_input_file_2 = self.dummy_project_dir + "/test_final_2.inp"
        input_handler.write_input(settings, final_input_file_2, self.orbital_file, alter="2; 1 4 5; 3 6 8")

        settings.post_cas_method = "NEVpt2"
        settings.method = "casci"
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
        molecule = Molecule(atoms=["H", "H"])
        settings = Molcas.Settings(molecules=[molecule])
        settings.xyz_file = "blubblub.xyz"

        # write input for initial orbitals
        initial_input_file = self.dummy_project_dir + "/test_initial.inp"
        input_handler.write_input(settings, initial_input_file)

        # write dmrg input file
        settings.initial_orbitals = False
        settings.method = "casscf"
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

    def test_dumper(self):
        dumper = Dumper()
        dumper.project_dir = self.dummy_project_dir
        dumper.project_name = "xyz"
        dumper.create_project_dir()
        dumper.setup_sub_dir()
        dummy_orbital_file = dumper.current_dir + "/xyz.scf.h5"
        pathlib.Path(dummy_orbital_file).touch()
        dumper.setup_sub_dir()
        dumper.setup_sub_dir()
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/xyz/initial/xyz.scf.h5"))
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/xyz/dmrg/xyz.scf.h5_sel"))
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/xyz/final/xyz.scf.h5_sel"))

    def test_dumper_large_cas(self):
        dumper = Dumper()
        dumper.project_dir = self.dummy_project_dir
        dumper.project_name = "abc"
        dumper.large_cas = True
        dumper.create_project_dir()
        dumper.setup_sub_dir()

        dummy_orbital_file = dumper.current_dir + "/abc.scf.h5"
        pathlib.Path(dummy_orbital_file).touch()

        dumper.setup_sub_dir()
        dumper.setup_sub_dir()
        dumper.large_cas = False
        dumper.setup_sub_dir()
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/abc/initial/abc.scf.h5"))
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/abc/dmrg_1/abc.scf.h5_sel"))
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/abc/dmrg_2/abc.scf.h5_sel"))
        self.assertTrue(os.path.isfile(self.dummy_project_dir + "/abc/final/abc.scf.h5_sel"))
        dumper.orbital_file = self.dummy_project_dir + "/abc/initial/non_existing_file"
        self.assertRaises(FileNotFoundError, lambda: dumper.copy_orbital_file())

    def test_molcas(self):
        molecule = Molecule(self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])
        molcas.project_name = "test123"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # if molcas is installed
        try:
            # fmt: off
            molcas_results = molcas.calculate([2, 2, 0, 0], [3, 4, 5, 6])
            # fmt: on
        except OSError:
            print("No MOLCAS binary found.")
        else:
            energy = molcas_results[0]
            s1_entropy = molcas_results[1]
            self.assertTrue(abs(energy - -104.24521324932175) < 1e-6)
            # fmt: off
            self.assertTrue(
                np.allclose(
                    s1_entropy,
                    [0.0252543, 0.04813958, 0.03711964, 0.03711973]
                )
            )
            # fmt: on

    def test_charge(self):
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])
        molcas.project_name = "test_charge"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.settings.spin_multiplicity = 4
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.settings.charge = -1

        # if molcas is installed
        try:
            # fmt: off
            molcas_results = molcas.calculate([2, 2, 2, 2, 1, 1, 1, 0], [2, 3, 4, 5, 6, 7, 8, 9])
            # fmt: on
        except OSError:
            print("No MOLCAS binary found.")
        else:
            energy = molcas_results[0]
            s1_entropy = molcas_results[1]
            self.assertTrue(abs(energy - -108.5830821627763) < 1e-6)
            # fmt: off
            self.assertTrue(
                np.allclose(
                    s1_entropy,
                    [0.00211653, 0.00265149, 0.00423185, 0.09447688, 0.00215388, 0.00243114, 0.04983321, 0.09760716]
                )
            )
            # fmt: on

    def test_excited_states(self):
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])
        molcas.settings.n_excited_states = 3
        molcas.project_name = "test_excited_states"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.make_ci_root()

        try:
            # fmt: off
            energy, s1_entropy, s2_entropy, _ = molcas.calculate([2, 2, 2, 2, 2, 0, 0, 0], [2, 3, 4, 5, 6, 7, 8, 9])
            # fmt: on
            _ = s2_entropy  # get rid of warning
        except OSError:
            print("No MOLCAS binary found.")
        else:
            self.assertEqual(len(s1_entropy), 3)
            # fmt: off
            self.assertTrue(abs(energy - -108.74985111) < 1e-7)
            ci_string = molcas.make_ci_root(3, 5, [2, 4, 5], [1, 1, 3])
            self.assertEqual(ci_string, "3 5; 2 4 5 ; 1 1 3 ")
            # fmt: on

    def test_analyze(self):
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])
        molcas.settings.n_excited_states = 3
        molcas.project_name = "test_analyze"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
        molcas.make_ci_root()
        # check if molcas is available
        try:
            # fmt: off
            _ = molcas.calculate([2, 2, 2, 2, 2, 0, 0, 0], [2, 3, 4, 5, 6, 7, 8, 9])
            # fmt: on
        except OSError:
            print("No MOLCAS binary found.")
        else:
            molecule1 = Molecule(xyz_file=self.path + "/files/n2.xyz")
            molcas1 = Molcas(molecules=[molecule1])
            hdf5_results = molcas1.analyze(
                molcas_hdf5_file=self.dummy_project_dir
                + "/test_analyze/dmrg/test_analyze.scf.h5_sel"
            )
            qcmaquis_results = molcas1.analyze(
                qcmaquis_result_file=self.dummy_project_dir
                + "/test_analyze/dmrg/test_analyze.results_state.0.h5"
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
            # fmt: on

    def test_symmetry(self):
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])

        # problem with cholesky together with symmetry
        molcas.settings.cholesky = False
        molcas.settings.fiedler = False  # fiedler reorders orbitals
        molcas.project_name = "test_for_symmetry_1_no_sym"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # fmt: off
        n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]
        # fmt: on

        # check if molcas is installed, otherwise cancel test
        try:
            (
                energy_init_no_sym,
                s1_init_no_sym,
                s2_init_no_sym,
                mut_inf_init_no_sym,
            ) = molcas.calculate(n2_occ_initial, n2_index_initial)
        except OSError:
            print("No MOLCAS binary found.")
        else:

            # fmt: off
            n2_occ = [2, 2, 2, 0, 0, 0]
            n2_index = [4, 5, 6, 7, 8, 9]
            # fmt: on

            energy_no_sym, s1_no_sym, s2_no_sym, mut_inf_no_sym = molcas.calculate(n2_occ, n2_index)

            # same with symmetry
            molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
            molcas1 = Molcas(molecules=[molecule])

            # problem with cholesky
            molcas1.settings.cholesky = False
            molcas1.settings.fiedler = False  # fiedler reorders orbitals
            molcas1.settings.work_dir = self.dummy_project_dir
            molcas1.settings.xyz_file = self.path + "/files/n2.xyz"
            molcas1.settings.point_group = "X Y"
            molcas1.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"
            molcas1.project_name = "test_for_symmetry_1_sym"

            # n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
            # n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]

            energy_init_sym, s1_init_sym, s2_init_sym, mut_inf_init_sym = molcas1.calculate(
                n2_occ_initial, n2_index_initial
            )

            # n2_occ = [2, 2, 2, 0, 0, 0]
            # n2_index = [4, 5, 6, 7, 8, 9]

            energy_sym, s1_sym, s2_sym, mut_inf_sym = molcas1.calculate(n2_occ, n2_index)

            np.testing.assert_allclose(s1_init_no_sym, s1_init_sym, atol=1e-5)
            np.testing.assert_allclose(s1_no_sym, s1_sym, atol=1e-5)

            # zero these blocks, since a comparison does not work with symmetry otherwise
            s2_init_sym[4:7, 4:7] = 0
            s2_init_sym[5:7, 3:4] = 0
            s2_init_sym[3:4, 5:7] = 0
            s2_init_no_sym[4:7, 4:7] = 0
            s2_init_no_sym[5:7, 3:4] = 0
            s2_init_no_sym[3:4, 5:7] = 0
            s2_sym[2:5, 2:5] = 0
            s2_sym[3:5, 1:2] = 0
            s2_sym[1:2, 3:5] = 0
            s2_no_sym[2:5, 2:5] = 0
            s2_no_sym[3:5, 1:2] = 0
            s2_no_sym[1:2, 3:5] = 0
            mut_inf_init_sym[4:7, 4:7] = 0
            mut_inf_init_sym[5:7, 3:4] = 0
            mut_inf_init_sym[3:4, 5:7] = 0
            mut_inf_init_no_sym[4:7, 4:7] = 0
            mut_inf_init_no_sym[5:7, 3:4] = 0
            mut_inf_init_no_sym[3:4, 5:7] = 0
            mut_inf_sym[2:5, 2:5] = 0
            mut_inf_sym[3:5, 1:2] = 0
            mut_inf_sym[1:2, 3:5] = 0
            mut_inf_no_sym[2:5, 2:5] = 0
            mut_inf_no_sym[3:5, 1:2] = 0
            mut_inf_no_sym[1:2, 3:5] = 0
            np.testing.assert_allclose(s2_init_no_sym, s2_init_sym, atol=1e-6)
            np.testing.assert_allclose(s2_no_sym, s2_sym, atol=1e-6)
            np.testing.assert_allclose(mut_inf_init_no_sym, mut_inf_init_sym, atol=1e-6)
            np.testing.assert_allclose(mut_inf_no_sym, mut_inf_sym, atol=1e-6)
            self.assertAlmostEqual(energy_init_no_sym, energy_init_sym, 6)
            self.assertAlmostEqual(energy_no_sym, energy_sym, 6)

    def test_fiedler(self):
        molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
        molcas = Molcas(molecules=[molecule])
        molcas.settings.fiedler = False
        molcas.settings.dmrg_bond_dimension = 3000
        molcas.settings.dmrg_sweeps = 20
        molcas.project_name = "test_for_fiedler_1_no_fiedler"
        molcas.settings.work_dir = self.dummy_project_dir
        molcas.settings.xyz_file = self.path + "/files/n2.xyz"
        molcas.environment.molcas_scratch_dir = self.dummy_project_dir + "/scratch"

        # fmt: off
        n2_occ_initial = [2, 2, 2, 2, 2, 0, 0, 0]
        n2_index_initial = [2, 3, 4, 5, 6, 7, 8, 9]
        # fmt: on

        # check if molcas is installed, otherwise cancel test
        try:
            (
                energy_init_no_sym,
                s1_init_no_sym,
                s2_init_no_sym,
                mut_inf_init_no_sym,
            ) = molcas.calculate(n2_occ_initial, n2_index_initial)
        except OSError:
            print("No MOLCAS binary found.")
        else:

            # fmt: off
            n2_occ = [2, 2, 2, 0, 0, 0]
            n2_index = [4, 5, 6, 7, 8, 9]
            # fmt: on

            energy_no_sym, s1_no_sym, s2_no_sym, mut_inf_no_sym = molcas.calculate(n2_occ, n2_index)

            # same with symmetry
            molecule = Molecule(xyz_file=self.path + "/files/n2.xyz")
            molcas1 = Molcas(molecules=[molecule])
            molcas1.settings.fiedler = True
            molcas1.settings.dmrg_bond_dimension = 3000
            molcas1.settings.dmrg_sweeps = 20
            molcas1.settings.work_dir = self.dummy_project_dir
            molcas1.settings.xyz_file = self.path + "/files/n2.xyz"
            molcas1.project_name = "test_for_fiedler_1_fiedler"
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
            self.assertAlmostEqual(energy_init_no_sym, energy_init_sym, 5)
            self.assertAlmostEqual(energy_no_sym, energy_sym, 5)


if __name__ == "__main__":
    np.set_printoptions(edgeitems=10, linewidth=180)
    unittest.main()
