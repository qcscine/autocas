# -*- coding: utf-8 -*-
# pylint: disable=C0114, C0115, C0116
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule


class TestAutocasClasses(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))

    def test_double_shell(self):
        molecule = Molecule(xyz_file=(self.path + "/files/fe.xyz"), charge=2, spin_multiplicity=1)
        autocas = Autocas(molecule)
        fe_occ, fe_index = autocas.make_initial_active_space()
        # fmt: off
        self.assertEqual(fe_occ, [2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(fe_index, [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
        # fmt: on
        #
        molecule = Molecule(
            xyz_file=(self.path + "/files/fe.xyz"),
            charge=2,
            spin_multiplicity=1,
            double_d_shell=False,
        )
        autocas = Autocas(molecule)
        fe_occ, fe_index = autocas.make_initial_active_space()
        # fmt: off
        self.assertEqual(fe_occ, [2, 2, 2, 0, 0, 0, 0, 0, 0])
        self.assertEqual(fe_index, [9, 10, 11, 12, 13, 14, 15, 16, 17])
        # fmt: on

    def test_spin_multiplicity_and_charge(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"), charge=1, spin_multiplicity=2)
        autocas = Autocas(molecule)
        n2_occ, n2_index = autocas.make_initial_active_space()

        # fmt: off
        self.assertEqual(autocas.molecule.core_orbitals, 2)
        self.assertEqual(autocas.molecule.valence_orbitals, 8)
        self.assertEqual(autocas.molecule.electrons, 13)
        self.assertEqual(autocas.molecule.occupation, [2, 2, 2, 2, 2, 2, 1, 0, 0, 0])
        self.assertEqual(autocas.cas.n_orbitals, 8)
        self.assertEqual(autocas.cas.n_electrons, 9)
        self.assertEqual(autocas.cas.orbital_indices, [2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(autocas.cas.occupation, [2, 2, 2, 2, 1, 0, 0, 0])
        self.assertEqual(n2_occ, [2, 2, 2, 2, 1, 0, 0, 0])
        self.assertEqual(n2_index, [2, 3, 4, 5, 6, 7, 8, 9])
        # fmt: on

    def test_spin_multiplicity(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"), charge=-1, spin_multiplicity=4)
        autocas = Autocas(molecule)
        n2_occ, n2_index = autocas.make_initial_active_space()

        # fmt: off
        self.assertEqual(autocas.molecule.core_orbitals, 2)
        self.assertEqual(autocas.molecule.valence_orbitals, 8)
        self.assertEqual(autocas.molecule.electrons, 15)
        self.assertEqual(autocas.molecule.occupation, [2, 2, 2, 2, 2, 2, 1, 1, 1, 0])
        self.assertEqual(autocas.cas.n_orbitals, 8)
        self.assertEqual(autocas.cas.n_electrons, 11)
        self.assertEqual(autocas.cas.orbital_indices, [2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(autocas.cas.occupation, [2, 2, 2, 2, 1, 1, 1, 0])
        self.assertEqual(n2_occ, [2, 2, 2, 2, 1, 1, 1, 0])
        self.assertEqual(n2_index, [2, 3, 4, 5, 6, 7, 8, 9])
        # fmt: on

    def test_valence_space(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule)
        n2_occ, n2_index = autocas.make_initial_active_space()

        # fmt: off
        self.assertEqual(autocas.molecule.core_orbitals, 2)
        self.assertEqual(autocas.molecule.valence_orbitals, 8)
        self.assertEqual(autocas.molecule.electrons, 14)
        self.assertEqual(autocas.molecule.occupation, [2, 2, 2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(autocas.cas.n_orbitals, 8)
        self.assertEqual(autocas.cas.n_electrons, 10)
        self.assertEqual(autocas.cas.orbital_indices, [2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(autocas.cas.occupation, [2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(n2_occ, [2, 2, 2, 2, 2, 0, 0, 0])
        self.assertEqual(n2_index, [2, 3, 4, 5, 6, 7, 8, 9])
        # fmt: on

        molecule = Molecule(xyz_file=(self.path + "/files/h2o.xyz"))
        autocas = Autocas(molecule)
        h2o_occ, h2o_index = autocas.make_initial_active_space()

        # fmt: off
        self.assertEqual(autocas.molecule.core_orbitals, 1)
        self.assertEqual(autocas.molecule.valence_orbitals, 6)
        self.assertEqual(autocas.molecule.electrons, 10)
        self.assertEqual(autocas.molecule.occupation, [2, 2, 2, 2, 2, 0, 0])
        self.assertEqual(autocas.cas.n_orbitals, 6)
        self.assertEqual(autocas.cas.n_electrons, 8)
        self.assertEqual(autocas.cas.occupation, [2, 2, 2, 2, 0, 0])
        self.assertEqual(autocas.cas.orbital_indices, [1, 2, 3, 4, 5, 6])
        self.assertEqual(h2o_occ, [2, 2, 2, 2, 0, 0])
        self.assertEqual(h2o_index, [1, 2, 3, 4, 5, 6])
        # fmt: on

    def test_autocas_space(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule)
        n2_occ, _ = autocas.make_initial_active_space()
        # fmt: off
        s1_entropy = np.array(
            [0.04520418, 0.03036112, 1.22331744, 1.22357522, 0.87789977, 1.20122905, 0.89020811, 1.20159583]
        )
        # fmt: on
        n2_cas_occ, n2_cas_index = autocas.get_active_space(n2_occ, s1_entropy, force_cas=True)

        self.assertEqual(n2_cas_occ, [2, 2, 2, 0, 0, 0])
        self.assertEqual(n2_cas_index, [4, 5, 6, 7, 8, 9])

    def test_large_active_spaces(self):
        molecule = Molecule(xyz_file=(self.path + "/files/3n2.xyz"))
        autocas = Autocas(molecule)
        # to allow predictions in the test. Otherwise it is None, so that it uses data from /dev/urandom or clock
        autocas.large_spaces.seed = 1
        # TODO: make test pass for max entanglement
        autocas.large_spaces.average_entanglement = True
        autocas.large_spaces.max_orbitals = 8
        _ = autocas.make_initial_active_space()
        large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()

        # fmt: off
        self.assertEqual(large_cas_occupations[0],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[1],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[2],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[3],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[4],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[5],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[6],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[7],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[8],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[9],  [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[10], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[11], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[12], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[13], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[14], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[15], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[16], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[17], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[18], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[19], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[20], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[21], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[22], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[23], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[24], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[25], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[26], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[27], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[28], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[29], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[30], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[31], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[32], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[33], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[34], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[35], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[36], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[37], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[38], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_occupations[39], [2, 2, 2, 2, 0, 0, 0, 0])
        self.assertEqual(large_cas_indices[0], [12, 13, 14, 15, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[1], [12, 13, 14, 15, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[2], [12, 13, 14, 15, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[3], [12, 13, 14, 15, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[4], [12, 13, 14, 15, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[5], [16, 17, 18, 19, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[6], [16, 17, 18, 19, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[7], [16, 17, 18, 19, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[8], [16, 17, 18, 19, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[9], [16, 17, 18, 19, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[10], [20, 21, 22, 23, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[11], [20, 21, 22, 23, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[12], [20, 21, 22, 23, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[13], [20, 21, 22, 23, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[14], [20, 21, 22, 23, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[15], [24, 25, 26, 27, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[16], [24, 25, 26, 27, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[17], [24, 25, 26, 27, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[18], [24, 25, 26, 27, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[19], [24, 25, 26, 27, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[20], [28, 29, 30, 31, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[21], [28, 29, 30, 31, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[22], [28, 29, 30, 31, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[23], [28, 29, 30, 31, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[24], [28, 29, 30, 31, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[25], [32, 33, 34, 35, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[26], [32, 33, 34, 35, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[27], [32, 33, 34, 35, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[28], [32, 33, 34, 35, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[29], [32, 33, 34, 35, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[30], [36, 37, 38, 39, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[31], [36, 37, 38, 39, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[32], [36, 37, 38, 39, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[33], [36, 37, 38, 39, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[34], [36, 37, 38, 39, 47, 53, 58, 59])
        self.assertEqual(large_cas_indices[35], [17, 23, 40, 41, 42, 43, 44, 45])
        self.assertEqual(large_cas_indices[36], [17, 23, 40, 41, 46, 47, 48, 49])
        self.assertEqual(large_cas_indices[37], [17, 23, 40, 41, 50, 51, 52, 53])
        self.assertEqual(large_cas_indices[38], [17, 23, 40, 41, 54, 55, 56, 57])
        self.assertEqual(large_cas_indices[39], [17, 23, 40, 41, 47, 53, 58, 59])
        # fmt: on

    def test_collect_entropies(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        molecule.core_orbitals = 2
        molecule.electrons = 8
        autocas = Autocas(molecule)
        autocas.cas.n_orbitals = 6
        # fmt : off
        s1_list = [
            np.array([1.0, 0.5, 0.5, 1.0]),
            np.array([0.0, 0.0, 2.0, 1.0]),
            np.array([0.5, 0.25, 0.5, 0.0]),
            np.array([0.0, 1.0, 1.0, 1.0]),
        ]

        index_list = [
            [2, 3, 4, 5],
            [3, 4, 6, 7],
            [2, 3, 4, 6],
            [4, 5, 6, 7]
        ]

        occupation_list = [
            [2, 2, 0, 0],
            [2, 0, 0, 0],
            [2, 2, 0, 0],
            [0, 0, 0, 0]
        ]

        s2_list = [
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
        ]
        mutual_information_list = s2_list
        # fmt: on

        autocas.large_spaces.average_entanglement = True
        occupation, s1_entropy, s2_entropy, _ = autocas.collect_entropies(
            index_list, occupation_list, s1_list, s2_list, mutual_information_list
        )
        _ = s2_entropy  # get rid of warning

        self.assertListEqual(list(map(int, list(occupation))), [2, 2, 0, 0, 0, 0])
        self.assertEqual(s1_entropy.tolist(), np.array([0.75, 0.25, 0.25, 1.0, 1.0, 1.0]).tolist())

    def test_get_cas_from_large_cas(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        molecule.core_orbitals = 2
        molecule.electrons = 8
        autocas = Autocas(molecule)
        autocas.cas.n_orbitals = 6

        # fmt: off
        s1_list = [
            np.array([1.0, 0.5, 0.5, 1.0]),
            np.array([0.0, 0.0, 2.0, 0.0]),
            np.array([0.5, 0.25, 0.5, 0.0]),
            np.array([0.0, 1.0, 1.0, 0.0]),
        ]

        index_list = [
            [2, 3, 4, 5],
            [3, 4, 6, 7],
            [2, 3, 4, 6],
            [4, 5, 6, 7]
        ]

        occupation_list = [
            [2, 2, 0, 0],
            [2, 0, 0, 0],
            [2, 2, 0, 0],
            [0, 0, 0, 0]
        ]

        s2_list = [
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.5, 0.5, 1.0],
                    [0.0, 0.0, 2.0, 1.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.0, 1.0, 1.0, 1.0],
                ]
            ),
        ]
        mutual_information_list = s2_list
        # fmt: on

        occupation, indices = autocas.get_cas_from_large_cas(
            index_list, occupation_list, s1_list, s2_list, mutual_information_list
        )

        self.assertListEqual(list(map(int, list(occupation))), [2, 2, 0, 0, 0])
        self.assertListEqual(indices, [2, 3, 4, 5, 6])

    def test_get_cas_from_large_cas_excited_states(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        molecule.core_orbitals = 2
        molecule.electrons = 8
        autocas = Autocas(molecule)
        autocas.cas.n_orbitals = 6
        autocas.large_spaces.max_orbitals = 4

        # fmt: off
        s1_list = [
            [
                np.array([1.0, 0.5, 0.5, 1.0]),
                np.array([0.5, 0.5, 0.0, 0.0]),
            ],
            [
                np.array([0.0, 0.5, 0.0, 0.0]),
                np.array([0.0, 0.0, 2.0, 0.0]),
            ],
            [
                np.array([1.0, 0.5, 0.5, 0.0]),
                np.array([0.0, 0.0, 2.0, 2.0]),
            ],
            [
                np.array([0.0, 0.5, 2.0, 0.0]),
                np.array([0.0, 0.0, 2.0, 0.0]),
            ],
        ]

        index_list = [
            [2, 3, 4, 5],
            [3, 4, 6, 7],
            [2, 3, 4, 6],
            [4, 5, 6, 7]
        ]

        occupation_list = [
            [2, 2, 0, 0],
            [2, 0, 0, 0],
            [2, 2, 0, 0],
            [0, 0, 0, 0]
        ]

        s2_list = [
            [
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
            ],
            [
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
            ],
            [
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
            ],
            [
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
                np.array(
                    [
                        [1.0, 0.5, 0.5, 1.0],
                        [0.0, 0.0, 2.0, 1.0],
                        [0.5, 0.0, 0.5, 0.0],
                        [0.0, 1.0, 1.0, 1.0],
                    ]
                ),
            ],
        ]
        mutual_information_list = s2_list
        # fmt: off

        occupation, indices = autocas.get_cas_from_large_cas_excited_states(
            index_list, occupation_list, s1_list, s2_list, mutual_information_list
        )

        self.assertListEqual(indices, [2, 3, 4, 5, 6])
        self.assertListEqual(list(map(int, list(occupation))), [2, 2, 0, 0, 0])

    def test_get_cas_from_excited_states(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        molecule.core_orbitals = 0
        autocas = Autocas(molecule)

        # fmt: off
        s1_list = [
            np.array([1.0, 0.5, 0.5, 1.0, 0.0, 1.0, 0.0, 0.4]),
            np.array([0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 4.0, 0.0]),
            np.array([0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 4.0, 0.0]),
        ]
        occupation = [2, 2, 2, 2, 2, 0, 0, 0]
        # fmt: on

        occupation, indices = autocas.get_cas_from_excited_states(occupation, s1_list)

        self.assertListEqual(list(map(int, list(occupation))), [2, 2, 2, 2, 0, 0, 0])
        self.assertListEqual(sorted(indices), [2, 3, 4, 5, 7, 8, 9])

    def test_exclude_orbitals(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        molecule.core_orbitals = 0
        autocas = Autocas(molecule)

        # fmt: off
        s1_entropy = np.array(
            [
                1.4,
                0.91 * 1.4,
                0.83 * 1.4,
                0.75 * 1.4,
                0.67 * 1.4,
                0.59 * 1.4,
                0.51 * 1.4,
                0.015 * 1.4,
                0.43 * 1.4,
                0.35 * 1.4,
                0.27 * 1.4,
                0.19 * 1.4,
                0.015 * 1.4,
                0.11 * 1.4,
            ]
        )
        occupation = [2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # fmt: on

        occupation, indices = autocas.get_active_space(occupation, s1_entropy)

        # fmt: off
        self.assertListEqual(list(map(int, list(occupation))), [2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0])
        self.assertListEqual(sorted(indices), [2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15])
        # fmt: on

    def test_custom_initial_space(self):
        molecule = Molecule(xyz_file=(self.path + "/files/n2.xyz"))
        autocas = Autocas(molecule)
        _ = autocas.make_initial_active_space()

        # fmt: off
        s1_entropy = np.array(
            [0.03036112, 1.22331744, 1.22357522, 0.87789977, 1.20122905, 0.00089020, 1.20159583]
        )
        # fmt: on

        n2_cas_occ, n2_cas_index = autocas.get_active_space(
            [2, 2, 2, 2, 0, 0, 0], s1_entropy, force_cas=True
        )

        self.assertEqual(n2_cas_index, [4, 5, 6, 7, 9])
        self.assertEqual(n2_cas_occ, [2, 2, 2, 0, 0])


if __name__ == "__main__":
    unittest.main()
