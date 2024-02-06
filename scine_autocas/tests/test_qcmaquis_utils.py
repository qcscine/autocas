# pylint: disable=C0114, C0115, C0116
# pylint: disable=protected-access
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np

from scine_autocas.interfaces.qcmaquis import Qcmaquis
from scine_autocas.interfaces.qcmaquis.qcmaquis_hdf5_utils import Hdf5Converter
from scine_autocas.interfaces.qcmaquis.qcmaquis_orbital_rdm_builder import OrbitalRDMBuilder


class TestQcMquisClasses(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))

    def test_hdf5_converter(self):
        hdf5_converter = Hdf5Converter()
        qcmaquis_result_file = self.path + "/files/n2.results_state.0.h5"
        hdf5_converter.read_hdf5(qcmaquis_result_file)
        self.assertEqual(hdf5_converter.L, 8)
        self.assertAlmostEqual(hdf5_converter.energy, -108.7499817)
        self.assertEqual(hdf5_converter.symmetry, "su2u1pg")
        # fmt: off
        self.assertTrue(
            np.array_equal(hdf5_converter.orbital_order, np.array([1, 2, 3, 4, 5, 6, 7, 8]))
        )

        self.assertRaises(
            AssertionError,
            lambda: hdf5_converter._Hdf5Converter__mat_measurement([(1, 1), (2, 2)], np.array([0, 1, 2]))
        )
        self.assertRaises(
            AssertionError,
            lambda: hdf5_converter._Hdf5Converter__mat_merge_transpose([(1, 1), (2, 2)], np.array([0, 1, 2]),
                                                                       [(1, 1), (2, 2)], np.array([0, 1, 2]))
        )
        # fmt: on
        new_hdf5_converter = Hdf5Converter()
        new_hdf5_converter.result_file = qcmaquis_result_file
        new_hdf5_converter.read_hdf5()
        self.assertEqual(new_hdf5_converter.L, 8)
        self.assertAlmostEqual(new_hdf5_converter.energy, -108.7499817)
        self.assertEqual(new_hdf5_converter.symmetry, "su2u1pg")

    def test_orbital_rdm_builder(self):
        hdf5_converter = Hdf5Converter()
        ordm_builder = OrbitalRDMBuilder()
        self.assertRaises(AttributeError, lambda: ordm_builder.make_one_ordm(hdf5_converter))
        self.assertRaises(AttributeError, lambda: ordm_builder.make_two_ordm(hdf5_converter))
        qcmaquis_result_file = self.path + "/files/n2.results_state.0.h5"
        hdf5_converter.read_hdf5(qcmaquis_result_file)
        one_ordm = ordm_builder.make_one_ordm(hdf5_converter)
        ordm_builder.make_two_ordm(hdf5_converter)
        # fmt: off
        test_one_ordm = np.array(
            [
                [2.53220e-03, 2.53220e-03, 2.52018e-04, 9.94683e-01],
                [2.40389e-03, 2.40389e-03, 1.35218e-03, 9.93840e-01],
                [7.87185e-02, 7.87185e-02, 1.24039e-01, 7.18523e-01],
                [1.24453e-01, 1.24453e-01, 2.49871e-01, 5.01222e-01],
                [1.24453e-01, 1.24453e-01, 2.49871e-01, 5.01221e-01],
                [1.24453e-01, 1.24453e-01, 5.00725e-01, 2.50368e-01],
                [1.24453e-01, 1.24453e-01, 5.00725e-01, 2.50368e-01],
                [8.04333e-02, 8.04333e-02, 7.11261e-01, 1.27871e-01],
            ]
        )
        # fmt: on
        self.assertTrue(np.allclose(one_ordm, test_one_ordm))

    def test_qcmaquis(self):
        qcmaquis_result_file = self.path + "/files/n2.results_state.0.h5"
        qc_maquis = Qcmaquis()
        # this lambda is necessary to check for the AttributeError
        # pylint: disable=unnecessary-lambda
        self.assertRaises(AttributeError, lambda: qc_maquis.make_diagnostics())
        # pylint: enable=unnecessary-lambda
        self.assertRaises(
            AttributeError,
            lambda: qc_maquis.make_s1(
                qc_maquis.orbital_rdm_builder.make_one_ordm(qc_maquis.hdf5_converter)
            ),
        )
        self.assertRaises(
            AttributeError,
            lambda: qc_maquis.make_s2(
                qc_maquis.orbital_rdm_builder.make_two_ordm(qc_maquis.hdf5_converter)
            ),
        )
        qc_maquis.read_hdf5(qcmaquis_result_file)
        qc_maquis.make_diagnostics()
        qc_maquis.make_s1(qc_maquis.orbital_rdm_builder.make_one_ordm(qc_maquis.hdf5_converter))
        qc_maquis.make_s2(qc_maquis.orbital_rdm_builder.make_two_ordm(qc_maquis.hdf5_converter))
        qc_maquis.make_mutual_information()
        qc_maquis.hdf5_converter.L = None
        self.assertRaises(
            AttributeError, lambda: qc_maquis.make_s1(qc_maquis.orbital_rdm_builder.one_ordm)
        )
        self.assertRaises(
            AttributeError, lambda: qc_maquis.make_s2(qc_maquis.orbital_rdm_builder.two_ordm)
        )
        self.assertRaises(
            AttributeError,
            lambda: qc_maquis.make_s1(qc_maquis.orbital_rdm_builder.make_one_ordm(qc_maquis.hdf5_converter))
        )
        qc_maquis.s1_entropy = np.array([])
        # pylint: disable=unnecessary-lambda
        self.assertRaises(AttributeError, lambda: qc_maquis.make_mutual_information())
        self.assertRaises(AttributeError, lambda: qc_maquis.make_diagnostics())
        # pylint: enable=unnecessary-lambda


if __name__ == "__main__":
    unittest.main()
