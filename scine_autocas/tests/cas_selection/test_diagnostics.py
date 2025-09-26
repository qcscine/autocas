# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest

import numpy as np

from scine_autocas.cas_selection.diagnostics import Diagnostics


class TestDiagnostics(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.s1 = np.array([0.1, 0.2, 0.12, 0.18, 0.65])

    def test_single_reference(self):
        diagnostics = Diagnostics()
        self.assertFalse(diagnostics.is_single_reference(self.s1))
        diagnostics.single_reference_threshold = 0.7
        self.assertTrue(diagnostics.is_single_reference(self.s1))

    def test_zs1(self):
        diagnostics = Diagnostics()
        zs1 = diagnostics.zs1(self.s1)
        self.assertTrue((zs1 - 0.18033688) < 1e-7)


if __name__ == "__main__":
    unittest.main()
