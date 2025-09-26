# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
import unittest

from scine_autocas.interfaces.qcmaquis_utils.environment import Environment


class TestQCMaquisEnvironment(unittest.TestCase):
    def setUp(self):
        self.environment_bck = os.environ.copy()

    def tearDown(self):
        os.environ = self.environment_bck

    def test_path_error(self):
        # make sure it is not set
        if "QCMAQUIS_BINARY_DIR" in os.environ:
            del os.environ["QCMAQUIS_BINARY_DIR"]
        # pylint: disable=unnecessary-lambda
        self.assertRaises(OSError, lambda: Environment())
        self.assertRaises(OSError, lambda: Environment("random/path"))
        # pylint: enable=unnecessary-lambda

    def test_getter(self):
        try:
            env = Environment()
        except OSError:
            pass
        else:
            self.assertEqual(
                os.environ["QCMAQUIS_BINARY_DIR"] + "/" + env.transform_binary_name, env.get_transform_binary()
            )


if __name__ == "__main__":
    unittest.main()
