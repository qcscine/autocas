# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
import shutil
import unittest

from scine_autocas.io.file_handler import FileHandler


class TestFileHandler(unittest.TestCase):
    # dont put in setup, it may change
    this_dir = os.getcwd()

    def setUp(self):
        os.chdir(TestFileHandler.this_dir)
        self.this_dir = os.getcwd()

        self.mock_dir = os.getcwd() + "/mock"
        try:
            os.mkdir(self.mock_dir)
        except FileExistsError:
            pass

    def tearDown(self):
        try:
            shutil.rmtree(self.mock_dir)
        except FileNotFoundError:
            pass
        pass

    def test_setters(self):
        FileHandler.set_current_dir()
        self.assertEqual(FileHandler.current_dir, self.this_dir)
        FileHandler.set_current_dir("abc")
        self.assertEqual(FileHandler.current_dir, self.this_dir + "/abc")
        FileHandler.set_current_dir("abc/ghd")
        self.assertEqual(FileHandler.current_dir, self.this_dir + "/abc/ghd")
        FileHandler.set_current_dir("abc/ghd/")
        self.assertEqual(FileHandler.current_dir, self.this_dir + "/abc/ghd")
        FileHandler.set_current_dir("/abc/ghd/")
        self.assertEqual(FileHandler.current_dir, "/abc/ghd")

        FileHandler.set_project_dir("abc")
        self.assertEqual(FileHandler.project_dir, self.this_dir + "/abc")
        FileHandler.set_project_dir("abc/ghd")
        self.assertEqual(FileHandler.project_dir, self.this_dir + "/abc/ghd")
        FileHandler.set_project_dir("abc/ghd/")
        self.assertEqual(FileHandler.project_dir, self.this_dir + "/abc/ghd")
        FileHandler.set_project_dir("/abc/ghd/")
        self.assertEqual(FileHandler.project_dir, "/abc/ghd")

        FileHandler.set_project_dir()
        self.assertEqual(FileHandler.project_dir, self.this_dir)

    def test_make_custom_project_dir(self):
        FileHandler.set_project_dir(self.mock_dir)
        FileHandler.make_project_dir()
        self.assertTrue(os.path.exists(FileHandler.project_dir))
        self.assertTrue(os.path.exists(FileHandler.get_project_path()))

    def test_makers(self):
        os.chdir(self.mock_dir)
        FileHandler.set_project_dir()
        FileHandler.make_project_dir()
        self.assertTrue(os.path.exists(FileHandler.project_dir))
        self.assertTrue(os.path.exists(FileHandler.get_project_path()))
        FileHandler.make_initial_orbital_dir()
        self.assertTrue(os.path.exists(FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.initial_orbs}"))
        FileHandler.make_initial_dmrg_dir()
        self.assertTrue(os.path.exists(FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.initial_dmrg}"))
        FileHandler.make_final_calc_dir()
        self.assertTrue(os.path.exists(FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.final_calc}"))

    def test_checkers(self):
        os.chdir(self.mock_dir)
        FileHandler.set_project_dir()
        self.assertFalse(FileHandler.check_project_dir_exists())
        FileHandler.make_project_dir()
        self.assertTrue(FileHandler.check_project_dir_exists())

        self.assertFalse(FileHandler.check_dir_exists(FileHandler.DirectoryNames.initial_orbs))
        FileHandler.make_initial_orbital_dir()
        self.assertTrue(FileHandler.check_dir_exists(FileHandler.DirectoryNames.initial_orbs))

    def test_changers(self):
        os.chdir(self.mock_dir)
        FileHandler.set_project_dir()
        FileHandler.make_project_dir()
        FileHandler.make_initial_orbital_dir()
        FileHandler.ch_to_initial_orbital_dir()
        self.assertEqual(os.getcwd(), FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.initial_orbs}")
        self.assertEqual(os.getcwd(), FileHandler.current_dir)

        FileHandler.make_initial_dmrg_dir()
        FileHandler.ch_to_initial_dmrg_dir()
        self.assertEqual(os.getcwd(), FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.initial_dmrg}")
        self.assertEqual(os.getcwd(), FileHandler.current_dir)

        FileHandler.make_final_calc_dir()
        self.assertEqual(os.getcwd(), FileHandler.get_project_path() + f"/{FileHandler.DirectoryNames.final_calc}")
        self.assertEqual(os.getcwd(), FileHandler.current_dir)


if __name__ == "__main__":
    unittest.main()
