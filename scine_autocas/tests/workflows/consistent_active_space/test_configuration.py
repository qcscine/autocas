# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


import pytest
import unittest
import os

from scine_autocas.tests.workflows.consistent_active_space import create_n2_dissociation


class TestConsistentActiveSpaceConfiguration(unittest.TestCase):
    def setUp(self):
        self.file_base_path = os.path.join("/", *__file__.split("/")[:-1], "files")
        self.xyz_files = []

    def tearDown(self):
        for f in self.xyz_files:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    def test_write_yaml(self):
        from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration

        self.xyz_files = create_n2_dissociation()
        config = ConsistentActiveSpaceConfiguration()
        config.xyz_files = self.xyz_files
        config.autocas_indices = [2]
        config.project_name = "test_write_yaml"
        file_name = config.write_yaml_file("test_write_yaml.yaml")
        assert os.path.isfile(file_name), "YAML file was not written."
        os.remove(file_name)

    def test_create_from_options(self):
        from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration

        class DummyOptions:
            def __init__(self):
                self.load_orbitals = False
                self.cas_method = "CASPT2"
                self.autocas_molecule_indices = "0"
                self.basis_set = "def2-SVP"
                self.large_active_space = False
                self.exclude_core = False
                self.always_include_unmapables = False

        options = DummyOptions()
        xyz_files = [os.path.join(self.file_base_path, file_name) for file_name in ["n2_0.xyz", "n2_1.xyz"]]
        ConsistentActiveSpaceConfiguration.from_options(options, xyz_files)

    def test_create_from_yaml(self):
        from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration

        self.xyz_files = create_n2_dissociation()

        valid_path = os.path.join(self.file_base_path, "valid_input.yaml")
        assert os.path.isfile(valid_path), "Configuration file does not exist."
        ConsistentActiveSpaceConfiguration.from_file(valid_path)
        broken_path = os.path.join(self.file_base_path, "does_not_exist.yaml")
        with pytest.raises(FileNotFoundError):
            ConsistentActiveSpaceConfiguration.from_file(broken_path)

        broken_input_files = os.path.join(self.file_base_path, "broken_input_1.yaml")
        with pytest.raises(ValueError):
            ConsistentActiveSpaceConfiguration.from_file(broken_input_files)
        broken_input_files = os.path.join(self.file_base_path, "broken_input_2.yaml")
        with pytest.raises(ValueError):
            ConsistentActiveSpaceConfiguration.from_file(broken_input_files)
