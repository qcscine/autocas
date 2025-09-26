# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


import pytest
import unittest
import os
import shutil
from scine_autocas.tests.workflows.consistent_active_space import create_n2_dissociation


class TestConsistentActiveSpaceProtocol(unittest.TestCase):
    def setUp(self):
        self.xyz_files = []
        self.file_base_path = os.path.join("/", *__file__.split("/")[:-1], "files")
        self.project_path = ""

    def tearDown(self):
        for f in self.xyz_files:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
        if self.project_path:
            self.__remove_project_dir(self.project_path)

    def __remove_project_dir(self, directory: str) -> None:
        if os.path.exists(directory):
            shutil.rmtree(directory)

    def test_n2_dissociation_from_xyz(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        try:
            import qcserenity  # type: ignore # noqa: F401
        except ImportError:
            pytest.skip("qcserenity not found")
        from scine_autocas.workflows.consistent_active_space.protocol import run_consistent_active_space_protocol
        from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration

        start_directory = os.getcwd()
        self.xyz_files = create_n2_dissociation()
        config = ConsistentActiveSpaceConfiguration()
        config.xyz_files = self.xyz_files
        config.autocas_indices = [100]
        config.project_name = "test_n2_dissociation_from_xyz"
        self.project_path = os.path.join(os.getcwd(), config.project_name)
        with pytest.raises(ValueError):
            ConsistentActiveSpaceConfiguration.input_sanity_checks(config)
        config.autocas_indices = [len(self.xyz_files) - 1]
        ConsistentActiveSpaceConfiguration.input_sanity_checks(config)

        occupations, indices, energies = run_consistent_active_space_protocol(config)
        reference_energies = [-109.25459589, -108.96763019, -108.94287279]
        for occ, energy, energy_reference in zip(occupations, energies, reference_energies):
            assert occ == [2, 2, 2, 0, 0, 0]
            assert abs(energy - energy_reference) < 1e-6
        self.__remove_project_dir(self.project_path)
        os.chdir(start_directory)

    def test_n2_dissociation_from_hdf5(self):
        try:
            os.environ["MOLCAS"]
        except KeyError:
            pytest.skip("Molcas not found")
        try:
            import qcserenity  # type: ignore # noqa: F401
        except ImportError:
            pytest.skip("qcserenity not found")
        from scine_autocas.workflows.consistent_active_space.protocol import run_consistent_active_space_protocol
        from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration

        start_directory = os.getcwd()
        self.xyz_files = create_n2_dissociation()
        config = ConsistentActiveSpaceConfiguration()
        config.base_load_path = self.file_base_path
        config.load_orbitals = True
        config.system_names = ["system_" + str(i) for i in range(3)]
        config.autocas_indices = [len(config.system_names) - 1]
        config.project_name = "test_n2_dissociation_from_hdf5"
        self.project_path = os.path.join(os.getcwd(), config.project_name)
        ConsistentActiveSpaceConfiguration.input_sanity_checks(config)

        occupations, indices, energies = run_consistent_active_space_protocol(config)
        reference_energies = [-109.25459589, -108.96763019, -108.94287279]
        for occ, energy, energy_reference in zip(occupations, energies, reference_energies):
            assert occ == [2, 2, 2, 0, 0, 0]
            assert abs(energy - energy_reference) < 1e-6
        self.__remove_project_dir(self.project_path)
        os.chdir(start_directory)
