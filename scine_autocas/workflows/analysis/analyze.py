# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import glob
import os
from typing import List, Optional, Tuple

from scine_autocas.cas_selection import Autocas
from scine_autocas.interfaces.qcmaquis_utils import QcmaquisUtils
from scine_autocas.workflows.workflow import Workflow


class Analyze(Workflow):
    """Base class for analyze workflows."""

    def __init__(self, autocas: Autocas):
        """Don't nee an interface to analyze stuff.

        Parameters
        ----------
        autocas : Autocas
            an initialized autocas object
        """
        super().__init__(autocas)

    def _initial_orbitals_impl(self):
        """This function is not needed."""

    def _initial_dmrg_impl(self):
        """This function is not needed."""

    def _final_calc_impl(self):
        """This function is not needed."""

    def _get_results_file_and_mps(self, all_files: List[str]) -> Tuple[List[str], List[str]]:
        """Get all result and mps files in autocas project.

        Parameters
        ----------
        all_files : List[str]
            all files from an autocas project
        """
        result_files = []
        mps = []
        for file in all_files:
            if os.path.isfile(file):
                result_files.append(file)
            else:
                mps.append(file)
        return mps, result_files

    def analyze_dmrg_dir(self, state: Optional[int] = None):
        """Analyze a dmrg dir for qcmaquis.

        Parameters
        ----------
        state : int, optional
            if set, only analyze files for given state, else analyze all
        """
        state_str = "0"
        if state:
            state_str = str(state)
        result_file = glob.glob(f"*{state_str}.h5")
        mps, result_files = self._get_results_file_and_mps(result_file)

        print(f"result_files {result_files}")
        print(f"mps {mps}")

        qcmaquis_autocas = QcmaquisUtils()
        qcmaquis_autocas.read_hdf5(result_files[0])
        qcmaquis_autocas.make_diagnostics()

        return (qcmaquis_autocas.s1_entropy,
                qcmaquis_autocas.s2_entropy,
                qcmaquis_autocas.mutual_information)

    def run(self):
        """run"""
        raise NotImplementedError("Overwrite this")
