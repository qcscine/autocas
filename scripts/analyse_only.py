"""Script to analyze autocas projects.

This script provides all required functions to plot an entanglement Diagram
from an autocas project dir or an qcmaquis result state file.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This file is part of SCINE AutoCAS.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details
"""

import argparse
import glob
import os
import re
from typing import Any, Dict, Optional, Tuple, cast

import h5py
import numpy as np

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.plots.entanglement_plot import EntanglementPlot
from scine_autocas.plots.threshold_plot import ThresholdPlot


class Analyzer:
    """Analyze autoCAS Projects.

    Class to provide functions to analyse a single qcmaquis *.results_state.h5
    file or a whole autocas project.

    Attributes
    ----------
    thing_to_analyze : str
        path to dir or qcmaquis file
    save_dir :
    molecule : Molecule
        calculated molecule
    molcas : Molcas
        molcas interface
    autocas : Autocas
        autocas object
    thresh_plot : ThresholdPlot
        plot threshold diagrams
    entang_plot : EntanglementPlot
        plot entanglement diagrams
    excited_state : int
        which excited state to analyze
    """

    def __init__(self, settings: Dict[str, Any]):
        """Initialize from settings.

        Parameters
        ----------
        settings : Dict[str, Any]
            holds all command line arguments
        """
        self.thing_to_analyze: str = settings["QcMaquis_output"]

        tmp_dir = None
        if settings["save_dir"]:
            tmp_dir = settings["save_dir"]
        self.save_dir: Optional[str] = tmp_dir

        tmp_molecule = Molecule(atoms=["H", "H"])
        if settings["molecule"]:
            tmp_molecule = Molecule(settings["molecule"])

        self.molecule: Molecule = tmp_molecule
        self.molcas: Molcas = Molcas([self.molecule])
        self.autocas: Autocas = Autocas(self.molecule)
        self.thresh_plot: ThresholdPlot = ThresholdPlot()
        self.entang_plot: EntanglementPlot = EntanglementPlot()

        tmp_excited_states = 0
        if settings["excited_states"]:
            tmp_excited_states = settings["excited_states"]

        self.excited_state: int = tmp_excited_states

    def analyze(self):
        """Main function.

        Function determines if single file or whole dir
        should be analyzed.

        Raises
        ------
        FileNotFoundError
            if no qcmaquis file can be found
        """
        # single file
        if os.path.isfile(self.thing_to_analyze):
            s1_init, mut_inf_init = self.analyze_single_qcmaquis_file(self.thing_to_analyze)
            if self.save_dir:
                plty = self.entang_plot.plot(s1_init, mut_inf_init)
                entanglement_diagram = self.save_dir + "/entang.pdf"
                plty.savefig(entanglement_diagram)
                plty.show()
                plty_2 = self.thresh_plot.plot(s1_init)
                threshold_diagram = self.save_dir + "/thresh.pdf"
                plty_2.savefig(threshold_diagram)
                plty_2.show()

        # whole project
        elif os.path.isdir(self.thing_to_analyze):
            self.analyze_full_autocas_project(self.thing_to_analyze)

        else:
            raise FileNotFoundError(
                f"The file {args['QcMaquis_output']} cannot be found"
            )

    def analyze_single_qcmaquis_file(self, qcmaquis_output: str):
        """Analyze single qcmaquis result file.

        Parameters
        ----------
        qcmaquis_output : str
            path to qcmaquis result file
        """
        full_path = os.path.abspath(qcmaquis_output)
        # read result state
        print(f"Analyzing file: {full_path}")
        analysis_results = self.molcas.analyze(full_path)
        # we know it is qcmaquis
        analysis_results = cast(Tuple[np.ndarray, np.ndarray, np.ndarray], analysis_results)

        s1_init = analysis_results[0]
        mut_inf_init = analysis_results[2]
        print(s1_init)
        return s1_init, mut_inf_init

    def analyze_full_autocas_project(self, project_path: str):
        """Analyze full autocas project dir.

        Parameters
        ----------
        project_path : str
            path to autocas project to analyze

        Raises
        ------
        ValueError
            if dmrg and dmrg_* dirs are found
        ValueError
            if no dmrg or dmrg_* dirs are found
        FileNotFoundError
            if no dmrg based directory found in path
        """
        print(f"Analyzing autoCAS Project: {os.path.abspath(project_path)}")
        dmrg_dir = os.path.abspath(project_path) + "/dmrg"
        dmrg_large_cas_dir = os.path.abspath(project_path) + "/dmrg_1"
        final_dir = os.path.abspath(project_path) + "/final"

        # too many different calculations in project dir
        if os.path.isdir(dmrg_dir) and os.path.isdir(dmrg_large_cas_dir):
            raise ValueError(
                """Found initial DMRG and initial dmrg_1 dir.
            Please modify your autoCAS Project, that either DMRG or DMRG_1, DMRG_2, ... is present in project dir"""
            )

        # no calculations in project dir
        if (
            not os.path.isdir(dmrg_dir)
            and not os.path.isdir(dmrg_large_cas_dir)
            and not os.path.isdir(final_dir)
        ):
            raise ValueError(
                "No Folder named dmrg or dmrg_1, ... found in project path."
            )

        # standard initial dmrg
        if os.path.isdir(dmrg_dir):
            (
                indices, initial_s1, initial_mut_inf, occupation, atoms, excited_states_flag
            ) = self.analyze_dmrg_dir(dmrg_dir)
            self.molecule = Molecule(atoms=atoms)
            self.autocas = Autocas(self.molecule)
            self.molcas = Molcas([self.molecule])

        elif os.path.isdir(dmrg_large_cas_dir):
            indices, initial_s1, initial_mut_inf = self.analyze_large_cas_dmrg_dir(dmrg_large_cas_dir)

        # elif os.path.isdir(final_dir):
        #    analyse_final_dir(final_dir)

        else:
            raise FileNotFoundError("No calculations found.")

        print()
        print(f"Used Active space: CAS(e,o): ({sum(occupation)}, {len(occupation)})")
        print(f"Orbital indices:     {indices}")
        print(f"Orbital occupations: {occupation}")

        # os.chdir(os.path.expanduser("~/Programs/autoCAS/scripts/"))
        if self.save_dir is not None:
            plty = self.entang_plot.plot(initial_s1, initial_mut_inf, indices)
            entanglement_diagram = self.save_dir + "/entang.pdf"
            plty.savefig(entanglement_diagram)  # type: ignore
            plty.show()  # type: ignore
            plty_2 = self.thresh_plot.plot(initial_s1)
            threshold_diagram = self.save_dir + "/thresh.pdf"
            plty_2.savefig(threshold_diagram)
            plty_2.show()

        occupation, indices = self.autocas.get_active_space(
            occupation=occupation, s1_entropy=initial_s1, mutual_information=initial_mut_inf
        )
        print()
        print(f"Best active space: CAS(e,o): ({sum(occupation)}, {len(occupation)})")
        print(f"Orbital indices:     {indices}")
        print(f"Orbital occupations: {occupation}")

    def analyze_dmrg_dir(self, dir_path: str):
        """Analyze full dmrg dir.

        Parameters
        ----------
        dir_path : str
            path to dmrg dir to analyze

        Returns
        -------
        indices : List[int]
            orbital indices
        s1_init : np.ndarray
            initial single orbital entrpies
        mut_inf_init : np.ndarray
            initial mutual information
        occupation : List[int]
            orbital read_occupations
        atoms : List[str]
            atoms found in result file
        bool
            true if excited states were found
        """

        print(f"Analyzing DMRG calculation in {dir_path}")
        os.chdir(dir_path)

        result_state_string = f"*.results_state.{self.excited_state}.h5"
        qcmaquis_output = sorted(glob.glob(result_state_string))
        print(f"Reading entropies from {qcmaquis_output}")
        molcas_output = glob.glob("*.scf.h5_sel")[0]
        # molcas_output = glob.glob("*.dmrgscf.h5")[0]

        # read result state
        s1_list = []
        mut_inf_list = []

        # for iterate over each state
        if not qcmaquis_output:
            raise FileNotFoundError("No result state found.")
        for output in qcmaquis_output:
            s1_init, mut_inf_init = self.analyze_single_qcmaquis_file(output)
            s1_list.append(s1_init)
            mut_inf_list.append(mut_inf_init)

        # get active space
        typeindices = self.read_type_indices(molcas_output)
        indices = self._get_orbital_numbers_from_type_indices(typeindices)
        # total occupation
        occupation_tmp = self.read_occupations(molcas_output)
        occupation = []
        for index in indices:
            occupation.append(occupation_tmp[index])

        atoms = self.read_atoms(molcas_output)

        if len(qcmaquis_output) > 1:
            print("Found excited states")
            return indices, s1_list, mut_inf_list, occupation, atoms, True
        return indices, s1_init, mut_inf_init, occupation, atoms, False

    def read_type_indices(self, hdf5_file: str):
        """Read typeindices.

        Parameters
        ----------
        hdf5_file : str
            path to hdf5 file to analyze

        Returns
        -------
        type_indices : np.ndarray
            the found type indices
        """
        h5_file = h5py.File(hdf5_file, "r")

        # get type indices
        type_indices = h5_file.get("MO_TYPEINDICES")
        type_indices = np.array(type_indices)
        h5_file.close()
        return type_indices

    def _get_orbital_numbers_from_type_indices(self, typeindices: np.ndarray):
        """Read orbital indices from molcas or qcmaquis output.

        Parameters
        ----------
        typeindices : np.ndarray
            the typeindices from a hdf5 file

        Returns
        -------
        indices : List[int]
            orbital indices found
        """
        indices = []
        for i, item in enumerate(typeindices):
            # Not occupations
            if (
                item.decode("UTF-8") == "2"
                or item.decode("UTF-8") == "1"
            ):
                indices.append(i)
        return indices

    def read_occupations(self, hdf5_file: str):
        """Read occupations.

        Parameters
        ----------
        hdf5_file : str
            path to hdf5 file to analyze

        Returns
        -------
        occupation : List[int]
            occupation for each orbital
        """
        h5_file = h5py.File(hdf5_file, "r")
        # spin mult > 1
        occupation_alpha = h5_file.get("MO_ALPHA_OCCUPATIONS")
        if occupation_alpha is not None:
            occupation_alpha = np.array(occupation_alpha)
            occupation_beta = h5_file.get("MO_BETA_OCCUPATIONS")
            occupation_beta = np.array(occupation_beta)
            occupation = occupation_alpha + occupation_beta
        else:
            occupation = h5_file.get("MO_OCCUPATIONS")
            occupation = np.array(occupation)
        occupation = occupation.tolist()
        for i, occ in enumerate(occupation):
            occupation[i] = int(occ)
        h5_file.close()
        return occupation

    def read_atoms(self, hdf5_file: str):
        """Read atoms.

        Parameters
        ----------
        hdf5_file : str
            path to hdf5 file to analyze

        Returns
        -------
        atoms: List[str]
            element string for each atom found
        """
        h5_file = h5py.File(hdf5_file, "r")
        atoms = h5_file.get("CENTER_LABELS")
        atoms = np.array(atoms)
        atoms = atoms.tolist()
        for i, atom in enumerate(atoms):
            atoms[i] = re.sub(r'\d+', '', str(atom.decode("UTF-8"))).strip()
            if len(atoms[i]) > 1:
                atoms[i] = atoms[i][0] + atoms[i][1].lower()
        h5_file.close()
        return atoms

    def analyze_large_cas_dmrg_dir(self, project_path: str):
        """Analyze large cas dmrg calculation.

        Parameters
        ----------
        project_path : str
            path to autocas project to analyze

        Returns
        -------
        indices : List[int]
            orbital indices found
        s1_entropy : np.ndarray
            single orbital entropy
        mutual_information : np.ndarray
            mutual_information
        """
        print("Analyzing all sub-CAS DMRG calculations")
        excited_states = False
        partial_s1_list = []
        partial_s2_list = []
        partial_mut_inf_list = []
        large_cas_indices = []
        large_cas_occupations = []
        indices = []

        for dmrg_dir in sorted(os.walk(os.path.abspath(self.thing_to_analyze))):
            tmp_dir = dmrg_dir[0].split("/")[-1]
            if tmp_dir[0:5] == "dmrg_":
                indices_tmp, s1_partial, mut_inf_partial, occupation_tmp, atoms, excited_states = self.analyze_dmrg_dir(
                    os.path.abspath(dmrg_dir[0]))
                partial_s1_list.append(s1_partial)
                partial_mut_inf_list.append(mut_inf_partial)
                large_cas_indices.append(indices_tmp)
                large_cas_occupations.append(occupation_tmp)
                # dummy list
                partial_s2_list.append(mut_inf_partial)
                for i in indices_tmp:
                    if i not in indices:
                        indices.append(i)

        indices.sort()
        try:
            molecule = Molecule(atoms=atoms)
        except ValueError:
            molecule = Molecule(atoms=atoms, spin_multiplicity=2)
        # print(molecule.occupation)
        self.autocas = Autocas(molecule)
        self.autocas.cas.n_orbitals = len(indices)

        if excited_states:
            print("====================================")
            final_occupation, final_orbital_indices = self.autocas.get_cas_from_large_cas_excited_states(
                large_cas_indices,
                large_cas_occupations,
                partial_s1_list,
                partial_s2_list,
                partial_mut_inf_list,
                force_cas=False
            )
        else:
            final_occupation, final_orbital_indices = self.autocas.get_cas_from_large_cas(
                large_cas_indices,
                large_cas_occupations,
                partial_s1_list,
                partial_s2_list,
                partial_mut_inf_list,
                force_cas=False
            )
        return indices, self.autocas.diagnostics.s1_entropy, self.autocas.diagnostics.mutual_information


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="SCINE autoCAS - Entanglement Diagram",
        epilog="""Since this is not the main of autocas, a way to make your life easier
        could be an alias in your .bashrc, e.g.
        alias entanglement='conda activate autocas_gui; python3 /path/to/autoCAS/scripts/analyse_only.py'""",
    )
    parser.add_argument(
        "QcMaquis_output",
        metavar="QcMaquis_output",
        type=str,
        help="The output to extracte entanglement diagram from. Can be a directory or a qcmaquis output",
    )
    parser.add_argument(
        "-s",
        "--save_dir",
        type=str,
        help="""Path to store the diagrams. If not specified, it will be saved in current dir. """,
        # Please specify a image compatible file suffix (.pdf, .png, ...), otherwise it will be a png.""",
    )
    parser.add_argument(
        "-e",
        "--excited_states",
        type=int,
        help="""For excited state calculations, define which state to analyse.""",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=str,
        help="""If an active space is required from analysis, specify an xyz file here.""",
    )
    args = vars(parser.parse_args())

    print("********************************************************************************")
    print("*                                                                              *")
    print("*                              AutoCAS - Analysis                              *")
    print("*                                                                              *")
    print("********************************************************************************")
    print("Settings:")
    for key in args:
        print(f"    {key}: {args[key]}")
    print()

    analyzer = Analyzer(args)
    analyzer.analyze()
