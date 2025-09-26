# -*- coding: utf-8 -*-
__copyright__ = """This file is part of SCINE AutoCAS.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details
"""

import os

from scine_autocas import Autocas, Molecule
from scine_autocas.interfaces import Molcas
from scine_autocas.io import FileHandler
from scine_autocas.plots import EntanglementPlot
from scine_autocas.workflows import ClassicWorkflow


def standard_autocas_procedure(path: str, xyz_file: str):
    """Do ground state cas"""
    # create a molecule
    molecule = Molecule(xyz_file)

    # initialize autoCAS and Molcas interface
    autocas = Autocas(molecule)
    molcas = Molcas(molecule)

    # setup interface
    molcas.project_name = "test_run"
    molcas.environment.molcas_scratch_dir = path + "/scratch"

    # cas and hyphen do not matter for method names
    molcas.set_cas_method("dmrgci")

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 250
    molcas.settings.dmrg_sweeps = 5

    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()

    # no input means HF calculation
    molcas.calculate()

    # do cas calculation
    cas_results = molcas.calculate(occ_initial, index_initial)

    # energy = cas_results[0]
    s1_entropy = cas_results[1]
    # s2_entropy = cas_results[2]
    mut_inf = cas_results[3]

    # plot entanglement diagram
    plot = EntanglementPlot()
    plt = plot.plot(s1_entropy, mut_inf)  # type: ignore
    plt.savefig("entang.pdf")  # type: ignore

    # make active space based on single orbital entropies
    cas_occ, cas_index = autocas.get_active_space(
        s1_entropy   # type: ignore
    )

    # cas and hyphen do not matter for method names
    molcas.set_cas_method("casscf")
    molcas.set_post_cas_method("caspt2")

    # Do a calculation with this CAS
    final_energy, final_s1, final_s2, final_mut_inf = molcas.calculate(cas_occ, cas_index)

    # use results
    n_electrons = sum(cas_occ)
    n_orbitals = len(cas_occ)
    print(f"final energy:      {final_energy}")
    print(f"final CAS(e, o):  ({n_electrons}, {n_orbitals})")
    print(f"final cas indices: {cas_index}")
    print(f"final occupation:  {cas_occ}")
    print(f"final s1:          {final_s1}")
    return cas_occ, cas_index, final_energy


if __name__ == "__main__":
    path_to_this_file = os.path.dirname(os.path.abspath(__file__))
    xyz = path_to_this_file + "/../scine_autocas/tests/interfaces/files/n2.xyz"
    occupation, orbitals, energy = standard_autocas_procedure(path_to_this_file, xyz)
    print(f"n2 \ncas: {occupation} \norbs: {orbitals} \nenergy: {energy}")

    # to verify workflow
    print("\n\nverification run")

    # change project name, to prevent calculation in same project
    FileHandler.project_dir = path_to_this_file
    FileHandler.DirectoryNames.project_name = "verification_run"
    FileHandler.setup_project()

    # setup autocas classes and interface
    test_molecule = Molecule(xyz)
    test_autocas = Autocas(test_molecule)
    test_interface = Molcas(test_molecule)

    # change project name and scratch dir
    test_interface.project_name = "verify_example"
    test_interface.environment.molcas_scratch_dir = path_to_this_file + "/scratch_verify"

    # set methods for production run
    test_interface.set_cas_method("casscf")
    test_interface.set_post_cas_method("caspt2")
    test_interface.settings.dmrg_bond_dimension = 250
    test_interface.settings.dmrg_sweeps = 5

    # run workflow
    workflow = ClassicWorkflow(test_autocas, test_interface)
    workflow.run()

    assert abs(workflow.results["final_energy"][0] - energy[0]) < 1e-9
    assert workflow.results["final_occupation"] == occupation
    assert workflow.results["final_orbital_indices"] == orbitals
