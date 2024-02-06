"""Basic example script.

This script is a modified version from the ground state cas calculation
in scine_autocas.main_functions
"""
# -*- coding: utf-8 -*-
__copyright__ = """This file is part of SCINE AutoCAS.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details
"""

import os

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.main_functions import MainFunctions
from scine_autocas.plots.entanglement_plot import EntanglementPlot


def standard_autocas_procedure(path: str, xyz_file: str):
    """Do ground state cas"""
    # create a molecule
    molecule = Molecule(xyz_file)

    # initialize autoCAS and Molcas interface
    autocas = Autocas(molecule)
    molcas = Molcas([molecule])

    # setup interface
    molcas.project_name = "example"
    molcas.settings.work_dir = path + "/../test/example"
    molcas.environment.molcas_scratch_dir = path + "/../test/scratch"
    molcas.settings.xyz_file = xyz_file

    # cas and hyphen do not matter for method names
    molcas.settings.method = "DMRGCI"

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
    plt.savefig(molcas.settings.work_dir + "/entang.pdf")  # type: ignore

    # make active space based on single orbital entropies
    cas_occ, cas_index = autocas.get_active_space(
        occ_initial, s1_entropy   # type: ignore
    )

    # cas and hyphen do not matter for method names
    molcas.settings.method = "dmrg-scf"

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 2000
    molcas.settings.dmrg_sweeps = 20

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
    print(f"final s2: \n{final_s2}")
    print(f"final mut_inf: \n{final_mut_inf}")
    return cas_occ, cas_index, final_energy


if __name__ == "__main__":
    path_to_this_file = os.path.dirname(os.path.abspath(__file__))
    xyz = path_to_this_file + "/../scine_autocas/tests/files/n2.xyz"
    occupation, orbitals, energy = standard_autocas_procedure(path_to_this_file, xyz)
    print(f"n2 \ncas: {occupation} \norbs: {orbitals} \nenergy: {energy}")

    print("\n\n")
    print("only verification run")
    # to verify that everything is still valid
    main_functions = MainFunctions()
    test_molecule = Molecule(xyz)
    test_autocas = Autocas(test_molecule)
    test_interface = Molcas([test_molecule])

    test_interface.project_name = "verify_example"
    test_interface.settings.work_dir = path_to_this_file + "/../test/verify_example"
    test_interface.settings.xyz_file = xyz
    test_interface.environment.molcas_scratch_dir = path_to_this_file + "/../test/scratch_tmp"

    test_interface.settings.method = "dmrg-ci"
    test_interface.settings.dmrg_bond_dimension = 250
    test_interface.settings.dmrg_sweeps = 5
    test_cas, test_indices = main_functions.conventional(test_autocas, test_interface)

    test_interface.settings.method = "dmrg-scf"
    test_interface.settings.dmrg_bond_dimension = 2000
    test_interface.settings.dmrg_sweeps = 20
    test_results = test_interface.calculate(test_cas, test_indices)

    assert abs(test_results[0] - energy) < 1e-9
    assert test_cas == occupation
    assert test_indices == orbitals
