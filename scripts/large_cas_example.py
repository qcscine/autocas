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
from scine_autocas.workflows import LargeCasWorkflow


def large_cas_autocas_example(path: str, xyz_file: str):
    """Do Large active space."""
    # create a molecule
    molecule = Molecule(xyz_file)

    # initialize autoCAS and Molcas interface
    autocas = Autocas(molecule)
    molcas = Molcas(molecule)

    # setup interface
    molcas.project_name = "large_cas_example"
    molcas.environment.molcas_scratch_dir = path + "/../test/scratch"

    # cas and hyphen do not matter for method names
    molcas.set_cas_method("dmrgci")

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 250
    molcas.settings.dmrg_sweeps = 5

    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()

    n_electrons = sum(occ_initial)
    n_orbitals = len(occ_initial)
    print(f"initial CAS(e, o):  ({n_electrons}, {n_orbitals})")
    print(f"initial cas indices: {index_initial}")
    print(f"initial occupation:  {occ_initial}")
    print()

    # no input means HF calculation
    # interface automatically chooses unrestricted HF for spin spin_multiplicity > 1
    molcas.calculate()

    # allow interface to use large cas protocol
    molcas.settings.large_cas = True

    # set maximum number of orbitals in a sub cas
    autocas.large_spaces.max_orbitals = 30

    # get a list of occupations and indices for each sub cas
    large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()

    # iterate over subspaces and evaluate them
    partial_s1_list = []
    partial_s2_list = []
    partial_mut_inf_list = []
    print(f"autocas will perform: {len(large_cas_indices)} dmrg calculations")
    print(f"with {autocas.large_spaces.max_orbitals} in each active space")

    for i, partial_occupation in enumerate(large_cas_occupations):
        # create dir for each dmrg calculation
        FileHandler.make_directory_from_project_root(f"dmrg_{i}")
        energy_partial, s1_partial, s2_partial, mut_inf_partial = molcas.calculate(
            partial_occupation, large_cas_indices[i]
        )
        # partial energy is useless
        _ = energy_partial
        partial_s1_list.append(s1_partial)
        partial_s2_list.append(s2_partial)
        partial_mut_inf_list.append(mut_inf_partial)

    molcas.initial_cas_prepared()

    # combine entropies and occupation to start autocas
    # useful for plotting stuff
    initial_occupation, initial_s1, initial_s2, initial_mut_inf = autocas.collect_entropies(  # type: ignore
        large_cas_indices,
        large_cas_occupations,
        partial_s1_list,  # type: ignore
        partial_s2_list,  # type: ignore
        partial_mut_inf_list,  # type: ignore
    )

    # plot stuff
    _ = initial_occupation  # not required for plotting
    _ = initial_s2  # not required for plotting
    plot = EntanglementPlot()
    plt = plot.plot(initial_s1, initial_mut_inf)
    plt.savefig("entang.pdf")

    # shortcut method for collect_enetropies and get_active_space
    cas_occ, cas_index = autocas.get_cas_from_large_cas(
        large_cas_indices,
        large_cas_occupations,
        partial_s1_list,  # type: ignore
        partial_s2_list,  # type: ignore
        partial_mut_inf_list,  # type: ignore
    )

    # cas and hyphen do not matter for method names
    molcas.set_cas_method("casscf")

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 2000
    molcas.settings.dmrg_sweeps = 20
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
    occupation, orbitals, energy = large_cas_autocas_example(path_to_this_file, xyz)
    print(f"no \ncas: {occupation} \norbs: {orbitals} \nenergy: {energy}")

    # to verify workflow
    print("\n\nonly verification run")

    # change project name, to prevent calculation in same project
    FileHandler.project_dir = path_to_this_file
    FileHandler.DirectoryNames.project_name = "verification_run"
    FileHandler.setup_project()
    test_molecule = Molecule(xyz)
    test_autocas = Autocas(test_molecule)
    test_interface = Molcas(test_molecule)

    test_interface.project_name = "verify_large_cas_example"
    test_interface.environment.molcas_scratch_dir = path_to_this_file + "/../test/scratch_tmp"

    test_interface.set_cas_method("casscf")
    test_interface.set_post_cas_method("caspt2")
    test_interface.settings.dmrg_bond_dimension = 250
    test_interface.settings.dmrg_sweeps = 5

    # run workflow
    workflow = LargeCasWorkflow(test_autocas, test_interface)
    workflow.run()

    assert abs(workflow.results["final_energy"][0] - energy[0]) < 1e-9
    assert workflow.results["final_occupation"] == occupation
    assert workflow.results["final_orbital_indices"] == orbitals
