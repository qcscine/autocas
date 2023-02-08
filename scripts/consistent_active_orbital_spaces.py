"""Consistent active orbital spaces for multiple structures on a reaction coordinate.

  The input options are described in `set_up_options()` or available through `python this_script.py -h`

  Example input: python this_script.py -m caspt2 -b cc-ptz path_to_first_structure.xyz path_to_second_structure.xyz
  Output: * A file with the combined active space indices.
          * A file with the orbital map.
          * A file with all final energies.
          * A directory called 'molcas' containing all Molcas files (including entranglement and threshold diagrams).
          * A directory called 'serenity' containing all Serenity files.
          * Energies, single orbital entropies, etc. are printed to the standard output.

  The occupied and virtual valence orbitals are first localized and then mapped between structures. The active space
  is then determined through autoCAS and made consistently through the orbital maps.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This file is part of SCINE AutoCAS.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
See LICENSE.txt for details
"""

import os
from typing import List, Tuple

from scine_autocas import Autocas
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.interfaces.serenity import Serenity
from scine_autocas.plots.entanglement_plot import EntanglementPlot
from scine_autocas.plots.threshold_plot import ThresholdPlot
from scine_autocas.autocas_utils.cas_combination import combine_active_spaces

# Requires the Serenity-Python wrapper.
import serenipy as spy


def set_up_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--skip_scf", dest="skip_scf", action="store_true", default=False,
                      help="If true, paths to the Serenity system directories are expected instead of the path to the XYZ files."
                           " The orbitals are then loaded from these directories.")
    parser.add_option("-m", "--cas_method", dest="cas_method", action="store", default="caspt2")
    parser.add_option("-i", "--autocas_indices", dest="autocas_molecule_indices", default="",
                      help="The molecule indices used for the autocas selection. By default the CAS is selected for all molecules.")
    parser.add_option("-b", "--basis", dest="basis_set", default="cc-pvdz",
                      help="The atomic basis set label.")
    parser.add_option("-l", "--large_active_space", dest="large_active_space", default=True, action="store_false",
                      help="If true, the large active space protocoll in autocas is used. By default true.")
    parser.add_option("-u", "--unmapable", dest="always_include_unmapables", default=False, action="store_true")
    return parser.parse_args()


def setup_molcas_and_molecule(path: str, xyz_file: str, name: str, options) -> Tuple[List[Molcas], List[Molecule]]:
    # create a molecule
    molecule = Molecule(xyz_file)

    # initialize autoCAS and Molcas interface
    molcas = Molcas([molecule])

    # setup interface
    molcas.project_name = name
    molcas.basis_set = options.basis_set
    molcas.settings.work_dir = path + "/molcas/" + name
    molcas.environment.molcas_scratch_dir = path + "/molcas/scratch"
    molcas.settings.xyz_file = xyz_file
    molcas.settings.fiedler = True
    molcas.settings.skip_scf = True  # We calculate the SCF with Serenity because it is faster.

    # Prepare settings for the DMRG calculation, generating the single-orbital entropies.
    molcas.settings.method = "DMRGCI"

    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = 250
    molcas.settings.dmrg_sweeps = 5

    # Run 1 SCF step to generate the x.scf.h5 files.
    molcas.calculate()
    return molcas, molecule


def run_small_active_space(autocas: Autocas, molcas: Molcas, name: str, occ_initial: List[float], index_initial: List[int]) -> Tuple[List[float], List[float]]:

    # do cas calculation
    energy, s1_entropy, s2_entropy, mut_inf = molcas.calculate(occ_initial, index_initial)

    # plot entanglement + threshold diagram
    plot = EntanglementPlot()
    plt = plot.plot(s1_entropy, mut_inf)  # type: ignore
    plt.savefig(molcas.settings.work_dir + "/entang" + name + ".pdf")  # type: ignore
    plot = ThresholdPlot()
    plt = plot.plot(s1_entropy)  # type: ignore
    plt.savefig(molcas.settings.work_dir + "/threshold" + name + ".pdf")  # type: ignore
    try:
        # make active space based on single orbital entropies
        cas_occ, cas_index = autocas.get_active_space(
            occ_initial, s1_entropy   # type: ignore
        )
    except Exception:
        return list(), list()
    return cas_occ, cas_index


def run_large_active_space(autocas: Autocas, molcas: Molcas, name: str) -> Tuple[List[float], List[float]]:
    large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()
    # iterate over subspaces and evaluate them
    partial_s1_list = []
    partial_s2_list = []
    partial_mut_inf_list = []
    print(f"autocas will perform: {len(large_cas_indices)} dmrg calculations")
    print(f"with {autocas.large_spaces.max_orbitals} in each active space")

    for partial_indices, partial_occupation in zip(large_cas_indices, large_cas_occupations):
        _, s1_partial, s2_partial, mut_inf_partial = molcas.calculate(
            partial_occupation, partial_indices
        )
        partial_s1_list.append(s1_partial)
        partial_s2_list.append(s2_partial)
        partial_mut_inf_list.append(mut_inf_partial)

    # combine entropies and occupation to start autocas
    # useful for plotting
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
    plt.savefig(molcas.settings.work_dir + "/entang" + name + ".pdf")
    plot = ThresholdPlot()
    plt = plot.plot(initial_s1)  # type: ignore
    plt.savefig(molcas.settings.work_dir + "/threshold" + name + ".pdf")  # type: ignore

    try:
        # shortcut method for collect_enetropies and get_active_space
        cas_occ, cas_index = autocas.get_cas_from_large_cas(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_list,  # type: ignore
            partial_s2_list,  # type: ignore
            partial_mut_inf_list,  # type: ignore
            force_cas=False
        )
    except Exception:
        return list(), list()

    # no large cas required anymore
    molcas.dumper.large_cas = False
    return cas_occ, cas_index


def run_autocas(molecule: Molecule, molcas: Molcas, name: str, options) -> Tuple[List[float], List[float]]:
    """Do ground state cas"""
    # initialize autoCAS
    autocas = Autocas(molecule, {"plateau_values": 12})
    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()

    if options.large_active_space and len(occ_initial) > autocas.large_spaces.max_orbitals:
        print("Large Active Space Selection")
        autocas.large_spaces.max_orbitals = min(len(occ_initial), autocas.large_spaces.max_orbitals)
        return run_large_active_space(autocas, molcas, name) 
    else:
        return run_small_active_space(autocas, molcas, name, occ_initial, index_initial)


def run_final_calculation(molcas: Molcas, cas_occ: List[float], cas_index: List[int], method: str):
    # cas and hyphen do not matter for method names
    molcas.settings.method = method
    if "PT2" in method:
        molcas.settings.post_cas_method = "CASPT2"
        molcas.settings.method = "CASSCF"

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


def print_orbital_map(orbital_map: List[List[List[int]]]):
    print("Orbital groups")
    for group in orbital_map:
        n_orbitals = len(group[0])
        to_out = ''
        for i_orb in range(n_orbitals):
            for i_sys in range(len(group)):
                orb = group[i_sys][i_orb]
                to_out += f'{orb:8}  '
            to_out += '\n'
        print(to_out)


def construct_molecules(xyz_files: List[str], path_to_this_file: str, options, names: str):
    interfaces = []
    molecules = []
    run_scf = not options.skip_scf
    for xyz, name in zip(xyz_files, names):
        molcas, molecule = setup_molcas_and_molecule(path_to_this_file, xyz, name, options)
        molecules.append(molecule)
        interfaces.append(molcas)
        expected_path = molcas.settings.work_dir + "/initial/" + name + ".scf.h5"
    return interfaces, molecules


def load_orbitals_from_Serenity(options, load_paths: List[str], settings: dict):
    serenity = Serenity([], settings, load_paths)
    molecules = serenity.molecules
    return serenity, molecules


def main(path_to_this_file: str):
    options, xyz_files = set_up_options()
    if not xyz_files:
        print("No XYZ files supplied for the calculation. Stopping!")
        return
    xyz_files = [path_to_this_file + "/" + file for file in xyz_files]

    settings = {
        "settings": {
            "uhf": False,
            "localisation_method": "IBO",
            "alignment": True,
            "localize_virtuals": True,
            "optimized_mapping": True,
            "work_dir": "serenity/",
            "basis_set": options.basis_set,
            "score_start": 0.5,
            "skip_localization": options.skip_scf
        }
    }
    print("Final energy evaulation method: " + options.cas_method)
    # We first clean up any leftovers from molcas.
    if os.path.exists("molcas"):
        import shutil
        shutil.rmtree("molcas")
    # Load orbitals if desired.
    serenity = None
    if options.skip_scf:
        settings["settings"]["system_names"] = [name.split("/")[-1] for name in xyz_files]
        load_paths = []
        for path, name in zip(xyz_files, settings["settings"]["system_names"]):
            reduced_path = path[:-len(name)]
            load_paths.append(reduced_path)
        serenity, molecules = load_orbitals_from_Serenity(options, load_paths, settings)
    else:
        settings["settings"]["system_names"] = [str(i) for i in range(len(xyz_files))]
        settings["settings"]["xyz_files"] = xyz_files
    # Generate the Molcas orbital files.
    interfaces, molecules = construct_molecules(settings["settings"]["xyz_files"], path_to_this_file, options, settings["settings"]["system_names"])
    settings["settings"]["molcas_orbital_files"] = [m.settings.work_dir + "/initial/" for m in interfaces]
    
    # If the orbitals are loaded from existing files, Serenity will be initialized. Otherwise, we run the SCF calculation with Serenity.
    if serenity is None:
        serenity = Serenity(molecules, settings)
        serenity.load_or_write_molcas_orbitals()
        serenity.calculate()
    else:
        serenity.settings.molcas_orbital_files = [m.settings.work_dir + "/initial/" for m in interfaces]
    names = serenity.settings.system_names
    # Write canonical orbitals back to Molcas
    serenity.load_or_write_molcas_orbitals(True)
    # Generate mapping
    orbital_map, unmappable_orbitals = serenity.get_orbital_map()
    # Write localized orbitals back to Molcas
    print_orbital_map(orbital_map)
    serenity.load_or_write_molcas_orbitals(True)

    autocas_system_indices = [int(str_i) for str_i in options.autocas_molecule_indices.split()]
    if not autocas_system_indices:
        autocas_system_indices = [i for i in range(len(names))]
    
    cas_occupations = [[] for _ in names]
    cas_indices = [[] for _ in names]
    # Run autoCAS
    for i in autocas_system_indices:
        molecule = molecules[i]
        molcas = interfaces[i]
        name = names[i]
        cas_occ, cas_index = run_autocas(molecule, molcas, name, options)
        cas_occupations[i] = cas_occ
        cas_indices[i] = cas_index

    print("*******************************************************************************************")
    print("*                                                                                         *")
    print("*                               Combined Active Spaces                                    *")
    print("*                                                                                         *")
    print("*******************************************************************************************")
    # Combine CAS
    combined_occupations, combined_indices = combine_active_spaces(cas_occupations, cas_indices, orbital_map)
    # Always include all unmappable orbitals if required.
    if options.always_include_unmapables:
        unmappable_occupied = unmappable_orbitals[0]
        unmappable_virtuals = []
        occ = 2 
        if len(unmappable_orbitals) > 1:
            unmappable_virtuals = unmappable_orbitals[1]
        for cas_index, cas_occ, u_occ, u_virt in zip(combined_indices, combined_occupations, unmappable_occupied, unmappable_virtuals):
            for i in u_occ:
                if i not in cas_index:
                    cas_index.append(i)
                    cas_occ.append(occ)
            for i in u_virt:
                if i not in cas_index:
                    cas_index.append(i)
                    cas_occ.append(0)

    combined_file = open("combined_cas_spaces", "w")
    for cas_index, cas_occ in zip(combined_indices, combined_occupations):
        print(f"combined cas indices: {cas_index}")
        print(f"combined occupation:  {cas_occ}")
        combined_file.write(f"combined cas indices: {cas_index}\n")
        combined_file.write(f"combined occupation:  {cas_occ}\n")
    combined_file.close()

    energies = []
    # Run DMRG-SCF or CAS-PT2
    for molcas, cas_occ, cas_index in zip(interfaces, combined_occupations, combined_indices):
        _, _, energy = run_final_calculation(molcas, cas_occ, cas_index, options.cas_method)
        energies.append(energy)

    f = open("energies.dat", "w")
    for e in energies:
        f.write(str(e)+"\n")
    f.close()


if __name__ == "__main__":
    main(os.path.dirname(os.path.abspath(__file__)))
