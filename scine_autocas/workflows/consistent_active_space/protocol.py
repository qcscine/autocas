# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Tuple
from optparse import OptionParser
from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration
from scine_autocas.cas_selection.cas_combination import combine_active_spaces
from scine_autocas.interfaces import Molcas, Serenity
from scine_autocas.utils.defaults import CasMethods, PostCasMethods
from scine_autocas.workflows.consistent_active_space.prepare_orbitals import (
    load_orbitals_from_serenity, construct_molecules, print_orbital_map
)
from scine_autocas.utils.defaults import Defaults
from scine_autocas.workflows.consistent_active_space.run_autocas import run_autocas
from scine_autocas.io import FileHandler


def run_from_command_line() -> None:
    """
    Run the consistent active space protocol from the command line.
    """
    parser = OptionParser(description="Run the consistent active space protocol for a set of molecules. The input can"
                                      " be either a set of XYZ files and command line options or a yaml file provided"
                                      " through the option -y. ")
    parser.add_option("-s", "--load_orbitals", dest="load_orbitals", action="store_true", default=False,
                      help="If true, paths to the Serenity system directories are expected instead of the"
                           " path to the XYZ files. The orbitals are then loaded from these directories.")
    parser.add_option("-m", "--cas_method", dest="cas_method", action="store", default="CASPT2")
    parser.add_option("-i", "--autocas_indices", dest="autocas_molecule_indices", default="",
                      help="The molecule indices used for the autocas selection. By default the CAS"
                           " is selected for all molecules.")
    parser.add_option("-b", "--basis", dest="basis_set", default=Defaults.Interface.basis_set,
                      help="The atomic basis set label.")
    parser.add_option("-l", "--large_active_space", dest="large_active_space", default=True, action="store_false",
                      help="If given, the large active space protocoll in autocas is not used. By default it is used.")
    parser.add_option("-c", "--exclude_core", dest="exclude_core", default=True, action="store_false",
                      help="If given, core orbitals can be selected for the active space. By default, all core"
                           " orbitals are removed.")
    parser.add_option("-u", "--unmapable", dest="always_include_unmapables", default=False, action="store_true")
    parser.add_option("-y", "--yaml", dest="yaml_file", default="", type="str",
                      help="The configuration yaml file to use. If given, xyz files/loading paths must not be set"
                           " and all other options provided through the command line are ignored.")
    parser.add_option("-n", "--create_no_yaml", dest="create_yaml", default=True, action="store_false",
                      help="If true, the configuration will be written to a yaml file. By default, a yaml"
                           " file is created.")
    (options, args) = parser.parse_args()
    if options.yaml_file:
        if len(args) > 0:
            raise ValueError("If a yaml file is given, no other arguments must be set.")
        configuration = ConsistentActiveSpaceConfiguration.from_file(options.yaml_file)
    else:
        configuration = ConsistentActiveSpaceConfiguration.from_options(options, args)
        if options.create_yaml:
            configuration.write_yaml_file()
    run_consistent_active_space_protocol(configuration)


def run_consistent_active_space_protocol(configuration: ConsistentActiveSpaceConfiguration)\
        -> Tuple[List[List[int]], List[List[int]], List[float]]:
    """
    Run the consistent active space protocol.

    Parameters
    ----------
    configuration: ConsistentActiveSpaceConfiguration
        The configuration for the protocol.

    Returns
    -------
    Tuple[List[List[int]], List[List[int]], List[float]]
        The combined active space occupations, indices, and energies for each structure.
    """
    initial_project_name = FileHandler.DirectoryNames.project_name
    FileHandler.DirectoryNames.project_name = configuration.project_name
    FileHandler.setup_project()
    print("Final energy evaluation method: " + configuration.cas_method)
    print("*******************************************************************************************")
    print("*                                                                                         *")
    print("*                               Orbital Preparation                                       *")
    print("*                                                                                         *")
    print("*******************************************************************************************")
    # Load orbitals if desired.
    serenity = None
    settings = configuration.get_serenity_interface_settings()
    if configuration.load_orbitals:
        load_paths = configuration.get_serenity_load_paths()
        serenity, molecules = load_orbitals_from_serenity(load_paths, settings)
        orbital_map, unmappable_orbitals = serenity.get_orbital_map()
    # Generate the Molcas orbital files.
    interfaces, molecules = construct_molecules(configuration)

    # If the orbitals are loaded from existing files, Serenity will be initialized.
    # Otherwise, we run the SCF calculation with Serenity.
    if serenity is None:
        serenity = Serenity(molecules, settings)
        serenity.load_or_write_molcas_orbitals()
        serenity.calculate()
        orbital_map, unmappable_orbitals = serenity.get_orbital_map()
    names = serenity.settings.system_names
    # Write canonical orbitals back to Molcas
    serenity.load_or_write_molcas_orbitals(True)
    # Write localized orbitals back to Molcas
    print_orbital_map(orbital_map)
    serenity.load_or_write_molcas_orbitals(True)

    cas_occupations: List[List[int]] = [[] for _ in names]
    cas_indices: List[List[int]] = [[] for _ in names]
    # Run autoCAS
    for i in configuration.autocas_indices:
        molecule = molecules[i]
        molcas = interfaces[i]
        name = names[i]
        # reset the interface
        molcas.set_initial_cas_state(False)
        cas_occup, cas_idx = run_autocas(molecule, molcas, name, configuration.large_active_space)
        cas_occupations[i] = cas_occup  # type: ignore
        cas_indices[i] = cas_idx

    print("*******************************************************************************************")
    print("*                                                                                         *")
    print("*                               Combined Active Spaces                                    *")
    print("*                                                                                         *")
    print("*******************************************************************************************")

    # Combine CAS
    combined_occupations, combined_indices = combine_active_spaces(cas_occupations, cas_indices, orbital_map)
    # Always include all unmappable orbitals if required.
    if configuration.unmappable:
        unmappable_occupied = unmappable_orbitals[0]
        unmappable_virtuals: List[List[int]] = []
        occ = 2
        if len(unmappable_orbitals) > 1:
            unmappable_virtuals = unmappable_orbitals[1]
        for cas_index, cas_occ, u_occ, u_virt in zip(combined_indices, combined_occupations, unmappable_occupied,
                                                     unmappable_virtuals):
            for i in u_occ:
                if i not in cas_index:
                    cas_index.append(i)
                    cas_occ.append(occ)
            for i in u_virt:
                if i not in cas_index:
                    cas_index.append(i)
                    cas_occ.append(0)

    if configuration.exclude_core:
        n_core_orbitals = molecules[0].core_orbitals
        for active_space, occupations in zip(combined_indices, combined_occupations):
            for i in range(n_core_orbitals):
                if i in active_space:
                    idx = active_space.index(i)
                    active_space.remove(i)
                    occupations.remove(occupations[idx])

    combined_file = open("combined_cas_spaces", "w")
    for cas_index, cas_occ in zip(combined_indices, combined_occupations):
        print(f"combined cas indices: {cas_index}")
        print(f"combined occupation:  {cas_occ}")
        combined_file.write(f"combined cas indices: {cas_index}\n")
        combined_file.write(f"combined occupation:  {cas_occ}\n")
    combined_file.close()

    print("*******************************************************************************************")
    print("*                                                                                         *")
    print("*                               Final Calculations                                        *")
    print("*                                                                                         *")
    print("*******************************************************************************************")

    energies: List[float] = []
    # Run DMRG-SCF or CAS-PT2
    for molcas, cas_occ, cas_index in zip(interfaces, combined_occupations, combined_indices):
        _, _, energy = run_final_calculation(molcas, cas_occ, cas_index, configuration.cas_method)
        energies.append(energy)

    f = open("energies.dat", "w")
    for e in energies:
        f.write(str(e) + "\n")
    f.close()
    FileHandler.DirectoryNames.project_name = initial_project_name
    return combined_occupations, combined_indices, energies


def run_final_calculation(molcas: Molcas, cas_occ: List[int], cas_index: List[int], method: str)\
        -> Tuple[List[int], List[int], float]:
    """
    Run the final calculation with the given method.

    Parameters
    ----------
    molcas: Molcas
        The molcas interface.
    cas_occ: List[int]
        The CAS occupations.
    cas_index: List[int]
        The orbital indices for the CAS.
    method: str
        The method to use for the final calculation, e.g., "CASSCF", "CASCI", "DMRGCI", "DMRGSCF", or "CASPT2".

    Returns
    -------
    Tuple[List[int], List[int], float]
        The CAS occupations, CAS indices, and the final energy.
    """
    # cas and hyphen do not matter for method names
    if "PT2" in method:
        print("Running CASPT2")
        molcas.settings.post_cas_method = PostCasMethods.CASPT2
        molcas.settings.cas_method = CasMethods.CASSCF
    else:
        if method == "CASSCF":
            molcas.settings.cas_method = CasMethods.CASSCF
        elif method == "CASCI":
            molcas.settings.cas_method = CasMethods.CASCI
        elif method == "DMRGCI":
            molcas.settings.cas_method = CasMethods.DMRGCI
        elif method == "DMRGSCF":
            molcas.settings.cas_method = CasMethods.DMRGSCF
        else:
            raise NotImplementedError(f"Method: {method} is not available")

    molcas.settings.dmrg_bond_dimension = Defaults.Interface.dmrg_bond_dimension
    molcas.settings.dmrg_sweeps = Defaults.Interface.dmrg_sweeps
    # ensure interface is set correctly
    molcas.set_initial_cas_state(True)
    molcas.set_orbital_state(True)

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
    if len(final_energy) < 1:
        raise RuntimeError("Final energy calculation did not return a valid energy value.")
    return cas_occ, cas_index, final_energy[0]
