# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Tuple
from scine_autocas.workflows.consistent_active_space.configuration import ConsistentActiveSpaceConfiguration
from scine_autocas import Molecule
from scine_autocas.interfaces import Molcas, Serenity
from scine_autocas.utils.defaults import Defaults


def setup_molcas_and_molecule(xyz_file: str, name: str, basis_set: str) -> Tuple[Molcas, Molecule]:
    """
    Setup Molcas interface and molecule from an XYZ file.

    Parameters
    ----------
    xyz_file: str
        The path to the XYZ file containing the molecule.
    name: str
        The system name.
    basis_set: str
        The AO basis set.

    Returns
    -------
    Tuple[Molcas, Molecule]
        The Molcas interface and the Molecule object.
    """
    # create a molecule
    molecule = Molecule(xyz_file)
    # initialize autoCAS and Molcas interface
    molcas = Molcas(molecule)
    # setup interface
    molcas.project_name = name
    molcas.settings.basis_set = basis_set
    # # Prepare settings for the DMRG calculation, generating the single-orbital entropies.
    molcas.set_cas_method("dmrgci")
    # manually set dmrg sweeps and bond dmrg_bond_dimension to low number
    molcas.settings.dmrg_bond_dimension = Defaults.Interface.init_dmrg_bond_dimension
    molcas.settings.dmrg_sweeps = Defaults.Interface.init_dmrg_sweeps
    molcas.calculate()
    return molcas, molecule


def print_orbital_map(orbital_map: List[List[List[int]]]) -> None:
    """
    Print the orbital map in a readable format.

    Parameters
    ----------
    orbital_map: List[List[List[int]]]
        The orbital map.
    """
    print("Orbital groups")
    for group in orbital_map:
        n_orbitals = len(group[0])
        to_out = ''
        for i_orb in range(n_orbitals):
            for idx in group:
                orb = idx[i_orb]
                to_out += f'{orb:8}  '
            to_out += '\n'
        print(to_out)


def construct_molecules(configuration: ConsistentActiveSpaceConfiguration) -> Tuple[List[Molcas], List[Molecule]]:
    """
    Construct the Molcas interfaces and Molecule objects from the configuration.

    Parameters
    ----------
    configuration: ConsistentActiveSpaceConfiguration
        The calculation configuration.

    Returns
    -------
    Tuple[List[Molcas], List[Molecule]]
        A tuple containing the list of Molcas interfaces and the list of Molecule objects.
    """
    interfaces = []
    molecules = []
    for xyz, name in zip(configuration.xyz_files, configuration.system_names):
        molcas, molecule = setup_molcas_and_molecule(xyz, name, configuration.basis_set)
        molecules.append(molecule)
        interfaces.append(molcas)
    return interfaces, molecules


def load_orbitals_from_serenity(load_paths: List[str], settings: dict):
    """
    Load orbitals from Serenity and return the Serenity object and the list of molecules.

    Parameters
    ----------
    load_paths: List[str]
        The paths to load the orbitals/systems from.
    settings: dict
        The serenity interface settings.

    Returns
    -------
    Tuple[Serenity, List[Molecule]]
        A tuple containing the Serenity object and the list of Molecule objects.
    """
    serenity = Serenity([], settings, load_paths)
    serenity.settings.molcas_orbital_files = settings["Interface"]["molcas_orbital_files"]
    molecules = serenity.molecules
    return serenity, molecules
