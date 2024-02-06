"""Modify Molcas Hdf5 files.

This module implements a class to read and modify molcas input files.
The modification of the orbital hdf5 file allows to change the order of orbitals,
required for active space calculations, without using MOLCASs 'alter' keyword.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Union

import h5py
import numpy as np


class MolcasHdf5Utils:
    """Class to handle the Molcas HDF5 files.

    Attributes
    ----------
    self.nbas: np.ndarray
        number of basis functions per symmetry
    self.irrep_labels: np.ndarray
        name of the symmetries
    self.natoms: int = 0
        number of atoms
    self.module: str = ""
        name of the Molcas module
    self.orbital_type: str = ""
        type of orbitals
    # Datasets
    self.type_indices: np.ndarray
        molcas typeindex for each orbital
    self.mo_energies: List[float] = []
        array with MO energies ordered first by symmetry, then by energy
    # extracted from hdf5 file
    self.occupations: List[int] = []
        contains the occupation of each orbital,
       e.g. 0: virtual, 1: singly occupied (uhf), 2: doubly occupied (rhf)
    self.symmetries: List[int] = []
        contains the symmetry index of each orbtial
    self.energy: Union[np.ndarray, float] = 0
        contains the energy for every state
    """

    __slots__ = [
        "nbas",
        "irrep_labels",
        "natoms",
        "module",
        "orbital_type",
        "type_indices",
        "mo_energies",
        "occupations",
        "symmetries",
        "energy",
    ]

    def __init__(self):
        """Construct class."""
        # Attributes
        self.nbas: np.ndarray
        """number of basis functions per symmetry"""
        self.irrep_labels: np.ndarray
        """name of the symmetries"""
        self.natoms: int = 0
        """number of atoms"""
        self.module: str = ""
        """name of the Molcas module"""
        self.orbital_type: str = ""
        """type of orbitals"""
        # Datasets
        self.type_indices: np.ndarray
        """molcas typeindex for each orbital"""
        self.mo_energies: List[float] = []
        """array with MO energies ordered first by symmetry, then by energy"""
        # extracted from hdf5 file
        self.occupations: List[int] = []
        """contains the occupation of each orbital,
           e.g. 0: virtual, 1: singly occupied (uhf), 2: doubly occupied (rhf)"""
        self.symmetries: List[int] = []
        """contains the symmetry index of each orbtial"""
        self.energy: Union[np.ndarray, float] = 0
        """contains the energy for every state"""

    def modify_hdf5(self, hdf5_file: str, cas_orbitals: List[int]):
        """Modify the Molcas orbital file, to enable typeindices in Molcas.

        The active space is encoded in the hdf5 file.

        Parameters
        ----------
        hdf5_file : str
            name of the orbital file to modify
        cas_orbitals : List[int]
            orbital indices corresponding to orbitals in active space

        Raises
        ------
        KeyError
            if no typeindices are found in the hdf5 file.
        """
        h5_file = h5py.File(hdf5_file, "r+")
        # check if mo or mo_alpha typeindices are present
        mo_alpha_bool = False
        try:
            # hdf5_type_indices = h5_file["MO_ALPHA_TYPEINDICES"]
            mo_alpha_bool = True
        except KeyError:
            try:
                hdf5_type_indices = h5_file.get("MO_TYPEINDICES")
            except KeyError as exc:
                raise KeyError(
                    "No Typeindices found in the orbital file, smth went wrong"
                ) from exc

        # get unset typeindices
        type_indices = self.type_indices
        type_indices = np.array(type_indices)

        # set new type indices
        mo_alpha_bool = False
        if not mo_alpha_bool:
            hdf5_type_indices = h5_file["MO_TYPEINDICES"]
        else:
            hdf5_type_indices = h5_file["MO_ALPHA_TYPEINDICES"]
        for i in cas_orbitals:
            if mo_alpha_bool:
                type_indices[i] = "2"
            else:
                type_indices[i] = "2"

        hdf5_type_indices[...] = type_indices
        h5_file.close()

    def read_hdf5(self, hdf5_file: str):
        """Read required parameters from an Molcas HDF5 file.

        Parameters
        ----------
        hdf5_file : str
            name of the HDF5 file_name
        """
        # get basic stuff from h5 file
        h5_file = h5py.File(hdf5_file, "r")
        self.nbas = np.array(h5_file.attrs["NBAS"])
        self.irrep_labels = np.array(h5_file.attrs["IRREP_LABELS"])
        self.natoms = h5_file.attrs["NATOMS_UNIQUE"]
        self.module = h5_file.attrs["MOLCAS_MODULE"].item().decode("UTF-8")

        # get orbital type
        orbital_type = ""
        if self.module == "RASSCF":
            orbital_type = "RASSCF"
        else:
            orbital_type = h5_file.attrs["ORBITAL_TYPE"].item().decode("UTF-8")
        self.orbital_type = orbital_type

        # get type indices
        alpha = True
        # type_indices = h5_file.get("MO_TYPEINDICES")
        type_indices = np.array(h5_file.get("MO_ALPHA_TYPEINDICES"))
        try:
            if not type_indices:
                type_indices = h5_file.get("MO_TYPEINDICES")
                alpha = False
        except ValueError:
            pass
        self.type_indices = np.array(type_indices)

        indices = []
        for _, index in enumerate(type_indices):
            indices.append(index.decode("UTF-8"))
        type_indices = np.array(indices)

        mo_energies = np.array(h5_file.get("MO_ENERGIES"))
        if mo_energies.all() is None:
            mo_energies = np.array(h5_file.get("MO_ALPHA_ENERGIES"))
        mo_energies = np.array(mo_energies)

        # because numpy 1.19
        mo_energies = mo_energies.astype(float, copy=False)
        self.mo_energies = list(mo_energies)

        last_indices = 0
        tmp_index = 0
        if len(self.occupations) > 0:
            self.occupations = []
        if len(self.symmetries) > 0:
            self.symmetries = []
        for i in range(self.nbas.shape[0]):
            for _ in range(self.nbas[i]):
                electrons = 0
                if type_indices[tmp_index] == "I":
                    if alpha:
                        electrons = 1
                    else:
                        electrons = 2
                self.symmetries.append(i)
                self.occupations.append(electrons)
                tmp_index += 1
            last_indices += self.nbas[i]
        # self.type_indices = type_indices
        h5_file.close()

    def get_energy(self, hdf5_file: str):
        """Read the DMRG energy from a Molcas HDF5 file.

        Parameters
        ----------
        hdf5_file : str
            name of the HDF5 file

        Returns
        -------
        self.energy : Union[float, np.ndarray]
            if only one state found, then just the energy for this state.
            else an ndarray with an energy for each root.
        """
        h5_file = h5py.File(hdf5_file, "r")
        try:
            self.energy = np.array(h5_file.get("ROOT_ENERGIES"))[0]
        except KeyError:
            self.energy = np.array(h5_file.get("ROOT_ENERGIES"))
        h5_file.close()
        return self.energy
