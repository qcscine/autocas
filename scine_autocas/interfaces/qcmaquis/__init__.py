"""Module to handle qcmaquis output files.

This module implements the Qcmaquis class which evaluates entropies and
mutual information from the qcmaquis hdf5 output.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import numpy as np

from scine_autocas.interfaces.qcmaquis.qcmaquis_hdf5_utils import Hdf5Converter
from scine_autocas.interfaces.qcmaquis.qcmaquis_orbital_rdm_builder import OrbitalRDMBuilder


class Qcmaquis:
    """Class to handle all qcmaquis output.

    This class purpose is to read qcmaquis output and evaluate entropies from
    the provided measurements.

    Attributes
    ----------
    hdf5_converter : Hdf5Converter
        reads the qcmaquis hdf5 output file
    orbital_rdm_builder : OrbitalRDMBuilder
        evaluates the one and two orbital rdm from data provided by the hdf5 converter
    s1_entropy : np.ndarray
        stores the sinlge orbital entropy
    s2_entropy : np.ndarray
        stores the two orbital two orbital entropy
    mutual_information : np.ndarray
        stores the mutual information
    n_orbitals : int
        the number of orbitals in the active space
    energy : float
        the energy of the qcmaquis calculation
    """

    __slots__ = (
        "hdf5_converter",
        "orbital_rdm_builder",
        "s1_entropy",
        "s2_entropy",
        "mutual_information",
        "n_orbitals",
        "energy",
    )

    def __init__(self):
        """Init."""
        self.hdf5_converter: Hdf5Converter = Hdf5Converter()
        """extracts one and two orbital RDM from QCMaquis results fil"""
        self.orbital_rdm_builder: OrbitalRDMBuilder = OrbitalRDMBuilder()
        """evaluates and stores the orbital RMDs"""
        self.s1_entropy: np.ndarray = np.array([])
        """single orbital entrop"""
        self.s2_entropy: np.ndarray = np.array([])
        """two orbital entrop"""
        self.mutual_information: np.ndarray = np.array([])
        """mutual information"""
        self.n_orbitals: int = 0
        """number of orbitals in CA"""
        self.energy: float = 0
        """DMRG energ"""

    def make_s1(self, one_ordm: np.ndarray) -> np.ndarray:
        """Evaluate single orbital entropy from one orbital RDM.

        Parameters
        ----------
        one_ordm : np.ndarray
            one orbital RDM

        Returns
        -------
        s1 : np.ndarray
            single orbital entropy

        Raises
        ------
        AttributeError
            if hdf5_converter has not read the output yet
        """
        if not self.hdf5_converter.L:
            raise AttributeError("Read hdf5 file before, doing measurments")
        self.s1_entropy = np.zeros((self.n_orbitals))
        for site in range(self.n_orbitals):
            s1_entropy = 0
            for alpha in range(len(one_ordm[site])):
                eigenvalue = one_ordm[site, alpha]
                if eigenvalue > 0:
                    s1_entropy = s1_entropy - eigenvalue * np.log(eigenvalue)
            self.s1_entropy[site] = s1_entropy
        return self.s1_entropy

    def make_s2(self, two_ordm: np.ndarray) -> np.ndarray:
        """Evaluate two orbital entropy from two orbital RDM.

        Parameters
        ----------
        two_ordm : np.ndarray
            two orbital RDM

        Returns
        -------
        s2 : np.ndarray
            two orbital entropy

        Raises
        ------
        AttributeError
            if hdf5_converter has not read the output yet
        """
        if self.hdf5_converter.L is None:
            raise AttributeError("Read hdf5 file before, doing measurments")
        self.s2_entropy = np.zeros((self.n_orbitals, self.n_orbitals))
        for site_1 in range(self.n_orbitals):
            for site_2 in range(site_1 + 1, self.n_orbitals):
                sub_matrix = two_ordm[site_1][site_2][:][:]
                eigenvalue, _ = np.linalg.eig(sub_matrix)  # type: ignore[attr-defined]
                s2_entropy = 0
                for alpha in range(16):
                    if eigenvalue[alpha] > 0:
                        s2_entropy = s2_entropy - (eigenvalue[alpha] * np.log(eigenvalue[alpha]))
                self.s2_entropy[site_1, site_2] = s2_entropy.real
                self.s2_entropy[site_2, site_1] = s2_entropy.real
        return self.s2_entropy

    def make_mutual_information(self) -> np.ndarray:
        """Evaluate the mutual information.

        Returns
        -------
        mutual_information : np.ndarray
            mutual information

        Raises
        ------
        AttributeError
            if s1 and s2 entropies are not evaluated
        """
        if not (self.s1_entropy.size and self.s2_entropy.size):
            raise AttributeError("Evaluate s2 and s2 entropy before evaluating mutual information")
        self.mutual_information = np.zeros((self.n_orbitals, self.n_orbitals))
        for site_1 in range(self.n_orbitals):
            for site_2 in range(site_1 + 1, self.n_orbitals):
                self.mutual_information[site_1, site_2] = 0.5 * (
                    self.s1_entropy[site_1] + self.s1_entropy[site_2] - self.s2_entropy[site_1, site_2]
                )
                self.mutual_information[site_2, site_1] = self.mutual_information[site_1, site_2]
        return self.mutual_information

    def make_diagnostics(self):
        """Generate orbital RDMs from QCMaquis output file and evaluate entropies.

        Raises
        ------
        AttributeError
            if hdf5 converter has not read the output yet
        """
        if self.hdf5_converter.L is None:
            raise AttributeError("Read hdf5 file before, doing measurments")

        self.n_orbitals = self.hdf5_converter.L
        self.energy = self.hdf5_converter.energy

        self.make_s1(self.orbital_rdm_builder.make_one_ordm(self.hdf5_converter))
        self.make_s2(self.orbital_rdm_builder.make_two_ordm(self.hdf5_converter))
        self.make_mutual_information()

        # handle fiedler ordering
        sort_key = self.hdf5_converter.orbital_order.argsort()
        self.s1_entropy = self.s1_entropy[sort_key]
        self.s2_entropy = self.s2_entropy[sort_key][:, sort_key]
        self.mutual_information = self.mutual_information[sort_key][:, sort_key]

    def read_hdf5(self, file_name: str):
        """Read the hdf5 ouput file from qcmaquis.

        Parameters
        ----------
        file_name : str
            path to the qcmaquis hdf5 output file
        """
        self.hdf5_converter.read_hdf5(file_name)
