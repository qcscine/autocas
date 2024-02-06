"""Module to provide functionalities to evaluate RDMs from qcmaquis calculations.

This module implements a class to build orbital RDM from qcmaquis output.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import numpy as np

from scine_autocas.interfaces.qcmaquis.qcmaquis_hdf5_utils import Hdf5Converter


class OrbitalRDMBuilder:
    """Calculate one and two orbital reduced density matrices from qcmaquis output.

    The orbital rdms are build in their completeness, hence the one orbital rdm is
    a tensor with the dimensions of (n_orbitals x 4) and the two orbital rdm a tensor
    with the dimensions of (n_orbitals x n_orbitals x 16 x 16).

    Attributes
    ----------
    one_rdm : np.ndarray
        stores the one orbital rdm
    two_rdm : np.ndarray
        stores the two orbital rdm
    """

    __slots__ = ("one_ordm", "two_ordm")

    def __init__(self):
        """Init."""
        self.one_ordm: np.ndarray
        """stores the one orbital rdm"""
        self.two_ordm: np.ndarray
        """stroes the two orbital rdm"""

    def _check_hdf5_converter(self, hdf5_converter: Hdf5Converter):
        if not hdf5_converter.data:
            raise AttributeError(
                "read qcmaquis result file before calculating one orbital reduced density matrix"
            )
        if not hdf5_converter.L:
            raise AttributeError(
                "read qcmaquis result file before calculating one orbital reduced density matrix"
            )

    def make_one_ordm(self, hdf5_converter: Hdf5Converter):
        """
        Evaluate one orbital RDM from qcmaquis measurements.

        Parameters
        ----------
        hdf5_converter : Hdf5Converter
            contains all measurements
        """
        self._check_hdf5_converter(hdf5_converter)
        # not implemented yet
        # if self.symmetry == "u1dg":
        #   self.one_ordm = np.zeros((self.L,4))
        #   for i in range(self.L):
        #     self.one_ordm[i, 0] = 1 - self.occ_num_vector[i]
        #     self.one_ordm[i, 1] = self.occ_num_vector[i]
        # else:
        self.one_ordm = np.zeros((hdf5_converter.L, 4))
        for i in range(hdf5_converter.L):
            self.one_ordm[i, 0] = (
                hdf5_converter.data.nup[i] - hdf5_converter.data.nupdown[i]
            )
            self.one_ordm[i, 1] = (
                hdf5_converter.data.ndown[i] - hdf5_converter.data.nupdown[i]
            )
            self.one_ordm[i, 2] = (
                1
                - hdf5_converter.data.nup[i]
                - hdf5_converter.data.ndown[i]
                + hdf5_converter.data.nupdown[i]
            )
            self.one_ordm[i, 3] = hdf5_converter.data.nupdown[i]
        return self.one_ordm

    def make_two_ordm(self, hdf5_converter: Hdf5Converter):
        """
        Evaluate two orbital RDM from qcmaquis measurements.

        Parameters
        ----------
        hdf5_converter : Hdf5Converter
            contains all measurements
        """
        self._check_hdf5_converter(hdf5_converter)
        self.two_ordm = np.zeros((hdf5_converter.L, hdf5_converter.L, 16, 16))

        # snake case convetion, since p, q are okay names here
        # p, q is the orbital index from whole space
        # pylint: disable=C0103
        for p in range(hdf5_converter.L):
            for q in range(p + 1, hdf5_converter.L):
                self.two_ordm[p, q, 0, 0] = (
                    1.0
                    + hdf5_converter.data.nupdown[p]
                    + hdf5_converter.data.nupdown[q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                    - hdf5_converter.data.ndown[p]
                    - hdf5_converter.data.ndownDoccAsym[p, q]
                    - hdf5_converter.data.ndown[q]
                    - hdf5_converter.data.doccNdownAsym[p, q]
                    + hdf5_converter.data.ndownNdownSym[p, q]
                    - hdf5_converter.data.nup[p]
                    - hdf5_converter.data.nupDoccAsym[p, q]
                    + hdf5_converter.data.nupNdownSym[p, q]
                    - hdf5_converter.data.nup[q]
                    - hdf5_converter.data.doccNupAsym[p, q]
                    + hdf5_converter.data.ndownNupSym[p, q]
                    + hdf5_converter.data.nupNupSym[p, q]
                )
                self.two_ordm[p, q, 1, 1] = (
                    -hdf5_converter.data.nupdown[p]
                    - hdf5_converter.data.doccDoccSym[p, q]
                    + hdf5_converter.data.ndown[p]
                    + hdf5_converter.data.ndownDoccAsym[p, q]
                    + hdf5_converter.data.doccNdownAsym[p, q]
                    - hdf5_converter.data.ndownNdownSym[p, q]
                    + hdf5_converter.data.doccNupAsym[p, q]
                    - hdf5_converter.data.ndownNupSym[p, q]
                )
                self.two_ordm[p, q, 2, 2] = (
                    -hdf5_converter.data.nupdown[p]
                    - hdf5_converter.data.doccDoccSym[p, q]
                    + hdf5_converter.data.doccNdownAsym[p, q]
                    + hdf5_converter.data.nup[p]
                    + hdf5_converter.data.nupDoccAsym[p, q]
                    - hdf5_converter.data.nupNdownSym[p, q]
                    + hdf5_converter.data.doccNupAsym[p, q]
                    - hdf5_converter.data.nupNupSym[p, q]
                )
                self.two_ordm[p, q, 3, 3] = (
                    hdf5_converter.data.nupdown[p]
                    - hdf5_converter.data.doccNdownAsym[p, q]
                    - hdf5_converter.data.doccNupAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 4, 4] = (
                    -hdf5_converter.data.nupdown[q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                    + hdf5_converter.data.ndownDoccAsym[p, q]
                    + hdf5_converter.data.ndown[q]
                    + hdf5_converter.data.doccNdownAsym[p, q]
                    - hdf5_converter.data.ndownNdownSym[p, q]
                    + hdf5_converter.data.nupDoccAsym[p, q]
                    - hdf5_converter.data.nupNdownSym[p, q]
                )
                self.two_ordm[p, q, 5, 5] = (
                    hdf5_converter.data.ndownNdownSym[p, q]
                    - hdf5_converter.data.ndownDoccAsym[p, q]
                    - hdf5_converter.data.doccNdownAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 6, 6] = (
                    hdf5_converter.data.nupNdownSym[p, q]
                    - hdf5_converter.data.doccNdownAsym[p, q]
                    - hdf5_converter.data.nupDoccAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 7, 7] = (
                    hdf5_converter.data.doccNdownAsym[p, q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 8, 8] = (
                    -hdf5_converter.data.nupdown[q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                    + hdf5_converter.data.ndownDoccAsym[p, q]
                    + hdf5_converter.data.nupDoccAsym[p, q]
                    + hdf5_converter.data.nup[q]
                    + hdf5_converter.data.doccNupAsym[p, q]
                    - hdf5_converter.data.ndownNupSym[p, q]
                    - hdf5_converter.data.nupNupSym[p, q]
                )
                self.two_ordm[p, q, 9, 9] = (
                    hdf5_converter.data.ndownNupSym[p, q]
                    - hdf5_converter.data.ndownDoccAsym[p, q]
                    - hdf5_converter.data.doccNupAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 10, 10] = (
                    hdf5_converter.data.nupNupSym[p, q]
                    - hdf5_converter.data.nupDoccAsym[p, q]
                    - hdf5_converter.data.doccNupAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 11, 11] = (
                    hdf5_converter.data.doccNupAsym[p, q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 12, 12] = (
                    hdf5_converter.data.nupdown[q]
                    - hdf5_converter.data.nupDoccAsym[p, q]
                    - hdf5_converter.data.ndownDoccAsym[p, q]
                    + hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 13, 13] = (
                    hdf5_converter.data.ndownDoccAsym[p, q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 14, 14] = (
                    hdf5_converter.data.nupDoccAsym[p, q]
                    - hdf5_converter.data.doccDoccSym[p, q]
                )
                self.two_ordm[p, q, 15, 15] = hdf5_converter.data.doccDoccSym[p, q]
                self.two_ordm[p, q, 1, 4] = (
                    hdf5_converter.data.dmDownSym[p, q]
                    - hdf5_converter.data.transDownUp1Asym[p, q]
                    - hdf5_converter.data.transDownUp2Asym[p, q]
                    + hdf5_converter.data.transferDownWhileUpSym[p, q]
                )
                self.two_ordm[p, q, 4, 1] = self.two_ordm[p, q, 1, 4]
                self.two_ordm[p, q, 2, 8] = (
                    hdf5_converter.data.dmUpSym[p, q]
                    - hdf5_converter.data.transUpDown1Asym[p, q]
                    - hdf5_converter.data.transUpDown2Asym[p, q]
                    + hdf5_converter.data.transferUpWhileDownSym[p, q]
                )
                self.two_ordm[p, q, 8, 2] = self.two_ordm[p, q, 2, 8]
                self.two_ordm[p, q, 3, 6] = (
                    hdf5_converter.data.transDownUp1Asym[p, q]
                    - hdf5_converter.data.transferDownWhileUpSym[p, q]
                )
                self.two_ordm[p, q, 6, 3] = self.two_ordm[p, q, 3, 6]
                self.two_ordm[p, q, 3, 9] = (
                    -hdf5_converter.data.transUpDown1Asym[p, q]
                    + hdf5_converter.data.transferUpWhileDownSym[p, q]
                )
                self.two_ordm[p, q, 9, 3] = self.two_ordm[p, q, 3, 9]
                self.two_ordm[p, q, 6, 9] = hdf5_converter.data.spinflipSym[p, q]
                self.two_ordm[p, q, 9, 6] = self.two_ordm[p, q, 6, 9]
                self.two_ordm[p, q, 3, 12] = hdf5_converter.data.transferPairSym[p, q]
                self.two_ordm[p, q, 12, 3] = self.two_ordm[p, q, 3, 12]
                self.two_ordm[p, q, 6, 12] = -(
                    -hdf5_converter.data.transUpDown2Asym[p, q]
                    + hdf5_converter.data.transferUpWhileDownSym[p, q]
                )
                self.two_ordm[p, q, 12, 6] = self.two_ordm[p, q, 6, 12]
                self.two_ordm[p, q, 9, 12] = -(
                    hdf5_converter.data.transDownUp2Asym[p, q]
                    - hdf5_converter.data.transferDownWhileUpSym[p, q]
                )
                self.two_ordm[p, q, 12, 9] = self.two_ordm[p, q, 9, 12]
                self.two_ordm[
                    p, q, 7, 13
                ] = -hdf5_converter.data.transferUpWhileDownSym[p, q]
                self.two_ordm[p, q, 13, 7] = self.two_ordm[p, q, 7, 13]
                self.two_ordm[
                    p, q, 11, 14
                ] = -hdf5_converter.data.transferDownWhileUpSym[p, q]
                self.two_ordm[p, q, 14, 11] = self.two_ordm[p, q, 11, 14]
        return self.two_ordm
        # pylint: enable=C0103
