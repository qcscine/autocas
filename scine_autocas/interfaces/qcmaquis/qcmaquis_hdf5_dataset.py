"""Module to related objects.

This module stores all required measurements from a qcmaquis calculation.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import numpy as np

# pylint: disable=too-many-instance-attributes
# pylint: disable=C0103
# This removes the snake case convention.
# The variable names here are the same as in qcmaquis and the corresponding ouput.


class Datasets:
    """
    Class to store the measurements from the QCMaquis HDF5 file.

    All class member names correspond are similar to the names in the HDF5 file.

    Attributes
    ----------
    self.nup : np.ndarray
        nup measurement for orbital RDM
    self.ndown : np.ndarray
        ndown measurement for orbital RDM
    self.nupdown : np.ndarray
        nupdown measurement for orbital RDM
    self.dmUpSym : np.ndarray
        dmUpSym measurement for orbital RDM
    self.dmDownSym : np.ndarray
        dmDownSym measurement for orbital RDM
    self.nupNupSym : np.ndarray
        nupNupSym measurement for orbital RDM
    self.nupNdownSym : np.ndarray
        nupNdownSym measurement for orbital RDM
    self.ndownNupSym : np.ndarray
        ndownNupSym measurement for orbital RDM
    self.ndownNdownSym : np.ndarray
        ndownNdownSym measurement for orbital RDM
    self.doccDoccSym : np.ndarray
        doccDoccSym measurement for orbital RDM
    self.transferUpWhileDownSym : np.ndarray
        transferUpWhileDownSym measurement for orbital RDM
    self.transferDownWhileUpSym : np.ndarray
        transferDownWhileUpSym measurement for orbital RDM
    self.transferPairSym : np.ndarray
        transferPairSym measurement for orbital RDM
    self.spinflipSym : np.ndarray
        spinflipSym measurement for orbital RDM
    self.transUpDown2Asym : np.ndarray
        transUpDown2Asym measurement for orbital RDM
    self.transUpDown1Asym : np.ndarray
        transUpDown1Asym measurement for orbital RDM
    self.transDownUp2Asym : np.ndarray
        transDownUp2Asym measurement for orbital RDM
    self.transDownUp1Asym : np.ndarray
        transDownUp1Asym measurement for orbital RDM
    self.nupDoccAsym : np.ndarray
        nupDoccAsym measurement for orbital RDM
    self.doccNupAsym : np.ndarray
        doccNupAsym measurement for orbital RDM
    self.ndownDoccAsym : np.ndarray
        ndownDoccAsym measurement for orbital RDM
    self.doccNdownAsym : np.ndarray
        doccNdownAsym measurement for orbital RDM
    """

    __slots__ = [
        "nup",
        "ndown",
        "nupdown",
        "dmUpSym",
        "dmDownSym",
        "nupNupSym",
        "nupNdownSym",
        "ndownNupSym",
        "ndownNdownSym",
        "doccDoccSym",
        "transferUpWhileDownSym",
        "transferDownWhileUpSym",
        "transferPairSym",
        "spinflipSym",
        "transUpDown2Asym",
        "transUpDown1Asym",
        "transDownUp2Asym",
        "transDownUp1Asym",
        "nupDoccAsym",
        "doccNupAsym",
        "ndownDoccAsym",
        "doccNdownAsym",
        "occ_num_vector",
        "spinorCdagCSymMat",
        "spinorDoccDoccSymMat",
    ]

    def __init__(self):
        """Init."""
        self.nup: np.ndarray
        """nup measurement for orbital RDM"""
        self.ndown: np.ndarray
        """ndown measurement for orbital RDM"""
        self.nupdown: np.ndarray
        """nupdown measurement for orbital RDM"""
        self.dmUpSym: np.ndarray
        """dmUpSym measurement for orbital RDM"""
        self.dmDownSym: np.ndarray
        """dmDownSym measurement for orbital RDM"""
        self.nupNupSym: np.ndarray
        """nupNupSym measurement for orbital RDM"""
        self.nupNdownSym: np.ndarray
        """nupNdownSym measurement for orbital RDM"""
        self.ndownNupSym: np.ndarray
        """ndownNupSym measurement for orbital RDM"""
        self.ndownNdownSym: np.ndarray
        """ndownNdownSym measurement for orbital RDM"""
        self.doccDoccSym: np.ndarray
        """doccDoccSym measurement for orbital RDM"""
        self.transferUpWhileDownSym: np.ndarray
        """transferUpWhileDownSym measurement for orbital RDM"""
        self.transferDownWhileUpSym: np.ndarray
        """transferDownWhileUpSym measurement for orbital RDM"""
        self.transferPairSym: np.ndarray
        """transferPairSym measurement for orbital RDM"""
        self.spinflipSym: np.ndarray
        """spinflipSym measurement for orbital RDM"""
        self.transUpDown2Asym: np.ndarray
        """transUpDown2Asym measurement for orbital RDM"""
        self.transUpDown1Asym: np.ndarray
        """transUpDown1Asym measurement for orbital RDM"""
        self.transDownUp2Asym: np.ndarray
        """transDownUp2Asym measurement for orbital RDM"""
        self.transDownUp1Asym: np.ndarray
        """transDownUp1Asym measurement for orbital RDM"""
        self.nupDoccAsym: np.ndarray
        """nupDoccAsym measurement for orbital RDM"""
        self.doccNupAsym: np.ndarray
        """doccNupAsym measurement for orbital RDM"""
        self.ndownDoccAsym: np.ndarray
        """ndownDoccAsym measurement for orbital RDM"""
        self.doccNdownAsym: np.ndarray
        """doccNdownAsym measurement for orbital RDM"""

# pylint: enable=too-many-instance-attributes
# pylint: enable=C0103
