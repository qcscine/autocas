"""Module to read and store qcmaquis measurements

This module implements a class to read measurements from the qcmaquis hdf5 output file
and stores these measurements in a Datasets object.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import re
from typing import List, Optional, Tuple

import h5py
import numpy as np

from scine_autocas.interfaces.qcmaquis.qcmaquis_hdf5_dataset import Datasets


class Hdf5Converter:
    """Class to handle QCMaquis HDF5 file.

    Reads the HDF5 file, generates one orbital RDMs from corresponding
    measurements.

    Attributes
    ----------
    L : int
        the number of orbitals in the active space
    energy : float
        energy of the dmrg calculation
    symmetry : str
        symmetry used in the dmrg calculation
    orbital_order : Optional[np.ndarray]
        order of the orbitals, if changed by fiedler ordering
    results_file : str
        path to the hdf5 file
    data : Datasets
        stores all measurements from qcmaquis
    """

    __slots__ = [
        "L",
        "energy",
        "symmetry",
        "orbital_order",
        "data",
        "result_file",
    ]

    def __init__(self):
        """Construct class."""
        # L is usually used in DMRG context
        self.L: int = 0  # pylint: disable=C0103
        """number of orbitals in CAS"""
        self.energy: float = 0
        """DMRG energy"""
        self.symmetry: str = ""
        """symmetry used in the DMRG calculation"""
        self.orbital_order: Optional[np.ndarray] = None
        """orbital order from the DMRG calculation"""
        self.result_file: str = ""
        """name of the HDF5 file"""
        self.data: Datasets = None

    # def run(self):
    #    """
    #    reads HDF5 file and generates one and two orbital reduced density matrices
    #    """
    #    data = self._read_hdf5()
    #    self.make_one_ordm(data)
    #    self.make_two_ordm(data)

    def __get_labels(self, hdf5_label_vector: List[str]) -> List[Tuple[int, int]]:
        """Convert labels from HDF5 file into correct format.

        example string ( 12 ) -- ( 19 )

        Parameters
        ----------
        hdf5_label_vector : List[str]
            contains the label strings from HDF5 file

        Returns
        -------
        label_vector : List[Tuple[int]]
            contains the labels as indices
        """
        label_vector = []
        for i in hdf5_label_vector:
            # search for both integer in string
            match = re.search(
                r"^.*?\([^\d]*(\d+)[^\d]*\).*?\([^\d]*(\d+)[^\d]*\).*$", str(i)
            )
            if match is not None:
                label_vector.append((int(match.group(1)), int(match.group(2))))
        return label_vector

    def __mat_measurement(self, label_vector: List[Tuple[int, int]], value_vector: np.ndarray) -> np.ndarray:
        """Build a matrix from a label_vector and a value_vector.

        Parameters
        ----------
        label_vector : List[Tuple[int]]
            contains the converted labels from a measurement
        value_vector : np.ndarray
            contains the values from a measurement

        Returns
        -------
        matrix : np.ndarray
            the full, ordered matrix of the corresponding measurement
        """
        if len(label_vector) != len(value_vector):
            raise AssertionError("problem with vector length in 'matMeasurement")
        matrix = np.zeros((self.L, self.L))
        for i, label in enumerate(label_vector):
            matrix[label[0], label[1]] = value_vector[i]
            matrix[label[1], label[0]] = value_vector[i]
        return matrix

    def __mat_merge_transpose(
            self, label_vector_1: List[Tuple[int, int]], value_vector_1: np.ndarray,
            label_vector_2: List[Tuple[int, int]], value_vector_2: np.ndarray
    ) -> np.ndarray:
        """Build merged matrix from two measurements and corresponding label
        vectors.

        Parameters
        ----------
        label_vector_1 : List[Tuple[int]]
            converted labels for the first measurement
        value_vector_1 : np.ndarray
            contains the values from the first measurement
        label_vector_2 : List[Tuple[int]]
            converted labels for the second measurement
        value_vector_2 : np.ndarray
            contains the values from the second measurement

        Returns
        -------
        matrix : np.ndarray
            the full, ordered, merged matrix of the measurements
        """
        if (len(label_vector_1) != len(value_vector_1)) or (
            len(label_vector_2) != len(value_vector_2)
        ):
            raise AssertionError("problem with vector length in 'matMergeTranspose")
        matrix = np.zeros((self.L, self.L))
        for i, label in enumerate(label_vector_1):
            matrix[label[0], label[1]] = value_vector_1[i]
        for i, label in enumerate(label_vector_2):
            matrix[label[1], label[0]] = value_vector_2[i]
        return matrix

    def read_hdf5(self, file_name: Optional[str] = None) -> Datasets:
        """Read measurements form the QCMaquis HDF5 file.

        Parameters
        ----------
        file_name : str, optional
            name of the HDF5 file

        Returns
        -------
        data : Datasets
            contains all measurements as ordered vectors/matrices
        """
        if file_name:
            hdf5_file = h5py.File(file_name, "r")
        else:
            hdf5_file = h5py.File(self.result_file, "r")
        # read following Datasets for 1o-RDM,2o-RDM from hdf5File
        meas_dataset_list = [
            "/parameters/L",
            "/spectrum/iteration/3/results/Energy/mean/value",
            "/spectrum/results/Nup/mean/value",
            "/spectrum/results/Ndown/mean/value",
            "/spectrum/results/Nupdown/mean/value",
            "/spectrum/results/dm_up/labels",
            "/spectrum/results/dm_up/mean/value",
            "/spectrum/results/dm_down/labels",
            "/spectrum/results/dm_down/mean/value",
            "/spectrum/results/nupnup/labels",
            "/spectrum/results/nupnup/mean/value",
            "/spectrum/results/nupndown/labels",
            "/spectrum/results/nupndown/mean/value",
            "/spectrum/results/ndownnup/labels",
            "/spectrum/results/ndownnup/mean/value",
            "/spectrum/results/ndownndown/labels",
            "/spectrum/results/ndownndown/mean/value",
            "/spectrum/results/doccdocc/labels",
            "/spectrum/results/doccdocc/mean/value",
            "/spectrum/results/transfer_up_while_down/labels",
            "/spectrum/results/transfer_up_while_down/mean/value",
            "/spectrum/results/transfer_down_while_up/labels",
            "/spectrum/results/transfer_down_while_up/mean/value",
            "/spectrum/results/transfer_pair/labels",
            "/spectrum/results/transfer_pair/mean/value",
            "/spectrum/results/spinflip/labels",
            "/spectrum/results/spinflip/mean/value",
            "/spectrum/results/transfer_up_while_down_at_2/labels",
            "/spectrum/results/transfer_up_while_down_at_2/mean/value",
            "/spectrum/results/transfer_up_while_down_at_1/labels",
            "/spectrum/results/transfer_up_while_down_at_1/mean/value",
            "/spectrum/results/transfer_down_while_up_at_2/labels",
            "/spectrum/results/transfer_down_while_up_at_2/mean/value",
            "/spectrum/results/transfer_down_while_up_at_1/labels",
            "/spectrum/results/transfer_down_while_up_at_1/mean/value",
            "/spectrum/results/nupdocc/labels",
            "/spectrum/results/nupdocc/mean/value",
            "/spectrum/results/doccnup/labels",
            "/spectrum/results/doccnup/mean/value",
            "/spectrum/results/ndowndocc/labels",
            "/spectrum/results/ndowndocc/mean/value",
            "/spectrum/results/doccndown/labels",
            "/spectrum/results/doccndown/mean/value",
            "/parameters/symmetry",
            "/parameters/orbital_order",
        ]
        # go through list and extract values
        tmp_meas_list = []
        for i in meas_dataset_list:
            tmp_var = hdf5_file.get(i)
            tmp_var = np.array(tmp_var)
            tmp_meas_list.append(tmp_var)
        # extract values from list
        self.data = Datasets()
        self.L = tmp_meas_list[0]
        try:
            self.energy = tmp_meas_list[1][-1]
        except IndexError:
            self.energy = tmp_meas_list[1]
        self.data.nup = tmp_meas_list[2][0]
        self.data.ndown = tmp_meas_list[3][0]
        self.data.nupdown = tmp_meas_list[4][0]
        dmup_labels = tmp_meas_list[5]
        dmup = tmp_meas_list[6][0]
        dmdown_labels = tmp_meas_list[7]
        dmdown = tmp_meas_list[8][0]
        nupnup_labels = tmp_meas_list[9]
        nupnup = tmp_meas_list[10][0]
        nupndown_labels = tmp_meas_list[11]
        nupndown = tmp_meas_list[12][0]
        ndownnup_labels = tmp_meas_list[13]
        ndownnup = tmp_meas_list[14][0]
        ndownndown_labels = tmp_meas_list[15]
        ndownndown = tmp_meas_list[16][0]
        doccdocc_labels = tmp_meas_list[17]
        doccdocc = tmp_meas_list[18][0]
        tuwd_labels = tmp_meas_list[19]
        tuwd = tmp_meas_list[20][0]
        tdwu_labels = tmp_meas_list[21]
        tdwu = tmp_meas_list[22][0]
        transferpair_labels = tmp_meas_list[23]
        transferpair = tmp_meas_list[24][0]
        spinflip_labels = tmp_meas_list[25]
        spinflip = tmp_meas_list[26][0]
        tuwd_at2_labels = tmp_meas_list[27]
        tuwd_at2 = tmp_meas_list[28][0]
        tuwd_at1_labels = tmp_meas_list[29]
        tuwd_at1 = tmp_meas_list[30][0]
        tdwu_at2_labels = tmp_meas_list[31]
        tdwu_at2 = tmp_meas_list[32][0]
        tdwu_at1_labels = tmp_meas_list[33]
        tdwu_at1 = tmp_meas_list[34][0]
        nupdocc_labels = tmp_meas_list[35]
        nupdocc = tmp_meas_list[36][0]
        doccnup_labels = tmp_meas_list[37]
        doccnup = tmp_meas_list[38][0]
        ndowndocc_labels = tmp_meas_list[39]
        ndowndocc = tmp_meas_list[40][0]
        doccndown_labels = tmp_meas_list[41]
        doccndown = tmp_meas_list[42][0]
        self.symmetry = tmp_meas_list[43].item().decode("UTF-8")
        self.orbital_order = np.array(list(map(int, tmp_meas_list[44].item().decode("UTF-8").split(","))))
        # print("Orbital Order", self.orbital_order)
        # # u1dg symmetry not implemented yet
        # if(self.symmetry == "u1dg"):
        #   print("not implemented yet")
        #   # additional values for u1dg
        #   measDatasetList_C = [
        #     '/spectrum/results/dm/labels',
        #     '/spectrum/results/dm/mean/value',
        #     '/spectrum/results/doccdocc/labels',
        #     '/spectrum/results/doccdocc/mean/value',
        #     '/spectrum/results/N/mean/value'
        #   ]
        #   tmpMeasList_C = []
        #   for i in measDatasetList:
        #     tmpVar = hdf5_file.get(i)
        #     tmpVar = np.array(tmpVar)
        #     tmpMeasList_C.append(tmpVar)
        #   cdagc_labels        = tmpMeasList_C[0]
        #   cdagc               = tmpMeasList_C[1][0]
        #   doccdoccc_labels    = tmpMeasList_C[2]
        #   doccdoccc           = tmpMeasList_C[3]
        #   self.occ_num_vector = tmpMeasList_C[4][0]
        #   # convert label vector
        #   cdagc_labels     = self.get_labels(cdagc_labels)
        #   doccdoccc_labels = self.get_labels(doccdoccc_labels)
        #   # build matrices from labels and corresponding vector
        #   # has to be complex
        #   self.spinorCdagCSymMat    = self.matMeasurementC(cdagc_labels, cdagc)
        #   self.spinorDoccDoccSymMat = self.matMeasurementC(doccdoccc_labels, doccdoccc)
        hdf5_file.close()

        # convert labels from hdf5 to new format
        dmup_labels = self.__get_labels(dmup_labels)
        dmdown_labels = self.__get_labels(dmdown_labels)
        nupnup_labels = self.__get_labels(nupnup_labels)
        nupndown_labels = self.__get_labels(nupndown_labels)
        ndownnup_labels = self.__get_labels(ndownnup_labels)
        ndownndown_labels = self.__get_labels(ndownndown_labels)
        doccdocc_labels = self.__get_labels(doccdocc_labels)
        tuwd_labels = self.__get_labels(tuwd_labels)
        tdwu_labels = self.__get_labels(tdwu_labels)
        transferpair_labels = self.__get_labels(transferpair_labels)
        spinflip_labels = self.__get_labels(spinflip_labels)
        tuwd_at2_labels = self.__get_labels(tuwd_at2_labels)
        tuwd_at1_labels = self.__get_labels(tuwd_at1_labels)
        tdwu_at2_labels = self.__get_labels(tdwu_at2_labels)
        tdwu_at1_labels = self.__get_labels(tdwu_at1_labels)
        nupdocc_labels = self.__get_labels(nupdocc_labels)
        doccnup_labels = self.__get_labels(doccnup_labels)
        ndowndocc_labels = self.__get_labels(ndowndocc_labels)
        doccndown_labels = self.__get_labels(doccndown_labels)
        # build matrices from labels and corresponding vectors
        self.data.dmUpSym = self.__mat_measurement(dmup_labels, dmup)
        self.data.dmDownSym = self.__mat_measurement(dmdown_labels, dmdown)
        self.data.nupNupSym = self.__mat_measurement(nupnup_labels, nupnup)
        self.data.nupNdownSym = self.__mat_measurement(nupndown_labels, nupndown)
        self.data.ndownNupSym = self.__mat_measurement(ndownnup_labels, ndownnup)
        self.data.ndownNdownSym = self.__mat_measurement(ndownndown_labels, ndownndown)
        self.data.doccDoccSym = self.__mat_measurement(doccdocc_labels, doccdocc)
        self.data.transferUpWhileDownSym = self.__mat_measurement(tuwd_labels, tuwd)
        self.data.transferDownWhileUpSym = self.__mat_measurement(tdwu_labels, tdwu)
        self.data.transferPairSym = self.__mat_measurement(transferpair_labels, transferpair)
        self.data.spinflipSym = self.__mat_measurement(spinflip_labels, spinflip)
        self.data.transUpDown2Asym = self.__mat_merge_transpose(tuwd_at2_labels, tuwd_at2, tuwd_at1_labels, tuwd_at1)
        self.data.transUpDown1Asym = self.__mat_merge_transpose(tuwd_at1_labels, tuwd_at1, tuwd_at2_labels, tuwd_at2)
        self.data.transDownUp2Asym = self.__mat_merge_transpose(tdwu_at2_labels, tdwu_at2, tdwu_at1_labels, tdwu_at1)
        self.data.transDownUp1Asym = self.__mat_merge_transpose(tdwu_at1_labels, tdwu_at1, tdwu_at2_labels, tdwu_at2)
        self.data.nupDoccAsym = self.__mat_merge_transpose(nupdocc_labels, nupdocc, doccnup_labels, doccnup)
        self.data.doccNupAsym = self.__mat_merge_transpose(doccnup_labels, doccnup, nupdocc_labels, nupdocc)
        self.data.ndownDoccAsym = self.__mat_merge_transpose(ndowndocc_labels, ndowndocc, doccndown_labels, doccndown)
        self.data.doccNdownAsym = self.__mat_merge_transpose(doccndown_labels, doccndown, ndowndocc_labels, ndowndocc)
        return self.data
