# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import numpy as np
import os
import h5py

class CQ_HDF5_Reader:

    def __init__(self, fName:str = None ):
        self.file_name: str = fName
        """
        str
            Name of the HDF5 file to read
        """
        self.energy: float = 0.
        """
        float
            Energy read from hdf5 file
        """
        self.num_mo: int = 0
        """
        int
            Number of molecular orbitals
        """

    def read_hdf5(self):
        # Open HDF5 File
        if not self.file_name:
            raise Exception("In the CQ hdf5 reader, the file name was never initialized\n")
        f = h5py.File(self.file_name,'r')

        # Read DMRG/CAS Energy from file
        try:
            data = f['/MCWFN/STATE_ENERGY']
            self.energy = data[0]
        except:
            try:
                data = f['/SCF/TOTAL_ENERGY']
                self.energy = data[0]
            except:
                raise RuntimeError("Could not find energy in ChronusQ binary file (" + self.file_name + ")")

        # Read number of  molecular orbitals
        try:
            data = f['/SCF/MO1']
            self.num_mo = int(data.shape[1])
        except:
            raise RuntimeError("Could not determine the number of MOs in ChronusQ binary file("+self.file_name+")")
        
