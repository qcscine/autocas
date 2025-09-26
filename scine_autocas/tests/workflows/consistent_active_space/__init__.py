# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import os


def create_n2_dissociation():
    xyz_names = [os.path.join(os.getcwd(), "n2_" + str(i) + ".xyz") for i in range(3)]
    bond_lengths = [1.1, 2.0, 4.0]
    for xyz_name, bond_length in zip(xyz_names, bond_lengths):
        n2_xyz = open(xyz_name, "w")
        n2_xyz.write("2\n\n")
        n2_xyz.write("N 0 0 0\n")
        n2_xyz.write(f"N 0 0 {bond_length}\n")
        n2_xyz.close()
    return xyz_names
