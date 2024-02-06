"""The threshold diagram class.

This module implements a class to plot a threshold diagram.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import Any

import matplotlib.pyplot as plt
import numpy as np


class ThresholdPlot:
    """A Class to plot the threshold diagram for s1 values.

    A threshold diagram provides a graphical way to check back the autocas
    plateau algorithm.

    Attributes
    ----------
    plateau_elements : int, default = 10
        number of steps required to form a plateau
    """

    __slots__ = ["plateau_elements"]

    def __init__(self):
        """Construct a treshholdplot object."""
        self.plateau_elements: int = 10
        """defines the number of elements necessary to count as a plateau"""

    def set_number_of_platueau_elements(self, n_elements: int):
        """
        Set number of steps required for a plateau.

        Parameters
        ----------
        n_elements : int
            new number of n_elements
        """
        self.plateau_elements = n_elements

    # currently not in use
    # def plot(self, s1=None, mutinf = None):
    #   if(s1 != None):
    #     plot_s1(s1)
    #   elif(mutinf != None):
    #     plot_mutual_information(mutinf)
    #   else:
    #     print("please provide either 's1', or 'mutinf' as argument")

    def plot(self, s1_entropy: np.ndarray) -> Any:
        """
        Plot the threshold diagram from the single orbital entropy.

        Parameters
        ----------
        s1 : np.ndarray
            single orbital entropy

        Returns
        -------
        plt : matplotlib.pyplot
            the matplotlib object
        """
        max_s1 = max(s1_entropy)
        number_of_orbitals = len(s1_entropy)
        orbitals_index = np.arange(1, number_of_orbitals + 1)
        thresholds_list = np.zeros((number_of_orbitals))
        for i in range(number_of_orbitals):
            thresholds_list[i] = s1_entropy[i] / max_s1
        # sort arrays decreasing
        sortkey = np.argsort(-thresholds_list)
        # because numpy 1.19
        sortkey = sortkey.astype(int, copy=False)
        thresholds_list = thresholds_list[sortkey]
        orbitals_index = orbitals_index[sortkey]
        x_values = np.arange(0, 1.01, 0.01)
        y_values = np.zeros((101))
        for i, x_value in enumerate(x_values):
            y_values[i] = sum(thresholds_list > x_value)
        # get plateau
        tmp_val = y_values[0]
        plateau_vector = []
        plateau_index = []
        thresh_count = 1
        for i, y_value in enumerate(y_values):
            if y_value == tmp_val:
                thresh_count = thresh_count + 1
            else:
                thresh_count = 1
                tmp_val = y_value
            if thresh_count >= 10:
                plateau_vector.append(y_value)
                plateau_index.append(i / 100)
        # plot threshold diagram
        plt.figure()
        plt.plot(x_values, y_values, marker="o", markersize=3, ls="", color="#203f78")
        # plt.plot(plateau_index, plateau_vector, ls="", marker="x")
        plt.ylabel("# selected Orbitals", fontsize=16)
        plt.xlabel("threshold in % of largest element", fontsize=16)
        if number_of_orbitals > 40:
            plt.yticks(np.arange(0, number_of_orbitals + 1, 5), fontsize=16)
        else:
            plt.yticks(np.arange(0, number_of_orbitals + 1, 1), fontsize=16)
        plt.xticks(np.arange(0, 1.1, 0.1), fontsize=16)
        return plt
