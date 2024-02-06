# pylint: disable=C0114, C0115, C0116
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

import unittest

import matplotlib.pyplot as plt
import numpy as np

from scine_autocas.plots.entanglement_plot import EntanglementPlot
from scine_autocas.plots.threshold_plot import ThresholdPlot


class TestPlots(unittest.TestCase):
    def test_entanglement_plot(self):
        plot = EntanglementPlot()
        # fmt: off
        s1_entropy = np.array(
            [0.8676371, 1.20708225, 1.20708227, 1.20712177, 1.20712175, 0.86743313]
        )
        mutual_information = np.array(
            [
                [0.00000000, 0.05103154, 0.05103153, 0.02192151, 0.02192151, 0.55965736],
                [0.05103154, 0.00000000, 0.14465414, 0.00000000, 0.00000000, 0.02185112],
                [0.05103153, 0.14465414, 0.00000000, 0.00000000, 0.00000000, 0.02185112],
                [0.02192151, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.05100015],
                [0.02192151, 0.00000000, 0.00550000, 0.00000000, 0.00000000, 0.05100015],
                [0.55965736, 0.02185112, 0.02185112, 0.05100015, 0.05100015, 0.00000100],
            ]
        )
        # fmt: on
        plt.close()
        plot.plot(s1_entropy, mutual_information)

        plot_in_plot = EntanglementPlot()
        plt.close()
        plt_test2 = plot_in_plot.plot_in_plot(s1_entropy, mutual_information, s1_entropy, mutual_information, order_1=[
                                              1, 4, 5, 6, 7, 8], order_2=[1, 2, 3, 4, 5, 6])
        _ = plt_test2
        # plt_test2.show()

    def test_plot_in_plot(self):
        plot = EntanglementPlot()
        # fmt: off
        s1_entropy = np.array(
            [0.8676371, 1.20708225, 1.20708227, 1.20712177, 1.20712175, 0.86743313]
        )
        mutual_information = np.array(
            [
                [0.00000000, 0.05103154, 0.05103153, 0.02192151, 0.02192151, 0.55965736],
                [0.05103154, 0.00000000, 0.14465414, 0.00550000, 0.00000000, 0.02185112],
                [0.05103153, 0.14465414, 0.00001000, 0.00001000, 0.00001000, 0.02185112],
                [0.02192151, 0.00001000, 0.00001000, 0.00001000, 0.00001000, 0.05100015],
                [0.02192151, 0.00001000, 0.00001000, 0.00001100, 0.00001000, 0.05100015],
                [0.55965736, 0.02185112, 0.02185112, 0.05100015, 0.05100015, 0.00000100],
            ]
        )
        # fmt: on
        plt.close()
        plot.plot(s1_entropy, mutual_information)
        plot_in_plot = EntanglementPlot()
        plt.close()
        plot_in_plot.plot_in_plot(
            s1_entropy,
            mutual_information,
            s1_entropy[2:],
            mutual_information[2:, 2:],
            order_1=[2, 3, 4, 5, 6, 7],
            order_2=[4, 5, 6, 7],
        )
        plot_in_plot.plot_in_plot(s1_entropy, mutual_information, s1_entropy, mutual_information)
        # new_plot.show()

    def test_threshold_plots(self):
        plot = ThresholdPlot()
        # fmt: off
        s1_entropy = np.array(
            [0.8676371, 1.20708225, 1.20708227, 1.20712177, 1.20712175, 0.86743313]
        )
        # fmt: on
        plt.close()
        plot.plot(s1_entropy)

        plot.set_number_of_platueau_elements(20)
        self.assertEqual(plot.plateau_elements, 20)
        more_s1 = np.random.rand(50)
        plot.plot(more_s1)


if __name__ == "__main__":
    np.set_printoptions(edgeitems=10, linewidth=180)
    unittest.main()
