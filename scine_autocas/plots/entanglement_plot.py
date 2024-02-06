"""The entanglement plot class.

This module implements a class to conveniently plot an entanglement diagram.
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import lines  # type: ignore[attr-defined]
from numpy import pi


class EntanglementPlot:
    """A class to plot the entanglement diagram from s1 entropies and mutual information.

    The size of the dots is scaled by the corresponding s1 value and the thickness of the lines
    by the corresponding mututal information.

    Attributes
    ----------
    s1_color : str
        color of the s1 circles
    mutual_information_color1 : str
        color for large mututal information
    mutual_information_color2 : str
        color for medium mututal information
    mutual_information_color3 : str
        color for small mututal information
    background_color : str
        color of the background
    s1_scaled : bool
        flage to enable scaling of the s1 entropies with respect to the maxium s1 values
    alpha_s1 : float
        opacity of the s1 circles, used if 2 diagrams are plotted
    alpha_mut_inf : float
        opacity of the mutual information, used if 2 diagrams are plotted
    s1_position : float
        position of the second cirlce of s1 blobds, used if 2 diagrams are plotted
    label_offset : float
        ofset for labels in entanglement plot
    """

    __slots__ = (
        "s1_color",
        "mutual_information_color1",
        "mutual_information_color2",
        "mutual_information_color3",
        "background_color",
        "s1_scaled",
        "alpha_s1",
        "alpha_mut_inf",
        "s1_position",
        "label_offset",
    )

    def __init__(self):
        """Construct an entanglement plot object."""
        self.s1_color: str = "firebrick"
        """color of the s1 circles"""
        self.mutual_information_color1: str = "royalblue"
        """color of the strong mutual information lines"""
        self.mutual_information_color2: str = "cornflowerblue"
        """color of the medium strong mutual information lines"""
        self.mutual_information_color3: str = "lightsteelblue"
        """color of the least strong mutual information lines"""
        self.background_color: str = "white"
        """color of the background"""
        self.s1_scaled: bool = False
        """flag to enable scaling of s1 with respect to max(s1), to get circles with relative size"""
        self.alpha_s1: float = 1.0
        """opacity of s1 circles and labels"""
        self.alpha_mut_inf: float = 1.00
        """opacity of mutual information lines"""
        self.s1_position: float = 1.15
        """distance from the center for s1 circles"""
        self.label_offset: float = 0.2
        """offset for the labels"""

    def _plot_s1_circles(
        self,
        axes: matplotlib.axes._axes.Axes,  # type: ignore[name-defined]
        theta: np.ndarray,
        radii: np.ndarray,
        labels: List[str],
        area: np.ndarray,
        alpha: float,
        number_of_orbitals: int,
    ):
        """Help function to plot circles with labels.

        Parameters
        ----------
        axes : matplotlib.axes._axes.Axes
            The matplotlib axes
        theta : np.ndarray
            Position of s1 circles
        radii : np.ndarray
            Distance form the center for s1 circles
        labels : List[str]
            Orbital labels
        area : np.ndarray
            Size of each circle
        alpha : float
            opacity of the circles and labels
        number_of_orbitals : int
            number of oritals
        """
        axes.scatter(
            theta,
            radii,
            color=self.s1_color,
            s=area,
            alpha=alpha,
            zorder=1,
            edgecolor="black",
            linewidth=0.75,
        )
        # plot s1 labels
        for i in range(number_of_orbitals):
            axes.text(
                theta[i],
                (float(radii[i]) + self.label_offset),
                labels[i],
                alpha=alpha,
                ha="center",
                va="center",
                fontsize=3 + 7 * (10 / number_of_orbitals),
            )

    def _plot_s1(
        self,
        axes: matplotlib.axes._axes.Axes,  # type: ignore[name-defined]
        s1_entropy: np.ndarray,
        orbitals_index: List[int],
    ) -> np.ndarray:
        """Help function to plot the s1 circles.

        Parameters
        ----------
        axes : matplotlib.axes._axes.Axes
            The diagram axes
        s1_entropy : np.ndarray
            Single orbital entropies
        orbitals_index : List[int]
            Orbital labels for outer diagram

        Returns
        -------
        theta : np.ndarray
            Position of S1 circles
        """
        number_of_orbitals = len(s1_entropy)
        theta = np.zeros(number_of_orbitals)
        labels = []
        # s1 circles
        radii = np.zeros(number_of_orbitals)
        area = np.zeros(number_of_orbitals)
        slice_ = -2 * pi / number_of_orbitals
        for i in range(number_of_orbitals):
            # degree
            theta[i] = i * slice_ + pi / 2 + slice_ / 2
            # distance from mid point
            radii[i] = self.s1_position
            # labels
            labels.append(str(orbitals_index[i]))
            # scale so that dot area is larger
            if not self.s1_scaled:
                area[i] = s1_entropy[i] * 10 + s1_entropy[i] * 300 * (10 / number_of_orbitals)
            else:
                area[i] = s1_entropy[i] * 300 / max(s1_entropy)
        # plot s1 dots
        self._plot_s1_circles(axes, theta, radii, labels, area, self.alpha_s1, number_of_orbitals)

        return theta

    def _plot_mutual_information(
        self,
        axes: matplotlib.axes._axes.Axes,  # type: ignore[name-defined]
        mut_inf: np.ndarray,
        theta: np.ndarray,
    ):
        """Help function to plot the mutual information lines.

        Parameters
        ----------
        axes : matplotlib.axes._axes.Axes
            The diagram axes
        mut_inf : np.ndarray
            mutual information for the outer diagram
        theta : np.ndarray
            the position of the s1 circle
        """
        number_of_orbitals = mut_inf.shape[0]
        legendlines = {}
        try:
            for i in range(number_of_orbitals):
                for j in range(i, number_of_orbitals):
                    x_range = [theta[i], theta[j]]
                    y_range = [self.s1_position, self.s1_position]
                    mut_inf_line = mut_inf[i, j]
                    if mut_inf_line > 0.1:
                        line = lines.Line2D(
                            x_range,
                            y_range,
                            linewidth=5 * mut_inf_line / number_of_orbitals,
                            color=self.mutual_information_color1,
                            linestyle="-",
                            alpha=self.alpha_mut_inf,
                            label="0.1",
                            zorder=-1,
                        )
                        legendlines["0.1"] = line
                        axes.add_line(line)
                    elif mut_inf_line > 0.01:
                        line = lines.Line2D(
                            x_range,
                            y_range,
                            linewidth=20 * mut_inf_line / number_of_orbitals,
                            color=self.mutual_information_color2,
                            linestyle="--",
                            alpha=self.alpha_mut_inf,
                            label="0.01",
                            zorder=-1,
                        )
                        legendlines["0.01"] = line
                        axes.add_line(line)
                    elif mut_inf_line > 0.001:
                        line = lines.Line2D(
                            x_range,
                            y_range,
                            linewidth=30 * mut_inf_line / number_of_orbitals,
                            color=self.mutual_information_color3,
                            linestyle=":",
                            alpha=self.alpha_mut_inf,
                            label="0.001",
                            zorder=-1,
                        )
                        legendlines["0.001"] = line
                        axes.add_line(line)
        except TypeError:
            pass

    def _entanglement_plot(
        self,
        axes: matplotlib.axes._axes.Axes,  # type: ignore[name-defined]
        s1_entropy: np.ndarray,
        mut_inf: np.ndarray,
        order1: Optional[List[int]] = None,
    ):
        """Help function for the entanglement diagram within an entanglement
        diagram.

        Parameters
        ----------
        axes : matplotlib.axes._axes.Axes
            The diagram axes
        s1_entropy : np.ndarray
            single orbital entropy for the outer diagram
        mut_inf : np.ndarray
            mutual information for the outer diagram
        order1 : List[Union[int, str, float]]
            orbital labels for the outer diagram
        """
        # try:
        number_of_orbitals = len(s1_entropy)
        axes.set_xticklabels([])
        axes.set_yticklabels([])
        axes.grid(visible=False)
        axes.scatter(1.5, 1.5, color="white", zorder=-1, alpha=self.alpha_s1)

        if order1 is not None:
            orbitals_index = np.array(order1) - 1
        else:
            orbitals_index = np.arange(0, number_of_orbitals)
        orbitals_index_np = orbitals_index.tolist()

        # because numpy 1.19
        orbitals_index_np = list(orbitals_index_np)

        theta = self._plot_s1(axes, s1_entropy, orbitals_index_np)
        self._plot_mutual_information(axes, mut_inf, theta)

    def plot(
        self,
        s1_entropy: np.ndarray,
        mut_inf_1: np.ndarray,
        order_1: Optional[List[int]] = None,
    ) -> matplotlib.pyplot:  # type: ignore[valid-type]
        """Generate entanglement diagram, within an entanglement diagram.

        S1 entropies with a label appearing in both diagrams are highlighted in the outer one.

        Parameters
        ----------
        s1_entropy : np.ndarray
            single orbital entropy for the outer diagram
        mut_inf_1 : np.ndarray
            mutual information for the outer diagram
        order_1 : Optional[List[int]]
            orbital labels for the outer diagram

        Returns
        -------
        plt : matplotlib.pyplot
            the plot object
        """
        fig = plt.figure()
        axes = fig.add_subplot(polar=True, position=[0.0, 0.0, 1, 1])
        self._entanglement_plot(axes, s1_entropy, mut_inf_1, order_1)
        return plt

    def plot_in_plot(
        self,
        s1_entropy_1: np.ndarray,
        mut_inf_1: np.ndarray,
        s1_entropy_2: np.ndarray,
        mut_inf_2: np.ndarray,
        order_1: Optional[List[int]] = None,
        order_2: Optional[List[int]] = None,
    ) -> matplotlib.pyplot:  # type: ignore[valid-type]
        """Generate entanglement diagram, within an entanglement diagram.

        S1 entropies with a label appearing in both diagrams are highlighted in the outer one.

        Parameters
        ----------
        s1_entropy_1 : np.ndarray
            single orbital entropy for the outer diagram
        mut_inf_1 : np.ndarray
            mutual information for the outer diagram
        s1_entropy_2 : np.ndarray
            single orbital entropy for the inner diagram
        mut_inf_2 : np.ndarray
            mutual information for the inner diagram
        order_1 : Optional[List[int]]
            orbital labels for the outer diagram
        order_2 : Optional[List[int]]
            orbital labels for the inner diagram

        Returns
        -------
        plt : matplotlib.pyplot
            the plot object
        """
        fig = plt.figure()
        axes = fig.add_subplot(polar=True, position=[0.0, 0.0, 1, 1])
        # outer circles\
        self.alpha_mut_inf = 0.5
        self.alpha_s1 = 0.5

        self._entanglement_plot(axes, s1_entropy_2, mut_inf_2, order_2)
        # make same orbitals withput opacity
        if order_2 is not None and order_1 is not None:
            self.alpha_s1 = 1
            tmp_s1 = np.zeros(len(s1_entropy_2))
            tmp_mut_inf = np.zeros(mut_inf_2.shape)
            for i, index in enumerate(order_2):
                if index in order_1:
                    tmp_s1[i] = s1_entropy_2[i]
            self._entanglement_plot(axes, tmp_s1, tmp_mut_inf, order_2)

        self.alpha_mut_inf = 1
        self.alpha_s1 = 1
        self.s1_position = 1.0
        self.label_offset = 0.35
        # inner circle
        new_ax = fig.add_subplot(polar=True, position=[0.25, 0.25, 0.5, 0.5])

        self._entanglement_plot(new_ax, s1_entropy_1, mut_inf_1, order_1)
        return plt
