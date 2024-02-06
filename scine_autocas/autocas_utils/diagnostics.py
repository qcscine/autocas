"""Provide diagnostics for autoCAS.

This module implements the Diagnostics class, which stores entropies and
threshold values for the determination of the active space.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


from typing import Any, Dict, Optional

import numpy as np


class Diagnostics:
    """A Class to store the entropies and threshold values to evaluate
    diagnostics from it.

    Currently there is only a check for single reference and zs1 implemented, hence the
    only values are thresholds which are stored.

    Attributes
    ----------
    s1_entropy : np.ndarray
        the single orbital entropy matrix
    s2_entropy : np.ndarray
        the two orbtial entropy tensor
    mutual_information : np.ndarray
        the mutual information tensor
    single_reference_threshold : float, default = 0,14
        if the largest element in s1 is smaller, the system is considered single reference
    weak_correlation threshold : float, default = 0.02
        if no active space is found, orbitals with a lower single orbital entropy are remove.


    Notes
    -----
    For further information see
    [1] C. J. Stein and M. Reiher J. Comput. Chem. 40, 2216-2226 (2019).
    [2] C. J. Stein and M. Reiher J. Chem. Theory Comput. 12, 1760-1771 (2016).
    [3] C. J. Stein, V. von Burg and M. Reiher J. Chem. Theory Comput. 12, 3764-3773 (2016).
    [4] C. J. Stein and M. Reiher Chimia 71, 170-176 (2017).
    """

    __slots__ = (
        "s1_entropy",
        "s2_entropy",
        "mutual_information",
        "single_reference_threshold",
        "weak_correlation_threshold",
    )

    def __init__(self, settings_dict: Optional[Dict[str, Any]] = None):
        """Contruct a Diagnostics object.

        If a settings_dict is provided is provided with attributes in this class, these
        attributes will be overwritten.

        Parameters
        ----------
        settings_dict : Dict[str, Any], optional
            a dict, usually provided by the input_handler, which stores attributes and corresponding values

        See Also
        --------
        settings_dict : InputHandler
        """
        self.s1_entropy: np.ndarray
        """single orbital entropy"""
        self.s2_entropy: np.ndarray
        """two orbital entropy"""
        self.mutual_information: np.ndarray
        """mutual information"""
        self.single_reference_threshold: float = 0.14
        """if max(s1) is lower then this threshold the system will be defined as single reference"""
        self.weak_correlation_threshold: float = 0.02
        """if an orbital has a lower s1 than max(s1) * weak_correlation_threshold, it is excluded"""
        if settings_dict is not None:
            for key in settings_dict:
                if hasattr(self, key):
                    setattr(self, key, settings_dict[key])

    def zs1(self, s1_entropy: np.ndarray) -> float:
        """Evaluate the Zs(1) diagnostics for a set of s1 values.

        Zs(1) indicates a MC character if Zs(1) > 0.2.
        If 0.1 < Zs(1) < 0.2 than a SR method can be considered, since it is usually cheaper to
        use it.

        Parameters
        ----------
        s1_entropy : np.ndarray
            set of s1 values to evaluate zs1 from.

        Returns
        -------
        zs1 : float
            The Zs(1) diagnostics with respect to the set of s1 values.
        """
        n_orbitals = len(s1_entropy)
        zs1 = 1 / (n_orbitals * np.log(4)) * sum(s1_entropy)
        return zs1

    def is_single_reference(self) -> bool:
        """Check if the system is weak correlated.

        Maning the maximum s1 value is lower than 10% of the theoretical maximum, e.g. 0.14.

        Returns
        -------
        bool
            True if the system is single reference, else False meaning multi reference
        """
        if max(self.s1_entropy) < self.single_reference_threshold:
            return True
        return False
