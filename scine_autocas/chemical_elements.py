# it is okay to have many lines in this module
# pylint: disable=C0302
"""Module to provide all information on chemical elements.

This module implements the Elements class, which provides all chemical
element specific information required by autocas. The ElementDict is a
type, to satisfy mypy and other linters.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from typing import List

from mypy_extensions import TypedDict

# new type for mypy
ElementDict: TypedDict = TypedDict(
    "ElementDict",
    {
        "name": str,
        "atomic number": int,
        "number of core orbitals": int,
        "number of valence orbitals": int,
        "core orbitals": List[str],
        "valence orbitals": List[str],
    },
)


class Elements:
    """List of dictionaries, which contains information for elements H - Xe.

    Is used to create inital active space for valence electrons and orbitals.

    Attributes
    ----------
    double_d_shell : bool
        flag to add 4d orbitals for 3rd row transition metals
    double_d_shell_elements : List[str]
        defines which elements can have 4d orbitals through double d shell effect
    element_list : List[ElementDict]
        defines name, atomic number, number of core orbitals, number of valence
        orbitals, core orbitals, and valence orbitals for each element.
    """

    __slots__ = (
        "double_d_shell",
        "double_d_shell_elements",
        "element_list",
    )

    def __init__(self, double_d_shell: bool = False):
        """Construct the element orbject.

        Parameters
        ----------
        double_d_shell : bool, optional
            set the double_d_shell attribute
        """
        self.double_d_shell: bool = double_d_shell
        """add 4d orbitals to valence space"""
        # fmt: off
        self.double_d_shell_elements: List[str] = [
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        ]
        # fmt: on
        """elements suited for double d shell option"""
        # fmt: off
        self.element_list: List[ElementDict] = [
            {
                "name": "H",
                "atomic number": 1,
                "number of core orbitals": 0,
                "number of valence orbitals": 1,
                "core orbitals": [],
                "valence orbitals": ["1s"],
            },
            {
                "name": "He",
                "atomic number": 2,
                "number of core orbitals": 0,
                "number of valence orbitals": 1,
                "core orbitals": [],
                "valence orbitals": ["1s"],
            },
            {
                "name": "Li",
                "atomic number": 3,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "Be",
                "atomic number": 4,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "B",
                "atomic number": 5,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "C",
                "atomic number": 6,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "N",
                "atomic number": 7,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "O",
                "atomic number": 8,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "F",
                "atomic number": 9,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "Ne",
                "atomic number": 10,
                "number of core orbitals": 1,
                "number of valence orbitals": 4,
                "core orbitals": ["1s"],
                "valence orbitals": ["2s", "2p"],
            },
            {
                "name": "Na",
                "atomic number": 11,
                "number of core orbitals": 5,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s"],
            },
            {
                "name": "Mg",
                "atomic number": 12,
                "number of core orbitals": 5,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s"],
            },
            {
                "name": "Al",
                "atomic number": 13,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "Si",
                "atomic number": 14,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "P",
                "atomic number": 15,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "S",
                "atomic number": 16,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "Cl",
                "atomic number": 17,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "Ar",
                "atomic number": 18,
                "number of core orbitals": 5,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p"],
                "valence orbitals": ["3s", "3p"],
            },
            {
                "name": "K",
                "atomic number": 19,
                "number of core orbitals": 9,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s"],
            },
            {
                "name": "Ca",
                "atomic number": 20,
                "number of core orbitals": 9,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s"],
            },
            {
                "name": "Sc",
                "atomic number": 21,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Ti",
                "atomic number": 22,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "V",
                "atomic number": 23,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Cr",
                "atomic number": 24,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Mn",
                "atomic number": 25,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Fe",
                "atomic number": 26,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            # example for double d-shell
            # {
            #    "name": "Fe",
            #    "atomic number": 26,
            #    "number of core orbitals": 9,
            #    "number of valence orbitals": 14,
            #    "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
            #    "valence orbitals": ["4s", "3d", "4p", "4d"],
            # },
            {
                "name": "Co",
                "atomic number": 27,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Ni",
                "atomic number": 28,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Cu",
                "atomic number": 29,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["3d", "4s", "4p"],
            },
            {
                "name": "Zn",
                "atomic number": 30,
                "number of core orbitals": 9,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p"],
                "valence orbitals": ["4s", "3d", "4p"],
            },
            {
                "name": "Ga",
                "atomic number": 31,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "Ge",
                "atomic number": 32,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "As",
                "atomic number": 33,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "Se",
                "atomic number": 34,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "Br",
                "atomic number": 35,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "Kr",
                "atomic number": 36,
                "number of core orbitals": 14,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d"],
                "valence orbitals": ["4s", "4p"],
            },
            {
                "name": "Rb",
                "atomic number": 37,
                "number of core orbitals": 18,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s"],
            },
            {
                "name": "Sr",
                "atomic number": 38,
                "number of core orbitals": 18,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s"],
            },
            {
                "name": "Y",
                "atomic number": 39,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Zr",
                "atomic number": 40,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Nb",
                "atomic number": 41,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Mo",
                "atomic number": 42,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Tc",
                "atomic number": 43,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Ru",
                "atomic number": 44,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Rh",
                "atomic number": 45,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "Pd",
                "atomic number": 46,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["4d", "5s", "5p"],
            },
            {
                "name": "Ag",
                "atomic number": 47,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["4d", "5s", "5p"],
            },
            {
                "name": "Cd",
                "atomic number": 48,
                "number of core orbitals": 18,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p"],
                "valence orbitals": ["5s", "4d", "5p"],
            },
            {
                "name": "In",
                "atomic number": 49,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "Sn",
                "atomic number": 50,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "Sb",
                "atomic number": 51,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "Te",
                "atomic number": 52,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "I",
                "atomic number": 53,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "Xe",
                "atomic number": 54,
                "number of core orbitals": 23,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                "valence orbitals": ["5s", "5p"],
            },
            {
                "name": "Cs",
                "atomic number": 55,
                "number of core orbitals": 27,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s"],
            },
            {
                "name": "Ba",
                "atomic number": 56,
                "number of core orbitals": 27,
                "number of valence orbitals": 1,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s"],
            },
            {
                "name": "La",
                "atomic number": 57,
                "number of core orbitals": 27,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Ce",
                "atomic number": 58,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Pr",
                "atomic number": 59,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Nd",
                "atomic number": 60,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Pm",
                "atomic number": 61,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Sm",
                "atomic number": 62,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Eu",
                "atomic number": 63,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Gd",
                "atomic number": 64,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Tb",
                "atomic number": 65,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Dy",
                "atomic number": 66,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Ho",
                "atomic number": 67,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Er",
                "atomic number": 68,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Tm",
                "atomic number": 69,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Yb",
                "atomic number": 70,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Lu",
                "atomic number": 71,
                "number of core orbitals": 27,
                "number of valence orbitals": 16,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "4f", "6p"],
            },
            {
                "name": "Hf",
                "atomic number": 72,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Ta",
                "atomic number": 73,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "W",
                "atomic number": 74,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Re",
                "atomic number": 75,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Os",
                "atomic number": 76,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Ir",
                "atomic number": 77,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Pt",
                "atomic number": 78,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Au",
                "atomic number": 79,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Hg",
                "atomic number": 80,
                "number of core orbitals": 34,
                "number of valence orbitals": 9,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p"],
                "valence orbitals": ["6s", "5d", "6p"],
            },
            {
                "name": "Tl",
                "atomic number": 81,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "Pb",
                "atomic number": 82,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "Bi",
                "atomic number": 83,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "Po",
                "atomic number": 84,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "At",
                "atomic number": 85,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "Rn",
                "atomic number": 86,
                "number of core orbitals": 39,
                "number of valence orbitals": 4,
                "core orbitals": ["1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d"],
                "valence orbitals": ["6s", "6p"],
            },
            {
                "name": "Fr",
                "atomic number": 87,
                "number of core orbitals": 43,
                "number of valence orbitals": 1,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s"],
            },
            {
                "name": "Ra",
                "atomic number": 88,
                "number of core orbitals": 43,
                "number of valence orbitals": 1,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s"],
            },
            {
                "name": "Ac",
                "atomic number": 89,
                "number of core orbitals": 43,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Th",
                "atomic number": 90,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Pa",
                "atomic number": 91,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "U",
                "atomic number": 92,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Np",
                "atomic number": 93,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Pu",
                "atomic number": 94,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Am",
                "atomic number": 95,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Cm",
                "atomic number": 96,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Bk",
                "atomic number": 97,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Cf",
                "atomic number": 98,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Es",
                "atomic number": 99,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Fm",
                "atomic number": 100,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Md",
                "atomic number": 101,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "No",
                "atomic number": 102,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Lr",
                "atomic number": 103,
                "number of core orbitals": 43,
                "number of valence orbitals": 16,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "5f", "7p"],
            },
            {
                "name": "Rf",
                "atomic number": 104,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Db",
                "atomic number": 105,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Sg",
                "atomic number": 106,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Bh",
                "atomic number": 107,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Hs",
                "atomic number": 108,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Mt",
                "atomic number": 109,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Ds",
                "atomic number": 110,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Rg",
                "atomic number": 111,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Cn",
                "atomic number": 112,
                "number of core orbitals": 50,
                "number of valence orbitals": 9,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p"],
                "valence orbitals": ["7s", "6d", "7p"],
            },
            {
                "name": "Nh",
                "atomic number": 113,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
            {
                "name": "Fl",
                "atomic number": 114,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
            {
                "name": "Mc",
                "atomic number": 115,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
            {
                "name": "Lv",
                "atomic number": 116,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
            {
                "name": "Ts",
                "atomic number": 117,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
            {
                "name": "Og",
                "atomic number": 118,
                "number of core orbitals": 55,
                "number of valence orbitals": 4,
                "core orbitals": [
                    "1s", "2s", "2p", "3s", "3p", "3d", "4s", "4p", "4d", "4f", "5s", "5p", "5d", "5f", "6s", "6p", "6d"
                ],
                "valence orbitals": ["7s", "7p"],
            },
        ]
        # fmt: on
        """contains information on elements H - Og"""
        if self.double_d_shell:
            for element in self.element_list:
                if element.get("name") in self.double_d_shell_elements:
                    element["number of valence orbitals"] += 5
                    element["valence orbitals"].append("4d")

    def get(self, element_string: str) -> ElementDict:
        """Return an dictionary of an element.

        Parameters
        ----------
        element_string : str
            name of an element

        Returns
        -------
        element : ElementDict
            conatins required information of an atom

        Raises
        ------
        NotImplementedError
            if element string is not found in element_list
        """
        for element in self.element_list:
            if element.get("name") == element_string:
                return element
        raise NotImplementedError(
            f"Element '{element_string}' is currently not implemented in chemical_elements.py"
        )

    def get_core_orb_labels(self, atom: str) -> List[str]:
        """Return orbital labels of core orbitals.

        Parameters
        ----------
        atom : str
            name of the element

        Returns
        -------
        labels : List[str]
            List with all core orbital labels
        """
        element = self.get(atom)
        labels = element["core orbitals"]
        return labels

    def get_valence_orb_labels(self, atom: str) -> List[str]:
        """Return orbital labels of valence orbitals.

        Parameters
        ----------
        atom : str
            name of the element

        Returns
        -------
        labels : List[str]
            List with all core orbital labels
        """
        element = self.get(atom)
        labels = element["valence orbitals"]
        return labels

    def get_core_orbitals(self, atom: str) -> int:
        """Get the number of core orbitals for an element.

        Parameters
        ----------
        atom : str
            name of the atom

        Returns
        -------
        n_core_orbitals : int
            number of core orbitals
        """
        element = self.get(atom)
        n_core_orbitals = element["number of core orbitals"]
        return n_core_orbitals

    def get_valence_orbitals(self, atom: str) -> int:
        """Get the number of valence orbitals for an element.

        Parameters
        ----------
        atom : str
            name of the atom

        Returns
        -------
        valence_orbitals : int
            number of valence orbitals
        """
        element = self.get(atom)
        valence_orbitals = element["number of valence orbitals"]
        return valence_orbitals

    def get_electrons(self, atom: str) -> int:
        """Get the number of electrons for an element.

        Parameters
        ----------
        atom : str
            name of the atom

        Returns
        -------
        atomic_number : int
            number of electrons
        """
        element = self.get(atom)
        atomic_number = element["atomic number"]
        return atomic_number
