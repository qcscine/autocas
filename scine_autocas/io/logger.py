# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
from typing import List

import numpy as np


def init_cas(init_occ: List[int], init_index: List[int]):
    """Convenience function for initial active space printing.

    Parameters
    ----------
    init_occ: List[int]
        initial orbital occupation
    init_index: List[int]
        orbital indices(0-based)
    """
    print(f"CAS(e, o) = CAS({sum(init_occ)}, {len(init_occ)})")
    print(f"Initial occupation:        {init_occ}")
    print(f"Initial indices (0-based): {init_index}")


def final_cas(final_occ: List[int], final_index: List[int]):
    """Convenience function for final active space printing.

    Parameters
    ----------
    final_occ: List[int]
        final orbital occupation
    final_index: List[int]
        orbital indices(0-based)
    """
    print(f"CAS(e, o) = CAS({sum(final_occ)}, {len(final_occ)})")
    print(f"Final occupation:        {final_occ}")
    print(f"Final indices (0-based): {final_index}")


def s1_entropies(entropies: np.ndarray):
    """Convenience function for single orbital entropy printing.

    Parameters
    ----------
    s1_entropies: np.ndarray
        The single orbital entropies
    """
    print(f"s1 entropies: {entropies}")
    print(f"max(s1):      {entropies}")


def banner():
    """Print the autocas banner."""
    print("                            _            _____             _____                ")
    print("                           | |          / ____|    /\\     / ____|               ")
    print("               __ _  _   _ | |_   ___  | |        /  \\   | (___                 ")
    print("              / _` || | | || __| / _ \\ | |       / /\\ \\   \\___ \\                ")
    print("             | (_| || |_| || |_ | (_) || |____  / ____ \\  ____) |               ")
    print("              \\__,_| \\__,_| \\__| \\___/  \\_____|/_/    \\_\\|_____/                ")
    print("                                                                                ")


def frame(message: str):
    """Print frame around message.

    Parameters
    ----------
    message: str
        The string to print.
    """
    white_space_length = 80 - len(message) - 2
    print("*"*80)
    print("*" + " "*78 + "*")
    print("*" + " " * int(((white_space_length+1)/2)) + message + " " * int((white_space_length/2)) + "*")
    print("*" + " "*78 + "*")
    print("*"*80)


def empty_line():
    """Print an empty line."""
    print("")
