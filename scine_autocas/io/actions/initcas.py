# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
from scine_autocas.cas_selection import Autocas
from scine_autocas.io import logger
from scine_autocas.io.input_handler import InputHandler
from scine_autocas.utils.molecule import molecule_setup_from_dict


def initcas_action(inp: InputHandler):
    """Run init cas function.

    Parameters
    ----------
    input : InputHandler
        inputhandler initialized by cli arguments
    """
    molecule = molecule_setup_from_dict(inp.get_molecule_options())
    autocas = Autocas(molecule)
    autocas.setup_from_dict(inp.get_autocas_options())

    logger.frame("Initial CAS")
    init_occ, init_index = autocas.make_initial_active_space()
    logger.init_cas(init_occ, init_index)

    if inp.get_workflow() in ("large_cas", "large_cas_excited_states"):
        logger.empty_line()
        print("Sub-CAS in large cas protocol")
        large_cas_spaces = autocas.get_large_active_spaces()
        large_cas_occupations = large_cas_spaces[0]
        large_cas_indices = large_cas_spaces[1]
        print(f"Number of sub-CASs: {len(large_cas_occupations)}")
        for i, partial_occupation in enumerate(large_cas_occupations):
            print(f" Cas space: {i+1}/{len(large_cas_occupations)}")
            print(f"  occupation {partial_occupation}")
            print(f"  indices    {large_cas_indices[i]}")
            print("")
