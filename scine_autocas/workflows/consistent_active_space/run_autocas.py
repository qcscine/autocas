# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """

from typing import List, Tuple
from matplotlib import pyplot
from scine_autocas import Autocas, Molecule
from scine_autocas.interfaces import Molcas
from scine_autocas.plots import EntanglementPlot, ThresholdPlot


def run_small_active_space(autocas: Autocas, molcas: Molcas, name: str,
                           occ_initial: List[int], index_initial: List[int])\
        -> Tuple[List[int], List[int]]:
    """
    Run autoCAS with the small active space protocol.

    Parameters
    ----------
    autocas: Autocas
        The autocas object.
    molcas: Molcas
        The molcas interface.
    name: str
        The system name.
    occ_initial: List[int]
        The initial occupations of the active space.
    index_initial: List[int]
        The initial orbital indices of the active space.

    Returns
    -------
    Tuple[List[int], List[int]]
        AutoCAS suggestion for the occupations and indices of the CAS.
    """
    # do cas calculation
    _, s1_entropy, __, mut_inf = molcas.calculate(occ_initial, index_initial)

    # plot entanglement + threshold diagram
    if s1_entropy.size > 0:
        plot = EntanglementPlot()
        plot.plot(s1_entropy, mut_inf)  # type: ignore
        pyplot.savefig("entang" + name + ".pdf")
        pyplot.close()
        plot2 = ThresholdPlot()
        plot2.plot(s1_entropy)  # type: ignore
        pyplot.savefig("threshold" + name + ".pdf")  # type: ignore
    # make active space based on single orbital entropies
    cas_occ, cas_index = autocas.get_active_space(
        s1_entropy   # type: ignore
    )
    return cas_occ, cas_index


def run_large_active_space(autocas: Autocas, molcas: Molcas, name: str) -> Tuple[List[int], List[int]]:
    """
    Run autoCAS with the large active space protocol.

    Parameters
    ----------
    autocas: Autocas
        The autocas object.
    molcas: Molcas
        The molcas interface.
    name: str
        The system name.

    Returns
    -------
    Tuple[List[int], List[int]]
        AutoCAS suggestion for the occupations and indices of the CAS.

    """
    large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()
    # iterate over subspaces and evaluate them
    partial_s1_list = []
    partial_s2_list = []
    partial_mut_inf_list = []
    print(f"autocas will perform: {len(large_cas_indices)} dmrg calculations")
    print(f"with {autocas.large_spaces.max_orbitals} in each active space")

    for partial_indices, partial_occupation in zip(large_cas_indices, large_cas_occupations):
        _, s1_partial, s2_partial, mut_inf_partial = molcas.calculate(
            partial_occupation, partial_indices
        )
        partial_s1_list.append(s1_partial)
        partial_s2_list.append(s2_partial)
        partial_mut_inf_list.append(mut_inf_partial)

    # combine entropies and occupation to start autocas
    # useful for plotting
    initial_occupation, initial_s1, initial_s2, initial_mut_inf = autocas.collect_entropies(  # type: ignore
        large_cas_indices,
        large_cas_occupations,
        partial_s1_list,  # type: ignore
        partial_s2_list,  # type: ignore
        partial_mut_inf_list,  # type: ignore
    )

    # plot stuff
    _ = initial_occupation  # not required for plotting
    _ = initial_s2  # not required for plotting
    entang_plot = EntanglementPlot()
    entang_plot.plot(initial_s1, initial_mut_inf)
    pyplot.savefig("entang" + name + ".pdf")
    thresh_plot = ThresholdPlot()
    plt = thresh_plot.plot(initial_s1)  # type: ignore
    plt.savefig("threshold" + name + ".pdf")  # type: ignore

    cas_occ, cas_index = autocas.get_cas_from_large_cas(
        large_cas_indices,
        large_cas_occupations,
        partial_s1_list,  # type: ignore
        partial_s2_list,  # type: ignore
        partial_mut_inf_list,  # type: ignore
        force_cas=False
    )

    return cas_occ, cas_index


def run_autocas(molecule: Molecule, molcas: Molcas, name: str, large_active_space: bool) -> Tuple[List[int], List[int]]:
    """
    Run the autoCAS workflow.

    Parameters
    ----------
    molecule: Molecule
        The molecule object containing the system information.
    molcas: Molcas
        The Molcas interface for performing calculations.
    name: str
        The name of the system.
    large_active_space: bool
        Flag indicating whether to use the large active space protocol. The large active space protocol is only used for
        large valence spaces.

    Returns
    -------
    Tuple[List[int], List[int]]
        The occupations and indices of the CAS suggested by autoCAS.
    """
    # initialize autoCAS
    autocas = Autocas(molecule)
    autocas.plateau_values = 12
    # make initial active space and evaluate initial DMRG calculation
    occ_initial, index_initial = autocas.make_initial_active_space()

    if large_active_space and len(occ_initial) > autocas.large_spaces.max_orbitals:
        print("Large Active Space Selection")
        autocas.large_spaces.max_orbitals = min(len(occ_initial), autocas.large_spaces.max_orbitals)
        return run_large_active_space(autocas, molcas, name)
    return run_small_active_space(autocas, molcas, name, occ_initial, index_initial)
