"""Module to provide main autocas functions.

All functions here, provide basic autocas workflows for groundstates and excited states
in combination with standard and large cas protocols.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


import os
from typing import Any, Dict, List, Tuple, cast

import numpy as np

from scine_autocas import Autocas
from scine_autocas.autocas_utils.input_handler import InputHandler
from scine_autocas.autocas_utils.molecule import Molecule
from scine_autocas.interfaces import Interface
from scine_autocas.interfaces.molcas import Molcas
from scine_autocas.plots.entanglement_plot import EntanglementPlot
from scine_autocas.plots.threshold_plot import ThresholdPlot


class MainFunctions:
    """Class for autocas main.

    Provides some utility functions.

    Attributes
    ----------
    results : Dict[str, Any]
        stores all results
    """

    def __init__(self):
        """Init."""
        self.results: Dict[str, Any] = {}

    def _large_cas_loop(self, interface: Interface):
        partial_s1_list = []
        partial_s2_list = []
        partial_mutual_information_list = []

        large_cas_occupations = self.results["partial_occupations"]
        large_cas_indices = self.results["partial_indices"]

        print("")
        print(f"Number of Sub-CAS: {len(large_cas_occupations)}")
        for i, occupation in enumerate(large_cas_occupations):
            print("")
            print(f"Evaluating CAS: {i} / {len(large_cas_occupations)}")
            print(f"Current CAS(e,o): ({sum(occupation)}, {len(occupation)})")
            print(f"Current orbital indices:    {large_cas_indices[i]}")
            print(f"Current orbital occupation: {occupation}")
            (
                partial_energy,
                s1_partial,
                s2_partial,
                mutual_informaion_partial,
            ) = interface.calculate(occupation, large_cas_indices[i])
            # energy is meaningless
            _ = partial_energy
            # we know we have excited states
            s1_partial = cast(List[np.ndarray], s1_partial)
            s2_partial = cast(List[np.ndarray], s2_partial)
            mutual_informaion_partial = cast(List[np.ndarray], mutual_informaion_partial)

            partial_s1_list.append(s1_partial)
            partial_s2_list.append(s2_partial)
            partial_mutual_information_list.append(mutual_informaion_partial)

            # FOR CAST
            # no excited states, so we know its a ndarray not List[ndarray]
            # s1_partial = cast(np.ndarray, s1_partial)
            # s2_partial = cast(np.ndarray, s2_partial)
            # mutual_information_partial = cast(np.ndarray, mutual_information_partial)

        self.results["partial_s1"] = partial_s1_list
        self.results["partial_s2"] = partial_s2_list
        self.results["partial_mutual_information_list"] = partial_mutual_information_list

    def large_cas_excited_states(self, autocas: Autocas, interface: Interface):
        """Evaluate initial cas with excited states and large cas protocol.

        Parameters
        ----------
        autocas : Autocas
            the autocas object to search for the active space
        interface : Interface
            an interface to a suitable electronic structure program
        """
        print("Running large active spaces with excited states algorithm.")
        # raise NotImplementedError("large cas with excited states")
        if interface.dumper:
            interface.dumper.large_cas = True

        try:
            initial_occupation = self.results["initial_occupation"]
            initial_orbital_indices = self.results["initial_orbital_indices"]
        except KeyError:
            initial_occupation, initial_orbital_indices = autocas.make_initial_active_space()
            self.results["initial_occupation"] = initial_occupation
            self.results["initial_orbital_indices"] = initial_orbital_indices

        large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()

        self.results["partial_occupations"] = large_cas_occupations
        self.results["partial_indices"] = large_cas_indices

        # interface.calculate()
        self._large_cas_loop(interface)
        # for i, occupation in enumerate(large_cas_occupations):
        #     (
        #         partial_energy,
        #         s1_partial,
        #         s2_partial,
        #         mutual_informaion_partial,
        #     ) = interface.calculate(occupation, large_cas_indices[i])
        #     # energy is meaningless
        #     _ = partial_energy
        #     # we know we have excited states
        #     s1_partial = cast(List[np.ndarray], s1_partial)
        #     s2_partial = cast(List[np.ndarray], s2_partial)
        #     mutual_informaion_partial = cast(List[np.ndarray], mutual_informaion_partial)

        #     partial_s1_list.append(s1_partial)
        #     partial_s2_list.append(s2_partial)
        #     partial_mutual_information_list.append(mutual_informaion_partial)

        partial_s1_list = self.results["partial_s1"]
        partial_s2_list = self.results["partial_s2"]
        partial_mutual_information_list = self.results["partial_mutual_information_list"]

        partial_s1_per_state: List[List[np.ndarray]] = []
        partial_s2_per_state: List[List[np.ndarray]] = []
        partial_mut_inf_per_state: List[List[np.ndarray]] = []

        for i in range(len(partial_s1_list[0])):
            partial_s1_per_state.append([])
            partial_s2_per_state.append([])
            partial_mut_inf_per_state.append([])

        for i, part_s1 in enumerate(partial_s1_list):
            for j, part_s1_j in enumerate(part_s1):
                partial_s1_per_state[j].append(part_s1_j)
                partial_s2_per_state[j].append(partial_s2_list[i][j])
                partial_mut_inf_per_state[j].append(partial_mutual_information_list[i][j])

        # only for storage not for autocas
        autocas_results = autocas.collect_entropies(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_per_state[0],
            partial_s2_per_state[0],
            partial_mut_inf_per_state[0],
        )

        if len(autocas_results) == 4:
            autocas_results = cast(Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], autocas_results)
            initial_occupation = autocas_results[0].tolist()
            initial_s1 = autocas_results[1]
            initial_s2 = autocas_results[2]
            initial_mut_inf = autocas_results[3]

        self.results["initial_occupation"] = initial_occupation
        self.results["initial_s1"] = initial_s1
        self.results["initial_s2"] = initial_s2
        self.results["initial_mutual_information"] = initial_mut_inf

        final_occupation, final_orbital_indices = autocas.get_cas_from_large_cas_excited_states(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_list,
            partial_s2_list,
            partial_mutual_information_list,
        )

        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices
        interface.dumper.large_cas = False

        return final_occupation, final_orbital_indices

    def large_cas(self, autocas: Autocas, interface: Interface):
        """Evaluate initial cas with large cas protocol.

        Parameters
        ----------
        autocas : Autocas
            the autocas object to search for the active space
        interface : Interface
            an interface to a suitable electronic structure program
        """
        print("Running large active spaces algorithm.")
        if interface.dumper:
            interface.dumper.large_cas = True

        try:
            initial_occupation = self.results["initial_occupation"]
            initial_orbital_indices = self.results["initial_orbital_indices"]
        except KeyError:
            initial_occupation, initial_orbital_indices = autocas.make_initial_active_space()
            self.results["initial_occupation"] = initial_occupation
            self.results["initial_orbital_indices"] = initial_orbital_indices

        large_cas_occupations, large_cas_indices = autocas.get_large_active_spaces()

        self.results["partial_occupations"] = large_cas_occupations
        self.results["partial_indices"] = large_cas_indices

        self._large_cas_loop(interface)

        # for i, occupation in enumerate(large_cas_occupations):
        #     (
        #         partial_energy,
        #         s1_partial,
        #         s2_partial,
        #         mutual_information_partial,
        #     ) = interface.calculate(occupation, large_cas_indices[i])
        #     # energy is meaningless
        #     _ = partial_energy
        #     # no excited states, so we know its a ndarray not List[ndarray]
        #     s1_partial = cast(np.ndarray, s1_partial)
        #     s2_partial = cast(np.ndarray, s2_partial)
        #     mutual_information_partial = cast(np.ndarray, mutual_information_partial)
        #     partial_s1_list.append(s1_partial)
        #     partial_s2_list.append(s2_partial)
        #     partial_mutual_information_list.append(mutual_information_partial)

        partial_s1_list = self.results["partial_s1"]
        partial_s2_list = self.results["partial_s2"]
        partial_mutual_information_list = self.results["partial_mutual_information_list"]

        autocas_results = autocas.collect_entropies(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_list,
            partial_s2_list,
            partial_mutual_information_list
        )

        if len(autocas_results) == 4:
            autocas_results = cast(Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], autocas_results)
            initial_occupation = autocas_results[0].tolist()
            initial_s1 = autocas_results[1]
            initial_s2 = autocas_results[2]
            initial_mut_inf = autocas_results[3]

        self.results["initial_occupation"] = initial_occupation
        self.results["initial_s1"] = initial_s1
        self.results["initial_s2"] = initial_s2
        self.results["initial_mutual_information"] = initial_mut_inf

        final_occupation, final_orbital_indices = autocas.get_cas_from_large_cas(
            large_cas_indices,
            large_cas_occupations,
            partial_s1_list,
            partial_s2_list,
            partial_mutual_information_list,
        )

        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices
        interface.dumper.large_cas = False

        return final_occupation, final_orbital_indices

    def excited_states(self, autocas: Autocas, interface: Interface):
        """Evaluate initial cas with excited states.

        Parameters
        ----------
        autocas : Autocas
            the autocas object to search for the active space
        interface : Interface
            an interface to a suitable electronic structure program
        """
        print("Running excited states algorithm.")
        try:
            initial_occupation = self.results["initial_occupation"]
            initial_orbital_indices = self.results["initial_orbital_indices"]
        except KeyError:
            initial_occupation, initial_orbital_indices = autocas.make_initial_active_space()
            self.results["initial_occupation"] = initial_occupation
            self.results["initial_orbital_indices"] = initial_orbital_indices

        energy, initial_s1, initial_s2, initial_mutual_information = interface.calculate(
            initial_occupation, initial_orbital_indices
        )
        # energy is meaningless
        _ = energy
        # we know we calculate excited states
        initial_s1 = cast(List[np.ndarray], initial_s1)
        initial_s2 = cast(List[np.ndarray], initial_s2)
        initial_mutual_information = cast(List[np.ndarray], initial_mutual_information)

        self.results["initial_s1"] = initial_s1[0]
        self.results["initial_s2"] = initial_s2[0]
        self.results["initial_mutual_information"] = initial_mutual_information[0]

        final_occupation, final_orbital_indices = autocas.get_cas_from_excited_states(
            initial_occupation,
            initial_s1,
            mutual_information_list=initial_mutual_information,
            force_cas=False,
        )

        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices

        return final_occupation, final_orbital_indices

    def conventional(self, autocas: Autocas, interface: Interface):
        """Evaluate initial cas.

        Parameters
        ----------
        autocas : Autocas
            the autocas object to search for the active space
        interface : Interface
            an interface to a suitable electronic structure program
        """
        try:
            initial_occupation = self.results["initial_occupation"]
            initial_orbital_indices = self.results["initial_orbital_indices"]
        except KeyError:
            initial_occupation, initial_orbital_indices = autocas.make_initial_active_space()
            self.results["initial_occupation"] = initial_occupation
            self.results["initial_orbital_indices"] = initial_orbital_indices

        # interface.calculate()
        energy, initial_s1, initial_s2, initial_mutual_information = interface.calculate(
            initial_occupation, initial_orbital_indices
        )
        # energy is meaningless
        _ = energy
        # we know we do not calculate excited states
        initial_s1 = cast(np.ndarray, initial_s1)
        initial_s2 = cast(np.ndarray, initial_s2)
        initial_mutual_information = cast(np.ndarray, initial_mutual_information)

        self.results["initial_s1"] = initial_s1
        self.results["initial_s2"] = initial_s2
        self.results["initial_mutual_information"] = initial_mutual_information

        final_occupation, final_orbital_indices = autocas.get_active_space(
            initial_occupation,
            initial_s1,
            mutual_information=initial_mutual_information,
            force_cas=False,
        )
        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices

        return final_occupation, final_orbital_indices

    # flake8: noqa: C901
    def main(self, settings_dict: Dict):
        """Provide main function of autocas.

        Parameters
        ----------
        settings_dict : Dict
            provide yaml settings
        """
        # looks better in output (escape \)
        print("                                  _            _____             _____                     ")
        print("                                 | |          / ____|    /\\     / ____|                    ")
        print("                     __ _  _   _ | |_   ___  | |        /  \\   | (___                      ")
        print("                    / _` || | | || __| / _ \\ | |       / /\\ \\   \\___ \\                     ")
        print("                   | (_| || |_| || |_ | (_) || |____  / ____ \\  ____) |                    ")
        print("                    \\__,_| \\__,_| \\__| \\___/  \\_____|/_/    \\_\\|_____/                     ")
        print("                                                                                           ")

        print("*******************************************************************************************")
        print("*                                                                                         *")
        print("*                                      Settings                                           *")
        print("*                                                                                         *")
        print("*******************************************************************************************")

        # if "yaml_input" in settings_dict:
        if settings_dict["yaml_input"] is not None:
            input_handler = InputHandler(settings_dict["yaml_input"])
            input_handler.print_settings()
            settings = input_handler.settings_dir
            molecule = input_handler.get_molecule()
            autocas = input_handler.get_autocas()
            interface = input_handler.get_interface()
        elif "xyz_file" in settings_dict:
            for key in settings_dict:
                print(f"    {key}: {settings_dict[key]}")
            settings = settings_dict
            molecule = Molecule(xyz_file=settings_dict["xyz_file"])
            autocas = Autocas(molecule=molecule)
            if settings_dict["xyz_file"][0] != "/":
                settings_dict["xyz_file"] = os.getcwd() + "/" + settings_dict["xyz_file"]
            interface = Molcas(molecules=[molecule])
            interface.settings.basis_set = settings_dict["basis_set"]
            interface.settings.xyz_file = settings_dict["xyz_file"]
            interface.environment.molcas_scratch_dir = os.getcwd() + "/molcas_scratch"
        else:
            raise ValueError("Provide yaml input or an xyz-file.")

        self.results["settings"] = settings
        self.results["autocas"] = autocas
        self.results["molecule"] = molecule
        self.results["interface"] = interface

        use_large_cas = False
        if "large_cas" in settings and settings["large_cas"] is True:
            use_large_cas = True
            print("The Large active space protocol is enabled!")

        interface.settings.method = "dmrg_ci"
        interface.settings.dmrg_sweeps = 5
        interface.settings.dmrg_bond_dimension = 250
        interface.settings.post_cas_method = ""

        print("")
        print("*******************************************************************************************")
        print("*                                                                                         *")
        print("*                                  Initial HF Calculation                                 *")
        print("*                                                                                         *")
        print("*******************************************************************************************")
        # make initial hf calculation
        interface.calculate()

        print("")
        print("*******************************************************************************************")
        print("*                                                                                         *")
        print("*                                    DMRG Calculation                                     *")
        print("*                                                                                         *")
        print("*******************************************************************************************")

        print("Settings for autoCAS DMRG calculation:")
        print("method:              'dmrg_ci'")
        print("dmrg_sweeps:           5")
        print("dmrg_bond_dimension: 250")
        print("Switching back to original settings for the final calculation.")
        print("")

        (
            initial_occupation,
            initial_orbital_indices,
        ) = autocas.make_initial_active_space()
        self.results["initial_occupation"] = initial_occupation
        self.results["initial_orbital_indices"] = initial_orbital_indices

        print(f"Initial active space CAS(e, o): ({sum(initial_occupation)}, {len(initial_occupation)})")
        print(f"Orbital indices:     {initial_orbital_indices}")
        print(f"Orbital occupations: {initial_occupation}")
        print("")

        # large cas and excited states
        if interface.settings.n_excited_states > 0 and use_large_cas:
            final_occupation, final_orbital_indices = self.large_cas_excited_states(autocas, interface)

        # large active space protocol groundstate only
        elif use_large_cas:
            final_occupation, final_orbital_indices = self.large_cas(autocas, interface)

        # conventional excited states
        elif interface.settings.n_excited_states > 0:
            final_occupation, final_orbital_indices = self.excited_states(autocas, interface)

        # conventional groundstate only
        else:
            final_occupation, final_orbital_indices = self.conventional(autocas, interface)

        if "entanglement_diagram" in settings or "plot" in settings:
            print("Creating entanglement diagram.")
            entanglement_plot = EntanglementPlot()
            plt = entanglement_plot.plot(
                self.results["initial_s1"],
                self.results["initial_mutual_information"],
                self.results["initial_orbital_indices"],
            )
            if "entanglement_diagram" in settings:
                plt.savefig(str(settings["entanglement_diagram"]))  # type: ignore[attr-defined]
            else:
                plt.show()  # type: ignore[attr-defined]

        # if "full_entanglement_diagram" in settings:
        #    entanglement_plot = EntanglementPlot()
        #    plt = entanglement_plot.plot_in_plot(
        #        self.results["initial_s1"],
        #        self.results["initial_mutual_information"],
        #        self.results["final_s1"],
        #        self.results["final_mutual_information"],
        #        self.results["initial_orbital_indices"],
        #        self.results["final_orbital_indices"],
        #    )
        #    plt.savefig(settings["full_entanglement_diagram"], dpi=300)

        if "threshold_diagram" in settings:
            print("Creating threshold diagram.")
            threshold_plot = ThresholdPlot()
            plt = threshold_plot.plot(self.results["initial_s1"])
            plt.savefig(str(settings["threshold_diagram"]))  # type: ignore[attr-defined]

        if settings_dict["yaml_input"] is not None:
            try:
                interface.settings.method = settings["interface"]["settings"]["method"]
            except KeyError:
                pass
            try:
                interface.settings.post_cas_method = settings["interface"]["settings"]["post_cas_method"]
            except KeyError:
                pass
            try:
                interface.settings.dmrg_sweeps = settings["interface"]["settings"]["dmrg_sweeps"]
            except KeyError:
                pass
            try:
                interface.settings.dmrg_bond_dimension = settings["interface"]["settings"]["dmrg_bond_dimension"]
            except KeyError:
                pass
        else:
            try:
                interface.settings.method = settings["method"]
            except KeyError:
                pass
            try:
                interface.settings.post_cas_method = settings["post_cas_method"]
            except KeyError:
                pass
            try:
                interface.settings.dmrg_bond_dimension = settings["dmrg_bond_dimension"]
            except KeyError:
                pass
            try:
                interface.settings.dmrg_sweeps = settings["dmrg_sweeps"]
            except KeyError:
                pass

        print("")
        print("*******************************************************************************************")
        print("*                                                                                         *")
        print("*                                    Final Calculation                                    *")
        print("*                                                                                         *")
        print("*******************************************************************************************")

        print("Starting final calculation with optimized CAS.")
        print(f"Final active space CAS(e, o): ({sum(final_occupation)}, {len(final_occupation)})")
        print(f"Final orbital indices:    {final_orbital_indices}")
        print(f"Final orbital occupation: {final_occupation}")
        print("")

        energy, s1_entropy, s2_entropy, mutual_information = interface.calculate(
            final_occupation, final_orbital_indices
        )
        self.results["energy"] = energy
        print(f"Final energy: {energy}")
        try:
            self.results["s1"] = s1_entropy
            self.results["s2"] = s2_entropy
            self.results["mutual_information"] = mutual_information
            if "full_entanglement_diagram" in settings:
                print("Creating full entanglement diagram.")
                entanglement_plot = EntanglementPlot()
                plt = entanglement_plot.plot_in_plot(
                    self.results["initial_s1"],
                    self.results["initial_mutual_information"],
                    self.results["final_s1"],
                    self.results["final_mutual_information"],
                    self.results["initial_orbital_indices"],
                    self.results["final_orbital_indices"],
                )
                plt.savefig(settings["full_entanglement_diagram"], dpi=600)  # type: ignore[attr-defined]
        except KeyError:
            pass
