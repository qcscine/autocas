# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
from scine_autocas.io import FileHandler, logger

from .workflow import Workflow


class ClassicWorkflow(Workflow):
    """The classic workflow"""

    def __init__(self, autocas, interface):
        super().__init__(autocas, interface)

    def _initial_dmrg_impl(self):
        """Initial DMRG calculation"""
        logger.frame("Initial DMRG calculation")
        FileHandler.make_initial_dmrg_dir()
        initial_energy, initial_s1, initial_s2, init_mut_inf = self.interface.calculate(
            self.results["initial_occupation"], self.results["initial_orbital_indices"]
        )
        # logging.debug(f"Initial CAS energy: {initial_energy}")
        print(f"Initial CAS energy: {initial_energy}")
        print(f"initial s1: {initial_s1}")
        print(f"max s1: {max(initial_s1)}")
        self.results["initial_energy"] = initial_energy
        self.results["initial_s1"] = initial_s1
        self.results["initial_s2"] = initial_s2
        self.results["initial_mutual_information"] = init_mut_inf
        try:
            # pylint: disable=W0212
            orbital_order = self.interface._pyscf_cas_solver.parameters.get('orbital_order')
            print(f"initial orbital_order: {orbital_order}")
            self.results["initial_orbital_order"] = orbital_order
        except AttributeError:
            pass

        logger.empty_line()

        # create entanglement and threshold diagrams
        self.plotting()
        logger.empty_line()

        logger.frame("Active space search")
        final_occupation, final_orbital_indices = self.autocas.get_active_space(
            initial_s1,
            mutual_information=init_mut_inf,
            occupation=self.results["initial_occupation"],
            indices=self.results["initial_orbital_indices"],
            force_cas=False,
        )
        logger.final_cas(final_occupation, final_orbital_indices)
        logger.empty_line()
        self.results["final_occupation"] = final_occupation
        self.results["final_orbital_indices"] = final_orbital_indices
