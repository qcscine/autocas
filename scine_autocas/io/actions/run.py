# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
from scine_autocas.cas_selection import Autocas
from scine_autocas.interfaces import Molcas, PyscfInterface
from scine_autocas.interfaces.interface import Interface
from scine_autocas.io.input_handler import InputHandler
from scine_autocas.utils import Molecule
from scine_autocas.utils.defaults import AvailableInteraces
from scine_autocas.utils.exceptions import UnkownWorkflowError
from scine_autocas.utils.molecule import molecule_setup_from_dict
from scine_autocas.workflows import (ClassicExcitedStatesWorkflow, ClassicWorkflow, LargeCasExcitedStateWorkflow,
                                     LargeCasWorkflow)
from scine_autocas.workflows.workflow import Workflow

# from .action import Action


def _get_interface(molecule: Molecule, inp: InputHandler) -> Interface:
    """Get correct interface based on settings.

    Parameters
    ----------
    molecule : Molecule
        the molecule object
    inp : InputHandler
        inputhandler initialized by cli arguments

    Returns
    -------
    interface : Interface
        an initialized interface object
    """
    interface: Interface
    if inp.interface == AvailableInteraces.PYSCF:
        # logging.info("Interface: PySCF")
        print("Interface: PySCF")
        interface = PyscfInterface(molecule)
    elif inp.interface in (AvailableInteraces.OPENMOLCAS, AvailableInteraces.MOLCAS):
        # logging.info("Interface: OpenMolcas")
        print("Interface: OpenMolcas")
        interface = Molcas(molecule)
    else:
        # logging.critical(f"Interface: {inp.interface} is not implemented")
        print(f"Interface: {inp.interface} is not implemented")
        raise NotImplementedError
    interface.setup_from_dict(inp.get_interface_options())
    return interface


def run_action(inp: InputHandler):
    """Run an autocas calculation.

    Parameters
    ----------
    inp : InputHandler
        inputhandler initialized by cli arguments
    """
    molecule = molecule_setup_from_dict(inp.get_molecule_options())
    print(f"Molecule: {molecule}")
    autocas = Autocas(molecule)  # , settings_dict=inp.get_autocas_options())
    autocas.setup_from_dict(inp.get_autocas_options())
    interface = _get_interface(molecule, inp)

    workflow_name = inp.get_workflow()
    workflow: Workflow
    if workflow_name == "conventional":
        print("Run conventional workflow\n")
        workflow = ClassicWorkflow(autocas, interface)
        workflow.run()
    elif workflow_name == "large_cas":
        print("Run large cas workflow\n")
        workflow = LargeCasWorkflow(autocas, interface)
        workflow.run()
    elif workflow_name == "excited_states":
        print("Run excited states workflow\n")
        workflow = ClassicExcitedStatesWorkflow(autocas, interface)
        workflow.run()
    elif workflow_name == "large_cas_excited_states":
        print("Run excited states workflow\n")
        workflow = LargeCasExcitedStateWorkflow(autocas, interface)
        workflow.run()
    else:
        raise UnkownWorkflowError(f"Workflow {workflow} is not Available")
