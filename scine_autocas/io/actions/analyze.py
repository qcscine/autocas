# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""
from typing import Union

from scine_autocas.cas_selection import Autocas
from scine_autocas.io import FileHandler
from scine_autocas.io.input_handler import InputHandler
from scine_autocas.utils.exceptions import ProjectError, UnkownWorkflowError
# from scine_autocas.utils import molecule
from scine_autocas.utils.molecule import molecule_setup_from_dict
from scine_autocas.workflows.analysis import AnalyzeConventional, AnalyzeLargeCAS


def _auto_check_workflow(workflow: Union[str, None]):
    """Automatically check folder structure to determine workflow."""
    if workflow is None:
        if FileHandler.check_dir_exists("dmrg"):
            workflow = "conventional"
        elif FileHandler.check_dir_exists("dmrg_1"):
            workflow = "large_cas"
        else:
            workflow = "single_folder"
    return workflow


def _get_autocas(inp: InputHandler):
    """Get correct autocas object."""
    mol = molecule_setup_from_dict(inp.get_molecule_options())
    print(f"mol: {mol}")
    autocas = Autocas(mol)
    autocas.setup_from_dict(inp.get_autocas_options())
    return autocas


def analyze_action(inp: InputHandler):
    """run action"""
    if not FileHandler.check_project_dir_exists():
        raise FileNotFoundError(
            f"""AutoCAS project dir: {FileHandler.project_dir} does not exist!
            You can only analyze existing AutoCAS dirs.""")

    workflow = inp.get_workflow()
    print(f"workflow from inp {workflow}")
    workflow = _auto_check_workflow(workflow)
    print(f"workflow after check {workflow}")

    autocas = _get_autocas(inp)

    print("")
    if workflow == "conventional":
        try:
            if not FileHandler.check_dir_exists("dmrg"):
                raise ProjectError("""AutoCAS project not well formated. No dmrg dir found.""")
            print("analyze conventional autocas dir")
            conv_analyzer = AnalyzeConventional(autocas)
            conv_analyzer.run()
        except ProjectError as exc:
            if not FileHandler.check_dir_exists("dmrg_1"):
                raise ProjectError("AutoCAS project not well formated. No dmrg_1 dir found.") from exc
            print("analyze large cas autocas dir")
            large_cas_analyzer = AnalyzeLargeCAS(autocas)
            large_cas_analyzer.run()

    elif workflow == "large_cas":
        try:
            if not FileHandler.check_dir_exists("dmrg_1"):
                raise ProjectError("AutoCAS project not well formated. No dmrg_1 dir found.")
            print("analyze large cas autocas dir")
            large_cas_analyzer = AnalyzeLargeCAS(autocas)
            large_cas_analyzer.run()
        except ProjectError as exc:
            if not FileHandler.check_dir_exists("dmrg"):
                raise ProjectError("""AutoCAS project not well formated. No dmrg dir found.""") from exc
            print("analyze conventional autocas dir")
            conv_analyzer = AnalyzeConventional(autocas)
            conv_analyzer.run()

    else:
        raise UnkownWorkflowError(f"Workflow {workflow} is not available for analysis")
