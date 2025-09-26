# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
from argparse import SUPPRESS, ArgumentParser, _ArgumentGroup, _SubParsersAction
from typing import Any

from scine_autocas.utils.defaults import AvailableInteraces, Defaults


def __help_with_default(help_str: str, default: Any) -> Any:
    """Add default arguments to help string, but don't set them.

    Because argparse is like it is, providing a default value results into a value set
    in args, which makes it hard to destinguish between user input and default inputs.

    Parameters
    ----------
    help_str : str
        What the help string should say
    default : Any
        Default argument for the cli argument.

    Returns
    -------
    str
        A modified version of the help_str.
    """
    return f"{help_str} (default: {default})"


def __make_create_parser(action_parsers: _SubParsersAction):
    """Create the argument parser for the 'create' action.

    Parameters
    ----------
    action_parsers: _SubParsersAction
        a sub_parser of the global argumentparser dedicated for all available actions.
    """

    create_parser = action_parsers.add_parser(
        "create",
        help="Create a autocas yaml file",
        argument_default=SUPPRESS,
        description="""Create an autocas input yaml with all defaults.""")

    create_group = create_parser.add_argument_group("create agrs")
    create_group.add_argument(
        "-o", "--output_name",
        type=str,
        help="Path for the autocas yaml. Per DEFAULT it will be created in the current directory.",
        default=os.getcwd() + "/autocas.yaml",
    )
    __add_yaml_argument(create_parser)


def __make_initial_cas_parser(action_parsers: _SubParsersAction):
    """Create the argument parser for the 'initcas' action.

    All CLI arguments for an autocas calculation are defined here.

    Parameters
    ----------
    action_parsers: _SubParsersAction
        a sub_parser of the global argumentparser dedicated for all available actions.
    """

    initcas_parser = action_parsers.add_parser(
        "initcas",
        help="Show the valence cas.",
        argument_default=SUPPRESS,
        description="""Get the initial active space.
        The initial active space is defined as an array of indices and occupations.
        NOTE: No calculation is done/required to select this CAS.
        """)

    initcas_group = initcas_parser.add_argument_group("initcas agrs")
    initcas_group.add_argument(
        "-d", "--double_d_shell",
        action="store_false",
        help="Disable the inclusion of 4d orbitals for 3d transition metals.")
    initcas_group.add_argument(
        "-l", "--large_cas",
        action="store_true",
        help="Show active space distribution for the large active space protocol.")
    __add_molecule_argument(initcas_parser)
    __add_yaml_argument(initcas_parser)


def __make_run_parser(action_parsers: _SubParsersAction):
    """Create the argument parser for the 'run' action.

    All CLI arguments for an autocas calculation are defined here.

    Parameters
    ----------
    action_parsers: _SubParsersAction
        a sub_parser of the global argumentparser dedicated for all available actions.
    """

    run_parser = action_parsers.add_parser(
        "run",
        help="Run an autocas workflow.",
        description="""Run a autocas calculation with the given parameters.\n
For a fully automated autocas calculation, e.g. automated selection of the large cas protocol, please use the\n
'automated' action instead of 'run'.""",
        argument_default=SUPPRESS,
    )

    run_group = run_parser.add_argument_group("optional run arguments")
    run_group.add_argument(
        "-i", "--interface",
        help=__help_with_default(
            f"Available interfaces are {[i.value for i in AvailableInteraces]}",
            Defaults.Interface.interface
        ),
        # default=Defaults.Interface.interface,
    )
    run_group.add_argument(
        "-a", "--init_dmrg_bond_dimension",
        help=__help_with_default(
            "DMRG bond dimension for the active space selection.",
            Defaults.Interface.init_dmrg_bond_dimension,
        ),
    )
    run_group.add_argument(
        "-I", "--init_dmrg_sweeps",
        help=__help_with_default(
            "DMRG sweeps for the active space selection.",
            Defaults.Interface.init_dmrg_sweeps,
        ),
    )
    run_group.add_argument(
        "-d", "--dump",
        action="store_true",
        help=__help_with_default(
            "Enable dumping of calculation results.",
            Defaults.Interface.dump,
        ),
    )
    run_group.add_argument(
        "-u", "--uhf",
        action="store_true",
        help=__help_with_default(
            "Force unrestricted Hartree-Fock for initial orbitals.",
            Defaults.Interface.uhf,
        ),
    )
    run_group.add_argument(
        "-l", "--large_cas",
        action="store_true",
        help="Enable active space distribution with the large active space protocol.")
    __add_molecule_argument(run_parser)
    __add_yaml_argument(run_parser)


def __make_analyze_parser(action_parsers: _SubParsersAction):
    """Create the argument parser for the 'analyze' action.

    All CLI arguments in order to analyze an autocas project are defined here.

    Parameters
    ----------
    action_parsers: _SubParsersAction
        a sub_parser of the global argumentparser dedicated for all available actions.
    """
    analyze_parser = action_parsers.add_parser(
        "analyze",
        help="Analyze an existing project",
        argument_default=SUPPRESS,
        description="""Analyze an existing autocas project.
        NOTE: Some properties may not be analyzed, if the autocas projected did not enable 'dump'.""")

    analyze_group = analyze_parser.add_argument_group("analyze agrs")
    analyze_group.add_argument(
        "-d",
        "--dir",
        help="autocas project dir to analyze.",
    )
    __add_molecule_argument(analyze_group)
    __add_yaml_argument(analyze_group)


def __add_molecule_argument(parser: _ArgumentGroup):
    """Create the argument parser for the general autocas inputs.

    Global or shared arguemnts can be specified here, in order to make them available to ALL
    other programs. So only arguments, required by all other actions (but dmrg) should be specified here.

    Parameters
    ----------
    action_parsers: _SubParsersAction
        a sub_parser of the global argumentparser dedicated for all available actions.
    """
    molecule_group = parser.add_argument_group("molecule")
    molecule_group.add_argument(
        "-x", "--xyz-file",
        help="""Path to the molecular structure file.""")

    molecule_group.add_argument(
        "-s", "--spin_multiplicity",
        help=__help_with_default(
            "Total spin multiplicity (2S+1) of the molecule.",
            Defaults.Molecule.spin_multiplicity
        )
    )
    molecule_group.add_argument(
        "-c", "--charge",
        help=__help_with_default(
            "Total charge of the molecule.",
            Defaults.Molecule.charge
        ),
    )
    molecule_group.add_argument(
        "-b", "--basis_set",
        help=__help_with_default(
            "Basis set for all calculations.",
            Defaults.Interface.basis_set
        )
    )
    molecule_group.add_argument(
        "-U", "--unit",
        choices=("bohr", "ang"),
        help=__help_with_default(
            "Length unit used in the molecular structure file.",
            Defaults.Molecule.unit
        ),
    )
    molecule_group.add_argument(
        "-E", "--ecp_electrons",
        help=__help_with_default(
            "Number of electrons in core potential",
            Defaults.Molecule.ecp_electrons
        ),
    )
    molecule_group.add_argument(
        "-m", "--molden_file",
        help="""Use mo coeffs from molden file.""")


def __add_yaml_argument(global_parser: ArgumentParser):
    """Create the argument parser for the general autocas inputs.

    Global or shared arguemnts can be specified here, in order to make them available to ALL
    other programs. So only arguments, required by all other actions (but dmrg) should be specified here.

    Parameters
    ----------
    action_parsers: ArgumentParser
        a sub_parser of the global argumentparser dedicated for all available actions.
    """

    global_parser.add_argument(
        "-y", "--yaml",
        type=str,
        help="""Use yaml input file.
        NOTE: Command line options overwrite settings set in the yaml file!""")


def make_parser() -> ArgumentParser:
    """Create the argument parser for the autocas CLI.

    The returned argument parser included also all subparsers for each defined action.
    Currently supported actions are: 'run', 'analyze' and 'dmrg'. See the corresponding
    <__make_'action'_parser> function for a detailed explanation.

    Returns
    -------
    action_parsers: ArgumentParser
        a parser for all CLI arguemnts and actions.
    """
    global_parser = ArgumentParser(
        prog="scine_autocas",
    )

    __add_yaml_argument(global_parser)

    action_parsers = global_parser.add_subparsers(
        title="Actions",
        help="Actions define the applied workflow in AutoCAS.",
        dest="action")

    __make_run_parser(action_parsers)
    __make_initial_cas_parser(action_parsers)
    __make_analyze_parser(action_parsers)
    __make_create_parser(action_parsers)

    return global_parser
