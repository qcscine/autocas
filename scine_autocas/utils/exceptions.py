# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """


class InterfaceError(Exception):
    """Raise if Interface is not found"""


class ProjectError(Exception):
    """Raise if Project is not well structured"""


class UnkownActionError(Exception):
    """Raise if action is not available"""


class UnkownWorkflowError(Exception):
    """Raise if workflow is not available"""


class InputError(Exception):
    """Raise if input is wrong."""


class SingleReferenceException(Exception):
    """
    Raised when autoCAS determines that the system is single reference
    instead of multireference. This allows the client to catch the exception
    in case they want to perform a single reference calculation afterward.
    """


class DumperError(Exception):
    """Raise if Dump file does not exists"""
