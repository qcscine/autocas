__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
from typing import Any, Dict


def get_all_params(obj: Any) -> Dict[str, Any]:
    """Get a Dict with all parameters from an object or a class."""
    return_dict = {}
    try:
        return_dict[obj.__name__] = _get_all_params_impl(obj)
    except AttributeError:
        return_dict[type(obj).__name__] = _get_all_params_impl(obj)
    return return_dict


def _get_all_params_impl(obj: Any) -> Dict[str, Any]:
    """Get a Dict from any object with all sub classes.

    Parameters
    ----------
    obj : Any
        Any python object or class.

    Returns
    -------
    return_dict : Dict[Any, str]
        nested dict with all subclasses as dict
    """
    return_dict: Dict[str, Any] = {}
    storage_string = ""

    if hasattr(obj, "__dict__"):
        loopy_doopy = obj.__dict__
        storage_string = "__dict__"
    else:
        loopy_doopy = obj.__slots__
        storage_string = "__slots__"

    # check for strings
    if isinstance(loopy_doopy, str):
        loopy_doopy = list(loopy_doopy)

    for i in loopy_doopy:
        if not i.startswith("__") and not i.endswith("__"):
            # the second check is to omitt enums
            if hasattr(loopy_doopy[i], storage_string) and hasattr(loopy_doopy[i], "__name__"):
                return_dict[i] = _get_all_params_impl(loopy_doopy[i])
            else:
                if hasattr(getattr(obj, i), "__dict__"):
                    return_dict[i] = get_all_params(loopy_doopy[i])
                else:
                    return_dict[i] = getattr(obj, i)
    return return_dict


def set_all_params(obj: Any, parameter_dict: Dict[str, Any]) -> Any:
    """Set all parameters to an object.

    Parameters
    ----------
    obj : Any
        any python object
    parameter_dict : Dict[str, Any]
        a dict with all parameters and parameters to subclasses
        the outest dict has only the name of the classes object

    Returns
    -------
    obj : Any
        the object with all parameters set.
    """
    try:
        object_name = obj.__name__
    except AttributeError:
        object_name = type(obj).__name__

    if hasattr(obj, "__dict__"):
        loopy_doopy = obj.__dict__
        storage_string = "__dict__"
    else:
        loopy_doopy = obj.__slots__
        storage_string = "__slots__"

    # check for strings
    if isinstance(loopy_doopy, str):
        loopy_doopy = [loopy_doopy]

    for i in loopy_doopy:
        if not i.startswith("__") and not i.endswith("__"):
            # the second check is to omitt enums
            if hasattr(loopy_doopy[i], storage_string):  # and hasattr(loopy_doopy[i], "__name__"):
                if hasattr(loopy_doopy[i], "__name__"):
                    setattr(obj, i, set_all_params(loopy_doopy[i], parameter_dict[object_name]))
                else:
                    setattr(obj, i, set_all_params(loopy_doopy[i], parameter_dict[object_name][i]))
            else:
                setattr(obj, i, parameter_dict[object_name][i])
    return obj
