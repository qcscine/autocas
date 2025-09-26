# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details. """
import os
from typing import Dict


class Environment:
    """Handle the environment for qcmaquis variables."""

    __slots__ = (
        "environment",
        "transform_binary_name",
        "qcmaquis_binary_path",
    )

    def __init__(self, qcmaquis_path: str = ""):
        """Initialize environment and read environment variable.

        Parameters
        ----------
        qcmaquis_path : str, optional
            The path to the qcmaquis binary

        Notes
        -----
        If not path provided, the path will be read from the environment variable
        QCMAQUIS_BINARY_DIR
        """
        self.environment: Dict[str, str]
        """Environment dict"""
        self.transform_binary_name: str = "mps_transform_pg"
        """name of the qcmaquis transform binary"""
        self.qcmaquis_binary_path: str = qcmaquis_path
        """path to all qcmaquis binaries"""

        self._make_environment()
        self._check_if_mps_transform_works()

    def _make_environment(self):
        """Make the environment for qcmaquis binaries.

        If not qcmaquis_path is provided, try to read the environment variable
        <QCMAQUIS_BINARY_DIR>, else raise Error.

        Raises
        ------
        OSError
            if not path to qcmaquis is provided, or environemnt variable is not set.
        """
        self.environment = os.environ.copy()
        if self.qcmaquis_binary_path != "":
            print(
                f"""Reading QCMaquis binary path from variable self.qcmaquis_binary_path as
                {self.qcmaquis_binary_path}"""
            )
            return

        if "QCMAQUIS_BINARY_DIR" in self.environment:
            self.qcmaquis_binary_path = self.environment["QCMAQUIS_BINARY_DIR"]
            print(
                f"""Reading QCMaquis binary path from environment variable QCMAQUIS_BINARY_DIR as
                {self.qcmaquis_binary_path}"""
            )
            return

        raise OSError(
            """QCMaquis is not installed or 'QCMAQUIS_BINARY_DIR' system variable is not set.
            Please use the following command: \n
            export QCMAQUIS_BINARY_DIR=/path/to/qcmaquis/dir \n
            and set the correct path to the qcmaquis binary directory.
            """
        )

    def _check_if_mps_transform_works(self):
        """Check if mps_transform_pq binary is available and works

        Raises
        ------
        OSError
            if not mps_transform_pg binary is not accessible
        """
        if os.path.isfile(self.get_transform_binary()) and os.access(self.get_transform_binary(), os.X_OK):
            print(f"""Found <mps_transform_pg> binary in {self.get_transform_binary()}""")
        else:
            print(f"qcmaquis_binary_path variable is set to {self.qcmaquis_binary_path}")
            raise OSError(
                """QCMaquis <mps_transform_pg> is not accessible.
                    Please check if QCMaquis is installed and compiled with the <mps_transform_pg> binary enabled and
                    compiled.
                    Otherwise, check if the environment variable <QCMAQUIS_BINARY_DIR> is set to the correct path,
                    or the qcmaquis_binary_path variable is set correctly."""
            )

    def get_transform_binary(self) -> str:
        """Get full path to mps_transform_pg binary

        Returns
        -------
        path : str
            full path to the binary
        """
        return self.qcmaquis_binary_path + "/" + self.transform_binary_name
