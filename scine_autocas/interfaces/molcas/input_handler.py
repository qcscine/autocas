"""Handle molcas input.

This module implements a class to write molcas input files from
a molcas settings object.
"""
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import io
import os
# rename class and corresponding file
from typing import Any, Optional

from scine_autocas.utils.defaults import CasMethods, PostCasMethods


class InputHandler:
    """Create Molcas input files."""

    def __initial_orbitals(self, settings: Any, input_file: io.IOBase):
        """Create Molcas input to evaluate initial orbitals.

        Parameters
        ----------
        settings : Molcas.Settings
            stores all molcas related settings
        input_file : io.IOBase
            The already opened input file
        """
        input_file.write("&GATEWAY\n")
        input_file.write(f"  COORD = {settings.get_molecule().xyz_file}\n")
        input_file.write(f"  BASIS = {settings.basis_set}\n")
        if settings.point_group != "":
            input_file.write(f"  GROUP = {settings.point_group}\n")
        input_file.write("\n&SEWARD\n")
        if settings.cholesky:
            input_file.write("  CHOLesky\n")
        input_file.write("\n&SCF\n")
        if settings.uhf or settings.get_molecule().spin_multiplicity != 1:
            input_file.write("  uhf\n")
        # if settings.skip_scf:
        #     input_file.write("  ITERations = 1\n")
        #     input_file.write("  THREshold  = 1e+9 1e+9 1e+9 1e+9\n")

        # multiplicity, e.g. 1 (singlet), 2 (doublet), ...
        input_file.write(f"  SPIN   = {settings.get_molecule().spin_multiplicity}\n")
        input_file.write(f"  CHARGE = {settings.get_molecule().charge}\n")
        if settings.cholesky:
            input_file.write("  CHOLesky\n")
        if settings.orbital_localisation:
            input_file.write("\n&LOCALISATION\n")
            input_file.write(f"  {settings.localisation_space}\n")
            input_file.write(f"  {settings.localisation_method}\n")

    def __casscf(self, settings: Any, input_file: io.IOBase, orbital_file: Optional[str], alter: Optional[str]):
        """Create Molcas input for casscf calculations.

        Parameters
        ----------
        settings : Settings
            the settings object from molcas
        file_name : str
            name of the input file
        orbital_file : str
           path to an modified orbital file. If specific orbitals are selected they are encoded in the
           typeindex in the orbital_file.
        alter : str
            a string to reorder orbitals, instead of doing that in the orbital file
        """
        input_file.write("&RASSCF\n")
        input_file.write(f"  SPIN    = {settings.get_molecule().spin_multiplicity}\n")
        input_file.write(f"  NACTEL  = {settings.active_electrons}\n")
        if orbital_file:
            input_file.write(f"  FILEORB = {os.path.basename(orbital_file)}\n")
            if not alter:
                input_file.write("  TYPEINDEX\n")
            else:
                input_file.write(f"  ALTEr={alter}\n")
        else:
            input_file.write(f"  RAS2    = {settings.active_orbitals}\n")
            if alter:
                input_file.write(f"  ALTEr={alter}\n")
        if settings.cas_method in (CasMethods.DMRGCI, CasMethods.CASCI):
            input_file.write("  CIONLY\n")
        if settings.ci_root_string != "":
            input_file.write(f"  CIRoot = {settings.ci_root_string}\n")

    def __dmrg(self, settings: Any, input_file: io.IOBase, orbital_file: Optional[str], alter: Optional[str]):
        """Create Molcas input for qcmaquis dmrg calculations.

        Parameters
        ----------
        settings : Settings
            the settings object from molcas
        file_name : str
            name of the input file
        orbital_file : str
           path to an modified orbital file. If specific orbitals are selected they are encoded in the
           typeindex in the orbital_file.
        alter : str
            a string to reorder orbitals, instead of doing that in the orbital file
        """
        input_file.write("&DMRGSCF\n")
        # Ensure formatted starting orbitals are used by Molcas DMRG
        input_file.write("LUMORB\n")
        input_file.write("StartOrb = INPORB\n")
        if settings.fiedler:
            input_file.write("  FIEDLER=ON\n")
        input_file.write("ActiveSpaceOptimizer = QCMaquis\n")
        input_file.write("DMRGSettings\n")
        input_file.write(f"  nsweeps = {settings.init_dmrg_sweeps}\n")
        input_file.write(f"  max_bond_dimension = {settings.init_dmrg_bond_dimension}\n")
        input_file.write("EndDMRGSettings\n")
        input_file.write("OOptimizationSettings\n")
        input_file.write(f"  SPIN    = {settings.get_molecule().spin_multiplicity}\n")
        input_file.write(f"  NACTEL  = {settings.active_electrons}\n")
        if orbital_file:
            input_file.write(f"  FILEORB = {os.path.basename(orbital_file)}\n")
            if not alter:
                input_file.write("  TYPEINDEX\n")
            else:
                input_file.write(f"  ALTEr={alter}\n")
        else:
            input_file.write(f"  RAS2    = {settings.active_orbitals}\n")
            if alter:
                input_file.write(f"  ALTEr={alter}\n")
        if settings.cas_method in (CasMethods.DMRGCI, CasMethods.CASCI):
            input_file.write("  CIONLY\n")
        if settings.ci_root_string != "" and settings.ci_root_string is not None:
            input_file.write(f"  CIRoot = {settings.ci_root_string}\n")
        if settings.get_molecule().spin_multiplicity > 1:
            pass
        input_file.write("EndOOptimizationSettings\n")

    def write_input(
        self, settings: Any, file_name: str,
        orbital_file: Optional[str] = None, alter: Optional[str] = None,
    ):
        """Write a basic Molcas input file, with respect of the strings in
        self.methods.

        Parameters
        ----------
        settings : Settings
            the settings object from molcas
        file_name : str
            name of the input file
        orbital_file : str
           path to an modified orbital file. If specific orbitals are selected they are encoded in the
           typeindex in the orbital_file.
        alter : str
            a string to reorder orbitals, instead of doing that in the orbital file
        """
        print()
        with open(file_name, "w", encoding="utf-8") as input_file:
            if settings.initial_orbitals:
                self.__initial_orbitals(settings, input_file)

            if not settings.initial_orbitals:
                # casscf
                if settings.cas_method in (CasMethods.CASSCF, CasMethods.CASCI):
                    self.__casscf(settings, input_file, orbital_file, alter)

                # dmrgscf
                if settings.cas_method in (CasMethods.DMRGSCF, CasMethods.DMRGCI):
                    self.__dmrg(settings, input_file, orbital_file, alter)

                try:
                    # caspt2
                    if settings.post_cas_method == PostCasMethods.CASPT2:
                        input_file.write("\n&CASPT2\n")
                        input_file.write(f"  IPEA = {settings.ipea}\n")
                    # nevpt2
                    if settings.post_cas_method == PostCasMethods.NEVPT2:
                        input_file.write("\n&NEVPT2\n")
                except AttributeError:
                    pass
