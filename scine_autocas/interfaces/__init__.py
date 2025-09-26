"""
Implemented Interfaces are:
    Molcas, PyscfInterface
    Experimental: ChronusQ, Serenity

Examples
--------
>>> from scine_autocas.interfaces import Molcas
>>> from scine_autocas import Molecule
>>> mol = Molecule("/path/to/mol.xyz")
>>> interface = Molcas(mol)
>>> # initial orbital generation
>>> e_hf = interface.calculate()
>>> cas_occ = [2, 2, 2, 0, 0, 0]
>>> cas_ind = [3, 4, 9, 10, 12, 15]
>>> # initial cas calculation
>>> e_cas, s1, s2, mut_inf = interface.calculate(cas_occ, cas_ind)
"""
# -*- coding: utf-8 -*-
__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""


from .chronusq import ChronusQ
from .molcas import Molcas
from .pyscf import PyscfInterface
from .serenity import Serenity
