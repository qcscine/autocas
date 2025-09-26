Interfaces
==========

This section explains basics on the implementation of an interface to an electronic structue program.

Interfaces
----------

An interface provides a way to start calculations from within the AutoCAS framework. In order to
do so each interface is required to inherit from the base interface base class. This ensures that
every interface provides the ``calculate`` function and is able to be initialized by a dictionary.

The general AutoCAS workflow requires an interface to evaluate Hartree-Fock as well as DMRG for
a given active space.

An active space is defined by a list of orbital indices and a list of occupations. Hence, a calculation
function requires an active space as input and has to return an energy, s1 and s2-entropies and mutual
information.

Supported Interfaces
--------------------

OpenMolcas
..........

The `OpenMolcas` program package is interfaces as electronic structure back-end for AutoCAS.
Since OpenMolcas does not provide any Python bindings, it has to be called by a `subprocess`,
requiring an input file.
All OpenMolcas specific environment variables can be set trough the environment-directive, which will
override already set variables from the current session.

PySCF
.....

The `PySCF` program package is an electronic structure back-end for autoCAS.
Since `PySCF` is a python-native program package, the interface can be further utilized.
The interface stores at every step the current pyscf wave function object, which can be 
either further utilized or manipulated in other ways, outside of the pre-implemented workflows.
For information on how to work with pyscf, see the pyscf manual.

