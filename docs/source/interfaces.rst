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

Required input settings:

   - basis_set (the basis set)
   - spin_multiplicity (the spin multiplicity of the system)
   - charge (the charge of the system)
   - dmrg_sweeps  (the number of DMRG sweeps)
   - dmrg_bond_dimension (the DMRG bond dimension)
   - xyz_file (the path to the molecular XYZ file)
   - method (the active space method)
   - post_cas_method (The post-CAS method for dynamic correlation)
   - work_dir (the path for a work directory to store output files in)

Supported Interfaces
--------------------

OpenMolcas
..........

The `OpenMolcas` program package is interfaces as electronic structure back-end for AutoCAS.
Since OpenMolcas does not provide any Python bindings, it has to be called by a `subprocess`,
requiring an input file.
All OpenMolcas specific environment variables can be set trough the environment-directive, which will
override already set variables from the current session.

Additionally the OpenMolcas interface comes with more settings, namely:

   - point_group (the point group of the molecule)
   - cholesky (enable Cholesky decomposition)
   - uhf (enable unrestricted calculations)
   - fiedler (enable Fiedler ordering)
   - ipea (IPEA shift for CASPT2)
   - only_hf (if only HF should be calculated)
   - orbital_localisation (if orbitals should be localized)
   - localisation_space (the localization space)
   - localisation_method (the localization method)
   - n_excited_states (the number of excited states)
   - ci_size (the Davidson size for excited states)
   - weights (the weight for each root)
   - roots (which roots should be evaluated)

For further information on point_group, orbital localization or excited states, see the OpenMolcas manual.

