Release History
===============

3.0.0
-----

- Updated the interface to Serenity to be compatible with Serenity 1.6.3.
- Added interface to PySCF in combination with QCMaquis. 
- Updated Large-CAS protocol to properly handle S>0 cases.
- Updated requirements
- Improved command line interface
- Added a requirements file for Serenity.
- The consistent active space protocol is now available directly on the command line
  by calling `scine_autocas_consistent_active_space`. Options for the protocol can be
  provided as command line arguments or as a yaml file.
- AutoCAS now stops if an error is detected in an underlying MOLCAS calculation.

2.3.1
-----

- Removed unused package from requirements file

2.3.0
-----
- Fixed a typo in the Serenity interface.
- Updated default orbital mapping approach in `consistent_active_orbital_spaces.py` and
  added an option to always exclude core orbitals from the active space if they cannot
  be mapped consistently.

2.2.0
-----

- Raise specific exception when single-reference case is encountered
- Update address in license
- Fix typo in Serenity interface

2.1.0
-----

- Added an interface to Serenity that provides a map, mapping orbitals between multiple structures.
  This map can be used to ensure a consistent CAS along reaction coordinates.

2.0.0
-----

Initial Features of the Python Versions:
   - The standard AutoCAS algorithm to select active spaces.
   - Selection of a general active space for excited states.
   - The large active space protocol for ground and excited states.

