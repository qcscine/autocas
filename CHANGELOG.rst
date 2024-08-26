Release History
===============

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

