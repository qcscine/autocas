Quickstart
----------

After installing AutoCAS it can be started from the command line. To show all possible options, please run:

.. code-block:: bash

   python3 -m scine_autocas -h

For example, AutoCAS can be started by passing a valid XYZ file to it, and running all calculations with the corresponding defaults.

.. code-block:: bash

   python3 -m scine_autocas -x <molecule.xyz>

To pass a basis set, a different interface or enable the creation of entanglement diagrams
the following directives can be passed:

.. code-block:: bash

   python3 -m scine_autocas --xyz_file <molecule.xyz> --basis_set cc-pvtz --plot --interface Molcas

However we would strongly recommend providing a ``.yml``-input file, to make calculations reproducible and
allowing higher customization of AutoCAS.

