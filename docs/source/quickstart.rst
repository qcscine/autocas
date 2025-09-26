Quickstart
----------

After installing AutoCAS it can be started from the command line. To show all possible options, please run:

.. code-block:: bash

   scine_autocas -h

For example, AutoCAS can be started by passing a valid XYZ file to it, and running all calculations with the corresponding defaults.

.. code-block:: bash

   scine_autocas run -x <molecule.xyz>

To pass a basis set, a different interface or enable the creation of entanglement diagrams
the following directives can be passed:

.. code-block:: bash

   scine_autocas run --xyz_file <molecule.xyz> --basis_set cc-pvtz 

However we would strongly recommend providing a ``.yml``-input file, to make calculations reproducible and
allowing higher customization of autoCAS. A ``.yml``-input including all default values, can be generated via:

.. code-block:: bash

   scine_autocas create

AutoCAS provides a functionality top re-plot the entanglement and threshold diagram for a finished, existing, valid
autoCAS project through:

.. code-block:: bash

   scine_autocas analyze -y <yaml_config.yml> -x <molecule.xyz>

Quickstart for the Consistent Active Space Protocol
---------------------------------------------------
To use the consistent active space protocol, you can run the following command:

.. code-block:: bash

   scine_autocas_consistent_active_space -y input.yaml

or

.. code-block:: bash

    scine_autocas_consistent_active_space -i 1 n2_0.xyz n2_1.xyz

Example yaml files are given in `scine_autocas/tests/workflows/consistent_active_space/files`. Furthermore,
when calling `scine_autocas_consistent_active_space` without a yaml file, the corresponding yaml file is
created.
