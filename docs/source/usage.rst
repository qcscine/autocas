How to use AutoCAS
==================

.. include:: quickstart.rst

Advanced Usage
--------------

YAML Input
..........

Instead of relying on the provided input information from the command line, AutoCAS can be fully
controlled by passing a ``.yml``-input file to it:

.. code-block:: bash

   python3 -m scine_autocas -y <input.yml>

A minimal ``.yml``-input requires to have at least a molecule section, with the location of the corresponding
XYZ file.

.. code-block:: yaml

   molecule:
     xyz_file: "/path/to/molecule.xyz"

Most options from the ``.yml``-input are populated through the whole AutoCAS framework, meaning
every variable can be set and further it can be used to start from any point in an AutoCAS calculation.
A comprehensive ``.yml``-file could look like:

.. literalinclude:: ../../scripts/full.yml
   :language: yaml

For more information on keywords, take a look at the API section.

Custom Scripts
..............

Until now, everything discussed is based on the front-end of AutoCAS implemented in the usual ``__main__.py`` in
``scine_autocas/``. However, AutoCAS comes as a Python3 library, which can be utilized on its own,
writing custom workflows incorporating the AutoCAS framework.
A basic script to set up an AutoCAS-based calculation could look like:

.. literalinclude:: ../../scripts/example.py
   :language: python

More scripts can be found in ``/path/to/autoCAS/scripts``

Analysing Completed Projects
............................

In order to analyze already finished calculations, AutoCAS provides a script to do so. The analysis
can be accomplished by:

.. code-block:: bash

    python3 /path/to/autoCAS/scripts/analyse_only.py

In order to make this task more convenient, we suggest adding an alias to your ``.bashrc``, for example:

.. code-block:: bash

   echo -e "alias entanglement='python3 /path/to/autoCAS/scripts/analyse_only.py'" >> ~/.bashrc
   source .bashrc

The analysis script requires a `QCMaquis` output (assuming the alias set before):

.. code-block:: bash

   entanglement qcmaquis_output.results_state.0.h5

and can optionally save the entanglement diagram, select which state to analyze or set the molecule to
get an active space suggestion from AutoCAS. Instead of directly analyzing only one output, the script
can be applied to the root of the project, e.g. a folder containing ``initial/``, ``dmrg/``, ``final/``.
Use the following command to save the entanglement diagram and analyze the ground state for a molecule:

.. code-block:: bash

   entanglement -s entanglement_diagram.pdf -e 0 -m <molecule.xyz> <autocas_dir>

Note: ensure either one ``dmrg/`` folder is provided, or a number of folders ``dmrg_1/``, ``dmrg_2/``, etc. from
a large active space calculation, but **not both**.
For more information run:

.. code-block:: bash

    python3 /path/to/autoCAS/scripts/analyse_only.py -h

