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

   scine_autocas -y <input.yml> -x <molecule.xyz>

.. note::
   A ``.yml``-file overwrites the defaults with the corresponding values. However, all provided command-line options overwrite 
   the set options in the ``.yml`` file.

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

