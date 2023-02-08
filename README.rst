SCINE - AutoCAS
===============

.. image:: docs/source/_static/autocas_logo.svg
   :align: center
   :alt: SCINE AutoCAS

.. inclusion-marker-do-not-remove

Introduction
------------

SCINE autoCAS automates the crucial active-orbital-space selection step in multi-configurational calculations. Based on orbital entanglement
measures derived from an approximate DMRG wave function, it identifies all strongly correlated orbitals to be included in the active space
of a final, converged calculation. All steps can be carried out in a fully automated fashion.

Installation
------------

Currently autoCAS can be installed via pip or manually with git and pip (see sections below).

Basic Requirements
..................

AutoCAS utilizes a couple of third-party packages, which are defined in the ``requirements.txt``.

All requirements are automatically installed by installing autoCAS over one of the following methods.

pip
...

Prerequisites
``````````````

    #. python3.6+

Install
```````
This methods allows you to install the package via pip.

.. code-block:: bash

   pip install scine-autocas

Git + pip
.........

Prerequisites
``````````````

    #. git
    #. python3.6+

Install
```````
This methods requires you to first clone the repository and install the package over pip.

.. code-block:: bash

   git clone <autocas-repo>
   cd autoCAS
   pip install -r requirements.txt
   pip install .

Set up OpenMolcas
-----------------
Up to this point `OpenMolcas` does not provide any way to call itself from a native Python interfaces.
Hence, autoCAS calls `OpenMolcas` directly over `pymolcas`. In order to do so, please set the environment
variable `MOLCAS` pointing to the ``build`` directory.

.. code-block:: bash

   export MOLCAS=/path/to/Molcas/build

Quickstart
----------

After installing autoCAS it can be started from the command line. To show all possible options, please run:

.. code-block:: bash

   python3 -m scine_autocas -h

For example, autoCAS can be started by passing a valid XYZ file to it, and running all calculations with the corresponding defaults.

.. code-block:: bash

   python3 -m scine_autocas -x <molecule.xyz>

To pass a basis set, a different interface or enable the creation of entanglement diagrams
the following directives can be passed:

.. code-block:: bash

   python3 -m scine_autocas --xyz_file <molecule.xyz> --basis_set cc-pvtz --plot --interface Molcas

However we would strongly recommend providing a ``.yml``-input file, to make calculations reproducible and
allowing higher customization of autoCAS.

License and Copyright Information
---------------------------------

AutoCAS is distributed under the BSD 3-clause "New" or "Revised" License. For more
license and copyright information, see the file ``LICENSE.txt`` in the repository.

How to Cite
-----------

When publishing results obtained with the autoCAS program, please cite the corresponding
release as archived on `Zenodo <https://zenodo.org>`_ (please use the DOI of
the respective release) and the following publications:

* Primary reference:
  C. J. Stein and M. Reiher, "autoCAS: A Program for Fully Automated Multiconfigurational Calculations", *J. Comput. Chem.*, **2019**, *40*, 2216-2226.

* Original presentation of the approach:
  C. J. Stein and M. Reiher, "Automated Selection of Active Orbital Spaces”", *J. Chem. Theory Comput.*, **2016**, *12*, 1760.

* Automated active space selection with multi-reference perturbation theory:
  C. J. Stein, V. von Burg and M. Reiher, "The Delicate Balance of Static and Dynamic Electron Correlation", *J. Chem. Theory Comput.*, **2016**, *12*, 3764.

* Multi-configurational diagnostic:
  C. J. Stein and M. Reiher, "Measuring Multi-Configurational Character by Orbital Entanglement", *Mol. Phys.*, **2017**, *115*, 2110.

* Excited states and reaction paths:
  C. J. Stein and M. Reiher, "Automated Identification of Relevant Frontier Orbitals for Chemical Compounds and Processes", *Chimia*, **2017**, *71*, 170.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message to
autocas@phys.chem.ethz.ch.

