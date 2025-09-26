Installation
------------

Currently AutoCAS can be installed via pip or manually with git and pip (see sections below).

Basic Requirements
..................

AutoCAS utilizes a couple of third-party packages, which are defined in the ``requirements.txt``.

All requirements are automatically installed by installing AutoCAS over one of the following methods.

pip
...

Prerequisites
``````````````

    #. python3.10

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
    #. python3.10 (other Python versions may also work, but are not tested and may require a manual update of the requirements file)

Install
```````
This methods requires you to first clone the repository and install the package over pip.

.. code-block:: bash

   git clone <autocas-repo>
   cd autoCAS
   pip install -r requirements.txt
   pip install .

Building the Documentation
--------------------------

The documentation can be found online, or it can be built using:

.. code-block:: bash

   pip install -r requirements-dev.txt
   make -C docs html

It is then available at `docs/build/html/index.html`.

Set up OpenMolcas
-----------------
Up to this point `OpenMolcas` does not provide any way to call itself from a native Python interfaces.
Hence, AutoCAS calls `OpenMolcas` directly over `pymolcas`. In order to do so, please set the environment
variable `MOLCAS` pointing to the corresponding directory.

.. note::
   We tested the OpenMolcas interface with the OpenMolcas version v23.10.

.. note::
   AutoCAS relies on `QCMaquis` to generate the required entropies. Hence, `OpenMolcas` has to be compiled
   with the `DMRG` flag enabled (this flag downloads and compiles `QCMaquis` and enables the interface in `OpenMolcas`).

.. code-block:: bash

   export MOLCAS=/path/to/pymolcas_binary

Set up PySCF
------------
`PySCF` is a python-native program package, which is utilized in the interface.
Hence, `PySCF` can be installed via pip:

.. code-block:: bash

   pip install pyscf

The interface utilizes the `QCMaquis` program package for DMRG calculations. The DMRG solver is connected to `PySCF` through Python bindings, which have to be enabled in the installation of `QCMaquis`. See the QCMaquis manual for more information.

.. note::
   We tested the PySCF interface with the PySCF version 2.4.0.

Set up Serenity
---------------

To use Serenity (e.g., for the consistent CAS protocol; version 1.6.3 is tested), please install the Python bindings
of Serenity. For example, you can do this with

    pip install qcserenity

or

    pip install -r requirements-serenity.txt

Alternatively, you can also compile Serenity from source. In this case, don't forget to source the `serenity.sh` script to add the Python bindings to your `PYTHONPATH`.

Note that for the consistent active space protocol, you also need to install OpenMolcas with the QCMaquis interface (see above).
