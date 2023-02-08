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
Hence, AutoCAS calls `OpenMolcas` directly over `pymolcas`. In order to do so, please set the environment
variable `MOLCAS` pointing to the ``build`` directory.

.. code-block:: bash

   export MOLCAS=/path/to/Molcas/build
