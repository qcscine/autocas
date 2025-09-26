.. set_up

Setting up AutoCAS
==================

AutoCAS itself does not provide any utilities to perform quantum chemical calculations. It relies on
interfaces to quantum chemistry packages. Currently,  the packages `OpenMolcas`, `PySCF` and `Serenity`
in combination with `QCMaquis` are supported (and it is easy to set up a new interface - see section
"New Interfaces").

Therefore, in order to use autoCAS, you need to install the autoCAS package itself, at least one of the
quantum chemistry codes mentioned above, as well as QCMaquis. We provide information how to install
these packages in the following sections.

.. include:: installation.rst
