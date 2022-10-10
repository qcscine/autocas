Plotting
========

Even though AutoCAS can select active spaces without any required user interaction and provides scripts
to run CAS calculations in a black-box fashion, it is still important to monitor calculations and
verify results.
In order to do so, AutoCAS comes with ways to generate entanglement diagrams as well as threshold diagrams
to verify the automatically selected active space.

For further information, see the corresponding articles.

Entanglement Diagrams
---------------------

An entanglement diagram shows the single orbital entropy from an DMRG calculation as scaled dots, which are
connected with each other through lines. These lines correspond to the mutual information and are scaled by
the corresponding values.

Advanced Diagrams
.................

In case also single orbital entropies for the final calculation are provided, a custom entanglement diagram
can be created. This diagram shows the entropies and mutual information in an outer circle, while the final
single orbital entropies are shown in an inner diagram.

Threshold Diagrams
------------------

The main way to verify AutoCAS algorithms is a threshold diagram, since internally AutoCAS active
space search utilizes an threshold based algorithm.
A threshold diagram plots the amount of orbitals against the percentage with respect to the maximal s1-entropy.
