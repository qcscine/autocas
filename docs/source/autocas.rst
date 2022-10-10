AutoCAS Design and Terminologies
================================

The following sections describe the most important workflows. For further information see the articles
provided on the main page.

Workflow
--------

Ground State Active Space Search
................................

AutoCAS is completely separated from any electronic structure program and only provides
methods to search for active spaces. An active space is defined as a list of orbital indices and
occupations.
AutoCAS' default algorithm to search for an active space for the ground state consists of two steps.
First an initial active space is selected, the valence space, based on the molecule of interest.
After successfully evaluating the corresponding s1-entropies, autoCAS provides a method to select
the final active space on these. Hence the only required input for AutoCAS is a molecule (either a
list with atom labels, or a XYZ file) and for the active space selection the single orbital
entropies are required.

Advanced Workflows
------------------

Large Active Space Protocol
...........................

The single orbital entropies for AutoCAS are evaluated based on an unconverged DMRG wave function
with a low bond dimension. To allow large initial active spaces, e.g. valence spaces with more
than 100 orbitals, the large active space protocol can be utilized. In this protocol, the active
space is first divided into occupied and virtual space. Subsequently, these spaces are further
divided into small sub-spaces, consisting of a feasible number of orbitals.
Afterwards the occupied and virtual sub-spaces are recombined into many small active spaces,
which are evaluated as usual by DMRG. The resulting s1-entropies from each small active space
are then recombined, to get an approximate s1-entropy for the large active space, which can
be used as usual by AutoCAS to select the active space.

For a comprehensive explanation of the algorithm, see the corresponding article:
  C. J. Stein and M. Reiher, "autoCAS: A Program for Fully Automated Multiconfigurational Calculations", *J. Comput. Chem.*, **2019**, *40*, 2216-2226.

Excited States
..............

One issue for excited states is that each state comes with different requirements on an active
space. AutoCAS handles active spaces by selecting an active space for each root and recombining the
active space later, to have one master active space, which is suited for each root.

Combination of Both
...................

AutoCAS can combine the large active space protocol with the excited state active space search, in
order to handle even excited state calculations for valence spaces with more than 100 orbitals.
