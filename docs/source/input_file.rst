Input File
==========

A plasmonX `input_file.yaml` supports the following sections.

.. code-block:: yaml

   what:
     # defines the type of calculation [energy, static/dynamic response, restart].

   forcefield:
     # defines the static forcefield for atomistic simulations

   field:
     # defines the field parameters [static/dynamic, frequencies,...]

   algorithm:
     # defines the numerical algorithm to be exploited to solve response equations.

   output:
     # defines the parameters of the output printing

   bem:
     # defines the parameters for implicit BEM calculations

   control:
     # defines the parameters to control geometry/creation of files

   atom_types:
     # defines the parameters of the atomtypes for atomistic simulations

   parameters:
     # defines the parameters of the ωFQ and ωFQFμ models

   input_geometry:
     # defines the input geometry [in Angstrom]

-------------------------------------------------------------------------------

**Section Overview**

The following sections describe the full list of available keywords, expected types, and allowed values:

- :doc:`input_sections/what`
- :doc:`input_sections/algorithm`
- :doc:`input_sections/control`
- :doc:`input_sections/output`
- :doc:`input_sections/field`
- :doc:`input_sections/forcefield`
- :doc:`input_sections/atom_types`
- :doc:`input_sections/parameters`
- :doc:`input_sections/bem`
- :doc:`input_sections/input_geometry`

-------------------------------------------------------------------------------

**Minimal Input**


**plasmonX** provides the possibility to use default values to give the most user friendly experience. 

Here are some minimal input examples to get started quickly with **plasmonX**.
Click on any of the following to jump directly to the corresponding configuration:

- :ref:`Minimal Example -- Sodium <example_wfq_sodium>` (:math:`\omega\mathrm{FQ}`)
- :ref:`Minimal Example -- Graphene <example_wfq_graphene>` (:math:`\omega\mathrm{FQ}`)
- :ref:`Minimal Example -- Silver <example_wfqfmu_silver>` (:math:`\omega\mathrm{FQF}\mu`)
- :ref:`Minimal BEM Example <example_bem>`

----

.. _example_wfq_sodium:

Minimal :math:`\omega\text{FQ}` -- Sodium

.. code-block:: yaml

   what:
     - dynamic response

   forcefield:
     dynamic: wFQ

   field:
     min freq: 0.01
     max freq: 3.00
     step freq: 0.01

   input_geometry:
     shape: sphere
     atomtype: na
     radius: 10.0

----


.. _example_wfq_graphene:

Minimal :math:`\omega\text{FQ}` -- Graphene

.. code-block:: yaml

   what:
     - dynamic response

   forcefield:
     dynamic: wFQ

   field:
     nfreq: 10
     min freq: 0.1
     max freq: 1.0

   input_geometry:
     shape: disk
     atomtype: C
     radius: 20.0

----

.. _example_wfqfmu_silver:

Minimal :math:`\omega\text{FQF}\mu` -- Silver

.. code-block:: yaml

   what:
     - dynamic response

   forcefield:
     dynamic: wFQFMu

   field:
     external freq: 3.0

   input_geometry:
     shape: decahedron
     atomtype: Ag
     radius: 14.0

----

.. _example_bem:

Minimal BEM Example

.. code-block:: yaml

   what:
     - dynamic response

   field:
     nfreq: 10
     min freq: 3.00
     max freq: 4.00

   bem:
     permittivity: silver jc

   input_geometry:
     type: implicit
     shape: rod
     main_axis: Z
     length: 100.0
     width: 30.0


.. toctree::
   :maxdepth: 2
   :hidden:

   input_sections/what
   input_sections/algorithm
   input_sections/control
   input_sections/output
   input_sections/field
   input_sections/forcefield
   input_sections/atom_types
   input_sections/parameters
   input_sections/bem
   input_sections/bem_solvents
   input_sections/input_geometry
   input_sections/geom_interface
