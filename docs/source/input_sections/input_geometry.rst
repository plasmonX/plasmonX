.. _input_geometry:

input_geometry section
======================

The ``input_geometry`` section defines the atomic structure used in the simulation.

The geometry is assumed in all cases to be in Angstrom. 

There are **three supported modes** for declaring geometries:

1. :ref:`Direct YAML string input <yaml_geometry>`
2. :ref:`External .xyz file <xyz_geometry>`
3. :ref:`Generation with GEOM interface <geom_geometry>`

The ``input_geometry`` section is **mandatory** for atomistic simulations. For BEM, it can be used to construct new geometries (see :ref:`Generation with GEOM interface <geom_geometry>`)

-------------------------------------------------------------------------------

.. _yaml_geometry:

1. YAML inline geometry
-----------------------

Geometry can be written directly in the YAML file as a multiline string, using the ``|`` character. 

*Important* : Remember the first blank space at the beginning of each new line for YAML format.

Example
-------

.. code-block:: yaml

   input_geometry: |
     Ag 0.0 0.0 0.0
     Ag 1.0 0.0 0.0
     Ag 0.0 1.0 0.0
     Ag 0.0 0.0 1.0

If you are running a FQ or FQFMu calculations, you can specify the molecule number using the `[IMol=X]` tag:

.. code-block:: yaml

   input_geometry: |
     O-OW[IMol= 1]      0.09543100      0.65652700      1.21555100
     H-HW[IMol= 1]      0.67989369      0.65652700      1.97965006
     H-HW[IMol= 1]      0.67989417      0.65652700      0.45145232
     O-OW[IMol= 2]     -0.09543100     -0.65652700     -1.21555100
     H-HW[IMol= 2]     -0.67989417     -0.65652700     -0.45145232
     H-HW[IMol= 2]     -0.67989369     -0.65652700     -1.97965006

- Atom types must be defined in the `atom_types` section or defaulted.
- IMol must be a positive integer.
- Coordinates must be three floating-point numbers.

-------------------------------------------------------------------------------

.. _xyz_geometry:

2. External .xyz file
---------------------

You can reference an external geometry file in XYZ format using the keyword ``external xyz file``.

*Important* : This is not recommended for FQ or FQFMu calculations, because all the atoms are assumed to constitute one molecule.

Example
-------

.. code-block:: yaml

   input_geometry:
     external xyz file: my_structure.xyz

- The file must be in standard XYZ format.
- Atom types must match those used in `atom_types`.

-------------------------------------------------------------------------------

.. _geom_geometry:

3. Generation with GEOM
-----------------------

You can use **plasmonX** interface to `GEOM <https://github.com/pgrobasillobre/geom/tree/branch-v1.0.0>`_ to construct complex nanostructures directly.

The allowed shapes and related keywords are given in: :doc:`geom_interface <geom_interface>`.
