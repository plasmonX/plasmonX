.. _bem:

bem section
===========

The ``bem`` section defines parameters for boundary element method (BEM) simulations based on implicit solvation models.
This section is **mandatory** when using the BEM approach.

Valid keywords
--------------

+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| Keyword                  | Type              | Allowed Values                                                         | Default            |
+==========================+===================+========================================================================+====================+
| ``mesh file``            | string            | Path to mesh in `.msh` format (generally obtained from GMSH)           | ``input_file.msh`` |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``green function``       | string            | ``accurate``, ``approximate``                                          | ``approximate``    |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``variant``              | string            | ``dpcm``, ``iefpcm``                                                   | ``dpcm``           |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``sphere radius``        | float or string   | --                                                                     | --                 |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``normal scalar factor`` | float             | ``1.0``, ``-1.0`` (sign of surface normals)                            | ``1.0``            |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``permittivity file``    | string            | Path to CSV table with frequency-dependent complex dielectric function | --                 |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``permittivity``         | string            | Predefined permittivity functions (see Notes)                          | --                 |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``solvent``              | string            | --                                                                     | ``vacuum``         |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+
| ``solvent epsilon``      | float             | --                                                                     | --                 |
+--------------------------+-------------------+------------------------------------------------------------------------+--------------------+

Notes
-----

- ``mesh file`` : 3D molecular surface in `.msh` format with vertex positions, normals, etc. This is generally obtained from `GMSH <https://gmsh.info/>`_ software. Assumed in Angstrom.
- ``green function``: specifies whether accurate or approximate Green's function integrals must be computed

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``accurate``: full numerical evaluation of Green’s function integrals.
   - ``approximate``: uses simplified expressions (faster, less accurate, centroid-based).
- ``variant`` : specifies the right-hand-side formalism.

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``dpcm``: Dielectric Polarizable Continuum Model (field, faster).
   - ``iefpcm``: Integral Equation Formalism Polarizable Continuum Model (potential, more stable numerically) .
- ``sphere radius`` : radius of enclosing cavity (required for ``green function: approximate``). If it is not given, it is calculated from ``mesh file``.
- ``normal scalar factor`` : Orientation of surface normals (``-1.0`` must only be used for compatibility with other softwares).
- ``permittivity`` : Built-in frequency-dependent permittivities . Available models are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``silver etchegoin`` or ``silver johnson-christy`` or ``silver jc`` : Silver permittivity taken from `Phys. Rev. B 1972, 6, 4370 <https://doi.org/10.1103/physrevb.6.4370>`_
   - ``silver brendel-bormann`` or ``silver bb`` : Silver permittivity taken from `Appl. Opt. 1998, 37, 5271 <https://doi.org/10.1364/ao.37.005271>`_
   - ``silver palik`` : Silver permittivity taken from `Handbook of Optical Constants of Solids (Elsevier, 1997) <https://www.sciencedirect.com/book/9780125444156/handbook-of-optical-constants-of-solids>`_
   - ``gold etchegoin`` or ``gold johnson-christy`` or ``gold jc`` : Gold permittivity taken from `Phys. Rev. B 1972, 6, 4370 <https://doi.org/10.1103/physrevb.6.4370>`_
   - ``gold brendel-bormann`` or ``gold bb`` : Gold permittivity taken from `Appl. Opt. 1998, 37, 5271 <https://doi.org/10.1364/ao.37.005271>`_
   - ``gold palik`` : Gold permittivity taken from `Handbook of Optical Constants of Solids (Elsevier, 1997) <https://www.sciencedirect.com/book/9780125444156/handbook-of-optical-constants-of-solids>`_
- ``permittivity file`` : points to a CSV file with frequency-dependent complex dielectric function of the material. The file must contains 3 columns: frequencies (eV) Re(Eps) Im(Eps)
- ``solvent`` : specifies the solvent embedding the nanostructure. Tabulated solvents are given in :doc:`bem_solvents <bem_solvents>`.
- ``solvent epsilon`` : specifies the solvent optical dielectric constant.

Example
-------

.. code-block:: yaml

   bem:
     mesh file: ag_nanocube.msh
     normal scalar factor: 1.0
     permittivity: silver etchegoin
     green function: approximate
     variant: iefpcm
     sphere radius: 5.0
     solvent epsilon: 1.77768

.. toctree::
   :maxdepth: 1
   :hidden:

   bem_solvents
