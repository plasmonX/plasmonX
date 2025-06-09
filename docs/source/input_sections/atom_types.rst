.. _atom_types:

atom_types section
==================

The ``atom_types`` section defines the physical parameters associated with each atom type used in atomistic simulations.

This section is **optional**, and must be omitted if you want to rely on built-in parameters. 

See example for how to structure the input file.

Valid keywords
--------------

+-----------+---------------------+---------------+-------------+
| Keyword   | Type                | Allowed Values| Default     |
+===========+=====================+===============+=============+
| ``chi``   | float or string     | >= 0.0        | units: a.u. |
+-----------+---------------------+---------------+-------------+
| ``eta``   | float or string     | > 0.0         | units: a.u. |
+-----------+---------------------+---------------+-------------+
| ``alpha`` | float or string     | > 0.0         | units: a.u. |
+-----------+---------------------+---------------+-------------+
| ``rq``    | float or string     | > 0.0         | units: a.u. |
+-----------+---------------------+---------------+-------------+
| ``rmu``   | float or string     | > 0.0         | units: a.u. |
+-----------+---------------------+---------------+-------------+

Notes
-----


- Each atom type used in the ``input_geometry`` must have a corresponding entry here, unless default values are exploited.
- Each atom type defined in ``input_geometry`` can define its own ``atom_type`` block.
- These parameters define fundamental physical constants for each atomtype:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``chi``: electronegativity [required for fq/fqfmu forcefield]
   - ``eta``: chemical hardness [required for fq/fqfmu forcefield]
   - ``alpha``: polarizability  [required for fqfmu forcefield]
   - ``rq``: charge damping radius for gaussian kernel  [optional for fq/fqfmu forcefield]
   - ``rmu``: dipole damping radius for gaussian kernel [optional for fq/fqfmu forcefield]
- All values can be provided either as **pure float** (e.g. `3.2`) or **with units** as strings (e.g. `"3.2 a.u."`)

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``chi`` and ``eta`` can be expressed in:

      .. raw:: html

         <div style="margin: 0; padding: 0;"></div>

      - ``a.u.`` or `au`
      - ``eV``
   - ``alpha`` can be expressed in: 

      .. raw:: html

         <div style="margin: 0; padding: 0;"></div>

      - ``a.u.`` or `au`
      - ``nm^3`` or ``nm3``
      - ``ang^3`` or ``ang3`` or ``angstrom3`` or ``angstrom^3``
   - ``rq`` and ``rmu`` can be expressed in: 

      .. raw:: html

         <div style="margin: 0; padding: 0;"></div>

      - ``a.u.`` or `au`
      - ``nm`` or ``nm``
      - ``ang`` or ``angstrom``

Example
-------

.. code-block:: yaml

   atom_types:
     Ag:
       chi: 0.0
       eta: 0.483
       alpha: 49.9843
       rq: 1.78756061721763 ang
       rmu: 1.25362081177874 ang
