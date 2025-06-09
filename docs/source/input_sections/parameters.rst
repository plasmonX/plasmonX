.. _parameters:

parameters section
==================

The ``parameters`` section defines additional physical and numerical parameters for atomistic models for nanoplasmonics (:math:`omegatext{FQ}` and :math:`omegatext{FQF}\mu` ).

This section is **optional**, and must be omitted if you want to rely entirely on default parameters [currently defined for: sodium, carbon, silver, gold].  

Each atom type is given as a key, and parameters are defined in nested blocks (see Example)

Valid keywords
--------------

+-----------------------------+---------------------+---------------------------------------------+-------------+
| Keyword                     | Type                | Allowed Values                              | Default     |
+=============================+=====================+=============================================+=============+
| ``tau``                     | float or string     | > 0.0                                       | units: a.u. |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``sigma0``                  | float or string     | > 0.0                                       | units: S/m  |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``w_p``                     | float or string     | > 0.0                                       | units: eV   |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``gamma``                   | float or string     | > 0.0                                       | units: a.u. |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``scaling sigma0-tau``      | float or string     | > 0.0                                       | --          |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``a_ij``                    | float or string     | > 0.0                                       | units: a.u. |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``fermi_function.d``        | float or string     | > 0.0                                       | --          |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``fermi_function.s``        | float or string     | > 0.0                                       | --          |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``fermi energy``            | float or string     | > 0.0                                       | units: eV   |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``wfqfmu file``             | string              | valid path to .csv                          | --          |
+-----------------------------+---------------------+---------------------------------------------+-------------+
| ``permittivity``            | string              | ``silver etchegoin``, ``gold etchegoin``    | --          |
+-----------------------------+---------------------+---------------------------------------------+-------------+


Notes
-----

- Each atom type used in the ``input_geometry`` must have a corresponding entry here, unless default values are exploited.
- Each atom type defined in ``input_geometry`` can define its own ``parameters`` block.
- ``tau`` : Electron relaxation/scattering time (Drude model). ``tau`` can be expressed in:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``a.u.`` or ``au``
   - ``s``
   - ``fs``
- ``sigma0`` : Static electrical conductivity (Drude model). ``sigma0`` can be expressed in:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``s/m``
   - ``a.u.`` or ``au``
- ``w_p`` : Plasma frequency (Drude model). ``w_p`` can be expressed in:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``eV``
   - ``nm`` 
   - ``micron`` or ``um`` 
   - ``cm-1`` or ``cm^-1``
   - ``thz``
   - ``a.u.`` or ``au``
- ``gamma`` : Drude damping constant. ``gamma`` can be expressed in:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``au`` or ``a.u.``
   - ``fs^-1`` or ``fs-1``
   - ``s^-1`` or ``s-1``

- To define the Drude part of the response, any combination between (``sigma0``, ``w_p``) and (``tau``, ``gamma``) must be provided.
- ``scaling sigma0-tau`` : Scaling factor applied to both sigma0 and tau, used to adjust the conductivity model for specific materials (e.g. sodium for which fefault: ``10``)
- ``a_ij``: effective area between the atoms. ``a_ij`` can be expressed in:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``au`` or ``a.u.``
   - ``nm^2`` or ``nm2``
   - ``ang^2`` or ``ang2`` or ``angstrom^2`` or ``angstrom2``

- ``fermi function``: defines the Fermi-like function parameters for tunneling regime. 

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``fermi_function.d``: damping width
   - ``fermi_function.s``: damping slope

- ``fermi energy`` : Fermi energy. Used only for graphene-based materials
- ``wfqfmu file`` : points to a CSV file with frequency-dependent parameters. The file must contains 3 columns: frequencies (eV) Re(alpha) Im(alpha)
- ``permittivity`` refers to a predefined dielectric function, from which the interband polarizability is computed.

- When more than one atom type is defined in the ``input_geometry``, it is possible to define interaction-specific damping functions using the ``interaction AtomA->AtomB`` syntax (see Example).

   This subsection is **optional**, and must be omitted if you want to rely entirely on default parameters [currently defined for: silver-gold bimetallic systems].  

   +-----------------------------+---------------------+-------------------------------+
   | Keyword                     | Type                | Description                   |
   +=============================+=====================+===============================+
   | ``fermi_function.d``        | float or string     | Width of the damping function | 
   +-----------------------------+---------------------+-------------------------------+
   | ``fermi_function.s``        | float or string     | Slope of the damping function | 
   +-----------------------------+---------------------+-------------------------------+

Example
-------

.. code-block:: yaml

   parameters:
     Au:
       tau: 318.018
       sigma0: 11834849.81
       A_ij: 9.61
       fermi_function:
         d: 12.0
         s: 1.1
       wfqfmu file: jc_gold_polar.csv
     Ag:
       tau: 1633.608
       sigma0: 65313000.0
       A_ij: 9.61
       fermi_function:
         d: 12.0
         s: 1.1
       wfqfmu file: jc_silver_polar.csv
     interaction Au->Ag:
       fermi_function:
         d: 12.0
         s: 0.82555
     interaction Ag->Au:
       fermi_function:
         d: 12.0
         s: 0.82555
