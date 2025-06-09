.. _forcefield:

forcefield section
==================

The ``forcefield`` section selects the static and dynamic force field models used in atomistic simulations.
This section is **mandatory** when running atomistic models (i.e. not BEM).

Valid keywords
--------------

+-------------+--------+----------------------------------------+---------+
| Keyword     | Type   | Allowed Values                         | Default |
+=============+========+========================================+=========+
| ``static``  | string | ``fq``, ``fqfmu``                      | --      |
+-------------+--------+----------------------------------------+---------+
| ``dynamic`` | string | ``wfq``, ``wfqfmu``                    | --      |
+-------------+--------+----------------------------------------+---------+
| ``kernel``  | string | ``ohno``, ``gaussian``, ``coulomb``    | --      |
+-------------+--------+----------------------------------------+---------+

Notes
-----

In the ``forcefield`` section, the default values are flexible and depends on the user selected keywords.

- ``static`` selects the model for the static force field:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``fq`` : fluctuating charges force field [default if ``dynamic: wfq``]
   - ``fqfmu`` : fluctuating charges and fluctuating dipoles force field [default if ``dynamic: wfqfmu``]
- ``dynamic`` selects the model for the frequency-dependent response:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``wfq``: frequency-dependent fluctuating charges [currently available for sodium and carbon-based materials]
   - ``wfqfmu``: frequency-dependent fluctuating charges and fluctuating dipoles [currently available for silver and gold]
- ``kernel`` defines the electrostatic interaction model between atomistic sites:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``coulomb``: pure Coulomb interaction [not recommended]
   - ``ohno``: Ohno kernel [default if ``static: fq``]
   - ``gaussian``: Gaussian kernel [default if ``static: fq`` or ``dynamic: wfq`` or ``dynamic: wfqfmu``]

Example
-------

.. code-block:: yaml

   forcefield:
     static: fqfmu
     dynamic: wfqfmu
     kernel: gaussian
