.. _output:

output section
==============

The ``output`` section defines parameters related to the verbosity and analysis of the printed results.

This section is **optional**.

Valid keywords
--------------

+---------------------+--------+------------------------------------------------+-----------------+
| Keyword             | Type   | Allowed Values                                 | Default         |
+=====================+========+================================================+=================+
| ``verbose``         | int    | --                                             | ``0``           |
+---------------------+--------+------------------------------------------------+-----------------+
| ``maxima analysis`` | string | ``absorption``, ``scattering``, ``extinction`` | ``absorption``  |
+---------------------+--------+------------------------------------------------+-----------------+

Notes
-----

- ``verbose`` sets the verbosity level of the printed output. Be carefull when running in parallel.

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``0`` : minimal output
   - ``1`` : print minimal additional information
   - ``2`` : print most relevant quantities 
   - ``3`` : print all matrices/rhs [very large output]
   - ``4`` : print almost all computed quantities [very large output]

- ``maxima analysis`` find the maximum/maxima of selected properties. Allowed values:
   - ``absorption``
   - ``scattering``
   - ``extinction``

Example
-------

.. code-block:: yaml

   output:
     verbose: 1
     maxima analysis: absorption
