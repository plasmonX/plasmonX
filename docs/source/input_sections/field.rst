.. _field:

field section
=============

The ``field`` section defines the parameters of the applied external electric field.

The section is **mandatory** for any simulation involving a response calculation (static or dynamic).

Valid keywords
--------------

+---------------------+-------------------+----------------------------------------+-----------------+
| Keyword             | Type              | Allowed Values                         |  Default        |
+=====================+===================+========================================+=================+
| ``type``            | string            | ``static``, ``dynamic``                |  --             |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``rhs type``        | string            | ``field``, ``potential``               |  ``field``      |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``field intensity`` | float             | –                                      | ``1.0e-4 a.u.`` |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``nfreq``           | integer           | –                                      |  --             |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``min freq``        | float or string   | –                                      |  units: eV      |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``max freq``        | float or string   | –                                      |  units: eV      |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``step freq``       | float or string   | –                                      |  units: eV      |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``external freq``   | float or string   | –                                      |  units: eV      |
+---------------------+-------------------+----------------------------------------+-----------------+
| ``polarization``    | string            | ``none``, ``x``, ``y``, ``z``, ``all`` |  ``all``        |
+---------------------+-------------------+----------------------------------------+-----------------+

Some additional notes on the keywords are given below.

Notes
-----

- ``type`` selects the nature of the field:

   .. raw:: html
 
      <div style="margin: 0; padding: 0;"></div>

   - ``static`` : static field.
   - ``dynamic`` : frequency-dependent field (linear absorption spectra, etc.).
- ``rhs type`` determines the formalism exploited for the right-hand side source term of the response equation (they are equivalent):
 
   .. raw:: html
 
      <div style="margin: 0; padding: 0;"></div>

   - ``field`` : electric field vector.
   - ``potential`` : scalar potential.
- ``field intensity`` sets the amplitude of the external field in atomic units.
- ``nfreq`` is number of external frequencies [optional]
- ``min freq``, ``max freq``, ``step freq`` define a frequency scan. The user must provide at least 3 between ``nfreq`` and these keywords.
- ``external freq`` can be used to specify specific custom frequencies instead of a range.
- All frequency values can be provided as either a float (e.g. `3.2` default units: eV) or a string with unit (e.g. `3.2 eV`). Allowed units are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``eV``
   - ``nm`` [in this case ``min freq`` and ``max freq`` are exchanged by plasmonX]
   - ``micron`` or ``um`` [in this case ``min freq`` and ``max freq`` are exchanged by plasmonX]
   - ``cm-1`` or ``cm^-1``
   - ``thz``
   - ``au`` or ``a.u.``
- ``polarization`` defines the direction of the field vector:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``x``, ``y``, ``z`` : linear polarization along specific axis
   - ``all`` : linear polarization along all axis
   - ``none`` : disables polarization (``energy`` calculation)

Example
-------

.. code-block:: yaml

   field:
     type: dynamic
     rhs type: field
     field intensity: 1.0
     min freq: 1.0
     max freq: 5.0
     step freq: 0.05
     polarization: all
