.. _control:

control section
===============

The ``control`` section provides optional parameters that influence the generation of internal files or the handling of geometries.

This section is **optional**.

Valid keywords
--------------

+------------------------+-----------------------------------------------------------------------------------+
| Keyword                | Description                                                                       |
+========================+===================================================================================+
| ``- no info file``     | disable the creation of info files used for post-processing                       |
+------------------------+-----------------------------------------------------------------------------------+
| ``- principal axes``   | forces alignment of the input geometry along its principal axes of inertia.       |
+------------------------+-----------------------------------------------------------------------------------+

Example
--------

.. code-block:: yaml

   control:
     - no info file
     - principal axes

