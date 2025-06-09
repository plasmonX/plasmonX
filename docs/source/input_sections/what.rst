.. _what:

what section
============

The ``what`` section specifies the type of calculation that **plasmonX** should perform.

This section is **optional**. Default: ``energy``

Valid keywords
--------------

The following keywords are supported:

+------------------------+-----------------------------------------------------------------------------------+
| Keyword                | Description                                                                       |
+========================+===================================================================================+
| ``- energy``           | Performs a single energy evaluation                                               |
+------------------------+-----------------------------------------------------------------------------------+
| ``- static response``  | Computes the linear response to static field                                      |
+------------------------+-----------------------------------------------------------------------------------+
| ``- dynamic response`` | Computes the linear response to dynamic field                                     |
+------------------------+-----------------------------------------------------------------------------------+
| ``- restart``          | Restarts a previous calculation. Only allowed combined with ``dynamic response``. |
+------------------------+-----------------------------------------------------------------------------------+

Examples
--------

.. code-block:: yaml

   what:
     - energy


.. code-block:: yaml

   what:
     - dynamic response

.. code-block:: yaml

   what:
     - dynamic response
     - restart
