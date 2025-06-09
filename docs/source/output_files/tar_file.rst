.. _tar_file:

Tar File
========

The `.tar.gz` file contains the following files:

- :ref:`.info <targz_info>`
- :ref:`.freq <freq_info>`
- :ref:`.csv <csv_info>`

----

.. _targz_info:

**.info File**

The `.info` file collects all relevant information about the **plasmonX** execution.

In particular, it contains:

- The number of variables
- Type of field and forcefield
- Input geometry
- Atomtype parameters
- Frequencies 

----

.. _freq_info:

**.freq File**

The `.freq` file contains the variable (charges/dipoles) values at the specific frequencies.

This is only present for ``what: - dynamic response`` calculations.

----

.. _csv_info:

**.csv File**

The `.csv` file contains the frequency dependent results (if ``what: dynamic response``)

See :doc:`csv_file`
