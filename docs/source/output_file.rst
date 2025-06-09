.. _output_file:

.. raw:: latex

   \newpage

Output File(s)
==============

By default, **PlasmonX** produces two output files:

1. A `.log` file containing the main results.
2. A `.tar.gz` folder containing info file (if not disable in :doc:`input_sections/control`).
3. A `.csv` file containing frequency dependent results (if ``what: dynamic response``)

The following sections describe the structures of all output files:

- :doc:`output_files/log_file`
- :doc:`output_files/tar_file`
- :doc:`output_files/csv_file`

.. toctree::
   :maxdepth: 2
   :hidden:

   output_files/log_file
   output_files/tar_file
   output_files/csv_file
