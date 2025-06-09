.. _log_file:

Log File
========

The `.log` file collects all relevant information about the simulation in a structured format which is characterized by the following sections. 
Click on any of the following to jump directly to the corresponding section:

- :ref:`System and Build Information <log_system_info>`
- :ref:`Parsed Input Summary <log_parse_input>`
- :ref:`Input Geometry <log_input_geom>`
- :ref:`Calculation Cycles <log_calc_cycle>`
- :ref:`Backup for Restart <log_backup_restart>`
- :ref:`Results <log_results>`
- :ref:`Maxima Analysis <log_max_analysis>`
- :ref:`Memory Usage <log_mem_usage>`
- :ref:`Final Messages & Exit Status <log_final_message>`

----

.. _log_system_info:

System and Build Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Includes metadata on the environment and code configuration:
   
   - Date and time of the run
   - Git branch used for compilation
   - Compiler version (e.g., GNU 9.4.0)
   - Math libraries used (e.g., MKL, OpenBLAS)
   - Whether integers are 64-bit
   - OpenMP usage
   - Total memory available

----

.. _log_parse_input:

Parsed Input Summary
^^^^^^^^^^^^^^^^^^^^

   Echoes the key inputs used during the run:
   
   - Input/output file names
   - Memory and OMP thread allocation
   - Type of calculation (`what`)
   - Field parameters: intensity, polarization, frequency range
   - Algorithm and forcefield settings
   - Principal axes / Restart / Verbose flags

----

.. _log_input_geom:

Input Geometry
^^^^^^^^^^^^^^

   For an atomistic calculation, the coordinates of all atoms, printed in tabular format under the header:`

   .. code-block:: console
   
      --------------------------------------------------------------------------------
                                       Input Geometry (Ang)
      --------------------------------------------------------------------------------

   If ``control: - principal axes`` is required, the rotated geometry is found under the header:

   .. code-block:: console
   
      --------------------------------------------------------------------------------
                                    Rotated Geometry (Ang)
      --------------------------------------------------------------------------------
      
   For a BEM calculation, the tesserae centroids are found under the header:

   .. code-block:: console
   
      --------------------------------------------------------------------------------
                         Input Tesserae Centroids (Ang)
      --------------------------------------------------------------------------------

----

.. _log_calc_cycle:

Calculation Cycles
^^^^^^^^^^^^^^^^^^

   For ``what: - dynamic response`` calculation, each macroiteration on the requested frequencies is printed with a header like:

   .. code-block:: console
   
      --------------------------------------------------------------------------------
      Cycle     1 out of    38
      --------------------------------------------------------------------------------

   This helps track progress of multi-frequency calculations.

----

.. _log_backup_restart:

Backup for restart
^^^^^^^^^^^^^^^^^^

   For ``what: - dynamic response`` calculation, during the execution of **plasmonX**, several plasmonx.bk files are created, which are used in case the calculation does not terminate and used for the restart. If the calculation, correctly ends, all .bk files are removed. This is indicated by the header:

   .. code-block:: console

      --------------------------------------------------------------------------------
      I am cleaning up the backup (*plasmonX.bk)
      --------------------------------------------------------------------------------

----

.. _log_results:

Results
^^^^^^^

   For ``what: -energy`` calculation, the energy is printed at the end of the calculation:

   .. code-block:: console

      --------------------------------------------------------------------------------
      Energy =    -0.06107855 a.u.
      --------------------------------------------------------------------------------

   For ``what: -static response`` calculation, the static polarizability is printed at the end of the calculation under the header:

   .. code-block:: console

      --------------------------------------------------------------------------------
                               polarizability tensor (a.u.)
      --------------------------------------------------------------------------------

   For ``what: -dynamic response`` calculation, the dynamic results are printed for each frequency. Example:

   .. code-block:: console

      --------------------------------------------------------------------------------
      Results for w = 0.220E-02 a.u.   0.207E+05 nm   0.600E-01 eV

      Isotr. Real Polar.       =         0.256586E+05 a.u.
      Isotr. Imag Polar.       =         0.988576E+01 a.u.
      Long. Real Polar. X      =         0.565633E+05 a.u.
      Long. Real Polar. Y      =         0.102063E+05 a.u.
      Long. Real Polar. Z      =         0.102063E+05 a.u.
      Long. Imag Polar. X      =         0.275545E+02 a.u.
      Long. Imag Polar. Y      =         0.105137E+01 a.u.
      Long. Imag Polar. Z      =         0.105137E+01 a.u.
      Iso. Abs. Cross. Sect.   =         0.199889E-02 a.u.
      Long Abs. Cross. Sect. X =         0.557149E-02 a.u.
      Long Abs. Cross. Sect. Y =         0.212585E-03 a.u.
      Long Abs. Cross. Sect. Z =         0.212585E-03 a.u.
      Iso. Sca. Cross. Sect.   =         0.369708E-09 a.u.
      Long Sca. Cross. Sect. X =         0.179664E-08 a.u.
      Long Sca. Cross. Sect. Y =         0.584964E-10 a.u.
      Long Sca. Cross. Sect. Z =         0.584964E-10 a.u.
      Iso. Ext. Cross. Sect.   =         0.199889E-02 a.u.
      Long Ext. Cross. Sect. X =         0.557149E-02 a.u.
      Long Ext. Cross. Sect. Y =         0.212585E-03 a.u.
      Long Ext. Cross. Sect. Z =         0.212585E-03 a.u.
      --------------------------------------------------------------------------------

   - The frequency value is printed in various units (a.u., nm, eV)
   - Real and imaginary components of polarizability (isotropic and along X, Y, Z)
   - Absorption, scattering, and extinction cross sections (isotropic and longitudinal components)

----

.. _log_max_analysis:

Maxima Analysis
^^^^^^^^^^^^^^^

   For ``what: -dynamic response`` calculation, the results for each frequency is followed by the maxima analysis section, which identifies peaks in the spectrum. Example:
   
   .. code-block:: console

      --------------------------------------------------------------------------------
      Maxima Analysis: NumExFreq =    300

               NState      Freq(eV)     Isotr. Abs. Cross. Sec.  (a.u.)
                 1          1.580              0.72424E+04
                 2          2.330              0.70713E+02

----

.. _log_mem_usage:

Memory Usage
^^^^^^^^^^^^

   At the end of the simulation, PlasmonX reports the peak memory usage during the run, e.g.:
   
   .. code-block:: console

      --------------------------------------------------------------------------------
      Peak memory used  :     13.768 MB       
      --------------------------------------------------------------------------------

----

.. _log_final_message:

Final Messages & Exit Status
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   A set of rock-inspired messages signals the simulation's outcome:
   
   - 🎸 **Ob-La-Di, Ob-La-Done!** – Normal Termination
   - 💀 **Ob-La-Di, Ob-La-Doom!** – Error Termination
   - 🚪 **Knock, knock, knockin’ on debug’s door** – Segmentation fault or memory issue
   
   These messages are printed together with the CPU and elapsed time, and a final statement of success or failure.

   Reference songs:

   - `Ob-La-Di, Ob-La-Da, The Beatles. <https://www.youtube.com/watch?v=9x5WY_jmsko>`_
   - `Knockin' On Heaven's Door, Bob Dylan. <https://www.youtube.com/watch?v=rm9coqlk8fY>`_
