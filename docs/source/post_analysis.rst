.. _post_analysis:

.. raw:: latex

   \newpage

Post-Processing
===============

**plasmonX** provides a command-line tool called ``plasmonX_analysis`` for post-processing output data generated from simulations.

This tool extracts and visualizes information from the compressed `.tar.gz` output file.

Basic usage:

.. code-block:: bash

   plasmonX_analysis.py -i input_file.tar.gz -w WHAT [options]

After execution, a `post_process_input_file` folder is created.

Required arguments
------------------

- ``-i`` or ``--input_file``: path to the `.tar.gz` archive produced by plasmonX.
- ``-w`` or ``--what``: type of analysis to perform. Valid values:

  - ``xyz``: export atomic/tesserae positions in XYZ format
  - ``pqr``: export atomic charges for a given frequency (available for :math:`\omega\text{FQ}` and BEM)
  - ``field``: export induced electric field (3D or 2D plots)
  - ``density``: export induced charge density (available :math:`\omega\text{FQ}` and :math:`\omega\text{FQ}F\mu`)

Optional arguments
------------------

+---------------------------+---------------------------------------+------------------------------+------------------+
| Option                    | Description                           | Allowed Values               | Default          |
+===========================+=======================================+==============================+==================+
| ``-n``, ``--num_ex_freq`` | Number of frequencies to process      | Integer                      | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-freq``                 | Explicit list of frequencies (in eV)  | Comma-separated floats       | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-min_freq``             | Start of frequency range (eV)         | Float                        | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-max_freq``             | End of frequency range (eV)           | Float                        | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-step``                 | Frequency step (eV)                   | Float                        | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-plane``                | Plane to slice data (2D plot)         | ``xy``, ``xz``, ``yz``       | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-n_plane``              | Number of planes                      | Integer                      | ``10``           |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-step_plane``           | Distance between planes (Ang)         | Float                        | ``1.0``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-start``                | Start coordinate for slicing (Ang)    | Float                        | ``0.0``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-field_dir``            | Field direction for plotting          | ``x``, ``y``, ``z``, ``all`` | ``all``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-volume``               | Enable volume integration             | Flag                         | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-nx_points``            | Number of grid points in X direction  | Integer                      | ``200``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-ny_points``            | Number of grid points in Y direction  | Integer                      | ``200``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-nz_points``            | Number of grid points in Z direction  | Integer                      | ``200``          |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-min_grid``             | Min corner of 2D/3D grid (x,y,z)      | String (comma-separated)     | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-max_grid``             | Max corner of 2D/3D grid (x,y,z)      | String (comma-separated)     | --               |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-offset_grid``          | Padding around structure (Ang)        | Float                        | ``10.0``         |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-format_grid``          | Format of output files                | ``cube``, ``plt``            | ``cube``         |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-omp``                  | Number of OpenMP threads              | Integer                      | Max available    |
+---------------------------+---------------------------------------+------------------------------+------------------+
| ``-mem``                  | Available memory (GB)                 | Float                        | Max available    |
+---------------------------+---------------------------------------+------------------------------+------------------+

**Field**

- ``-n`` or ``--num_ex_freq``: Number of external frequencies to analyze. Required when extracting field or density data
- ``-freq`` or ``--frequencies``: Manually specify a comma-separated list of frequencies (in eV)
- ``-min_freq``, ``-max_freq``, ``-step``: Define a uniform frequency grid. These are used as an alternative to ``-freq``
- ``-field_dir``: Direction of the external electric field inducing the density or the field

**2D plot**

- ``-plane``: Select the slicing plane for visualization of induced density/field. When requested a folder `planes/` will be created containing `.csv` files, `png` plots, and `py` python scripts to produce the plots
- ``-n_plane``: Number of parallel planes to generate
- ``-step_plane``: Distance between planes in Angstorm
- ``-start``: Starting coordinate for slicing the system

**Grid**

- ``-nx_points`` : number of points in X direction
- ``-ny_points`` : number of points in Y direction
- ``-nz_points`` : number of points in Y direction
- ``-offset_grid`` : Extends the grid around the nanostructured system (default: ``±10.0`` Angstrom in all directions)
- ``-format_grid`` : Defines the file format for output grid data
- ``-min_grid``, ``-max_grid``: User-defined grid by setting the minimum and maximum coordinates in all directions. For 2D plots, define the min. and max. coordinates for the 2 required axis in ``-plane``

**Effective Volume & Area**

- ``-volume``: Enables the calculation of effective volume and area in a thin slab of volume V and thicknes h to analyze the induced field localization. To be combined with a user-defined grid, possibly with a large number of points. 

   The code performs the following integrals according to `Nanolett 2015, 15, 3410 <https://doi.org/10.1021/acs.nanolett.5b00759>`_ :

   .. math::
   
      V^{\text{eff}} = \int_V \frac{|\mathbf{E}_{\text{ind}}(x, y, z)|^2}{|\mathbf{E}^{\text{max}}_{\text{ind}}|^2} \, dV \\
      A^{\text{eff}} = \frac{1}{h} \int_V \frac{|\mathbf{E}_{\text{ind}}(x, y, z)|^2}{|\mathbf{E}^{\text{max}}_{\text{ind}}|^2} \, dV

**Resource management:**

- ``--omp``: Number of OpenMP threads to use for parallelization. If not specified, uses all available threads.
- ``--mem``: Amount of memory (in GB) available for the analysis. If not set, uses all system-available memory.

Examples 
--------

- XYZ extraction

   .. code-block:: bash
   
      plasmonX_analysis.py -i results.tar.gz -w xyz

   This will create a `results.xyz` file

- Field plots on XY plane

   .. code-block:: bash
   
      plasmonX_analysis.py -i results.tar.gz -w field -n 1 -freq 1.58 --plane xy --start 0.0 --step_plane 0.5 -field_dir x

   This will create a folder `planes/xy` containing 10 `.csv`, `.png`, `.py` for the selected planes (from 0.0 to 4.5 Angstrom)

- Extract charge densities for specific frequencies

   .. code-block:: bash
   
      plasmonX_analysis.py -i results.tar.gz -w density -n 3 -freq 1.2,2.4,3.6 -nx_points 100 -ny_points 100 -nz_points
   
   This will produce six `.cube` files with real and imaginary charge distributions at the selected frequencies.

- Field plots on XY plane

   .. code-block:: bash
   
      plasmonX_analysis.py -i results.tar.gz -w field -n 1 -freq 1.58 -volume -min_grid=-5.0,-5.0,-1.0 -max_grid=5.0,5.0,1.0

   This will produce a `.cube` file for the field calculated in the selected volume, and will calculate the effective volume and area (considering h = 2.0 Ang in the Z direction)
