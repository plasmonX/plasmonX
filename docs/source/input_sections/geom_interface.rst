.. _geom_interface:

GEOM interface 
==============

plasmonX provides an internal engine for building complex nanostructures (both atomistic and implicit) through an interface to `GEOM <https://github.com/pgrobasillobre/geom/tree/branch-v1.0.0>`_

The geometry is defined by specifying shape-specific keywords in the ``input_geometry`` section.

GEOM supports two major geometry types:

- :ref:`Implicit geometries <implicit_geom>`
- :ref:`Atomistic geometries <atomistic_geom>`

In both cases, all the length parameters can be specified in:

- ``Ang`` or ``Angstrom`` : **Default**
- ``nm`` 
- ``a.u.`` or ``au``

-------------------------------------------------------------------------------

.. _implicit_geom:

Implicit Geometries
--------------------

Implicit (continuum) shapes are used with BEM models. These are defined via simple parameters and meshed for surface integration.

The implicit shapes are required by: ``type: implicit`` (mandatory)

The only supported shapes in implicit mode are:

+---------+-------------------------------+-------------------------------+
| Shape   | Required Parameters           | Description                   |
+=========+===============================+===============================+
| sphere  | ``radius``                    | Spherical nanoparticle        |
+---------+-------------------------------+-------------------------------+
| rod     | ``main_axis``, ``length``,    | Cylindrical nanorod aligned   |
|         | ``width``                     | along x/y/z                   |
+---------+-------------------------------+-------------------------------+

Notes
-----

Supported shapes for ``type: implicit`` geometries are:

- ``sphere`` : the allowed keywords are:

  - ``radius``: outer radius of the nanoparticle.
  - ``mesh_size`` *(optional)*: resolution of the mesh used to discretize the surface (default: ``10.0`` Angstrom).

- ``rod`` : the allowed keywords are:

  - ``main_axis``: orientation of the rod, one of ``x``, ``y``, or ``z``.
  - ``length``: total length of the rod.
  - ``width``: base diameter of the rod.
  - ``mesh_size`` *(optional)*: resolution of the mesh used to discretize the surface (default: ``10.0`` Angstrom).

Example
-------

.. code-block:: yaml

   input_geometry:
     type: implicit
     shape: sphere
     radius: 15.0
     mesh_size: 12.5


-------------------------------------------------------------------------------

.. _atomistic_geom:

Atomistic Geometries
---------------------

Atomistic geometries can be created by specifying ``type: atomistic``, which is optional (default). 

GEOM supports many shapes with both **single-material** and **core-shell** configurations. Dimers can also be generated.

For single-material configurations, the keyword ``atomtype`` is mandatory.

For core-shell configurations (available only for ``shape: sphere`` or ``shape: rod``) the keywords ``core_atomtype`` and ``shell_atomtype`` are mandatory.

+-------------+----------------------------------+--------------------------------------------------------------+
| Shape       | Required Parameters              | Description                                                  |
+=============+==================================+==============================================================+
| sphere      | ``radius``                       | Spherical nanoparticle                                       |
+-------------+----------------------------------+--------------------------------------------------------------+
| rod         | ``main_axis``, ``length``,       | Cylindrical nanorod                                          |
|             | ``width``                        |                                                              |
+-------------+----------------------------------+--------------------------------------------------------------+
| pyramid     | ``base_length``, ``height``      | Tip-like pyramids with square base                           |
+-------------+----------------------------------+--------------------------------------------------------------+
| icosahedron | ``radius``                       | Icosahedral nanoparticle (``ico`` or ``ih`` are allowed)     |
+-------------+----------------------------------+--------------------------------------------------------------+
| cone        | ``radius``, ``height``           | Cone                                                         |
+-------------+----------------------------------+--------------------------------------------------------------+
| tip         | ``a``, ``b``, ``height``         | Paraboloid tip                                               |
+-------------+----------------------------------+--------------------------------------------------------------+
| microscope  | ``a``, ``b``, ``height``,        | This is constructed by using a tip + base + pyramid,         |
|             | ``base_length``, ``height_tip``  | and resembles a TERS configuration                           |
+-------------+----------------------------------+--------------------------------------------------------------+
| disk        | ``radius``                       | Graphene-based disk                                          |
+-------------+----------------------------------+--------------------------------------------------------------+
| ring        | ``radius``, ``radius_in``,       | Graphene-based ring                                          |
+-------------+----------------------------------+--------------------------------------------------------------+
| triangle    | ``side_length``, ``edge_type``   | Graphene-based triangle armchair or zigzag                   |
+-------------+----------------------------------+--------------------------------------------------------------+
| ribbon      | ``length``, ``width``,           | Graphene nanoribbon armchair or zigzag                       |
|             | ``edge_type``                    |                                                              |
+-------------+----------------------------------+--------------------------------------------------------------+

Notes
-----

- ``sphere`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``radius``: radius of the sphere 
  - ``atomtype``: atomtype
  - For core-shell spheres:
   
    .. raw:: html
   
       <div style="margin: 0; padding: 0;"></div>
   
    - ``core_atomtype`` : core atomtype
    - ``shell_atomtype`` : shell atomtype
    - ``core_radius`` : radius of the core
    - ``shell_radius`` : radius of the shell 

- ``rod`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``main_axis``: the main axis of the rod (``x``, ``y``, ``z``)
  - ``length``: rod length
  - ``width``: cross-sectional width of the rod
  - ``atomtype``: atomtype
  - For core-shell rods:
   
    .. raw:: html
   
       <div style="margin: 0; padding: 0;"></div>
   
    - ``core_atomtype`` : core atomtype
    - ``shell_atomtype`` : shell atomtype
    - ``core_length`` : core length
    - ``shell_length`` : shell length
    - ``core_width`` : core width
    - ``shell_width`` : shell width

- ``pyramid`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``base_length`` : length of the square base
  - ``height`` : total height of the pyramid
  - ``atomtype`` : atomtype

- ``icosahedron`` (or ``ico``, ``ih``) : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``radius``: radius
  - ``atomtype``: atomtype

- ``cone`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``radius``: base radius
  - ``height``: height of the cone
  - ``atomtype``: atomtype

- ``tip`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``a``, ``b``: shape parameters controlling the paraboloid curvature :math:`z = x^2 / a^2 + y^2 / b^2`
  - ``height``: height of the tip
  - ``atomtype``: atomtype

- ``microscope`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``a``, ``b``: shape parameters controlling the paraboloid curvature :math:`z = x^2 / a^2 + y^2 / b^2`
  - ``height``: paraboloid height
  - ``base_length``: length of square base
  - ``height_tip``: height of the small pyramid at the top of the paraboloid
  - ``atomtype``: atomtype

- ``disk`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``radius``: disk radius
  - ``atomtype``: must be ``C`` (graphene)

- ``ring`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``radius``: outer radius
  - ``radius_in`` or ``shell_radius``: internal radius or shell thickness
  - ``atomtype``: must be ``C`` (graphene)

- ``triangle`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``side_length``: triangle side length
  - ``edge_type``: edge type. Allowed keywords are: ``armchair`` or ``zigzag``
  - ``atomtype``: must be ``C`` (graphene)

- ``ribbon`` : the allowed keywords are:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

  - ``length``: ribbon length
  - ``width``: ribbon width
  - ``edge_type``: edge type. Allowed keywords are: ``armchair`` or ``zigzag``.
  - ``atomtype``: must be ``C`` (graphene)

Additional options:

- ``dimer`` : if provided, creates a dimer of two structures. Allowed values are : 

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``default`` : the dimer is constructed as face-tip-face-tip configuration [optional if ``dimer_axis`` and ``distance`` are provided]
   - ``bowtie`` : the dimer is constructeda as face-tip-tip-face configuration
- ``distance`` : minimum distance between the two nanostructures
- ``dimer_axis`` : specifying the direction in which the two objects are aligned:

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``x``, ``y``, ``z`` : axis

Example
-------

.. code-block:: yaml

   input_geometry:
     shape: rod
     atomtype: Au
     main_axis: z
     length: 20.0 nm
     width: 8.0 nm
     distance: 5.0
     dimer_axis: z

