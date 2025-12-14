# plasmonX test suite

This directory contains a collection of tests used to check the installation (ctest) and showcase the main features of **plasmonX**.

## Energy & Static Response

- **fq_energy**: FQ energy (Ohno Kernel)
- **fq_principal_axis_energy**: FQ energy and rotation of the structure in principal axis
- **fq_polarizability**: FQ static polarizability
- **fqfmu_energy**: FQFμ energy (Gaussian kernel)
- **fqfmu_rq_rmu_energy**: FQFμ energy (Gaussian kernel, varying Rq and RMu parameters)

## ωFQ (sodium nanoparticles)

- **wfq_sodium_inversion**: ωFQ sodium NP with explicit parameters (solved by inversion)
- **wfq_sodium_minimal_input**: ωFQ sodium NP with minimal input (solved by inversion)
- **wfq_sodium_restart**: ωFQ sodium NP requesting restart of the calculation
- **wfq_sodium_geom_sphere**: ωFQ sodium spherical NP constructed by using GEOM
- **wfq_sodium_geom_rod**: ωFQ sodium rod NP constructed by using GEOM
- **wfq_sodium_iterative_on_memory**: ωFQ sodium NP (solved by iterative algorithm on memory)
- **wfq_sodium_iterative_on_the_fly**: ωFQ sodium NP (solved by iterative algorithm on the fly)

## ωFQ (graphene-based nanostructures)

- **wfq_graphene_disk_inversion**: ωFQ graphene disk with explicit parameters(solved by inversion)
- **wfq_graphene_ribbon_inversion**: ωFQ graphene ribbon with explicit parameters (solved by inversion)
- **wfq_graphene_triangle_inversion**: ωFQ graphene triangle with explicit parameters (solved by inversion)
- **wfq_graphene_minimal_input**: ωFQ graphene disk with minimal input (solved by inversion)
- **wfq_graphene_geom_disk**: ωFQ graphene disk constructed by using GEOM
- **wfq_graphene_geom_ribbon_main_armchair**: ωFQ graphene ribbon with armchair main axis constructed by using GEOM
- **wfq_graphene_geom_ribbon_main_zigzag**: ωFQ graphene ribbon with zigzag main axis constructed by using GEOM
- **wfq_graphene_geom_ring**: ωFQ graphene ring constructed by using GEOM
- **wfq_graphene_geom_triangle_armchair**: ωFQ graphene triangle with armchair edges constructed by using GEOM
- **wfq_graphene_geom_triangle_zigzag**: ωFQ graphene disk with armchair edges constructed by using GEOM

## ωFQFμ (homogeneous nanoparticles)

- **wfqfmu_ag_inversion**: ωFQFμ Ag NP with explicit parameters (solved by inversion)
- **wfqfmu_ag_internal_permittivity**: ωFQFμ Ag NP reading internal permittivity (solved by inversion)
- **wfqfmu_ag_minimal_input**: ωFQFμ Ag NP with minimal input (solved by inversion)
- **wfqfmu_au_minimal_input**: ωFQFμ Au NP with minimal input (solved by inversion)
- **wfqfmu_ag_read_xyz**: ωFQFμ Ag NP reading external xyz file (solved by inversion)
- **wfqfmu_au_read_xyz**: ωFQFμ Au NP reading external xyz file (solved by inversion)
- **wfqfmu_ag_diff_units**: ωFQFμ Ag NP using different units (solved by inversion)
- **wfqfmu_ag_restart**: ωFQFμ Ag NP requesting restart of the calculation (solved by inversion)
- **wfqfmu_ag_geom_sphere**: ωFQFμ Ag spherical NP constructed by using GEOM(solved by inversion)
- **wfqfmu_ag_geom_sphere_units**: ωFQFμ Ag spherical NP constructed by using GEOM and different units (solved by inversion)
- **wfqfmu_ag_geom_decahedron**: ωFQFμ Ag decahedral NP constructed by using GEOM (solved by inversion)
- **wfqfmu_ag_geom_rod_x**: ωFQFμ Ag rod (main axis: X) NP constructed by using GEOM (solved by inversion)
- **wfqfmu_ag_geom_dimer_bowtie**: ωFQFμ Ag dimer in bowtie configuration constructed by using GEOM (solved by inversion)
- **wfqfmu_au_geom_icosahedron**: ωFQFμ Au icosahedral NP constructed by using GEOM (solved by inversion)
- **wfqfmu_au_geom_cuboctahedron**: ωFQFμ Au cuboctahedral NP constructed by using GEOM (solved by inversion)
- **wfqfmu_au_geom_paraboloid**: ωFQFμ Au paraboloid NP constructed by using GEOM (solved by inversion)
- **wfqfmu_au_geom_dimer_z**: ωFQFμ Au dimer along Z axis constructed by using GEOM (solved by inversion)
- **wfqfmu_ag_iterative_on_memory**: ωFQFμ Ag NP (solved by iterative algorithm on memory)
- **wfqfmu_ag_iterative_on_the_fly**: ωFQFμ Ag NP (solved by iterative algorithm on the fly)

## ωFQFμ (heterogeneous nanoparticles)

- **wfqfmu_hetero_agau_inversion**: ωFQFμ AgAu NP with explicit parameters (solved by inversion)
- **wfqfmu_hetero_agau_read_1param_inversion**: ωFQFμ AgAu NP reading 1 Ag-Au interaction parameter (solved by inversion)
- **wfqfmu_hetero_agau_default_int_param_inversion**: ωFQFμ AgAu NP using default interaction parameters (solved by inversion)
- **wfqfmu_hetero_agau_minimal_input**: ωFQFμ AgAu NP with minimal input (solved by inversion)
- **wfqfmu_hetero_agau_read_xyz**: ωFQFμ AgAu NP reading external xyz file (solved by inversion)
- **wfqfmu_hetero_agau_geom_sphere_core_shell**: ωFQFμ AgAu core-shell NP constructed by using GEOM (solved by inversion)
- **wfqfmu_hetero_agau_geom_rod_core_shell**: ωFQFμ AgAu core-shell rod NP constructed by using GEOM (solved by inversion)
- **wfqfmu_hetero_agau_iterative_on_memory**: ωFQFμ AgAu NP (solved by iterative algorithm on memory)
- **wfqfmu_hetero_agau_iterative_on_the_fly**: ωFQFμ AgAu NP (solved by iterative algorithm on the fly)

## BEM

- **bem_dpcm_accurate_inversion**: BEM spherical NP using DPCM and accurate integration (solved by inversion)
- **bem_ief_accurate_inversion**: BEM spherical NP using IEF-PCM and accurate integration (solved by inversion)
- **bem_ief_approximate_inversion**: BEM spherical NP using IEF-PCM and approximate integration (solved by inversion)
- **bem_dpcm_geom_rod**: BEM rof NP using DPCM, constructed by using GEOM (solved by inversion)
- **bem_dpcm_geom_sphere**: BEM spherical NP using DPCM, constructed by using GEOM (solved by inversion)
- **bem_dpcm_brendel_bormann_ag**: BEM spherical Ag NP using DPCM and Brendel-Bormann epsilon (solved by inversion)
- **bem_dpcm_brendel_bormann_au**: BEM spherical Au NP using DPCM and Brendel-Bormann epsilon (solved by inversion)
- **bem_dpcm_palik_ag**: BEM spherical Ag NP using DPCM and Palik epsilon (solved by inversion)
- **bem_dpcm_palik_au**: BEM spherical Au NP using DPCM and Palik epsilon (solved by inversion)
- **bem_dpcm_external_brendel_bormann_ag**: BEM spherical Ag NP using DPCM and reading Brendel-Bormann epsilon from file (solved by inversion)
- **bem_dpcm_accurate_iterative_on_memory**: BEM spherical NP using DPCM (solved by iterative algorithm on memory)

## Post-Process Analysis with **plasmonX_analysis**

- **analysis_xyz**: xyz file
- **analysis_xyz_bem**: xyz file for BEM
- **analysis_pqr_wfq_sodium**: pqr file for ωFQ calculation (sodium NP)
- **analysis_pqr_bem**: pqr file for BEM calculation
- **analysis_field_3d_x_wfq_sodium**: ωFQ induced field along X axis in CUBE format (sodium NP)
- **analysis_field_3d_x_wfq_sodium_plt**: ωFQ induced field along X axis in PLT format (sodium NP)
- **analysis_field_3d_x_y_z_wfqfmu_ag**: ωFQFμ induced field along all axes in CUBE format (Ag NP)
- **analysis_field_3d_x_y_z_wfqfmu_ag_plt**: ωFQFμ induced field along all axes in PLT format (Ag NP)
- **analysis_field_3d_y_bem**: BEM induced field along Y axis in CUBE format
- **analysis_density_3d_x_wfq_sodium**: ωFQ induced density along X axis in CUBE format (sodium NP)
- **analysis_density_3d_x_wfq_sodium_plt**: ωFQ induced density along X axis in PLT format (sodium NP)
- **analysis_density_3d_y_wfqfmu_au**: ωFQFμ induced field along all Y axis in CUBE format (Au NP)
- **analysis_density_3d_y_wfqfmu_au_plt**: ωFQFμ induced field along all Y axis in PLT format (Au NP)
- **analysis_density_2d_x_wfq_sodium**: ωFQ induced density along X axis in a specific 2D plane (xy, sodium NP)
- **analysis_density_2d_x_wfq_sodium_n_planes**: ωFQ induced density along X axis in various 2D planes (xy, sodium NP)
- **analysis_density_2d_y_xy_wfqfmu_ag**: ωFQFμ induced density along Y axis in a specific 2D plane (xy, Ag NP)
- **analysis_field_2d_x_yz_wfq_sodium**: ωFQ induced field along X axis in a specific 2D plane (yz, sodium NP)
- **analysis_field_2d_z_xz_wfqfmu_ag**: ωFQFμ induced field along Z axis in a specific 2D plane (xz, Ag NP)
- **analysis_field_2d_z_xz_wfqfmu_ag_reduced**: ωFQFμ induced field along Z axis in a specific reduced 2D plane (xz, Ag NP)
- **analysis_field_2d_y_xy_bem**: BEM induced field along Y axis in a specific 2D plane (xy)
- **analysis_field_volume_x_wfq_sodium**: ωFQ field localization volume and area along X axis (sodium NP)
- **analysis_field_volume_z_wfqfmu_au**: ωFQ field localization volume and area along Z axis (Au NP)
