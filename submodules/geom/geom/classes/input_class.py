from geom.functions import general, output

class input_class:
   """
   Manages user-defined parameters for geometric transformations and structure generation.

   This class stores user inputs related to translation, rotation, structure creation, 
   and other transformations for metal nanoparticles and graphene structures.

   Attributes:
       - Transformation Flags: Enable operations such as translation, rotation, 
         mirroring, and structure generation.
       - Geometry Files: Stores input filenames for molecular structures.
       - Translation & Rotation: Contains direction, distances, angles, and axis information.
       - Structure Generation: Flags and parameters for generating different nanostructures.
       - Miscellaneous: Handles alloy composition, verbosity settings, and file extensions.

   Notes:
       - The full list of attributes is initialized in `__init__()`.
       - Methods validate inputs, read data, and apply transformations.
   """

   def __init__(self):
      """
      Initializes input parameters for geometric transformations and structure generation.
  
      Attributes Initialized:
          - Flags for operations (translation, rotation, mirroring, geometry creation).
          - File handling (input geometry, output file names, temporary directories).
          - Geometric parameters (distances, angles, lattice properties).
          - Structure-specific settings (graphene edge type, core-shell configuration, alloy composition).
      """

      # -- Small tasks
      self.small_tasks = False
       
      # -- Translate
      self.translate = False 
      self.translate_controlled_distance = False
      self.translate_1 = False
     
      self.move_geom_to_000   = False
      self.move_geom_1_to_000 = False
      self.move_geom_2_to_000 = False

      self.distances  = []
      self.dir_factor = []

      self.file_geom2_translated = ''
      self.distances_input = ''
      self.geom_file       = ''
      self.geom1_file      = ''
      self.geom2_file      = ''
      self.dir_axis_input  = ''
      self.origin_CM       = ''
      self.origin_CM_1     = ''
      self.origin_CM_2     = '' 

      self.direction = 1.0
      self.dimer_distance = 0.0

      # -- Minimum distance
      self.min_dist = False
      
      # -- Rotate
      self.rotate = False
      self.rotate_angles = False
      self.rotate_1      = False

      self.angles = []

      self.angles_input = ''
      self.file_geom_rotated = ''

      self.angle = 0.0
      
      # -- Geometrical center
      self.geom_center = False

      # -- Specular geometry
      self.geom_specular = False

      # -- Generate structure geometry
      self.create_dimer = False
      self.create_bowtie = False
      self.create_geom = False
      self.create_ase_bulk = False
      self.gen_3d_mesh = False
      self.gen_3d_mesh_sphere = False
      self.gen_3d_mesh_rod = False
      self.gen_graphene = False
      self.gen_core_shell = False
      self.gen_sphere = False
      self.gen_sphere_core_shell = False
      self.gen_rod = False
      self.gen_rod_core_shell = False
      self.gen_tip = False
      self.gen_cone = False
      self.gen_icosahedra = False
      self.gen_cto = False
      self.gen_idh = False
      self.gen_pyramid = False
      self.gen_microscope = False
      self.alloy = False

      self.mesh_size = 1.0
      self.alloy_perc = 0.0
      self.elliptic_parabola_a = 1.0 # Modifies stepness along x
      self.elliptic_parabola_b = 1.0 # Modifies stepness along y
      self.elliptic_parabola_c = 0.0 # Fixed for xy parabolloid
      self.z_min = 0.0 
      self.z_max = 0.0 
      self.z_max_paraboloid = 0.0
      self.z_max_pyramid = 0.0 
      self.rod_length = 0.0
      self.rod_length_in = 0.0
      self.rod_length_out = 0.0
      self.rod_width = 0.0
      self.rod_width_in = 0.0
      self.rod_width_out = 0.0
      self.side_length = 0.0
      self.radius     = 0.0
      self.radius_in  = 0.0
      self.radius_out = 0.0
      self.X_length = 0.0
      self.Y_length = 0.0
      self.sphere_center = [0.0,0.0,0.0]

      self.xyz_output = ''
      self.mesh_output = ''
      self.tmp_folder = ''
      self.alloy_string = ''
      self.atomtype = ''
      self.atomtype_in = ''
      self.atomtype_out = ''
      self.atomtypes_alloys     = ['ag','au']
      self.atomtypes_core_shell = ['ag','au']
      self.axes = ["x","y","z"]
      self.graphene_structures  = ["rib","disk","ring","triangle"]
      self.graphene_structure   = ""
      self.graphene_edge_types  = ["armchair","zigzag"]
      self.graphene_edge_type   = ""

      # -- Merge geometries
      self.merge = False
      self.merge_cutoff = 0.0 

      # -- Verbose
      self.verbose = False
      self.verbose_inp = ''

   # -------------------------------- #
   # ------- Check input case ------- #
   
   def check_input_case(self):
      """
      Validates and checks input requirements for selected operations.

      Returns:
          None: Performs necessary checks and raises errors if conditions are not met.

      Notes:
          - Ensures input files exist before processing.
          - Checks file extensions for validity.
          - Validates direction axis input.
          - Ensures selected operations are correctly configured.
      """

      if (self.translate_controlled_distance):
         general.check_file_exists(self.distances_input)
         self.read_input(what='distances')

         general.check_file_exists(self.geom1_file)
         general.check_file_exists(self.geom2_file)

         general.check_file_extension(self.geom1_file,'.xyz')
         general.check_file_extension(self.geom2_file,'.xyz')

         general.check_dir_axis(self)

      elif (self.translate_1):
         general.check_file_exists(self.geom_file)
         general.check_file_extension(self.geom_file,'.xyz')

         general.check_dir_axis(self)

      elif (self.rotate_angles):
         general.check_file_exists(self.angles_input)
         self.read_input(what='angles')

         general.check_file_exists(self.geom_file)
         general.check_file_extension(self.geom_file,'.xyz')

         general.check_dir_axis(self)

      elif (self.rotate_1):
         general.check_file_exists(self.geom_file)
         general.check_file_extension(self.geom_file,'.xyz')

         general.check_dir_axis(self)

      elif (self.min_dist):
         general.check_file_exists(self.geom1_file)
         general.check_file_exists(self.geom2_file)

         general.check_file_extension(self.geom1_file,'.xyz')
         general.check_file_extension(self.geom2_file,'.xyz')

      elif (self.merge):
         general.check_file_exists(self.geom1_file)
         general.check_file_exists(self.geom2_file)

         general.check_file_extension(self.geom1_file,'.xyz')
         general.check_file_extension(self.geom2_file,'.xyz')

      elif (self.geom_center):
         general.check_file_exists(self.geom_file)
         general.check_file_extension(self.geom_file,'.xyz')

      elif (self.geom_specular):
         general.check_file_exists(self.geom_file)
         general.check_file_extension(self.geom_file,'.xyz')

      elif (self.create_geom):
         if (not self.gen_graphene          and
             not self.gen_sphere            and
             not self.gen_sphere_core_shell and
             not self.gen_rod               and
             not self.gen_rod_core_shell    and
             not self.gen_tip               and 
             not self.gen_pyramid           and
             not self.gen_cone              and
             not self.gen_microscope        and
             not self.gen_icosahedra        and
             not self.gen_cto               and
             not self.gen_idh               and 
             not self.gen_3d_mesh_sphere    and
             not self.gen_3d_mesh_rod): output.error("Create geom option not recognised.")

         if self.create_ase_bulk:
            general.check_file_exists(self.geom_file)
            general.check_file_extension(self.geom_file,'.xyz')


   # ------------------------------------------- #
   # ------- Read distances/angles input ------- #

   def read_input(self,what):
      """
      Reads an input file containing distances or angles.
  
      Args:
          what (str): Specifies whether to read 'distances' or 'angles'.
  
      Returns:
          None: Populates the respective list attribute (`self.distances` or `self.angles`).
  
      Raises:
          RuntimeError: If a blank line or multiple values are found in a single line.
  
      Notes:
          - Reads a file line by line and converts values to floats.
          - Ensures valid input formatting before processing.
      """

      if (what=='distances'):
         with open(self.distances_input) as infile:
            for line in infile:
               if len(line.split()) == 0: output.error(f'blank line found in distances input file "{self.distances_input}"')
               if len(line.split())  > 1: output.error(f'more than one distance found in a single line of input file "{self.distances_input}"')

               self.distances.append(float(line.split()[0])) 

      elif (what=='angles'):
         with open(self.angles_input) as infile:
            for line in infile:
               if len(line.split()) == 0: output.error(f'blank line found in angles input file "{self.angles_input}"')
               if len(line.split())  > 1: output.error(f'more than one angle found in a single line of input file "{self.angles_input}"')

               self.angles.append(float(line.split()[0])) 

   # ---------------------------------------- #
   # ------- Change translation sense ------- #
   
   def change_trans_sense(self):
      """
      Reverses the direction factor for translations.

      Returns:
          input_class: The modified `input_class` instance with updated `dir_factor`.

      Notes:
          - Multiplies all direction factors by -1.
          - Used to invert translation direction dynamically.
      """

      self.dir_factor = [-x for x in self.dir_factor] 

      return(self)



