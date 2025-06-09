import numpy as np
import random

from geom.classes import parameters

from geom.functions import output
from ase.cluster import Icosahedron, Octahedron, Decahedron

class molecule:
   """
   Represents a molecular system and provides methods for geometry transformations, filtering, and nanostructure creation.

   This class includes:
       - Reading and storing atomic coordinates from XYZ files.
       - Geometrical transformations (translation, centering, mirroring).
       - Filtering atoms based on predefined geometric shapes (sphere, cylinder, ribbon, disk, triangle, etc.).
       - Creating nanostructures (icosahedra, cuboctahedra, decahedra) using ASE.
       - Generating alloy structures by random atom substitutions.

   Dependencies:
       - NumPy: Used for numerical operations on atomic coordinates.
       - ASE: Used for nanostructure generation (clusters and nanoparticles).
       - `functions.output`: Handles error messaging and output operations.
       - `classes.parameters`: Provides lattice constants and structural parameters.

   Attributes:
       atoms (list[str]): List of atom types in the molecule.
       nAtoms (int): Number of atoms in the molecule.
       xyz_center (numpy.ndarray): Coordinates of the geometrical center.
       xyz_min (numpy.ndarray): Minimum coordinate values along each axis.
       xyz_max (numpy.ndarray): Maximum coordinate values along each axis.
   """

   def __init__(self):
      """
      Initializes a molecule object with default attributes.

      Attributes Initialized:
          - `atoms` (list[str]): List to store atom types.
          - `nAtoms` (int): Number of atoms in the molecule.
          - `xyz` (numpy.ndarray): 3×N array storing atomic coordinates.
          - `xyz_center` (numpy.ndarray): Coordinates of the geometrical center.
          - `xyz_min` (numpy.ndarray): Minimum coordinate values along each axis.
          - `xyz_max` (numpy.ndarray): Maximum coordinate values along each axis.

      Notes:
          - `xyz_center`, `xyz_min`, and `xyz_max` are initialized as zero arrays.
          - The molecule object is updated when geometry is read from an XYZ file.
      """
       
      self.atoms = []

      self.nAtoms = 0

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

   # ----------------------------------------- #
   # ------- Translate geometry to 000 ------- #
   
   def trans_geom_center_to_000(self):
      """
      Translates the molecular geometry such that its geometrical center moves to (0,0,0).
      
      Returns:
          None: Updates the molecule's atomic coordinates in place.
      """

      for i in range(self.nAtoms):
         self.xyz[0][i] = self.xyz[0][i] - self.xyz_center[0] 
         self.xyz[1][i] = self.xyz[1][i] - self.xyz_center[1]
         self.xyz[2][i] = self.xyz[2][i] - self.xyz_center[2]

   # ----------------------------- #
   # ------- Read geometry ------- #
   
   def read_geom(self, geom_file, translate_geom_to_000):
      """
      Reads an XYZ file and stores atomic coordinates.
  
      Args:
          geom_file (str): Path to the XYZ file containing the molecular geometry.
          translate_geom_to_000 (bool): If True, translates the molecular geometry to the origin.
  
      Returns:
          molecule: The molecule object with updated atomic coordinates.
  
      Notes:
          - Computes the geometrical center of the structure.
          - Stores the minimum and maximum coordinates for bounding box calculations.
          - Can output a translated structure if `translate_geom_to_000` is set to True.
      """

      with open(geom_file, 'r') as infile:
         self.nAtoms = int(infile.readline())
         infile.readline()
       
         if self.nAtoms <= 0: output.error('Corrupt geometry file "' + geom_file + '"')
    
         self.atoms = []
         self.xyz = np.zeros((3,self.nAtoms))
       
         i = 0
         for line in infile:
            line = line.split()

            self.atoms.append(line[0])
            self.xyz[0][i] = float(line[1])
            self.xyz[1][i] = float(line[2])
            self.xyz[2][i] = float(line[3])

            i += 1

            if i > self.nAtoms: output.error('Corrupt geometry file "' + geom_file + '"')

      # Calculate geometrical center
      self.xyz_center[0] = np.mean(self.xyz[0,:]) 
      self.xyz_center[1] = np.mean(self.xyz[1,:]) 
      self.xyz_center[2] = np.mean(self.xyz[2,:]) 

      # Translate geometrical center to 000 and save, if requested  
      if (translate_geom_to_000):
         self.trans_geom_center_to_000()
         output.print_geom(self,geom_file[:-4]+'_000')

      # Save maximun/minimum coordinates limits
      self.xyz_max[0] = np.max(self.xyz[0,:])
      self.xyz_max[1] = np.max(self.xyz[1,:])
      self.xyz_max[2] = np.max(self.xyz[2,:])

      self.xyz_min[0] = np.min(self.xyz[0,:])
      self.xyz_min[1] = np.min(self.xyz[1,:])
      self.xyz_min[2] = np.min(self.xyz[2,:])

      return(self)


   # -------------------------------------------- #
   # ------- Change atom types of geometry------- #
   
   def change_atomtype(self, new_atomtype):
      """
      Changes the atom type for all atoms in the molecule.
  
      Args:
          new_atomtype (str): The new atom type to assign.
  
      Returns:
          molecule: Updated molecule object with modified atom types.
  
      Notes:
          - This function is mainly used in core-shell structure creation.
      """

      self.atoms = [new_atomtype] * self.nAtoms
       
      return(self)


   # ------------------------------------------------- #
   # ------- Translate geometry along dir_axis ------- #
   
   def translate_geom(self, shift, dir_factor):
      """
      Translates the molecular geometry along a specified axis.
  
      Args:
          shift (float): The amount to shift the molecule in angstroms.
          dir_factor (list[float]): The direction factors along x, y, and z axes.
  
      Returns:
          molecule: The translated molecule object.
  
      Notes:
          - Positive and negative direction factors determine the translation axis.
      """

      for i in range(self.nAtoms):
         self.xyz[0,i] = self.xyz[0,i] + dir_factor[0] * shift  
         self.xyz[1,i] = self.xyz[1,i] + dir_factor[1] * shift 
         self.xyz[2,i] = self.xyz[2,i] + dir_factor[2] * shift 
   
      return(self)

   # ------------------------------------------------ #
   # ------- Filter XYZ within a ribbon shape ------- #
   
   def filter_xyz_graphene_to_ribbon(self, inp):
      """
      Filters the graphene sheet to create a ribbon of a specified length and width.
  
      Args:
          inp (input_class): The input parameters containing X and Y lengths.
  
      Returns:
          molecule: The filtered molecule containing only atoms within the ribbon.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
   
      # Define X and Y bounds based on inp parameters
      x_min, x_max = -inp.X_length / 2, inp.X_length / 2
      y_min, y_max = -inp.Y_length / 2, inp.Y_length / 2

      # Condition for points to be within the ribbon
      condition = (
          (x_min <= x) & (x <= x_max) &
          (y_min <= y) & (y <= y_max)
      )
   
      # Filter coordinates based on condition
      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]
   
      # Update the atom list and geometry based on filtered data
      self.nAtoms = len(x_filtered)
      self.atoms = [inp.atomtype] * self.nAtoms  # Atom types remain consistent
   
      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)
   
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))
   
      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)
   
      # Save maximum/minimum coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)
   
      return self

   # ---------------------------------------------- #
   # ------- Filter XYZ within a disk shape ------- #
    
   def filter_xyz_graphene_to_disk(self, inp):
      """
      Filters the graphene sheet to create a disk of a given radius.
   
      Args:
          inp (input_class): The input parameters containing the disk radius.
   
      Returns:
          molecule: The filtered molecule containing only atoms within the disk.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
   
      # Condition for points to be within the disk
      condition = (x**2 + y**2 <= inp.radius**2)
   
      # Filter coordinates based on condition
      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]
   
      # Update the atom list and geometry based on filtered data
      self.nAtoms = len(x_filtered)
      self.atoms = [inp.atomtype] * self.nAtoms  # Atom types remain consistent
   
      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)
   
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))
   
      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)
   
      # Save maximum/minimum coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)
   
      return self

   # ---------------------------------------------- #
   # ------- Filter XYZ within a ring shape ------- #

   def filter_xyz_graphene_to_ring(self, inp):
      """
      Filters the graphene sheet to create a ring structure with an inner and outer radius.

      Args:
          inp (input_class): The input parameters containing inner and outer radius values.

      Returns:
          molecule: The filtered molecule containing only atoms within the ring.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
   
      # Condition for points to be within the ring (between radius_in and radius_out)
      condition = (inp.radius_in**2 <= x**2 + y**2) & (x**2 + y**2 <= inp.radius_out**2)
   
      # Filter coordinates based on condition
      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]
   
      # Update the atom list and geometry based on filtered data
      self.nAtoms = len(x_filtered)
      self.atoms = [inp.atomtype] * self.nAtoms  # Atom types remain consistent
   
      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)
   
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))
   
      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)
   
      # Save maximum/minimum coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)
   
      return self 

   # -------------------------------------------------- #
   # ------- Filter XYZ within a triangle shape ------- #

   def filter_xyz_graphene_to_triangle(self, inp):
      """
      Filters the graphene sheet to create a triangular structure based on the edge type.
   
      Args:
          inp (input_class): The input parameters containing side length and edge type.
   
      Returns:
          molecule: The filtered molecule containing only atoms within the triangle.
      """
   
      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
   
      # Define triangular region based on edge type
      if inp.graphene_edge_type == 'zigzag':
          # Armchair triangle aligned with the armchair direction
          condition = (
              (y >= 0) &
              (y <= inp.side_length * np.sqrt(3) / 2) &
              (x >= -y / np.sqrt(3)) &
              (x <= y / np.sqrt(3))
          )
      elif inp.graphene_edge_type == 'armchair':
          # Zigzag triangle aligned with the zigzag direction
          condition = (
              (x >= 0) &
              (x <= inp.side_length) &
              (y >= -x / np.sqrt(3)) &
              (y <= x / np.sqrt(3))
          )
      else:
          raise ValueError("Invalid edge type. Must be 'armchair' or 'zigzag'.")
   
      # Filter coordinates based on condition
      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]
   
      # Update the atom list and geometry based on filtered data
      self.nAtoms = len(x_filtered)
      self.atoms = [inp.atomtype] * self.nAtoms  # Atom types remain consistent
   
      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)
   
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))
   
      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)
   
      # Save maximum/minimum coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)
   
      return self


   # ---------------------------------------------- #
   # ------- Remove graphene dangling bonds ------- #

   def remove_dangling_bonds_graphene(self, inp):
      """
      Removes dangling bonds from a generated graphene structure.
  
      Args:
          inp (input_class): The input parameters for graphene filtering.
  
      Returns:
          molecule: The molecule object with dangling atoms removed.
  
      Notes:
          - Performs up to 3 iterations to remove all dangling atoms.
          - If dangling atoms remain after 3 iterations, an error is raised.
      """

      #debug
      #print("  Removing dangling C atoms on generated structure...")
      #print("")
   
      def get_neighbors(atom_index, cutoff=1.5):
          """Find neighbors of an atom within a distance cutoff."""
          distances = np.linalg.norm(self.xyz.T - self.xyz.T[atom_index], axis=1)
          neighbors = np.where((distances > 0) & (distances <= cutoff))[0]
          return neighbors
   
      max_iterations = 3
      for iteration in range(max_iterations):
          dangling_atoms = []
          for i in range(self.nAtoms):
              neighbors = get_neighbors(i)
              if len(neighbors) < 2:  # Dangling bond if fewer than 2 neighbors
                  dangling_atoms.append(i)
   
          # If no dangling atoms are found, exit the loop
          if not dangling_atoms:
              #debug
              #print("")
              #print(f"  --> All dangling bonds removed after {iteration + 1} iteration(s).")
              return self
   
          # Remove dangling atoms
          keep_atoms = np.setdiff1d(np.arange(self.nAtoms), dangling_atoms)
          self.xyz = self.xyz[:, keep_atoms]
          self.nAtoms = len(keep_atoms)
          self.atoms = [self.atoms[i] for i in keep_atoms]
   
          # Recalculate geometry
          self.xyz_center = np.mean(self.xyz, axis=1)
          self.xyz_max = np.max(self.xyz, axis=1)
          self.xyz_min = np.min(self.xyz, axis=1)
   
          #debug
          #print(f"  - Iteration {iteration + 1}: Removed {len(dangling_atoms)} dangling atom(s).")
   
      # If dangling bonds remain after 3 iterations, raise an error
      dangling_atoms = []
      for i in range(self.nAtoms):
          neighbors = get_neighbors(i)
          if len(neighbors) < 2:
              dangling_atoms.append(i)
   
      if dangling_atoms:
          output.error(f"{len(dangling_atoms)} dangling bonds on generated graphene sheet could not be eliminated after {max_iterations} iterations.")
   
      return self


   # ---------------------------------------- #
   # ------- Filter XYZ within sphere ------- #
   
   def filter_xyz_in_sphere(self,inp):
      """
      Filters atoms in the molecular geometry, keeping only those inside a sphere of a given radius.
   
      Args:
          inp (input_class): The input parameters containing sphere radius and center coordinates.
   
      Returns:
          None: Updates the molecule by removing atoms outside the sphere.
   
      Notes:
          - The sphere is centered at `inp.sphere_center`.
          - Any atom with a distance greater than `inp.radius` from the sphere center is removed.
          - This function is commonly used to create spherical nanoparticles.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
      
      # Condition for points to be within the sphere 
      condition = ((x-inp.sphere_center[0])**2 + 
                   (y-inp.sphere_center[1])**2 + 
                   (z-inp.sphere_center[2])**2 <= (inp.radius)**2)

      x_filtered = self.xyz[0, condition]
      y_filtered = self.xyz[1, condition]
      z_filtered = self.xyz[2, condition]

      # Fill previous geometry with current structure
      self.nAtoms = len(x_filtered)

      self.atoms = []
      self.atoms = [inp.atomtype] * self.nAtoms

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

      self.xyz = np.zeros((3,self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save maximum/minimum coordinates limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return(self)


   # ------------------------------------------ #
   # ------- Filter XYZ within cylinder ------- #
   
   def filter_xyz_in_cylinder(self, inp):
      """
      Filters atoms in the molecular geometry, keeping only those inside a cylindrical region.
  
      Args:
          inp (input_class): The input parameters containing cylinder radius, height, and axis orientation.
  
      Returns:
          None: Updates the molecule by removing atoms outside the defined cylindrical region.
  
      Notes:
          - The cylinder is aligned along a specified axis (typically the z-axis).
          - Any atom outside the given radius (`inp.radius`) or height (`inp.rod_length`) is removed.
          - Used for generating rod-like nanostructures, such as metallic nanorods.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
       
      # Compute base radius from the rod width
      # and length subtracting the radius of the spheres at the extremes
      inp.sphere_center = [0.0,0.0,0.0]
      radius = inp.rod_width / 2.0
      length = inp.rod_length - inp.rod_width 

      # Condition for points to be within the cylinder 
      if inp.main_axis == "x":
        condition = (
            ((y - inp.sphere_center[1])**2 + (z - inp.sphere_center[2])**2 <= radius**2) &
            (x >= inp.sphere_center[0] - length/2.0) &
            (x <= inp.sphere_center[0] + length/2.0)
        )

      elif inp.main_axis == "y":
          condition = (
              ((x - inp.sphere_center[0])**2 + (z - inp.sphere_center[2])**2 <= radius**2) &
              (y >= inp.sphere_center[1] - length/2.0) &
              (y <= inp.sphere_center[1] + length/2.0)
          )

      elif inp.main_axis == "z":
          condition = (
              ((x - inp.sphere_center[0])**2 + (y - inp.sphere_center[1])**2 <= radius**2) &
              (z >= inp.sphere_center[2] - length/2.0) &
              (z <= inp.sphere_center[2] + length/2.0)
          )

      x_filtered = self.xyz[0, condition]
      y_filtered = self.xyz[1, condition]
      z_filtered = self.xyz[2, condition]

      # Fill previous geometry with current structure
      self.nAtoms = len(x_filtered)

      self.atoms = []
      self.atoms = [inp.atomtype] * self.nAtoms

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

      self.xyz = np.zeros((3,self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save maximum/minimum coordinates limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return(self)


   # ----------------------------------------------------- #
   # ------- Filter XYZ within elliptic paraboloid ------- #
   
   def filter_xyz_in_elliptic_paraboloid(self, inp):
      r"""
      Filters atoms in the molecular geometry, keeping only those inside an elliptic paraboloid.
  
      Args:
          inp (input_class): The input parameters containing coefficients (`a`, `b`, `c`) and height constraints.
  
      Returns:
          None: Updates the molecule by removing atoms outside the elliptic paraboloid.
  
      Notes:
          - The elliptic paraboloid is defined by the equation:
            \[
            z = c - \frac{x^2}{a} - \frac{y^2}{b}
            \]
          - Atoms outside this equation’s boundary or beyond `inp.z_max_paraboloid` are removed.
          - Used for modeling paraboloidal nanostructures, such as scanning probe tips.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
      
      # Calculate the paraboloid limit for each (x, y)
      paraboloid_limit = inp.elliptic_parabola_a * x**2 + inp.elliptic_parabola_b * y**2 + inp.elliptic_parabola_c

      # Condition for points to be within the paraboloid 
      condition = ((z >= inp.z_min) &
                   (z <= inp.z_max) &
                   (z >= paraboloid_limit))

      x_filtered = self.xyz[0, condition]
      y_filtered = self.xyz[1, condition]
      z_filtered = self.xyz[2, condition]

      # Fill previous geometry with current structure
      self.nAtoms = len(x_filtered)

      self.atoms = []
      self.atoms = [inp.atomtype] * self.nAtoms

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

      self.xyz = np.zeros((3,self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save maximum/minimum coordinates limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return(self)


   # ----------------------------------------------------- #
   # ------- Filter XYZ within square-base pyramid ------- #
   
   def filter_xyz_in_pyramid(self, inp, centers, planes):
      """
      Filters atoms in the molecular geometry, keeping only those inside a defined pyramid.
  
      Args:
          inp (input_class): The input parameters containing the atomic type.
          centers (dict): Dictionary with center coordinates of the pyramid base and apex.
                          Expected keys: `"center_1"`, `"center_3"`, `"center_4"`, `"center_5"`.
          planes (dict): Dictionary with normal vectors and offsets for the pyramid planes.
                         Expected keys: `"n_125"`, `"n_235"`, `"n_345"`, `"n_415"`.
  
      Returns:
          molecule: The molecule object with updated atomic coordinates.
  
      Notes:
          - The function first checks whether each atom lies within the bounding box of the pyramid.
          - Then, the atoms are filtered using the plane equations for each pyramid face.
          - The molecule is updated with only the atoms that meet these conditions.
          - Computes the new geometrical center, bounding box limits, and atom count.
          - The pyramid structure is defined using four planes and a set of center points.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]

      # Condition for points to be within the pyramid
      condition = (
          (centers["center_1"][2] <= z) & (z <= centers["center_5"][2]) &
          (centers["center_4"][0] <= x) & (x <= centers["center_1"][0]) &
          (centers["center_3"][1] <= y) & (y <= centers["center_4"][1]) &
          (planes["n_125"][0][0] * x + planes["n_125"][0][1] * y + planes["n_125"][0][2] * z >= -planes["n_125"][1]) &
          (planes["n_235"][0][0] * x + planes["n_235"][0][1] * y + planes["n_235"][0][2] * z >= -planes["n_235"][1]) &
          (planes["n_345"][0][0] * x + planes["n_345"][0][1] * y + planes["n_345"][0][2] * z >= -planes["n_345"][1]) &
          (planes["n_415"][0][0] * x + planes["n_415"][0][1] * y + planes["n_415"][0][2] * z >= -planes["n_415"][1])
      )

      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]

      # Fill previous geometry with current structure
      self.nAtoms = len(x_filtered)

      self.atoms = []
      self.atoms = [inp.atomtype] * self.nAtoms

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

      self.xyz = np.zeros((3,self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save maximum/minimum coordinates limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return(self)


   # -------------------------------------- #
   # ------- Filter XYZ within cone ------- #
   
   def filter_xyz_in_cone(self, inp):
      r"""
      Filters atoms in the molecular geometry, keeping only those inside a cone.
  
      Args:
          inp (input_class): The input parameters containing:
              - `radius` (float): Base radius of the cone.
              - `z_max` (float): Height of the cone.
              - `atomtype` (str): The atomic type for the filtered structure.
  
      Returns:
          molecule: The molecule object with updated atomic coordinates.
  
      Notes:
          - The cone is aligned along the z-axis with the apex at the origin `(0,0,0)`.
          - The equation defining the cone is:
            \[
            x^2 + y^2 \leq \left(\frac{\text{radius}}{\text{z_max}}\right)^2 \cdot z^2
            \]
          - Atoms are filtered based on the cone equation and must satisfy `0 ≤ z ≤ inp.z_max`.
          - The function recomputes:
              - The geometrical center of the new structure.
              - The bounding box limits (`xyz_min`, `xyz_max`).
              - The updated list of atoms inside the cone.
          - Typically used to shape nano-cones or probe-like structures.
      """

      x = self.xyz[0, :]
      y = self.xyz[1, :]
      z = self.xyz[2, :]
 
      # Condition for points inside the cone
      condition = (x**2 + y**2 <= (inp.radius / inp.z_max)**2 * z**2) & (z >= 0) & (z <= inp.z_max)
      
      # Select the points that satisfy the condition
      x_filtered = x[condition]
      y_filtered = y[condition]
      z_filtered = z[condition]
      
      # Fill previous geometry with current structure
      self.nAtoms = len(x_filtered)

      self.atoms = []
      self.atoms = [inp.atomtype] * self.nAtoms

      self.xyz_center = np.zeros(3)
      self.xyz_min    = np.zeros(3)
      self.xyz_max    = np.zeros(3)

      self.xyz = np.zeros((3,self.nAtoms))
      self.xyz = np.vstack((x_filtered, y_filtered, z_filtered))

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save maximum/minimum coordinates limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return(self)


   # --------------------------------- #
   # ------- Create icosahedra ------- #

   def create_icosahedra(self, inp):
       """
       Generates an atomically perfect icosahedral nanoparticle using ASE.
   
       Args:
           inp (input_class): The input parameters containing radius and lattice constants.
   
       Returns:
           molecule: The generated icosahedral structure.
       """

       # Extract lattice constant 
       param = parameters.parameters()
       lattice_constant = param.lattice_constant.get(inp.atomtype)

       # Convert radius to number of shells
       #noshells = round((2 * inp.radius) / (np.sqrt(2) * lattice_constant))  # Correct scaling for FCC growth
       noshells = round((2 * inp.radius) / (np.sqrt(2) * lattice_constant)) + 1  # Ensure full layer growth

       # Generate icosahedral cluster using ASE
       icosahedron = Icosahedron(symbol=inp.atomtype.capitalize(), noshells=noshells, latticeconstant=lattice_constant)

       # Extract atom coordinates and store
       positions = icosahedron.get_positions()

       self.nAtoms = len(positions) 
       self.atoms = [inp.atomtype] * self.nAtoms  # Assign atom type to all atoms

       # Store positions
       self.xyz = np.zeros((3, self.nAtoms))
       self.xyz = positions.T  # Transpose so shape matches (3, nAtoms)

       # Calculate geometrical center
       self.xyz_center = np.mean(self.xyz, axis=1)

       # Save max/min coordinate limits
       self.xyz_max = np.max(self.xyz, axis=1)
       self.xyz_min = np.min(self.xyz, axis=1)

       return self


   # ----------------------------------- #
   # ------- Create cuboctahedra ------- #

   def create_cuboctahedra(self, inp):
      """
      Generates an atomically perfect cuboctahedral nanoparticle using ASE.
  
      Args:
          inp (input_class): The input parameters containing radius and lattice constants.
  
      Returns:
          molecule: The generated cuboctahedral structure.
      """
   
      # Extract lattice constant
      param = parameters.parameters()
      lattice_constant = param.lattice_constant.get(inp.atomtype)
   
      # Calculate cutoff ang length based on radius
      cutoff = ((inp.radius * 2) / (np.sqrt(2) * lattice_constant)) + 1.00  # Add a small buffer
      length = int(2 * cutoff + 1)  # Convert to ASE-compatible parameter

      # Convert radius to ASE-compatible cutoff value
      max_cutoff = (length - 1) / 2
      cutoff = min(cutoff, max_cutoff)  # Ensure cutoff does not exceed the limit

      # Generate cuboctahedral cluster using ASE
      cuboctahedron = Octahedron(symbol=inp.atomtype.capitalize(), length=length, cutoff=cutoff, latticeconstant=lattice_constant)
   
      # Extract atom coordinates
      positions = cuboctahedron.get_positions()
   
      # Store data inside self
      self.nAtoms = len(positions)
      self.atoms = [inp.atomtype] * self.nAtoms  # Assign atom type to all atoms
   
      # Store positions
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = positions.T  # Transpose so shape matches (3, nAtoms)
   
      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)
   
      # Save max/min coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)
   
      return self  # Return updated object


   # -------------------------------- #
   # ------- Create decahedra ------- #

   def create_decahedra(self, inp):
      """
      Generates an atomically perfect decahedral nanoparticle using ASE.
  
      Args:
          inp (input_class): The input parameters containing radius and lattice constants.
  
      Returns:
          molecule: The generated decahedral structure.
      """

      # Extract lattice constant
      param = parameters.parameters()
      lattice_constant = param.lattice_constant.get(inp.atomtype)

      # Convert radius to ASE-compatible p, q, r values
      p = q = round((2 * inp.radius) / (np.sqrt(2) * lattice_constant))  + 1 # Add a small buffer 
      r = 0  # No Marks re-entrance (standard decahedron)

      # Generate decahedral cluster using ASE
      decahedron = Decahedron(symbol=inp.atomtype.capitalize(), p=p, q=q, r=r, latticeconstant=lattice_constant)

      # Extract atom coordinates
      positions = decahedron.get_positions()

      # Store data inside self
      self.nAtoms = len(positions)
      self.atoms = [inp.atomtype] * self.nAtoms  # Assign atom type to all atoms

      # Store positions
      self.xyz = np.zeros((3, self.nAtoms))
      self.xyz = positions.T  # Transpose so shape matches (3, nAtoms)

      # Calculate geometrical center
      self.xyz_center = np.mean(self.xyz, axis=1)

      # Save max/min coordinate limits
      self.xyz_max = np.max(self.xyz, axis=1)
      self.xyz_min = np.min(self.xyz, axis=1)

      return self  # Return updated object

   
   # ----------------------------------- #
   # ------- Create random alloy ------- #

   def create_alloy(self, inp):
      """
      Generates an alloy by randomly substituting a percentage of atoms.
  
      Args:
          inp (input_class): The input parameters containing alloy composition.
  
      Returns:
          molecule: The modified molecule with alloyed atoms.
  
      Notes:
          - Selects a fraction of atoms and replaces them with the alloy type.
          - The number of atoms replaced is based on `alloy_perc`.
      """

      # Number of atoms to replace 
      replace_indices = [i for i, atom in enumerate(self.atoms)]
      n_replace = int(self.nAtoms * (inp.alloy_perc / 100.0))

      if n_replace > 0:
         # Randomly select indices to replace
         selected_indices = random.sample(replace_indices, n_replace)

         # Replace the selected atoms with the alloy type
         for idx in selected_indices:
             self.atoms[idx] = inp.atomtype_alloy.lower()

         #debug
         #print(f"Replaced {n_replace} {inp.atomtype} atoms with {inp.atomtype_alloy}")

      else:
         output.error('number of atoms to replace in alloy creation n_replace{}')


   # ---------------------------------------------- #
   # ------- Remove metal dangling atoms ------- #

   def remove_dangling_atoms_metals(self, inp):
      """
      Removes dangling atoms from a generated metal structure.
  
      Args:
          inp (input_class): The input parameters for metal filtering.
  
      Returns:
          molecule: The molecule object with dangling atoms removed.
  
      Notes:
          - Performs up to 3 iterations to remove all dangling atoms.
          - If dangling atoms remain after 3 iterations, an error is raised.
      """

      #debug
      #print("  Removing dangling metal atoms on generated structure...")
      #print("")

      param = parameters.parameters()
      cutoff_distance = param.min_dist.get(inp.atomtype) + 0.1 # Add small buffer

      def get_neighbors_metal(atom_index, cutoff):
          """Find neighbors of an atom within a distance cutoff."""
          distances = np.linalg.norm(self.xyz.T - self.xyz.T[atom_index], axis=1)
          neighbors = np.where((distances > 0) & (distances <= cutoff))[0]
          return neighbors

      max_iterations = 3
      for iteration in range(max_iterations):
          dangling_atoms = []
          for i in range(self.nAtoms):
              neighbors = get_neighbors_metal(i, cutoff_distance)
              if len(neighbors) == 0:  # No neighbors found
                  dangling_atoms.append(i)

          # If no dangling atoms are found, exit the loop
          if not dangling_atoms:
              #debug
              #print("")
              #print(f"  --> All dangling bonds removed after {iteration + 1} iteration(s).")
              return self

          # Remove dangling atoms
          keep_atoms = np.setdiff1d(np.arange(self.nAtoms), dangling_atoms)
          self.xyz = self.xyz[:, keep_atoms]
          self.nAtoms = len(keep_atoms)
          self.atoms = [self.atoms[i] for i in keep_atoms]

          # Recalculate geometry
          self.xyz_center = np.mean(self.xyz, axis=1)
          self.xyz_max = np.max(self.xyz, axis=1)
          self.xyz_min = np.min(self.xyz, axis=1)

          #debug
          #print(f"  - Iteration {iteration + 1}: Removed {len(dangling_atoms)} dangling atom(s).")

      # If dangling bonds remain after 3 iterations, raise an error
      dangling_atoms = []
      for i in range(self.nAtoms):
          neighbors = get_neighbors_metal(i, cutoff_distance)
          if len(neighbors) < 1:
              dangling_atoms.append(i)

      if dangling_atoms:
          output.error(f"{len(dangling_atoms)} dangling atoms on generated metal structure could not be eliminated after {max_iterations} iterations.")

      return self


