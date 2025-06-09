import math
import copy

from geom.classes import molecule
from geom.functions import tools, general, output

# -------------------------------------------------------------------------------------
def select_case(inp):
   """
   Selects and executes a small task based on user input.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Calls the corresponding function.

   Notes:
       - Supports calculating the minimum distance, geometrical center,
         specular transformation, and merging geometries.
   """

   if (inp.min_dist):      min_dist(inp) 
   if (inp.geom_center):   geom_center(inp) 
   if (inp.geom_specular): geom_specular(inp)
   if (inp.merge):         merge_geoms(inp)
   if (inp.create_dimer):  create_dimer(inp)
# -------------------------------------------------------------------------------------
def min_dist(inp):
   """
   Computes the minimum distance between two molecular geometries.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Prints the calculated minimum distance.

   Notes:
       - Reads the geometries of two molecules from input files.
       - Uses `calc_min_distance` from `tools.py` to determine the minimum distance.
       - Outputs the result via `output.print_min_dist()`.
   """

   # Check input
   inp.check_input_case()   
 
   # Initialize molecules and read geometries
   mol_1 = molecule.molecule()
   mol_2 = molecule.molecule()

   mol_1.read_geom(inp.geom1_file,False)
   mol_2.read_geom(inp.geom2_file,False)
 
   # Calc min distance
   distance = tools.calc_min_distance(mol_1,mol_2)

   # Print calculated minimum distance
   output.print_min_dist(inp,distance)
# -------------------------------------------------------------------------------------
def geom_center(inp):
   """
   Calculates the geometrical center of a molecule.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Prints the geometrical center coordinates.

   Notes:
       - Reads the molecular geometry from an input file.
       - Computes the center using the molecule’s atomic coordinates.
       - Outputs the computed center via `output.print_geom_center()`.
   """

   # Check input
   inp.check_input_case()   
 
   # Initialize molecule and read geometry
   mol = molecule.molecule()
   mol.read_geom(inp.geom_file,False)
 
   output.print_geom_center(inp,mol.xyz_center)
# -------------------------------------------------------------------------------------
def geom_specular(inp):
   """
   Creates a specular (mirror image) geometry of a molecule.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the mirrored geometry.

   Notes:
       - Reads the molecular geometry from an input file.
       - Reflects the molecule across the x-axis.
       - Translates the mirrored molecule by a shift distance (5 Å + molecule width).
       - Saves the new geometry file.
   """

   # Check input
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()
 
   # Initialize molecule and read geometry
   mol = molecule.molecule()
   mol.read_geom(inp.geom_file,True)
 
   # Create specular geometry along x and move at 5 Å 
   shift = (mol.xyz_max[0] - mol.xyz_min[0]) + 5.0
   dir_factor = [1.0,0.0,0.0]

   mol.xyz[0,:] = -mol.xyz[0,:]

   mol.translate_geom(shift,dir_factor)
 
   # Save specular geometry
   output.print_geom(mol, inp.geom_file[:-4]+'_000_mirror')

   # Close and save logfile
   #output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
def merge_geoms(inp):
   """
   Merges two molecular geometries into a single structure.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the merged geometry.

   Notes:
       - Reads two molecular geometries from input files.
       - Uses `merge_geoms` from `tools.py` to combine the structures.
       - Saves the merged structure to an output file.
   """

   # Check input
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()
 
   # Initialize molecules and read geometries
   mol_1 = molecule.molecule()
   mol_2 = molecule.molecule()

   mol_1.read_geom(inp.geom1_file,False)
   mol_2.read_geom(inp.geom2_file,False)
 
   # Merge two geometries
   mol_3 = tools.merge_geoms(inp,mol_1,mol_2)

   # Save merged geometry
   file_geom_merged = f"{inp.geom1_file[:-4]}_MERGED_{inp.geom2_file[:-4]}"
   output.print_geom(mol_3, file_geom_merged)

   # Close and save logfile
   output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
