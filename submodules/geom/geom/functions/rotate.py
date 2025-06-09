import math
import copy

from geom.classes import molecule
from geom.functions import general, tools, output

# -------------------------------------------------------------------------------------
def select_case(inp):
   """ 
   Selects the appropriate rotation function based on user input.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Calls the corresponding rotation function.
   """

   if (inp.rotate_angles): rotate_angles(inp)
   if (inp.rotate_1):      rotate_1(inp)
# -------------------------------------------------------------------------------------
def rotate_angles(inp):
   """ 
   Rotates a molecule at a list of angles specified in the input.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the rotated geometry for each angle.
   
   Notes:
       - Reads the molecule from the specified geometry file.
       - Adjusts angles based on the specified axis direction.
       - Iterates through the list of angles and applies rotation.
       - Saves each rotated geometry as an output file.
   """

   # Check input, create results folder, initialize logfile
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()

   # Initialize molecules and read geometry
   mol = molecule.molecule()
   mol_rot = molecule.molecule()

   mol.read_geom(inp.geom_file,inp.move_geom_to_000)

   # Adjust angles depending on direction
   if inp.dir_axis_input[0] == '-': inp.angles = [360 - angle for angle in inp.angles]
 
   # Rotate over all angles
   for angle in inp.angles:

      mol_rot = tools.rotate(mol,angle,inp.dir_axis_input,mol_rot)

      # Save rotate geometry
      if inp.dir_axis_input[0] == '-': inp.file_geom_rotated = f"{inp.geom_file[:-4]}_{inp.dir_axis_input}_degree_{abs(angle - 360)}"
      if inp.dir_axis_input[0] == '+': inp.file_geom_rotated = f"{inp.geom_file[:-4]}_{inp.dir_axis_input}_degree_{angle}" 
       
      output.print_geom(mol_rot, inp.file_geom_rotated)

      # Close and save logfile
      #output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
def rotate_1(inp):
   """ 
   Rotates a molecule to a specified single angle.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the rotated geometry.

   Notes:
       - Reads the molecule from the specified geometry file.
       - Adjusts the rotation angle based on the specified axis direction.
       - Applies the rotation transformation.
       - Saves the rotated geometry as an output file.
   """

   # Check input, create results folder, initialize logfile
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()

   # Initialize molecules and read geometry
   mol = molecule.molecule()
   mol_rot = molecule.molecule()

   mol.read_geom(inp.geom_file,inp.move_geom_to_000)

   # Adjust angles depending on direction
   if inp.dir_axis_input[0] == '-': inp.angle = 360.0 - inp.angle
 
   # Rotate 
   mol_rot = tools.rotate(mol,inp.angle,inp.dir_axis_input,mol_rot)

   # Save rotate geometry
   if inp.dir_axis_input[0] == '-': inp.file_geom_rotated = f"{inp.geom_file[:-4]}_{inp.dir_axis_input}_degree_{abs(inp.angle - 360)}"
   if inp.dir_axis_input[0] == '+': inp.file_geom_rotated = f"{inp.geom_file[:-4]}_{inp.dir_axis_input}_degree_{inp.angle}" 
       
   output.print_geom(mol_rot, inp.file_geom_rotated)

   # Close and save logfile
   #output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
