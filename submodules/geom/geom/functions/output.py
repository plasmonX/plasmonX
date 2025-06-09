import sys

# -------------------------------------------------------------------------------------
def error(error_message):
   """
   Prints an error message and terminates execution.

   Args:
       error_message (str): The error message to be displayed.

   Returns:
       None: The function exits the program.
   """

   print("")
   print("")
   print("   ERROR: " + error_message)
   print("")
   print("")
   sys.exit()
# -------------------------------------------------------------------------------------
def error_dir_axis(dir_axis_input):
   """
   Prints an error message related to an unsupported axis input and exits.

   Args:
       dir_axis_input (str): The invalid direction axis input.

   Returns:
       None: The function exits the program.

   Notes:
       - Accepted values: `+x`, `+y`, `+z`, `-x`, `-y`, `-z`.
   """

   print(' ')
   print(' ERROR: Sense or direction axis "' + dir_axis_input + '" not supported')
   print(' ')
   print('    Options:')
   print('    --------')
   print('      +x, +y, +z')
   print('      -x, -y, -z')
   print(' ')
   sys.exit()
# -------------------------------------------------------------------------------------
def logfile_init():
   """
   Initializes and opens a logfile for writing.

   Returns:
       file object: An open logfile in write mode (`logfile.txt`).
   """

   out_log = open('results_geom/logfile.txt','w')

   return(out_log)
# -------------------------------------------------------------------------------------
def logfile_close(out_log):
   """
   Closes the logfile.

   Args:
       out_log (file object): The logfile object to close.

   Returns:
       None
   """

   out_log.close()

   return(out_log)
# -------------------------------------------------------------------------------------
def print_geom(molecule,output_file):
   """
   Saves molecular geometry to an XYZ file.

   Args:
       molecule (molecule): The molecule object containing atomic data.
       output_file (str): The name of the output XYZ file.

   Returns:
       None: The function writes the geometry to `results_geom/{output_file}.xyz`.

   Notes:
       - The first line of the XYZ file contains the number of atoms.
       - The second line contains a header.
       - The atomic coordinates are printed with 8 decimal places.
   """

   with open(f'results_geom/{output_file}.xyz', 'w') as out_f:
       out_f.write(f"{molecule.nAtoms}\n")
       out_f.write('Generated with GEOM code\n')

       for i in range(molecule.nAtoms):
           atom, x, y, z = molecule.atoms[i], *molecule.xyz[:, i]
           out_f.write(f'{atom.capitalize():2} {x:20.8f} {y:20.8f} {z:20.8f}\n')
# -------------------------------------------------------------------------------------
def print_optimization_starts():
   """
   Prints a banner indicating the start of the distance optimization process.

   Returns:
       None
   """

   print(' ')
   print(' ')
   print(' ===================================== ')
   print(' DISTANCE OPTIMIZATION PROCESS STARTED')
   print(' ===================================== ')
   print(' ')
# -------------------------------------------------------------------------------------
def print_optimizing_distance(distance):
   """
   Prints a message indicating that distance optimization is in progress.

   Args:
       distance (float): The target distance for optimization.

   Returns:
       None
   """

   print(f'  ------ Optimizing d = {distance} Å ------ ')  
   print('\n')
# -------------------------------------------------------------------------------------
def print_computed_distance(dist):
   """
   Prints the computed distance after translation or optimization.

   Args:
       dist (float): The computed distance.

   Returns:
       None
   """

   print('  Computed distance = ' + str(round(dist,4)) + ' Å')
   print('\n') 
# -------------------------------------------------------------------------------------
def print_convergence_achieved(dist):
   """
   Prints a message indicating that convergence has been achieved.

   Args:
       dist (float): The final optimized distance.

   Returns:
       None
   """

   print('  Convergence achieved to distance ' +  str(round(dist,4)) + ' Å') 
   print('')
# -------------------------------------------------------------------------------------
def save_distance_opt(out_log,distance,dist_new,dir_axis_input):
   """
   Saves the optimized distance information in the logfile.

   Args:
       out_log (file object): The logfile object.
       distance (float): The initial distance target.
       dist_new (float): The final achieved distance.
       dir_axis_input (str): The translation or rotation axis.

   Returns:
       None
   """

   out_log.write(f"\n"
                 f" {' ------ Optimizing d =':>22} {distance:20.8f} {'Å ------ ':>12}\n\n")
   
   out_log.write(f" {'  Convergence achieved to distance':>34} {dist_new:20.8f} {'Å':>5}\n\n\n")
# -------------------------------------------------------------------------------------
def print_normal_termination(inp):
   """
   Prints a normal termination banner if the process completes successfully.

   Args:
       inp (input_class): The input class instance containing execution parameters.

   Returns:
       None

   Notes:
       - This function does not print if `min_dist` or `geom_center` calculations are performed.
   """

   if (not inp.min_dist and not inp.geom_center):
      print(' ')
      print(' ===================================== ')
      print('           NORMAL TERMINATION          ')
      print(' ===================================== ')
      print(' ')
# -------------------------------------------------------------------------------------
def print_min_dist(inp,distance):
   """
   Prints the minimum distance between two geometries.

   Args:
       inp (input_class): The input class instance containing file information.
       distance (float): The computed minimum distance.

   Returns:
       None
   """

   distance = round(distance,4)

   print('')
   print('  -----------------------------------------------')
   print(f'    Geometry 1: {inp.geom1_file}')
   print(f'    Geometry 2: {inp.geom2_file}')
   print('')
   print(f'    Minimum\n    distance  : {distance} Å')
   print('  -----------------------------------------------')
   print('')
# -------------------------------------------------------------------------------------
def print_geom_center(inp,xyz_c):
   """
   Prints the geometrical center of a molecule.

   Args:
       inp (input_class): The input class instance containing file information.
       xyz_c (list[float]): The computed (x, y, z) coordinates of the geometrical center.

   Returns:
       None
   """

   x = round(xyz_c[0],4)
   y = round(xyz_c[1],4)
   z = round(xyz_c[2],4)

   print('')
   print('  ------------------------------------------')
   print(f'    Geometry : {inp.geom_file}')
   print('')
   print(f'    Geometrical \n    Center (xyz,Å)  : {x} {y} {z} ')
   print('  ------------------------------------------')
   print('')
# -------------------------------------------------------------------------------------
