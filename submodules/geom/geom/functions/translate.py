import math

from geom.classes import molecule, parameters
from geom.functions import general, tools, output

param = parameters.parameters()

# -------------------------------------------------------------------------------------
def select_case(inp):
   """ 
   Selects the appropriate translation function based on user input.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Calls the corresponding translation function.
   """

   if (inp.translate_controlled_distance): translate_controlled_distance(inp)
   if (inp.translate_1):                   translate_1(inp)
# -------------------------------------------------------------------------------------
def translate_controlled_distance(inp):
   """ 
   Translates a second molecule to a controlled distance from a fixed first molecule.

   This function moves `geom2` while keeping `geom1` fixed, ensuring that the minimum 
   distance between the two structures meets the user-defined criteria.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the translated geometry at the optimized distances.

   Notes:
       - Reads the molecular geometries from input files.
       - Performs an initial shift to move `geom2` far from `geom1`.
       - Adjusts the translation direction based on calculated minimum distance.
       - Iteratively optimizes the distance between `geom1` and `geom2`.
       - Saves the final translated geometry for each desired distance.
       - If optimization fails, raises an error.
   """

   # Check input, create results folder, initialize logfile
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()

   # Initialize molecules and read geometries
   mol_1 = molecule.molecule()
   mol_2 = molecule.molecule()

   mol_1.read_geom(inp.geom1_file,inp.move_geom_1_to_000)
   mol_2.read_geom(inp.geom2_file,inp.move_geom_2_to_000)
   
   # =============================================================
   # STEP 1: Perform an initial shift to move geom2 below min_dist
   # =============================================================
   
   # 1.1 Move geom to a big distance
   mol_2.translate_geom(99999999.9,inp.dir_factor)

   # 1.2 Check if dist_opt < min_dist to evaluate dir_axis sense (-->  or <--)
   dist_ini = tools.calc_min_distance(mol_1,mol_2)
   
   shift = dist_ini - param.min_dist_translate + 0.1
   mol_2.translate_geom(shift, inp.dir_factor)

   dist_opt = tools.calc_min_distance(mol_1,mol_2)

   # -- If not, perform the same shift changing dir_axis sense
   if (dist_opt > dist_ini):

      # Modify translation sense
      inp.change_trans_sense()

      shift = dist_opt - param.min_dist_translate + 0.5
      mol_2.translate_geom(shift, inp.dir_factor)

      # Calculate new distance
      dist_opt_2 = tools.calc_min_distance(mol_1,mol_2)

      if(dist_opt_2 > dist_opt): output.error('something wrong happened. Please, check the geometries and try again.')
   
      dist_opt = dist_opt_2
   
   # 1.3 Now that we are below min_dist, change translation sense.
   inp.change_trans_sense()

   # Save the final minimum distance
   dist_pre = dist_opt


   # ===========================================
   # STEP 2: Optimize geom2 to minimum distances
   # ===========================================

   if (inp.verbose): output.print_optimization_starts()

   for distance in inp.distances:
      
      if (inp.verbose): output.print_optimizing_distance(distance)

      # Calculate difference between previous distance and desired one
      shift = distance - dist_pre
       
      # Initial shift
      mol_2.translate_geom(shift - 0.1, inp.dir_factor) # -0.1 as small shift to ensure translation sense

      # Calculate distance
      dist_new = tools.calc_min_distance(mol_1,mol_2) 

      if (dist_new > distance):
   
         mol_2.translate_geom(shift - 0.5, inp.dir_factor) # Another shift with less tolerance to avoid problems

         # Calculate distance
         dist_new = tools.calc_min_distance(mol_1,mol_2)

         if (dist_new > distance):

            mol_2.translate_geom(shift - 0.5, inp.dir_factor) # Try again #1

            # Calculate distance
            dist_new = tools.calc_min_distance(mol_1,mol_2)

            if (dist_new > distance):

               mol_2.translate_geom(shift - 0.5, inp.dir_factor) # Try again #2

               # Calculate distance
               dist_new = tools.calc_min_distance(mol_1,mol_2)

               if (dist_new > distance): output.error(f'optimization error 1. Distance could not be optimized: dist_new = {dist_new} ; min_dist = {distance}')
   
      # Move slowly towards desired distance
      diff_dist = abs(dist_new - distance)
   
      # If after first translation we are quite far from the desire distance
      # apply a new shift to be closer   
      if (diff_dist > param.convergence and diff_dist > 0.02):
         mol_2.translate_geom(diff_dist - param.convergence, inp.dir_factor)   

         # Calculate distance
         dist_new = tools.calc_min_distance(mol_1,mol_2)
        
         if (dist_new > distance): output.error('optimization error 2')
   

      while (diff_dist > param.convergence):
         mol_2.translate_geom(param.convergence_step, inp.dir_factor)

         # Recalculate minimum distance and check difference with desired distance
         dist_new  = tools.calc_min_distance(mol_1,mol_2)

         diff_dist = abs(dist_new - distance)

         if (inp.verbose): output.print_computed_distance(dist_new)
         if (dist_new > distance): output.error('reduce convergence criteria for distance optimization')
      
         if(diff_dist < param.convergence): 
            if (inp.verbose): output.print_convergence_achieved(dist_new)

            #output.save_distance_opt(out_log,distance,dist_new,inp.dir_axis_input) # Save to logfile
   
            # Save distance-optimized geometry
            dist_new_rounded = math.ceil(dist_new * 100) / 100
            inp.file_geom2_translated = f"{inp.geom2_file[:-4]}_{inp.dir_axis_input}_d_{dist_new_rounded:.2f}"
            
            output.print_geom(mol_2, inp.file_geom2_translated)

      dist_pre = dist_new
   
   # Close and save logfile
   #output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
def translate_1(inp):
   """ 
   Translates a single molecule by a given shift distance.

   Args:
       inp (input_class): An instance containing input parameters.

   Returns:
       None: Saves the translated geometry.

   Notes:
       - Reads the molecular geometry from an input file.
       - Translates the molecule by the specified shift along a defined axis.
       - Saves the new geometry after translation.
   """

   # Check input, create results folder, initialize logfile
   inp.check_input_case()   
   general.create_results_geom()
   #out_log = output.logfile_init()

   # Initialize molecule and read geometry
   mol = molecule.molecule()
   mol.read_geom(inp.geom_file,inp.move_geom_to_000)
 
   # Translate 
   mol.translate_geom(inp.shift_t1,inp.dir_factor)

   # Save shifted geometry
   shift_rounded = math.ceil(inp.shift_t1 * 100) / 100
   file_geom_translated = f"{inp.geom_file[:-4]}_{inp.dir_axis_input}_d_{shift_rounded:.2f}"

   output.print_geom(mol, file_geom_translated)

   # Close and save logfile
   #output.logfile_close(out_log)
# -------------------------------------------------------------------------------------
