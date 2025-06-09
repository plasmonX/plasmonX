import numpy as np
import pytest
import sys
import os
import filecmp

# Add the project root to sys.path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".")))
test_folder_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))

from geom.classes import input_class
from geom.functions import general, various, translate, rotate, create_geom

# -------------------------------------------------------------------------------------
def move_input_geom(folder, xyz_file, optional_file=None):
   """
   Moves the input geometry file into the scratch folder.

   Args:
       folder (str): The folder where the input file is located.
       xyz_file (str): The name of the input XYZ file to be moved.
       optional_file (str, optional): An additional file to move if provided.

   Returns:
       None: The function executes shell commands to move files.

   Notes:
       - Uses `os.system` to copy files.
       - If an optional file is provided, it is also copied.
   """

   file_path_xyz = os.path.join(os.path.dirname(__file__), folder, xyz_file)

   os.system(f'cp {file_path_xyz} .')

   if optional_file:
      file_path_optional = os.path.join(os.path.dirname(__file__), folder, optional_file)
      os.system(f'cp {file_path_optional} .')
# -------------------------------------------------------------------------------------
def move_managed_geom(folder, remove_optional_file=None):
   """
   Moves managed geometry files to the specified folder and removes temporary files.

   Args:
       folder (str): The destination folder for the created geometry files.
       remove_optional_file (str, optional): A specific file to be removed if provided.

   Returns:
       None: The function executes shell commands to move and remove files.

   Notes:
       - Moves generated `.xyz` files from `results_geom` to the specified folder.
       - Deletes temporary geometry files and the `results_geom` directory.
   """

   file_path = os.path.join(os.path.dirname(__file__))
   file_path_results = os.path.join(file_path, "results_geom")
   file_path_results_xyz  = os.path.join(file_path_results, "*xyz")

   os.system(f'mv {file_path_results_xyz} {folder}')
   os.system(f'rm -rf {file_path}/*xyz {file_path_results} ')

   if remove_optional_file: os.system(f'rm -rf {remove_optional_file}')
# -------------------------------------------------------------------------------------
def move_created_geom(folder):
   """
   Moves the generated geometry files to the specified folder and removes temporary files.

   Args:
       folder (str): The destination folder for the created geometry files.

   Returns:
       None: The function executes shell commands to move and remove files.

   Notes:
       - Moves all generated files from `results_geom` to the specified folder.
       - Deletes the `results_geom` directory afterward.
   """


   file_path_results = os.path.join(os.path.dirname(__file__), "results_geom")
   file_path_results = os.path.join(file_path_results, "*")

   os.system(f'mv {file_path_results} {folder}')
   os.system(f'rm -rf {file_path_results}')
# -------------------------------------------------------------------------------------
def compare_mesh_files(file1: str, file2: str, tol: float = 1e-10) -> bool:
   """
   Compare two .msh files while allowing small floating-point differences.

   This function reads two mesh files and compares them line by line. It ignores 
   non-numeric headers (lines starting with "$") and allows numerical differences 
   within a specified tolerance (`tol`). Non-numeric values (e.g., node indices) 
   are compared strictly.

   Args:
       file1 (str): Path to the first .msh file (generated mesh).
       file2 (str): Path to the second .msh file (reference mesh).
       tol (float, optional): Absolute tolerance for floating-point comparisons. 
           Defaults to 1e-10.

   Returns:
       bool: True if the files match within the tolerance, False otherwise.

   Raises:
       AssertionError: If the files differ in length, structure, or values 
       outside the specified tolerance.
   """
   
   with open(file1, 'r') as f1, open(file2, 'r') as f2:
       lines1 = f1.readlines()
       lines2 = f2.readlines()

   assert len(lines1) == len(lines2), f"Files have different number of lines: {len(lines1)} vs {len(lines2)}"

   for i, (line1, line2) in enumerate(zip(lines1, lines2)):
       # Ignore non-numeric headers that start with "$"
       if line1.startswith("$") or line2.startswith("$"):
           assert line1.strip() == line2.strip(), f"Header mismatch at line {i+1}: {line1.strip()} != {line2.strip()}"
           continue

       tokens1 = line1.strip().split()
       tokens2 = line2.strip().split()

       assert len(tokens1) == len(tokens2), f"Line {i+1}: Different number of values"

       for j, (token1, token2) in enumerate(zip(tokens1, tokens2)):
           try:
               num1 = float(token1)
               num2 = float(token2)
               assert np.isclose(num1, num2, atol=tol), \
                   f"Line {i+1}, Value {j+1}: {num1} != {num2} (tolerance {tol})"
           except ValueError:
               # If not a float, compare as string (e.g., node indices)
               assert token1 == token2, f"Line {i+1}: {token1} != {token2}"

   return True  # Files match within tolerance
# -------------------------------------------------------------------------------------
def test_create_sphere(monkeypatch):
   """
   Tests the generation of a spherical geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of generated output.

   Notes:
       - Mocks command-line arguments for sphere creation.
       - Runs the geometry creation process.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'sphere'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-sphere", "Ag", "20.0", "0.0", "0.0", "0.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "sphere_ag_r_20.0.xyz")
   generated_file = f"{test_folder}/sphere_{inp.atomtype}_r_{inp.radius}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
##def test_create_sphere_continuum(monkeypatch):
##   """
##   Tests the generation of a sphere 3D continuum mesh and compares it with a reference file.
##
##   Args:
##       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.
##
##   Returns:
##       None: Uses assertions to validate the generated `.msh` file.
##
##   Notes:
##       - Mocks command-line arguments for sphere continuum creation.
##       - Runs the geometry creation process.
##       - Compares the generated `.msh` file with an expected reference file.
##   """
##
##   # Test folder
##   test_folder = 'sphere_continuum'
##   
##   # Mock sys.argv to simulate the command line input
##   mock_args = ["dummy", "-create", "-sphere", "-continuum", "50.0", "10.0"]
##   monkeypatch.setattr(sys, "argv", mock_args)
##   
##   # Manually create and populate the input class
##   inp = input_class.input_class()
##   general.read_command_line(sys.argv, inp)
##   
##   # Run the geometry creation
##   create_geom.select_case(inp)
##   
##   # Define the expected and actual output files
##   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "sphere_r_50.0_mesh_size_10.0.msh")
##   generated_file = f"{test_folder}/sphere_r_{inp.radius}_mesh_size_{inp.mesh_size}.msh"
##
##   move_created_geom(test_folder)
##   
##   # Compare the generated file with the reference
##   assert compare_mesh_files(generated_file, expected_file), "Generated Mesh file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_sphere_core_shell(monkeypatch):
   """
   Tests the generation of a sphere core-shell geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for core-shell sphere creation.
       - Runs the geometry creation process.
       - Compares the generated `.xyz` file with an expected reference file.
   """


   # Test folder
   test_folder = 'sphere_core_shell'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-sphere", "-core", "au", "20.0", "-shell", "ag", "30"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "sphere_core_au_r_20.0_shell_ag_r_30.0.xyz")
   generated_file = f"{test_folder}/sphere_core_{inp.atomtype_in}_r_{inp.radius_in}_shell_{inp.atomtype_out}_r_{inp.radius_out}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_rod(monkeypatch):
   """
   Tests the generation of a rod geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for rod creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'rod'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-rod', 'Ag', 'X', '50.0', '20.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "rod_ag_X_l_50.0_w_20.0.xyz")
   generated_file = f"{test_folder}/rod_{inp.atomtype}_{inp.main_axis.upper()}_l_{inp.rod_length}_w_{inp.rod_width}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
##def test_create_rod_continuum(monkeypatch):
##   """
##   Tests the generation of a rod 3D continuum mesh and compares it with a reference file.
##
##   Args:
##       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.
##
##   Returns:
##       None: Uses assertions to validate the generated `.msh` file.
##
##   Notes:
##       - Mocks command-line arguments for rod continuum creation.
##       - Runs the geometry creation process.
##       - Compares the generated `.msh` file with an expected reference file.
##   """
##
##   # Test folder
##   test_folder = 'rod_continuum'
##   
##   # Mock sys.argv to simulate the command line input
##   mock_args = ["dummy", "-create", "-rod", "-continuum", "Z", "100.0", "30.0", "5.0"]
##   monkeypatch.setattr(sys, "argv", mock_args)
##   
##   # Manually create and populate the input class
##   inp = input_class.input_class()
##   general.read_command_line(sys.argv, inp)
##   
##   # Run the geometry creation
##   create_geom.select_case(inp)
##   
##   # Define the expected and actual output files
##   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "rod_Z_l_100.0_w_30.0_mesh_size_5.0.msh")
##   generated_file = f"{test_folder}/rod_{inp.main_axis.upper()}_l_{inp.rod_length}_w_{inp.rod_width}_mesh_size_{inp.mesh_size}.msh"
##
##   move_created_geom(test_folder)
##   
##   # Compare the generated file with the reference
##   assert compare_mesh_files(generated_file, expected_file), "Generated Mesh file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_rod_core_shell(monkeypatch):
   """
   Tests the generation of a rod core-shell geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for core-shell rod creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'rod_core_shell'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-rod", 'x', '-core', 'au', '20.0', '10.0', '-shell', 'ag', '50.0', '20.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "rod_X_core_au_l_20.0_r_10.0_shell_ag_l_50.0_r_20.0_shell_ag.xyz")
   generated_file = f"{test_folder}/rod_{inp.main_axis.upper()}_core_{inp.atomtype_in}_l_{inp.rod_length_in}_r_{inp.rod_width_in}_shell_{inp.atomtype_out}_l_{inp.rod_length_out}_r_{inp.rod_width_out}_shell_{inp.atomtype_out}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_tip(monkeypatch):
   """
   Tests the generation of a tip geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for tip creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'tip'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-tip', 'Ag', '50.0', '0.02', '0.02']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "tip_ag_elliptic_paraboloid_a-0.02_b-0.02_zmin-0.0_zmax-50.0.xyz")
   generated_file = f"{test_folder}/tip_{inp.atomtype}_elliptic_paraboloid_a-{inp.elliptic_parabola_a}_b-{inp.elliptic_parabola_b}_zmin-{inp.z_min}_zmax-{inp.z_max}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_pyramid(monkeypatch):
   """
   Tests the generation of a pyramid geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for pyramid creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'pyramid'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-pyramid', 'Ag', '50.0', '30.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "pyramid_ag_length-30.0_zmin-0.0_zmax-50.0.xyz")
   generated_file = f"{test_folder}/pyramid_{inp.atomtype}_length-{inp.side_length}_zmin-{inp.z_min}_zmax-{inp.z_max}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_cone(monkeypatch):
   """
   Tests the generation of a cone geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for cone creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'cone'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-cone', 'Ag', '50.0', '30.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "cone_ag_radius-30.0_zmin-0.0_zmax-50.0.xyz")
   generated_file = f"{test_folder}/cone_{inp.atomtype}_radius-{inp.radius}_zmin-{inp.z_min}_zmax-{inp.z_max}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_icosahedra(monkeypatch):
   """
   Tests the generation of an icosahedral geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for icosahedral geometry creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'icosahedron'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-ico', 'au', '50.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "icosahedron_au_r_50.0.xyz")
   generated_file = f"{test_folder}/icosahedron_{inp.atomtype}_r_{inp.radius}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_cuboctahedra(monkeypatch):
   """
   Tests the generation of a cuboctahedral geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for cuboctahedral geometry creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'cuboctahedron'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-cto', 'au', '50.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "cuboctahedron_au_r_50.0.xyz")
   generated_file = f"{test_folder}/cuboctahedron_{inp.atomtype}_r_{inp.radius}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_decahedra(monkeypatch):
   """
   Tests the generation of a decahedral geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for decahedral geometry creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'decahedron'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-idh', 'ag', '50.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "decahedron_ag_r_50.0.xyz")
   generated_file = f"{test_folder}/decahedron_{inp.atomtype}_r_{inp.radius}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_microscope(monkeypatch):
   """
   Tests the generation of a microscope geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for microscope geometry creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'microscope'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", '-create', '-microscope', 'ag', '40.0', '0.02', '0.02', '26.0', '33.0']
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "microscope_ag_parabola_40.0_0.02_0.02_pyramid_26.0_33.0.xyz")
   generated_file = f"{test_folder}/microscope_{inp.atomtype}_parabola_{inp.z_max_paraboloid}_{inp.elliptic_parabola_a}_{inp.elliptic_parabola_b}_pyramid_{inp.z_max_pyramid}_{inp.side_length}{inp.alloy_string}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_graphene_ribbon(monkeypatch):
   """
   Tests the generation of a graphene ribbon geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for graphene ribbon creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'graphene_ribbon'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-graphene", "rib", "40.0", "20.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "graphene_ribbon_40.0_20.0.xyz")
   generated_file = f"{test_folder}/graphene_ribbon_{inp.X_length}_{inp.Y_length}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_graphene_disk(monkeypatch):
   """
   Tests the generation of a graphene disk geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for graphene disk creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'graphene_disk'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-graphene", "disk", "30.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "graphene_disk_30.0.xyz")
   generated_file = f"{test_folder}/graphene_disk_{inp.radius}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_graphene_ring(monkeypatch):
   """
   Tests the generation of a graphene ring geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for graphene ring creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'graphene_ring'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-graphene", "ring", "60.0", "30.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "graphene_ring_out_60.0_in_30.0.xyz")
   generated_file = f"{test_folder}/graphene_ring_out_{inp.radius_out}_in_{inp.radius_in}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_graphene_triangle_armchair(monkeypatch):
   """
   Tests the generation of a graphene triangle (armchair type) geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for graphene armchair triangle creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'graphene_triangle_armchair'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-graphene", "triangle", "armchair", "50.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "graphene_triangle_armchair_50.0.xyz")
   generated_file = f"{test_folder}/graphene_triangle_{inp.graphene_edge_type}_{inp.side_length}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_graphene_triangle_zigzag(monkeypatch):
   """
   Tests the generation of a graphene triangle (zigzag type) geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for graphene zigzag triangle creation.
       - Runs the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder = 'graphene_triangle_zigzag'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-graphene", "triangle", "zigzag", "50.0"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "graphene_triangle_zigzag_50.0.xyz")
   generated_file = f"{test_folder}/graphene_triangle_{inp.graphene_edge_type}_{inp.side_length}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_specular_geometry(monkeypatch):
   """
   Tests the creation of a specular (mirrored) geometry and compares it with a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for mirroring a geometry.
       - Moves the input file temporarily.
       - Runs the mirroring process.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder    = 'specular'
   xyz_input_file = 'doxorubicin.xyz'

   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-mirror", xyz_input_file]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)

   # Temporaly move input file
   move_input_geom(test_folder,xyz_input_file)

   # Create specular geometry
   various.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "doxorubicin_000_mirror.xyz")
   generated_file = f"{test_folder}/{inp.geom_file[:-4]}_000_mirror.xyz"

   move_managed_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_controlled_distance(monkeypatch):
   """
   Tests the controlled translation of a molecule to maintain a specific distance.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for controlled distance translation.
       - Moves the input file temporarily.
       - Runs the translation process.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder      = 'control_distance'
   xyz_input_file_1 = 'doxorubicin.xyz'
   xyz_input_file_2 = 'sphere_r_10.0_center_0.0_0.0_0.0.xyz'
   distances_input  = 'distances_input'

   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-t", distances_input, "doxorubicin.xyz", "no", "sphere_r_10.0_center_0.0_0.0_0.0.xyz", "no", "+x", "verbose_no"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)

   # Temporaly move input file
   move_input_geom(test_folder,xyz_input_file_1, optional_file = distances_input)
   move_input_geom(test_folder,xyz_input_file_2)

   # Translate controlled distance
   translate.select_case(inp)

   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "sphere_r_10.0_center_0.0_0.0_0.0_+x_d_10.00.xyz")
   generated_file = f"{test_folder}/sphere_r_10.0_center_0.0_0.0_0.0_+x_d_10.00.xyz"

   move_managed_geom(test_folder, remove_optional_file = distances_input)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_controlled_rotation(monkeypatch):
   """
   Tests the controlled rotation of a molecule along a specified axis.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv`.

   Returns:
       None: Uses assertions to validate the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for controlled rotation.
       - Moves the input file temporarily.
       - Runs the rotation process.
       - Compares the generated `.xyz` file with an expected reference file.
   """

   # Test folder
   test_folder    = 'control_rotation'
   xyz_input_file = 'doxorubicin.xyz'
   angles_input = 'angles_input'

   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-r", angles_input, "doxorubicin.xyz", "no", "+x"]
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)

   # Temporaly move input file
   move_input_geom(test_folder,xyz_input_file, optional_file = angles_input)

   # Translate controlled distance
   rotate.select_case(inp)

   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "doxorubicin_+x_degree_90.0.xyz")
   generated_file = f"{test_folder}/doxorubicin_+x_degree_90.0.xyz"

   move_managed_geom(test_folder, remove_optional_file = angles_input)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_dimer(monkeypatch):
   """
   Tests the creation of a dimer geometry and validates it against a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv` for command-line argument simulation.

   Returns:
       None: Uses assertions to verify the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for dimer creation with a conical base structure.
       - Initializes and populates the input class.
       - Executes the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file against an expected reference file.
   """

   # Test folder
   test_folder = 'dimer'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-cone", "ag", "30.0", "40.0", "-dimer", "10.0", "+z"] 
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "dimer_cone_ag_radius-40.0_zmin-0.0_zmax-30.0_+z_d_10.0.xyz")
   generated_file = f"{test_folder}/dimer_{inp.xyz_output}_{inp.dir_axis_input}_d_{inp.distances[0]}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
def test_create_bowtie(monkeypatch):
   """
   Tests the creation of a bowtie geometry and validates it against a reference file.

   Args:
       monkeypatch (pytest.MonkeyPatch): A fixture to modify `sys.argv` for command-line argument simulation.

   Returns:
       None: Uses assertions to verify the correctness of the generated `.xyz` file.

   Notes:
       - Mocks command-line arguments for bowtie creation with a conical base structure.
       - Initializes and populates the input class.
       - Executes the geometry creation process.
       - Moves the created geometry files to the test folder.
       - Compares the generated `.xyz` file against an expected reference file.
   """

   # Test folder
   test_folder = 'bowtie'
   
   # Mock sys.argv to simulate the command line input
   mock_args = ["dummy", "-create", "-cone", "ag", "30.0", "40.0", "-bowtie", "10.0"] 
   monkeypatch.setattr(sys, "argv", mock_args)
   
   # Manually create and populate the input class
   inp = input_class.input_class()
   general.read_command_line(sys.argv, inp)
   
   # Run the geometry creation
   create_geom.select_case(inp)
   
   # Define the expected and actual output files
   expected_file = os.path.join(os.path.dirname(__file__), test_folder, "reference", "bowtie_cone_ag_radius-40.0_zmin-0.0_zmax-30.0_-z_d_10.0.xyz")
   generated_file = f"{test_folder}/bowtie_{inp.xyz_output}_{inp.dir_axis_input}_d_{inp.distances[0]}.xyz"

   move_created_geom(test_folder)
   
   # Compare the generated file with the reference
   assert filecmp.cmp(generated_file, expected_file, shallow=False), "Generated XYZ file does not match the expected output"
# -------------------------------------------------------------------------------------
